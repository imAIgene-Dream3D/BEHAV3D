"""
BEHAV3D napari plugin – Processing Queue.

Provides a collapsible bottom panel for queueing Train → Segment → Track
steps with preset workflows and sequential execution.
"""
import sys
import time
import traceback
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QComboBox, QFrame, QSizePolicy, QMessageBox, QScrollArea,
)
from qtpy.QtCore import Qt, Signal, QObject, QTimer
from qtpy.QtGui import QFont


# ═══════════════════════════════════════════════════════════════════════════
# Data model
# ═══════════════════════════════════════════════════════════════════════════

class StepType(Enum):
    TRAIN = "train"
    SEGMENT = "segment"
    TRACK = "track"

    @property
    def label(self):
        return {
            StepType.TRAIN: "🧠 Train Classifier",
            StepType.SEGMENT: "🦠 Batch Segmentation",
            StepType.TRACK: "📍 Batch Tracking",
        }[self]

    @property
    def order(self):
        return {StepType.TRAIN: 0, StepType.SEGMENT: 1, StepType.TRACK: 2}[self]


class StepStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    DONE = "done"
    ERROR = "error"


@dataclass
class QueueStep:
    step_type: StepType
    status: StepStatus = StepStatus.PENDING
    error_msg: str = ""
    elapsed: float = 0.0


# ═══════════════════════════════════════════════════════════════════════════
# Queue step row widget
# ═══════════════════════════════════════════════════════════════════════════

class QueueStepRow(QFrame):
    """One row in the queue list, showing status icon + name + delete button."""

    remove_requested = Signal()

    _STATUS_ICONS = {
        StepStatus.PENDING: "⏸️",
        StepStatus.RUNNING: "⏳",
        StepStatus.DONE: "✅",
        StepStatus.ERROR: "❌",
    }

    def __init__(self, step: QueueStep, parent=None):
        super().__init__(parent)
        self.step = step
        self.setFixedHeight(24)
        self.setStyleSheet(
            "QFrame { background: #2b2b2b; border: 1px solid #444; "
            "border-radius: 2px; padding: 0; margin: 0; }"
        )
        row = QHBoxLayout(self)
        row.setContentsMargins(4, 0, 4, 0)
        row.setSpacing(4)

        self.status_label = QLabel(self._STATUS_ICONS[step.status])
        self.status_label.setFixedWidth(28)
        row.addWidget(self.status_label)

        self.name_label = QLabel(step.step_type.label)
        self.name_label.setStyleSheet("color: #ddd; font-size: 10px;")
        row.addWidget(self.name_label, stretch=1)

        self.time_label = QLabel("")
        self.time_label.setStyleSheet("color: #888; font-size: 9px;")
        self.time_label.setFixedWidth(45)
        row.addWidget(self.time_label)

        self.btn_remove = QPushButton("x")
        self.btn_remove.setFixedSize(16, 16)
        self.btn_remove.setStyleSheet(
            "QPushButton { background: transparent; border: none; font-size: 10px; color: #888; }"
            "QPushButton:hover { background: #dc3545; color: white; border-radius: 2px; }"
        )
        self.btn_remove.setToolTip("Remove from queue")
        self.btn_remove.clicked.connect(self.remove_requested.emit)
        row.addWidget(self.btn_remove)

    def refresh(self):
        """Update visual state from the step data."""
        self.status_label.setText(self._STATUS_ICONS[self.step.status])
        if self.step.status == StepStatus.DONE:
            secs = self.step.elapsed
            if secs >= 60:
                self.time_label.setText(f"{secs/60:.1f}m")
            else:
                self.time_label.setText(f"{secs:.0f}s")
        elif self.step.status == StepStatus.ERROR:
            self.time_label.setText("failed")
            self.time_label.setStyleSheet("color: #dc3545; font-size: 10px;")
        # Disable delete while running/done
        self.btn_remove.setEnabled(self.step.status == StepStatus.PENDING)


# ═══════════════════════════════════════════════════════════════════════════
# Main Queue Panel
# ═══════════════════════════════════════════════════════════════════════════

class ProcessingQueuePanel(QWidget):
    """Collapsible bottom panel for managing the processing queue."""

    # Emitted when queue finishes (success or with errors)
    queue_finished = Signal()

    PRESETS = {
        "Train + Segment + Track": [StepType.TRAIN, StepType.SEGMENT, StepType.TRACK],
        "Segment + Track": [StepType.SEGMENT, StepType.TRACK],
    }

    def __init__(self, segmentation_tab, tracking_tab, metadata_loader, parent=None):
        super().__init__(parent)
        self.segmentation_tab = segmentation_tab
        self.tracking_tab = tracking_tab
        self.metadata_loader = metadata_loader

        self._steps: List[QueueStep] = []
        self._step_rows: List[QueueStepRow] = []
        self._is_running = False

        self._init_ui()

    # ── UI ─────────────────────────────────────────────────────────────────

    def _init_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # ── Toggle bar ───────────────────────────────────────────────
        self.toggle_bar = QFrame()
        self.toggle_bar.setStyleSheet(
            "QFrame { background: #1a1a2e; border-top: 2px solid #16213e; }"
        )
        self.toggle_bar.setCursor(Qt.PointingHandCursor)
        bar_layout = QHBoxLayout(self.toggle_bar)
        bar_layout.setContentsMargins(10, 5, 10, 5)

        self.badge_label = QLabel("🛒 Processing Queue (0)")
        self.badge_label.setStyleSheet(
            "color: #e0e0e0; font-size: 12px; font-weight: bold;"
        )
        bar_layout.addWidget(self.badge_label)
        bar_layout.addStretch()

        self.toggle_btn = QPushButton("▲")
        self.toggle_btn.setFixedSize(24, 24)
        self.toggle_btn.setStyleSheet(
            "QPushButton { background: #16213e; color: #aaa; border: none; "
            "border-radius: 4px; font-size: 12px; }"
            "QPushButton:hover { background: #0f3460; color: white; }"
        )
        bar_layout.addWidget(self.toggle_btn)

        main_layout.addWidget(self.toggle_bar)

        # ── Expandable body ──────────────────────────────────────────
        self.body = QWidget()
        body_layout = QVBoxLayout(self.body)
        body_layout.setContentsMargins(8, 2, 8, 4)
        body_layout.setSpacing(4)

        # Preset row
        preset_row = QHBoxLayout()
        preset_row.setSpacing(6)
        preset_label = QLabel("Preset:")
        preset_label.setStyleSheet("color: #aaa; font-size: 11px;")
        preset_row.addWidget(preset_label)

        self.preset_combo = QComboBox()
        self.preset_combo.addItems(["— Select preset —"] + list(self.PRESETS.keys()))
        self.preset_combo.setStyleSheet("font-size: 11px;")
        self.preset_combo.currentTextChanged.connect(self._on_preset_selected)
        preset_row.addWidget(self.preset_combo, stretch=1)

        self.btn_clear = QPushButton("Clear")
        self.btn_clear.setFixedHeight(24)
        self.btn_clear.setStyleSheet(
            "QPushButton { background: #6c757d; color: white; border-radius: 3px; "
            "font-size: 10px; padding: 2px 8px; }"
            "QPushButton:hover { background: #5a6268; }"
        )
        self.btn_clear.clicked.connect(self.clear_queue)
        preset_row.addWidget(self.btn_clear)
        body_layout.addLayout(preset_row)

        # Steps list
        self.steps_container = QWidget()
        self.steps_layout = QVBoxLayout(self.steps_container)
        self.steps_layout.setContentsMargins(0, 0, 0, 0)
        self.steps_layout.setSpacing(1)

        self.empty_label = QLabel("Queue is empty. Add steps using the 🛒+ buttons or select a preset.")
        self.empty_label.setStyleSheet("color: #666; font-style: italic; font-size: 10px;")
        self.empty_label.setAlignment(Qt.AlignCenter)
        self.steps_layout.addWidget(self.empty_label)

        body_layout.addWidget(self.steps_container)

        # Run button
        self.btn_run = QPushButton("▶  Run Queue")
        self.btn_run.setStyleSheet(
            "QPushButton { background: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 4px; font-size: 11px; }"
            "QPushButton:hover { background: #218838; }"
            "QPushButton:disabled { background: #555; color: #999; }"
        )
        self.btn_run.clicked.connect(self._on_run_queue)
        body_layout.addWidget(self.btn_run)

        main_layout.addWidget(self.body)
        self.body.setVisible(False)  # Start collapsed
        self.body.setMaximumHeight(250)  # Prevent inflating the Napari window

        # Don't let the queue panel push the window bigger when collapsed
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        # Toggle click
        self.toggle_btn.clicked.connect(self._toggle_body)
        self.toggle_bar.mousePressEvent = lambda e: self._toggle_body()

    # ── Toggle ─────────────────────────────────────────────────────────

    def _toggle_body(self):
        visible = not self.body.isVisible()
        self.body.setVisible(visible)
        self.toggle_btn.setText("▼" if visible else "▲")

    # ── Badge ──────────────────────────────────────────────────────────

    def _update_badge(self):
        n = len(self._steps)
        self.badge_label.setText(f"🛒 Processing Queue ({n})")
        if n > 0:
            self.badge_label.setStyleSheet(
                "color: #ffc107; font-size: 12px; font-weight: bold;"
            )
        else:
            self.badge_label.setStyleSheet(
                "color: #e0e0e0; font-size: 12px; font-weight: bold;"
            )

    # ── Adding / Removing steps ────────────────────────────────────────

    def add_step(self, step_type: StepType):
        """Add a step to the queue (maintains order, avoids duplicates)."""
        if self._is_running:
            return
        # No duplicates
        if any(s.step_type == step_type for s in self._steps):
            return
        step = QueueStep(step_type=step_type)
        self._steps.append(step)
        # Train requires Segment — auto-add it
        if step_type == StepType.TRAIN:
            if not any(s.step_type == StepType.SEGMENT for s in self._steps):
                self._steps.append(QueueStep(step_type=StepType.SEGMENT))
        self._steps.sort(key=lambda s: s.step_type.order)
        self._rebuild_list()
        # Auto-expand when adding
        if not self.body.isVisible():
            self._toggle_body()
        # Reset preset combo to custom if it doesn't match
        self.preset_combo.blockSignals(True)
        self.preset_combo.setCurrentIndex(0)
        self.preset_combo.blockSignals(False)

    def remove_step(self, step_type: StepType):
        """Remove a step from the queue (deferred to avoid crash)."""
        if self._is_running:
            return
        # Defer so the signal-emitting row widget fully returns before we touch it
        QTimer.singleShot(0, lambda: self._do_remove_step(step_type))

    def _do_remove_step(self, step_type: StepType):
        """Actually remove the step — runs after the signal handler has returned."""
        # If removing Segment, also remove Train (Train without Segment is invalid)
        types_to_remove = {step_type}
        if step_type == StepType.SEGMENT:
            types_to_remove.add(StepType.TRAIN)

        indices_to_remove = []
        for i, step in enumerate(self._steps):
            if step.step_type in types_to_remove:
                indices_to_remove.append(i)
        # Remove in reverse order to preserve indices
        for i in reversed(indices_to_remove):
            self._steps.pop(i)
            row = self._step_rows.pop(i)
            self.steps_layout.removeWidget(row)
            row.hide()
            row.setParent(None)
        self.empty_label.setVisible(len(self._steps) == 0)
        self._update_badge()
        self.btn_run.setEnabled(len(self._steps) > 0 and not self._is_running)

    def clear_queue(self):
        """Remove all steps."""
        if self._is_running:
            return
        self._steps.clear()
        self._detach_all_rows()
        self.preset_combo.blockSignals(True)
        self.preset_combo.setCurrentIndex(0)
        self.preset_combo.blockSignals(False)

    def has_step(self, step_type: StepType) -> bool:
        return any(s.step_type == step_type for s in self._steps)

    def _detach_all_rows(self):
        """Safely remove all row widgets without deleteLater."""
        for row in self._step_rows:
            self.steps_layout.removeWidget(row)
            row.hide()
            row.setParent(None)
        self._step_rows.clear()
        self.empty_label.setVisible(True)
        self._update_badge()
        self.btn_run.setEnabled(False)

    def _rebuild_list(self):
        """Rebuild the visual step rows from the data model."""
        # Safely detach old rows
        self._detach_all_rows()

        self.empty_label.setVisible(len(self._steps) == 0)

        for step in self._steps:
            row = QueueStepRow(step)
            row.remove_requested.connect(lambda st=step.step_type: self.remove_step(st))
            self.steps_layout.addWidget(row)
            self._step_rows.append(row)

        self._update_badge()
        self.btn_run.setEnabled(len(self._steps) > 0 and not self._is_running)

    # ── Presets ────────────────────────────────────────────────────────

    def _on_preset_selected(self, text):
        if text in self.PRESETS:
            self._steps.clear()
            for st in self.PRESETS[text]:
                self._steps.append(QueueStep(step_type=st))
            self._rebuild_list()
            if not self.body.isVisible():
                self._toggle_body()

    # ── Pre-run validation ─────────────────────────────────────────────

    def _check_existing_data(self) -> List[str]:
        """Return list of descriptions of data that will be overwritten."""
        warnings = []
        md = self.metadata_loader.metadata
        if md is None:
            return warnings
        out_dir = Path(self.metadata_loader.output_dir)
        pc_dir = out_dir / "images" / "PixelClassification"

        for step in self._steps:
            if step.step_type == StepType.TRAIN:
                # Check for existing .joblib classifiers
                if pc_dir.exists():
                    jobllibs = list(pc_dir.glob("PixelClassifier_*.joblib"))
                    if jobllibs:
                        names = ", ".join(j.stem for j in jobllibs)
                        warnings.append(f"Trained classifiers: {names}")

            elif step.step_type == StepType.SEGMENT:
                # Check for existing segment zarrs
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    sample_img_dir = out_dir / "images" / sn
                    if sample_img_dir.exists():
                        seg_files = list(sample_img_dir.glob("*_segments.zarr"))
                        if seg_files:
                            warnings.append(f"Segmentation data for {sn}")
                            break

            elif step.step_type == StepType.TRACK:
                # Check for existing track files
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    sample_img_dir = out_dir / "images" / sn
                    if sample_img_dir.exists():
                        track_files = list(sample_img_dir.glob("*_tracked.zarr"))
                        if track_files:
                            warnings.append(f"Tracking data for {sn}")
                            break

        return warnings

    # ── Run queue ──────────────────────────────────────────────────────

    def _on_run_queue(self):
        if not self._steps or self._is_running:
            return

        # Validate metadata
        if self.metadata_loader.metadata is None:
            QMessageBox.warning(self, "No Metadata", "Please load metadata first.")
            return

        # Pre-run overwrite check
        existing = self._check_existing_data()
        if existing:
            details = "\n".join(f"  • {w}" for w in existing)
            res = QMessageBox.warning(
                self, "Overwrite Existing Data?",
                f"The following data will be overwritten:\n\n{details}\n\n"
                "Do you want to continue?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if res != QMessageBox.Yes:
                return

        # ── Execute ──
        self._is_running = True
        self.btn_run.setEnabled(False)
        self.btn_clear.setEnabled(False)
        self.preset_combo.setEnabled(False)
        # Disable delete buttons
        for row in self._step_rows:
            row.btn_remove.setEnabled(False)

        print("", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        print("  🛒 BEHAV3D Processing Queue — Starting", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

        from qtpy.QtWidgets import QApplication
        total_start = time.time()

        for i, step in enumerate(self._steps):
            step.status = StepStatus.RUNNING
            self._step_rows[i].refresh()
            QApplication.processEvents()

            step_start = time.time()
            print(f"\n▶ [{i+1}/{len(self._steps)}] {step.step_type.label}...", file=sys.stderr)

            try:
                self._execute_step(step)
                step.elapsed = time.time() - step_start
                step.status = StepStatus.DONE
                print(f"✅ {step.step_type.label} completed in {step.elapsed:.1f}s", file=sys.stderr)
            except Exception as e:
                step.elapsed = time.time() - step_start
                step.status = StepStatus.ERROR
                step.error_msg = str(e)
                traceback.print_exc()
                print(f"❌ {step.step_type.label} failed: {e}", file=sys.stderr)
                
                # Show error dialog and stop queue
                QMessageBox.critical(
                    self, 
                    "Queue Error", 
                    f"An error occurred during {step.step_type.label}:\n\n{e}\n\nQueue execution stopped."
                )
                
                # Stop on first error
                self._step_rows[i].refresh()
                break

            self._step_rows[i].refresh()
            QApplication.processEvents()

        total_elapsed = time.time() - total_start
        print(f"\n{'=' * 60}", file=sys.stderr)
        print(f"  🛒 Queue finished in {total_elapsed:.1f}s", file=sys.stderr)
        print(f"{'=' * 60}\n", file=sys.stderr)

        self._is_running = False
        self.btn_run.setEnabled(True)
        self.btn_clear.setEnabled(True)
        self.preset_combo.setEnabled(True)
        self.queue_finished.emit()

    def _execute_step(self, step: QueueStep):
        """Run a single step, calling the backend directly."""
        if step.step_type == StepType.TRAIN:
            self._run_train()
        elif step.step_type == StepType.SEGMENT:
            self._run_segment()
        elif step.step_type == StepType.TRACK:
            self._run_track()

    # ── Step runners (non-interactive) ─────────────────────────────────

    def _run_train(self):
        """Run classifier training (non-interactive)."""
        # Switch to segmentation tab first (important for internal widget state)
        parent = self.parent()
        while parent and not hasattr(parent, 'tabs'):
            parent = parent.parent()
        if parent and hasattr(parent, 'tabs'):
            parent.tabs.setCurrentIndex(2) # Segmentation tab

        pc_widget = self.segmentation_tab.pixel_classifier_page
        pc_widget.run_train(interactive=False)

    def _run_segment(self):
        """Run batch segmentation (non-interactive)."""
        pc_widget = self.segmentation_tab.pixel_classifier_page
        pc_widget.run_batch_segmentation(interactive=False)

    def _run_track(self):
        """Run batch tracking (non-interactive)."""
        self.tracking_tab.run_batch_tracking(interactive=False)
