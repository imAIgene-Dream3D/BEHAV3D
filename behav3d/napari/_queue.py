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
from qtpy.QtCore import Qt, Signal, Slot, QObject, QTimer
from qtpy.QtGui import QFont


# ═══════════════════════════════════════════════════════════════════════════
# Data model
# ═══════════════════════════════════════════════════════════════════════════

class StepType(Enum):
    TRAIN = "train"
    SEGMENT = "segment"
    CELLPOSE_SEGMENT = "cellpose_segment"
    DEAD_MASK = "dead_mask"
    APOC_SEGMENT = "apoc_segment"
    TRACK = "track"
    FEATURE_EXTRACT = "feature_extract"
    FILTER = "filter"

    @property
    def label(self):
        return {
            StepType.TRAIN: "🧠 Train Classifier",
            StepType.SEGMENT: "🦠 Batch Segmentation",
            StepType.CELLPOSE_SEGMENT: "🔬 Cellpose Segmentation",
            StepType.DEAD_MASK: "☠ Dead Mask (Otsu)",
            StepType.APOC_SEGMENT: "⚡ APOC",
            StepType.TRACK: "📍 Batch Tracking",
            StepType.FEATURE_EXTRACT: "🧪 Feature Extraction",
            StepType.FILTER: "🧹 Filtering",
        }[self]

    @property
    def order(self):
        return {
            StepType.TRAIN: 0,
            StepType.SEGMENT: 1,
            StepType.CELLPOSE_SEGMENT: 1.5,
            StepType.DEAD_MASK: 1.6,
            StepType.APOC_SEGMENT: 1.5,
            StepType.TRACK: 2,
            StepType.FEATURE_EXTRACT: 3,
            StepType.FILTER: 4,
        }[self]


class StepStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    DONE = "done"
    ERROR = "error"


@dataclass
class QueueStep:
    step_type: StepType
    params: dict = field(default_factory=dict)
    status: StepStatus = StepStatus.PENDING
    error_msg: str = ""
    elapsed: float = 0.0

    @property
    def display_label(self):
        """Label for queue row display. Includes cell type/strategy where relevant."""
        cell_type = self.params.get("cell_type")
        if cell_type and self.step_type == StepType.CELLPOSE_SEGMENT:
            return f"🔬 Cellpose — {cell_type}"
        if self.step_type == StepType.APOC_SEGMENT:
            strat = self.params.get("strategy_name", "")
            return f"⚡ APOC — {strat}" if strat else self.step_type.label
        return self.step_type.label


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

        self.name_label = QLabel(step.display_label)
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
        "Segment → Filter": [StepType.SEGMENT, StepType.TRACK, StepType.FEATURE_EXTRACT, StepType.FILTER],
    }

    def __init__(self, segmentation_tab, tracking_tab, metadata_loader,
                 feature_extraction_tab=None, filtering_tab=None, parent=None):
        super().__init__(parent)
        self.segmentation_tab = segmentation_tab
        self.tracking_tab = tracking_tab
        self.metadata_loader = metadata_loader
        self.feature_extraction_tab = feature_extraction_tab
        self.filtering_tab = filtering_tab

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

    @Slot(object)
    def add_step(self, step_type: StepType, params: dict = None):
        """Add a step to the queue, with conditional dependency prompting.

        For CELLPOSE_SEGMENT steps, *params* must contain at least
        ``{"cell_type": ..., "model_path": ...}``.  Multiple cellpose steps
        with different cell types are allowed; duplicates (same cell type)
        are silently ignored.
        """
        if params is None:
            params = {}

        # Dedup: for cellpose, dedup by (step_type, cell_type); others by step_type
        if step_type == StepType.CELLPOSE_SEGMENT:
            ct = params.get("cell_type", "")
            if any(
                s.step_type == StepType.CELLPOSE_SEGMENT
                and s.params.get("cell_type") == ct
                for s in self._steps
            ):
                return  # same cell type already queued
        else:
            if any(s.step_type == step_type for s in self._steps):
                return

        from qtpy.QtWidgets import QMessageBox

        # Dependency chain logic
        added_already = [s.step_type for s in self._steps]

        # 1. TRACK needs segmentation data
        if step_type == StepType.TRACK and StepType.SEGMENT not in added_already:
            missing = self._check_input_data_missing(StepType.TRACK)
            if missing:
                # Check if all missing cell types are covered by cellpose steps
                if not self._cellpose_covers_all_cell_types():
                    res = QMessageBox.question(
                        self, "Missing Segmentation Data",
                        "Segmentation data is missing for some samples/cell types.\n\n"
                        "Would you like to add Batch Segmentation to the pipeline?",
                        QMessageBox.Yes | QMessageBox.No
                    )
                    if res == QMessageBox.Yes:
                        self.add_step(StepType.SEGMENT)
                    else:
                        return  # Cancel adding Track

        # 2. FEATURE_EXTRACT needs TRACK
        elif step_type == StepType.FEATURE_EXTRACT and StepType.TRACK not in added_already:
            missing = self._check_input_data_missing(StepType.FEATURE_EXTRACT)
            if missing:
                res = QMessageBox.question(
                    self, "Missing Tracking Data",
                    "Tracking data is missing for some samples/cell types.\n\n"
                    "Would you like to add Batch Tracking to the pipeline?",
                    QMessageBox.Yes | QMessageBox.No
                )
                if res == QMessageBox.Yes:
                    self.add_step(StepType.TRACK)
                else:
                    return  # Cancel adding Feature Extract

        # 3. FILTER needs FEATURE_EXTRACT
        elif step_type == StepType.FILTER and StepType.FEATURE_EXTRACT not in added_already:
            missing = self._check_input_data_missing(StepType.FILTER)
            if missing:
                res = QMessageBox.question(
                    self, "Missing Feature Data",
                    "Feature extraction data is missing for some cell types.\n\n"
                    "Would you like to add Feature Extraction to the pipeline?",
                    QMessageBox.Yes | QMessageBox.No
                )
                if res == QMessageBox.Yes:
                    self.add_step(StepType.FEATURE_EXTRACT)
                else:
                    return  # Cancel adding Filter

        # Finally add the requested step
        step = QueueStep(step_type=step_type, params=params)
        self._steps.append(step)

        # Train still auto-adds Segment
        if step_type == StepType.TRAIN:
            if not any(s.step_type == StepType.SEGMENT for s in self._steps):
                self._steps.append(QueueStep(step_type=StepType.SEGMENT))

        self._steps.sort(key=lambda s: s.step_type.order)
        self._rebuild_list()

        # Auto-expand when adding
        if not self.body.isVisible():
            self._toggle_body()

        # Reset preset combo to custom (0)
        self.preset_combo.blockSignals(True)
        self.preset_combo.setCurrentIndex(0)
        self.preset_combo.blockSignals(False)

    def _cellpose_covers_all_cell_types(self) -> bool:
        """Return True if every cell type that needs segmentation is covered
        by a CELLPOSE_SEGMENT step already in the queue (or has data on disk)."""
        md = self.metadata_loader.metadata
        if md is None:
            return False
        out_dir = Path(self.metadata_loader.output_dir)

        all_cts = []
        if self.tracking_tab:
            all_cts = list(self.tracking_tab.panels.keys())
        if not all_cts:
            return False

        queued_cellpose_cts = {
            s.params.get("cell_type")
            for s in self._steps
            if s.step_type == StepType.CELLPOSE_SEGMENT
        }

        for ct in all_cts:
            # Check if data already exists on disk
            all_exist = True
            for _, sample in md.iterrows():
                sn = sample.get("sample_name", "unknown")
                seg_zarr = out_dir / "images" / sn / f"{sn}_{ct}_segments.zarr"
                if not seg_zarr.exists():
                    all_exist = False
                    break
            if all_exist:
                continue  # this cell type is covered by existing data
            # Not on disk — check if a cellpose step covers it
            if ct not in queued_cellpose_cts:
                return False
        return True

    def remove_step(self, step):
        """Remove a step from the queue (deferred to avoid crash).

        *step* can be a StepType (removes first matching) or a QueueStep
        instance (removes that exact step).  For SEGMENT, also removes TRAIN.
        """
        if self._is_running:
            return
        QTimer.singleShot(0, lambda: self._do_remove_step(step))

    def _do_remove_step(self, step):
        """Actually remove the step — runs after the signal handler has returned."""
        if isinstance(step, QueueStep):
            # Identity-based removal (used for cellpose with per-cell-type steps)
            indices_to_remove = []
            for i, s in enumerate(self._steps):
                if s is step:
                    indices_to_remove.append(i)
                    # If removing Segment, also remove Train
                    if s.step_type == StepType.SEGMENT:
                        for j, s2 in enumerate(self._steps):
                            if s2.step_type == StepType.TRAIN and j not in indices_to_remove:
                                indices_to_remove.append(j)
                    break
        else:
            # Legacy: step is a StepType enum
            step_type = step
            types_to_remove = {step_type}
            if step_type == StepType.SEGMENT:
                types_to_remove.add(StepType.TRAIN)
            indices_to_remove = [
                i for i, s in enumerate(self._steps)
                if s.step_type in types_to_remove
            ]

        for i in reversed(sorted(indices_to_remove)):
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
            row.remove_requested.connect(lambda _s=step: self.remove_step(_s))
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

    def _check_input_data_missing(self, step_type: StepType) -> bool:
        """Check if input data required for a step is missing."""
        md = self.metadata_loader.metadata
        if md is None:
            return True
        out_dir = Path(self.metadata_loader.output_dir)

        if step_type == StepType.TRACK:
            # Requires segmented zarrs for all types & samples
            all_cts = []
            if self.tracking_tab:
                all_cts = list(self.tracking_tab.panels.keys())
            if not all_cts:
                return False # No types configured yet
            
            for ct in all_cts:
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    seg_zarr = out_dir / "images" / sn / f"{sn}_{ct}_segments.zarr"
                    if not seg_zarr.exists():
                        return True
            return False

        elif step_type == StepType.FEATURE_EXTRACT:
            # Requires tracked zarrs/csvs
            all_cts = []
            if self.feature_extraction_tab:
                all_cts = list(self.feature_extraction_tab.panels.keys())
            if not all_cts:
                return False
                
            for ct in all_cts:
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    track_zarr = out_dir / "images" / sn / f"{sn}_{ct}_tracked.zarr"
                    # Also check csv dir
                    csv_dir = out_dir / "trackdata" / sn / ct
                    if not track_zarr.exists() and not (csv_dir.exists() and list(csv_dir.glob("*.csv"))):
                        return True
            return False

        elif step_type == StepType.FILTER:
            # Requires combined feature CSV
            all_cts = []
            if self.filtering_tab:
                all_cts = list(self.filtering_tab.panels.keys())
            if not all_cts:
                return False
                
            for ct in all_cts:
                feat_csv = out_dir / "analysis" / ct / "track_features" / f"BEHAV3D_{ct}_combined_track_features.csv"
                if not feat_csv.exists():
                    return True
            return False

        return False

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
                if pc_dir.exists():
                    jobllibs = list(pc_dir.glob("PixelClassifier_*.joblib"))
                    if jobllibs:
                        names = ", ".join(j.stem for j in jobllibs)
                        warnings.append(f"Trained classifiers: {names}")

            elif step.step_type == StepType.SEGMENT:
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    sample_img_dir = out_dir / "images" / sn
                    if sample_img_dir.exists():
                        seg_files = list(sample_img_dir.glob("*_segments.zarr"))
                        if seg_files:
                            warnings.append(f"Segmentation data for {sn}")
                            break

            elif step.step_type == StepType.CELLPOSE_SEGMENT:
                ct = step.params.get("cell_type", "")
                if ct:
                    for _, sample in md.iterrows():
                        sn = sample.get("sample_name", "unknown")
                        seg_zarr = out_dir / "images" / sn / f"{sn}_{ct}_segments.zarr"
                        if seg_zarr.exists():
                            warnings.append(f"Cellpose segmentation for {ct} ({sn})")
                            break

            elif step.step_type == StepType.DEAD_MASK:
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    dead_zarr = out_dir / "images" / sn / f"{sn}_dead_segments.zarr"
                    if dead_zarr.exists():
                        warnings.append(f"Dead mask for {sn}")
                        break

            elif step.step_type == StepType.TRACK:
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    sample_img_dir = out_dir / "images" / sn
                    if sample_img_dir.exists():
                        track_files = list(sample_img_dir.glob("*_tracked.zarr"))
                        if track_files:
                            warnings.append(f"Tracking data for {sn}")
                            break

            elif step.step_type == StepType.FEATURE_EXTRACT:
                analysis_dir = out_dir / "analysis"
                if analysis_dir.exists():
                    for ct_dir in analysis_dir.iterdir():
                        feat_csv = ct_dir / "track_features" / f"BEHAV3D_{ct_dir.name}_combined_track_features.csv"
                        if feat_csv.exists():
                            warnings.append(f"Feature data for {ct_dir.name}")
                            break

            elif step.step_type == StepType.FILTER:
                analysis_dir = out_dir / "analysis"
                if analysis_dir.exists():
                    for ct_dir in analysis_dir.iterdir():
                        filt_csv = ct_dir / "track_features" / f"BEHAV3D_{ct_dir.name}_filtered_track_features.csv"
                        if filt_csv.exists():
                            warnings.append(f"Filtered data for {ct_dir.name}")
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
        skip_existing = False
        if existing:
            details = "\n".join(f"  • {w}" for w in existing)
            box = QMessageBox(self)
            box.setWindowTitle("Overwrite Existing Data?")
            box.setText(
                f"The following data already exists:\n\n{details}\n\n"
                "What do you want to do?"
            )
            btn_overwrite = box.addButton("Overwrite All", QMessageBox.DestructiveRole)
            btn_skip = box.addButton("Skip Existing", QMessageBox.AcceptRole)
            btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(btn_cancel)
            box.exec_()
            clicked = box.clickedButton()
            if clicked == btn_cancel:
                return
            skip_existing = (clicked == btn_skip)

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
            print(f"\n▶ [{i+1}/{len(self._steps)}] {step.display_label}...", file=sys.stderr)

            try:
                self._execute_step(step, skip_existing=skip_existing)
                step.elapsed = time.time() - step_start
                step.status = StepStatus.DONE
                print(f"✅ {step.display_label} completed in {step.elapsed:.1f}s", file=sys.stderr)
            except Exception as e:
                step.elapsed = time.time() - step_start
                step.status = StepStatus.ERROR
                step.error_msg = str(e)
                traceback.print_exc()
                print(f"❌ {step.display_label} failed: {e}", file=sys.stderr)
                
                # Show error dialog and stop queue
                QMessageBox.critical(
                    self, 
                    "Queue Error", 
                    f"An error occurred during {step.display_label}:\n\n{e}\n\nQueue execution stopped."
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

    def _execute_step(self, step: QueueStep, skip_existing: bool = False):
        """Run a single step, calling the backend directly."""
        if step.step_type == StepType.TRAIN:
            self._run_train()
        elif step.step_type == StepType.SEGMENT:
            self._run_segment(skip_existing=skip_existing)
        elif step.step_type == StepType.CELLPOSE_SEGMENT:
            self._run_cellpose_segment(step, skip_existing=skip_existing)
        elif step.step_type == StepType.DEAD_MASK:
            self._run_dead_mask(skip_existing=skip_existing)
        elif step.step_type == StepType.TRACK:
            self._run_track(skip_existing=skip_existing)
        elif step.step_type == StepType.APOC_SEGMENT:
            self._run_apoc_segment(step, skip_existing=skip_existing)
        elif step.step_type == StepType.FEATURE_EXTRACT:
            self._run_feature_extract(skip_existing=skip_existing)
        elif step.step_type == StepType.FILTER:
            self._run_filter(skip_existing=skip_existing)

    # ── Step runners (non-interactive) ─────────────────────────────────

    def _run_train(self):
        """Run classifier training (non-interactive)."""
        parent = self.parent()
        while parent and not hasattr(parent, 'tabs'):
            parent = parent.parent()
        if parent and hasattr(parent, 'tabs'):
            parent.tabs.setCurrentIndex(2) # Segmentation tab

        pc_widget = self.segmentation_tab.pixel_classifier_page
        pc_widget.run_train(interactive=False)

    def _run_segment(self, skip_existing: bool = False):
        """Run batch segmentation (non-interactive)."""
        pc_widget = self.segmentation_tab.pixel_classifier_page
        pc_widget.run_batch_segmentation(interactive=False, skip_existing=skip_existing)

    def _run_track(self, skip_existing: bool = False):
        """Run batch tracking (non-interactive)."""
        self.tracking_tab.run_batch_tracking(interactive=False, skip_existing=skip_existing)

    def _run_feature_extract(self, skip_existing: bool = False):
        """Run batch feature extraction (non-interactive)."""
        if self.feature_extraction_tab is not None:
            self.feature_extraction_tab.run_batch_feature_extraction(interactive=False, skip_existing=skip_existing)
        else:
            raise RuntimeError("Feature Extraction tab not wired to queue.")

    def _run_cellpose_segment(self, step: QueueStep, skip_existing: bool = False):
        """Run cellpose segmentation for the cell type stored in step.params."""
        cp_widget = self.segmentation_tab.cellpose_page
        cp_widget.run_batch_cellpose(
            interactive=False,
            skip_existing=skip_existing,
            cell_type_override=step.params.get("cell_type"),
            model_path_override=step.params.get("model_path"),
        )

    def _run_dead_mask(self, skip_existing: bool = False):
        """Run Otsu dead mask segmentation (non-interactive)."""
        cp_widget = self.segmentation_tab.cellpose_page
        cp_widget.run_otsu_threshold(interactive=False)

    def _run_filter(self, skip_existing: bool = False):
        """Run batch filtering (non-interactive)."""
        if self.filtering_tab is not None:
            self.filtering_tab.run_batch_filtering(interactive=False, skip_existing=skip_existing)
        else:
            raise RuntimeError("Filtering tab not wired to queue.")

    def _run_apoc_segment(self, step: QueueStep, skip_existing: bool = False):
        """Run APOC GPU batch segmentation, restoring snapshotted params first."""
        apoc_widget = self.segmentation_tab.apoc_page
        p = step.params
        # Restore snapshot into widgets so _on_run_segmentation reads the right values
        if "gpu_device_name" in p and hasattr(apoc_widget, "combo_gpu_device"):
            idx = apoc_widget.combo_gpu_device.findText(str(p["gpu_device_name"]))
            if idx >= 0:
                apoc_widget.combo_gpu_device.setCurrentIndex(idx)
        if "strategy_index" in p:
            apoc_widget.combo_strategy.setCurrentIndex(p["strategy_index"])
        if "edt_threshold" in p:
            apoc_widget.spin_edt_threshold.setValue(p["edt_threshold"])
        if "edt_min_size" in p:
            apoc_widget.spin_edt_min_size.setValue(p["edt_min_size"])
        if "prob_threshold" in p:
            apoc_widget.spin_prob_threshold.setValue(p["prob_threshold"])
        if "prob_min_size" in p:
            apoc_widget.spin_prob_min_size.setValue(p["prob_min_size"])
        # Restore per-CT size filter thresholds (Direct APOC strategy)
        if "min_sizes" in p:
            for ct, val in p["min_sizes"].items():
                spin = apoc_widget._size_filter_spins.get(ct)
                if spin is not None:
                    spin.setValue(val)
        if "workers" in p:
            apoc_widget.spin_workers.setValue(p["workers"])
        # Timepoint range
        if p.get("process_all", True):
            apoc_widget.check_process_all.setChecked(True)
        else:
            apoc_widget.check_process_all.setChecked(False)
            apoc_widget.spin_t_start.setValue(p.get("t_start", 0))
            apoc_widget.spin_t_end.setValue(p.get("t_end", 100))
        # skip_existing overrides the stored overwrite flag
        apoc_widget.check_overwrite.setChecked(not skip_existing)
        apoc_widget._on_run_segmentation(interactive=False)
