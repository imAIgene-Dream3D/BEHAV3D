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
    QProgressBar, QStackedWidget,
)
from qtpy.QtCore import Qt, Signal, Slot, QObject, QTimer
from qtpy.QtGui import QFont

from behav3d.napari._background_runner import ProgressBarRow


# ═══════════════════════════════════════════════════════════════════════════
# Data model
# ═══════════════════════════════════════════════════════════════════════════

class StepType(Enum):
    TRAIN = "train"
    SEGMENT = "segment"
    CELLPOSE_SEGMENT = "cellpose_segment"
    CELLPOSE_SAM_SEGMENT = "cellpose_sam_segment"
    DEAD_MASK = "dead_mask"
    APOC_SEGMENT = "apoc_segment"
    CONVPAINT_SEGMENT = "convpaint_segment"
    TRACK = "track"
    FEATURE_EXTRACT = "feature_extract"
    ACTIVE_KILLING = "active_killing"
    FILTER = "filter"
    DEATH_DYNAMICS = "death_dynamics"
    MULTI_ORG_DEATH = "multi_org_death"
    INTERACTION = "interaction"
    MULTI_ORG_INTERACTION = "multi_org_interaction"
    INVASIVENESS = "invasiveness"
    SC_STATE_CLUSTER = "sc_state_cluster"
    SC_TRAIN_STATE = "sc_train_state"
    SC_APPLY_STATE = "sc_apply_state"
    SC_TRACK_CLUSTER = "sc_track_cluster"

    @property
    def label(self):
        return {
            StepType.TRAIN: "🧠 Train Classifier",
            StepType.SEGMENT: "🦠 Segmentation",
            StepType.CELLPOSE_SEGMENT: "🦠 Segmentation",
            StepType.CELLPOSE_SAM_SEGMENT: "🦠 Segmentation",
            StepType.DEAD_MASK: "☠ Dead Mask (Otsu)",
            StepType.APOC_SEGMENT: "🦠 Segmentation",
            StepType.CONVPAINT_SEGMENT: "🦠 Segmentation",
            StepType.TRACK: "📍 Batch Tracking",
            StepType.FEATURE_EXTRACT: "🧪 Feature Extraction",
            StepType.ACTIVE_KILLING: "🔥 Active Killing",
            StepType.FILTER: "🧹 Filtering",
            StepType.DEATH_DYNAMICS: "💀 Death Dynamics",
            StepType.MULTI_ORG_DEATH: "💀 Combined Death Dynamics",
            StepType.INTERACTION: "🤝 Interaction Analysis",
            StepType.MULTI_ORG_INTERACTION: "🤝 Interaction Overview",
            StepType.INVASIVENESS: "🧗 Invasiveness Analysis",
            StepType.SC_STATE_CLUSTER: "🔬 State Clustering",
            StepType.SC_TRAIN_STATE: "🔬 Train State Classifier",
            StepType.SC_APPLY_STATE: "🔬 Apply State Classifier",
            StepType.SC_TRACK_CLUSTER: "🛤️ Track Clustering",
        }[self]

    @property
    def order(self):
        return {
            StepType.TRAIN: 0,
            StepType.SEGMENT: 1,
            StepType.CELLPOSE_SEGMENT: 1.5,
            StepType.CELLPOSE_SAM_SEGMENT: 1.5,
            StepType.DEAD_MASK: 1.6,
            StepType.APOC_SEGMENT: 1.5,
            StepType.CONVPAINT_SEGMENT: 1.5,
            StepType.TRACK: 2,
            StepType.FEATURE_EXTRACT: 3,
            StepType.ACTIVE_KILLING: 3.5,
            StepType.FILTER: 4,
            StepType.DEATH_DYNAMICS: 5,
            StepType.MULTI_ORG_DEATH: 5.5,
            StepType.INTERACTION: 6,
            StepType.MULTI_ORG_INTERACTION: 6.5,
            StepType.INVASIVENESS: 6.75,
            StepType.SC_STATE_CLUSTER: 7,
            StepType.SC_TRAIN_STATE: 7.5,
            StepType.SC_APPLY_STATE: 8,
            StepType.SC_TRACK_CLUSTER: 8.5,
        }[self]


#: Segmentation methods that run for one cell type per step. APOC and ConvPaint
#: instead segment every cell type in a single step, so they are excluded.
_PER_CELL_TYPE_SEGMENT_STEPS = frozenset({
    StepType.CELLPOSE_SEGMENT,
    StepType.CELLPOSE_SAM_SEGMENT,
})


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
        if self.step_type == StepType.SEGMENT:
            return "🦠 Segmentation — CPU Pixel Classifier"
        if self.step_type == StepType.APOC_SEGMENT:
            strat = self.params.get("strategy_name", "")
            base = "🦠 Segmentation — APOC (GPU)"
            return f"{base} — {strat}" if strat else base
        if self.step_type == StepType.CONVPAINT_SEGMENT:
            strat = self.params.get("strategy_name", "")
            base = "🦠 Segmentation — ConvPaint"
            return f"{base} — {strat}" if strat else base
        if self.step_type == StepType.CELLPOSE_SEGMENT:
            base = "🦠 Segmentation — Cellpose"
            return f"{base} — {cell_type}" if cell_type else base
        if self.step_type == StepType.CELLPOSE_SAM_SEGMENT:
            base = "🦠 Segmentation — Cellpose-SAM"
            if self.params.get("all_cell_types"):
                return f"{base} — all cell types"
            return f"{base} — {cell_type}" if cell_type else base
        if self.step_type == StepType.ACTIVE_KILLING:
            immune = self.params.get("immune_type", "")
            return f"🔥 Active Killing — {immune}" if immune else self.step_type.label
        if self.step_type in (
            StepType.DEATH_DYNAMICS,
            StepType.INTERACTION,
        ):
            cts = self.params.get("cell_types", [])
            if cts:
                return f"{self.step_type.label} — {', '.join(cts)}"
        if self.step_type in (
            StepType.MULTI_ORG_DEATH,
            StepType.MULTI_ORG_INTERACTION,
        ):
            cts = self.params.get("cell_types", [])
            if cts:
                return f"{self.step_type.label} — {', '.join(cts)}"
        if self.step_type == StepType.INVASIVENESS:
            immune = self.params.get("immune_cell_type", "")
            tgts = self.params.get("targets", [])
            if immune:
                suffix = f"{immune} → {', '.join(tgts)}" if tgts else immune
                return f"{self.step_type.label} — {suffix}"
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

        # Right-side stack: shows the elapsed-time label by default and
        # swaps to an inline mini progress bar while the step is RUNNING.
        # Both widgets share the same fixed footprint so the row height
        # stays at 24 px regardless of state.
        self.right_stack = QStackedWidget()
        self.right_stack.setFixedSize(60, 14)

        self.time_label = QLabel("")
        self.time_label.setStyleSheet("color: #888; font-size: 9px;")
        self.time_label.setAlignment(Qt.AlignCenter)
        self.right_stack.addWidget(self.time_label)

        self.mini_bar = QProgressBar()
        self.mini_bar.setRange(0, 0)  # busy by default; switches to determinate via update_progress
        self.mini_bar.setTextVisible(False)
        self.mini_bar.setFixedHeight(8)
        self.mini_bar.setStyleSheet(
            "QProgressBar { border: 1px solid #444; border-radius: 2px; "
            "background: #1a1a1a; }"
            "QProgressBar::chunk { background: #28a745; border-radius: 2px; }"
        )
        self.right_stack.addWidget(self.mini_bar)

        self.right_stack.setCurrentWidget(self.time_label)
        row.addWidget(self.right_stack)

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
        if self.step.status == StepStatus.RUNNING:
            # Show mini progress bar (busy by default; ``update_progress``
            # may switch it to determinate once sample-level ticks arrive).
            self.right_stack.setCurrentWidget(self.mini_bar)
            self.time_label.setStyleSheet("color: #888; font-size: 9px;")
        elif self.step.status == StepStatus.DONE:
            self.right_stack.setCurrentWidget(self.time_label)
            secs = self.step.elapsed
            if secs >= 60:
                self.time_label.setText(f"{secs/60:.1f}m")
            else:
                self.time_label.setText(f"{secs:.0f}s")
            self.time_label.setStyleSheet("color: #888; font-size: 9px;")
        elif self.step.status == StepStatus.ERROR:
            self.right_stack.setCurrentWidget(self.time_label)
            self.time_label.setText("failed")
            self.time_label.setStyleSheet("color: #dc3545; font-size: 10px;")
        else:
            self.right_stack.setCurrentWidget(self.time_label)
            self.time_label.setText("")
            self.time_label.setStyleSheet("color: #888; font-size: 9px;")
        # Disable delete while running/done
        self.btn_remove.setEnabled(self.step.status == StepStatus.PENDING)

    def update_progress(self, current: int, total: int) -> None:
        """Drive the inline mini progress bar (only meaningful while RUNNING).

        Called by :class:`ProcessingQueuePanel` from the active step's
        ``BackgroundOperation.progress`` signal.  ``total <= 0`` keeps the
        bar in busy/indeterminate mode; otherwise we switch to a
        determinate range and clamp the current value.
        """
        if total > 0:
            if self.mini_bar.maximum() != total:
                self.mini_bar.setRange(0, total)
            self.mini_bar.setValue(max(0, min(current, total)))
        else:
            self.mini_bar.setRange(0, 0)

    def reset_progress(self) -> None:
        """Reset the mini bar to busy/indeterminate state."""
        self.mini_bar.setRange(0, 0)
        self.mini_bar.setValue(0)


# ═══════════════════════════════════════════════════════════════════════════
# Main Queue Panel
# ═══════════════════════════════════════════════════════════════════════════

class ProcessingQueuePanel(QWidget):
    """Collapsible bottom panel for managing the processing queue."""

    # Emitted when queue finishes (success or with errors)
    queue_finished = Signal()

    # StepType.SEGMENT is a placeholder meaning "the active segmentation method".
    # _on_preset_selected replaces it at load-time with the real step type
    # (APOC_SEGMENT, CONVPAINT_SEGMENT, SEGMENT, CELLPOSE_SEGMENT) that
    # corresponds to whatever is currently selected in the method combo.
    PRESETS = {
        "Segment + Track": [StepType.SEGMENT, StepType.TRACK],
        "Segment → Filter": [StepType.SEGMENT, StepType.TRACK, StepType.FEATURE_EXTRACT, StepType.FILTER],
    }

    def __init__(self, segmentation_tab, tracking_tab, metadata_loader,
                 feature_extraction_tab=None, filtering_tab=None,
                 analysis_tab=None, parent=None):
        super().__init__(parent)
        self.segmentation_tab = segmentation_tab
        self.tracking_tab = tracking_tab
        self.metadata_loader = metadata_loader
        self.feature_extraction_tab = feature_extraction_tab
        self.filtering_tab = filtering_tab
        self.analysis_tab = analysis_tab

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

        # Always-visible thin progress bar inside the toggle bar; hidden
        # when the queue is idle and shown for the entire queue run.
        self.toggle_progress = QProgressBar()
        self.toggle_progress.setFixedSize(80, 8)
        self.toggle_progress.setRange(0, 0)
        self.toggle_progress.setTextVisible(False)
        self.toggle_progress.setStyleSheet(
            "QProgressBar { border: 1px solid #16213e; border-radius: 2px; "
            "background: #0f1320; }"
            "QProgressBar::chunk { background: #28a745; border-radius: 2px; }"
        )
        self.toggle_progress.setVisible(False)
        bar_layout.addSpacing(8)
        bar_layout.addWidget(self.toggle_progress)

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

        # Active segmentation indicator — updated whenever the method combo changes.
        self.active_seg_label = QLabel()
        self.active_seg_label.setStyleSheet(
            "color: #80c8ff; font-size: 10px; font-style: italic; padding-left: 4px;"
        )
        body_layout.addWidget(self.active_seg_label)
        self._refresh_active_seg_label()

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

        # Detailed queue-level progress row inserted just above the Run
        # Queue button.  Hidden when idle.  Driven by the same
        # _on_step_progress slot as the toggle-bar mini bar and per-row
        # mini bar, so all three stay in sync.
        self.body_progress_row = ProgressBarRow()
        self.body_progress_row.finish()
        self.body_progress_row.setVisible(False)
        body_layout.addWidget(self.body_progress_row)

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
        self.body.setMaximumHeight(290)  # Prevent inflating the Napari window
                                         # (raised slightly for the body progress row)

        # Don't let the queue panel push the window bigger when collapsed
        self.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)

        # Toggle click
        self.toggle_btn.clicked.connect(self._toggle_body)
        self.toggle_bar.mousePressEvent = lambda e: self._toggle_body()

        # Keep the active-segmentation label in sync when the user changes method.
        try:
            self.segmentation_tab.method_combo.currentIndexChanged.connect(
                lambda _: self._refresh_active_seg_label()
            )
        except Exception:
            pass

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
        are silently ignored. CELLPOSE_SAM_SEGMENT behaves the same way, except
        that an "all cell types" step covers every type in one go.
        """
        if params is None:
            params = {}

        # Dedup: per-cell-type methods dedup by (step_type, cell_type); others by step_type
        if step_type in _PER_CELL_TYPE_SEGMENT_STEPS:
            ct = params.get("cell_type", "")
            all_cts = bool(params.get("all_cell_types"))
            for s in self._steps:
                if s.step_type != step_type:
                    continue
                # An all-cell-types step subsumes any single-cell-type step, and
                # vice versa, so either combination is a duplicate.
                if s.params.get("all_cell_types") or all_cts:
                    return
                if s.params.get("cell_type") == ct:
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
                # Check if all missing cell types are covered by another
                # segmentation step already in the queue (APOC, ConvPaint,
                # or per-cell-type Cellpose).
                if not self._segmentation_covers_all_cell_types():
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

        # 4. ACTIVE_KILLING needs FEATURE_EXTRACT for the selected immune type
        elif step_type == StepType.ACTIVE_KILLING and StepType.FEATURE_EXTRACT not in added_already:
            missing = self._check_input_data_missing(StepType.ACTIVE_KILLING, params=params)
            if missing:
                res = QMessageBox.question(
                    self, "Missing Feature Data",
                    "Active Killing data is missing for the selected immune type.\n\n"
                    "Would you like to add Feature Extraction to the pipeline?",
                    QMessageBox.Yes | QMessageBox.No
                )
                if res == QMessageBox.Yes:
                    self.add_step(StepType.FEATURE_EXTRACT)
                else:
                    return  # Cancel adding Active Killing

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

    def _segmentation_covers_all_cell_types(self) -> bool:
        """Return True if every cell type that needs segmentation is covered.

        A cell type is considered covered when:

        * a non-per-cell-type segmentation step (``APOC_SEGMENT`` or
          ``CONVPAINT_SEGMENT``) is already queued — those run for every
          cell type at once, so a single step covers all of them, or
        * a ``CELLPOSE_SEGMENT`` / ``CELLPOSE_SAM_SEGMENT`` step with that
          ``cell_type`` is queued (or a Cellpose-SAM step set to all cell
          types), or
        * the ``*_segments.zarr`` already exists on disk for every sample.
        """
        md = self.metadata_loader.metadata
        if md is None:
            return False
        out_dir = Path(self.metadata_loader.output_dir)

        all_cts = []
        if self.tracking_tab:
            all_cts = list(self.tracking_tab.panels.keys())
        if not all_cts:
            return False

        # APOC / ConvPaint segment every cell type in a single batch run, and a
        # Cellpose-SAM step in all-cell-types mode does the same.
        if any(
            s.step_type in (StepType.APOC_SEGMENT, StepType.CONVPAINT_SEGMENT)
            or (s.step_type == StepType.CELLPOSE_SAM_SEGMENT and s.params.get("all_cell_types"))
            for s in self._steps
        ):
            return True

        queued_cellpose_cts = {
            s.params.get("cell_type")
            for s in self._steps
            if s.step_type in _PER_CELL_TYPE_SEGMENT_STEPS
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
            # Not on disk — check if a per-cell-type segmentation step covers it
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

    # ── Segmentation method resolution ─────────────────────────────────────

    # Maps method_combo index → (StepType, page_attr_name).
    # Positionally coupled to the `methods` list in _segmentation.py — keep both
    # in the same order. _assert_seg_method_alignment() checks this at runtime.
    _SEG_METHOD_MAP = [
        (StepType.APOC_SEGMENT,         "apoc_page"),
        (StepType.CONVPAINT_SEGMENT,    "convpaint_page"),
        (StepType.SEGMENT,              "pixel_classifier_page"),
        (StepType.CELLPOSE_SEGMENT,     "cellpose_page"),
        (StepType.CELLPOSE_SAM_SEGMENT, "cellpose_sam_page"),
        (None,                          None),   # Import — not queueable
    ]

    _SEG_METHOD_LABELS = [
        "APOC (GPU)",
        "ConvPaint",
        "Pixel Classifier (CPU)",
        "Cellpose",
        "Cellpose-SAM",
        "Import (not queueable)",
    ]

    def _assert_seg_method_alignment(self):
        """Warn if the method combo and _SEG_METHOD_MAP have drifted apart.

        These two lists live in different files and are coupled only by index, so
        a method added to one and not the other silently queues the wrong step.
        """
        try:
            n_combo = self.segmentation_tab.method_combo.count()
        except Exception:
            return
        if n_combo != len(self._SEG_METHOD_MAP):
            print(
                f"⚠️ Segmentation method lists are out of sync: method_combo has "
                f"{n_combo} entries, _SEG_METHOD_MAP has {len(self._SEG_METHOD_MAP)}. "
                "Queued segmentation steps may map to the wrong method."
            )

    def _active_segment_step(self) -> "tuple[StepType | None, dict]":
        """Return (StepType, params) for the currently selected segmentation method."""
        seg_tab = self.segmentation_tab
        idx = 0
        try:
            idx = seg_tab.method_combo.currentIndex()
        except Exception:
            pass
        if idx < 0 or idx >= len(self._SEG_METHOD_MAP):
            return None, {}
        step_type, page_attr = self._SEG_METHOD_MAP[idx]
        if step_type is None or page_attr is None:
            return None, {}
        params = {}
        try:
            page = getattr(seg_tab, page_attr, None)
            if page is not None and hasattr(page, "get_queue_params"):
                params = page.get_queue_params() or {}
        except Exception:
            pass
        return step_type, params

    def _refresh_active_seg_label(self):
        """Update the indicator label to show the current segmentation method."""
        try:
            idx = self.segmentation_tab.method_combo.currentIndex()
        except Exception:
            idx = -1
        if 0 <= idx < len(self._SEG_METHOD_LABELS):
            name = self._SEG_METHOD_LABELS[idx]
            step_type, _ = self._active_segment_step()
            if step_type is None:
                self.active_seg_label.setText(f"📍 Segmentation: {name} — cannot be queued")
            else:
                self.active_seg_label.setText(f"📍 Segmentation preset will use: {name}")
        else:
            self.active_seg_label.setText("")

    def _on_preset_selected(self, text):
        if text not in self.PRESETS:
            return
        seg_type, seg_params = self._active_segment_step()
        self._steps.clear()
        for st in self.PRESETS[text]:
            if st == StepType.SEGMENT:
                # Resolve to the currently active segmentation method.
                if seg_type is None:
                    from qtpy.QtWidgets import QMessageBox
                    QMessageBox.warning(
                        self,
                        "Cannot Apply Preset",
                        "The currently selected segmentation method cannot be queued.\n"
                        "Please select a different method in the Segmentation tab first.",
                    )
                    self._steps.clear()
                    self.preset_combo.blockSignals(True)
                    self.preset_combo.setCurrentIndex(0)
                    self.preset_combo.blockSignals(False)
                    self._rebuild_list()
                    return
                self._steps.append(QueueStep(step_type=seg_type, params=seg_params))
            else:
                self._steps.append(QueueStep(step_type=st))
        self._rebuild_list()
        if not self.body.isVisible():
            self._toggle_body()

    # ── Pre-run validation ─────────────────────────────────────────────

    def _check_input_data_missing(self, step_type: StepType, params: dict | None = None) -> bool:
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

        elif step_type == StepType.ACTIVE_KILLING:
            immune_type = None
            if params:
                immune_type = params.get("immune_type")
            if not immune_type and hasattr(self, "_pending_active_killing_params"):
                immune_type = self._pending_active_killing_params.get("immune_type")
            if not immune_type:
                return True
            feat_csv = (
                out_dir / "analysis" / immune_type / "track_features" /
                f"BEHAV3D_{immune_type}_combined_track_features.csv"
            )
            filtered_csv = (
                out_dir / "analysis" / immune_type / "track_features" /
                f"BEHAV3D_{immune_type}_combined_track_features_filtered.csv"
            )
            return not (feat_csv.exists() or filtered_csv.exists())

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

            elif step.step_type in _PER_CELL_TYPE_SEGMENT_STEPS:
                is_sam = step.step_type == StepType.CELLPOSE_SAM_SEGMENT
                method_label = "Cellpose-SAM" if is_sam else "Cellpose"
                if is_sam and step.params.get("all_cell_types"):
                    cts = [
                        ct for ct in (
                            list(self.tracking_tab.panels.keys()) if self.tracking_tab else []
                        )
                    ]
                else:
                    cts = [step.params.get("cell_type", "")]
                for ct in [c for c in cts if c]:
                    for _, sample in md.iterrows():
                        sn = sample.get("sample_name", "unknown")
                        seg_zarr = out_dir / "images" / sn / f"{sn}_{ct}_segments.zarr"
                        if seg_zarr.exists():
                            warnings.append(f"{method_label} segmentation for {ct} ({sn})")
                            break

            elif step.step_type == StepType.DEAD_MASK:
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    dead_zarr = out_dir / "images" / sn / f"{sn}_dead_segments.zarr"
                    if dead_zarr.exists():
                        warnings.append(f"Dead mask for {sn}")
                        break

            elif step.step_type in (StepType.APOC_SEGMENT, StepType.CONVPAINT_SEGMENT):
                # APOC and ConvPaint produce *_segments.zarr per cell type.
                method_label = (
                    "APOC" if step.step_type == StepType.APOC_SEGMENT else "ConvPaint"
                )
                cell_types = []
                if self.tracking_tab:
                    cell_types = list(self.tracking_tab.panels.keys())
                found = False
                for _, sample in md.iterrows():
                    sn = sample.get("sample_name", "unknown")
                    sample_img_dir = out_dir / "images" / sn
                    if not sample_img_dir.exists():
                        continue
                    if cell_types:
                        for ct in cell_types:
                            seg_zarr = sample_img_dir / f"{sn}_{ct}_segments.zarr"
                            if seg_zarr.exists():
                                warnings.append(f"{method_label} segments for {sn}")
                                found = True
                                break
                    else:
                        seg_files = list(sample_img_dir.glob("*_segments.zarr"))
                        if seg_files:
                            warnings.append(f"{method_label} segments for {sn}")
                            found = True
                    if found:
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

            elif step.step_type == StepType.ACTIVE_KILLING:
                immune_type = step.params.get("immune_type", "")
                if immune_type:
                    ak_dir = out_dir / "analysis" / immune_type / "active_killing"
                    if ak_dir.exists() and any(ak_dir.rglob("*.csv")):
                        warnings.append(f"Active killing data for {immune_type}")

            elif step.step_type == StepType.DEATH_DYNAMICS:
                cts = step.params.get("cell_types", []) or []
                for ct in cts:
                    res_csv = (
                        out_dir
                        / "analysis"
                        / ct
                        / "results"
                        / f"combined_general_{ct}_dynamics_analysis.csv"
                    )
                    if res_csv.exists():
                        warnings.append(f"Death dynamics results for {ct}")

            elif step.step_type == StepType.MULTI_ORG_DEATH:
                comp_dir = out_dir / "analysis" / "multi_organoid_comparison"
                if comp_dir.exists() and any(comp_dir.iterdir()):
                    warnings.append("Multi-organoid death dynamics outputs")

            elif step.step_type == StepType.INTERACTION:
                cts = step.params.get("cell_types", []) or []
                for ct in cts:
                    ia_dir = out_dir / "analysis" / ct / "interaction_analysis"
                    if ia_dir.exists() and any(ia_dir.iterdir()):
                        warnings.append(f"Interaction analysis for {ct}")

            elif step.step_type == StepType.MULTI_ORG_INTERACTION:
                comp_dir = out_dir / "analysis" / "multi_organoid_comparison"
                if comp_dir.exists() and any(comp_dir.iterdir()):
                    warnings.append("Interaction Overview outputs")

            elif step.step_type == StepType.INVASIVENESS:
                immune_list = step.params.get("immune_cell_types")
                if not immune_list:
                    single = step.params.get("immune_cell_type", "")
                    immune_list = [single] if single else []
                check_dirs = [
                    out_dir / "analysis" / im / "invasiveness_analysis"
                    for im in immune_list
                ]
                if len(immune_list) > 1:
                    check_dirs.append(out_dir / "analysis" / "invasiveness_analysis")
                for inv_dir in check_dirs:
                    if inv_dir.exists() and any(inv_dir.iterdir()):
                        warnings.append(
                            f"Invasiveness analysis for {', '.join(immune_list)}"
                        )
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

    # ── Run queue (event-driven state machine) ─────────────────────────
    #
    # The queue used to be a synchronous ``for step in self._steps`` loop
    # that called the backend directly and froze napari for the entire
    # run.  It now mirrors the per-tab background-execution pattern:
    # each step is dispatched via ``run_batch_*(block=False,
    # extra_callbacks=...)`` and the queue advances when the underlying
    # :class:`BackgroundOperation` fires its on-done / on-failed slot.
    # The state machine has four entry points:
    #   - :meth:`_on_run_queue`       — pre-flight, lock UI, kick off step 0.
    #   - :meth:`_start_next_step`    — mark RUNNING, dispatch, wire progress.
    #   - :meth:`_step_done`          — mark DONE, advance to idx+1.
    #   - :meth:`_step_failed`        — mark ERROR, dialog, finalize.
    # All CLI prints live verbatim inside these slots.

    def _on_run_queue(self):
        if not self._steps or self._is_running:
            return

        if self.metadata_loader.metadata is None:
            QMessageBox.warning(self, "No Metadata", "Please load metadata first.")
            return

        # Block conflicting segmentation strategies in the same queue run.
        step_types = {s.step_type for s in self._steps}
        if StepType.SEGMENT in step_types and StepType.APOC_SEGMENT in step_types:
            QMessageBox.critical(
                self,
                "Queue Conflict: Segmentation Steps",
                "❌ The queue contains both:\n"
                "  • 🦠 Batch Segmentation (CPU Pixel Classifier)\n"
                "  • ⚡ APOC\n\n"
                "These are two different segmentation strategies and cannot run together.\n\n"
                "Please remove one of the two segmentation steps before starting the queue.",
            )
            return

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

        # ── UI lock ──
        self._is_running = True
        self._skip_existing = skip_existing
        self._queue_total_start = time.time()
        self._active_step_idx = -1
        self._active_widget = None
        self._active_progress_conn = None
        self._step_done_signaled = False  # de-dup on_done/on_failed bursts
        self.btn_run.setEnabled(False)
        self.btn_clear.setEnabled(False)
        self.preset_combo.setEnabled(False)
        for row in self._step_rows:
            row.btn_remove.setEnabled(False)

        # ── Surface the queue-level progress widgets ──
        self.toggle_progress.setVisible(True)
        self.toggle_progress.setRange(0, 0)  # busy until first sample-level tick
        self.body_progress_row.setVisible(True)
        self.body_progress_row.start(f"Queue 0/{len(self._steps)}: starting…")

        print("", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        print("  🛒 BEHAV3D Processing Queue — Starting", file=sys.stderr)
        print("=" * 60, file=sys.stderr)

        # Kick off the state machine after the current event loop tick
        # so the UI lock paints before the first heavy step begins.
        QTimer.singleShot(0, lambda: self._start_next_step(0))

    def _start_next_step(self, idx: int) -> None:
        """Mark step ``idx`` as RUNNING and dispatch it asynchronously.

        Connects the active widget's :class:`BackgroundOperation`
        progress signal to :meth:`_on_step_progress` so the toggle-bar
        bar, body bar, and per-row mini bar all update in lock-step.
        """
        if idx >= len(self._steps):
            self._finalize_queue(success=True)
            return

        step = self._steps[idx]
        self._active_step_idx = idx
        self._step_done_signaled = False
        step.status = StepStatus.RUNNING
        step.elapsed = 0.0
        step.error_msg = ""
        self._step_start = time.time()
        row = self._step_rows[idx]
        row.refresh()
        row.reset_progress()

        total = len(self._steps)
        label = f"Queue {idx + 1}/{total}: {step.display_label}"
        # Body bar starts in busy mode; switches to determinate as soon
        # as a sample-level tick comes in via _on_step_progress.
        self.body_progress_row.start(label)
        self.toggle_progress.setRange(0, 0)

        print(f"\n▶ [{idx + 1}/{total}] {step.display_label}...", file=sys.stderr)

        try:
            widget = self._dispatch_step(step, self._skip_existing)
        except Exception as e:  # dispatch itself failed (sync exception)
            traceback.print_exc()
            QTimer.singleShot(0, lambda err=str(e): self._step_failed(idx, err))
            return

        self._connect_active_progress(widget)

    def _connect_active_progress(self, widget) -> None:
        """Subscribe to the dispatched widget's progress signal (if any).

        Synchronous helpers (analysis tab, training) return ``None``; the
        body bar simply stays in busy mode for those steps until the
        on_done / on_failed callback fires.
        """
        # Disconnect any previous connection.
        self._disconnect_active_progress()
        if widget is None:
            return
        bg = getattr(widget, "_bg", None)
        if bg is None:
            return
        try:
            bg.progress.connect(self._on_step_progress)
            self._active_widget = widget
            self._active_progress_conn = bg.progress
        except Exception:
            traceback.print_exc()

    def _disconnect_active_progress(self) -> None:
        """Detach from the previous step's progress signal, if any."""
        sig = self._active_progress_conn
        if sig is not None:
            try:
                sig.disconnect(self._on_step_progress)
            except Exception:
                pass
        self._active_progress_conn = None
        self._active_widget = None

    @Slot(int, int, str)
    def _on_step_progress(self, current: int, total: int, label: str) -> None:
        """Drive all three queue-level progress widgets in unison.

        ``current/total`` are the active step's sample-level
        coordinates; we project them into the global queue coordinate
        space (each step contributes one slice of size ``SCALE``).
        """
        idx = self._active_step_idx
        n_steps = len(self._steps)
        if idx < 0 or n_steps == 0:
            return

        SCALE = 1000  # sub-units per step (smooth ETA + overall percentage)
        if total > 0:
            frac = max(0.0, min(current / total, 1.0))
            global_curr = idx * SCALE + int(frac * SCALE)
            global_total = n_steps * SCALE
            # Toggle-bar mini bar — global percentage.
            if self.toggle_progress.maximum() != global_total:
                self.toggle_progress.setRange(0, global_total)
            self.toggle_progress.setValue(global_curr)

            # Body bar — global percentage with the active-step label.
            step_label = self._steps[idx].display_label if 0 <= idx < n_steps else ""
            full_label = (
                f"Queue {idx + 1}/{n_steps}: {step_label} — {label}"
                if label else f"Queue {idx + 1}/{n_steps}: {step_label}"
            )
            self.body_progress_row.update(global_curr, global_total, full_label)

            # Per-row mini bar — local sample percentage.
            row = self._step_rows[idx]
            row.update_progress(int(current), int(total))

    def _step_done_cb(self, _result=None) -> None:
        """``extra_callbacks['on_done']`` entry point (Qt thread)."""
        if self._step_done_signaled:
            return
        self._step_done_signaled = True
        idx = self._active_step_idx
        if idx < 0:
            return
        QTimer.singleShot(0, lambda i=idx: self._step_done(i))

    def _step_failed_cb(self, err: str) -> None:
        """``extra_callbacks['on_failed']`` entry point (Qt thread)."""
        if self._step_done_signaled:
            return
        self._step_done_signaled = True
        idx = self._active_step_idx
        if idx < 0:
            return
        QTimer.singleShot(0, lambda i=idx, e=str(err): self._step_failed(i, e))

    def _step_done(self, idx: int) -> None:
        """Mark step ``idx`` as DONE and queue the next step."""
        if idx < 0 or idx >= len(self._steps):
            return
        step = self._steps[idx]
        step.elapsed = time.time() - self._step_start
        step.status = StepStatus.DONE
        self._step_rows[idx].refresh()
        self._disconnect_active_progress()
        print(f"✅ {step.display_label} completed in {step.elapsed:.1f}s",
              file=sys.stderr)
        # Snap progress widgets to "step complete" before moving on.
        n_steps = len(self._steps)
        SCALE = 1000
        full = (idx + 1) * SCALE
        grand = n_steps * SCALE
        if self.toggle_progress.maximum() != grand:
            self.toggle_progress.setRange(0, grand)
        self.toggle_progress.setValue(full)
        QTimer.singleShot(0, lambda: self._start_next_step(idx + 1))

    def _step_failed(self, idx: int, err: str) -> None:
        """Mark step ``idx`` as ERROR, show dialog, finalize."""
        if idx < 0 or idx >= len(self._steps):
            return
        step = self._steps[idx]
        step.elapsed = time.time() - self._step_start
        step.status = StepStatus.ERROR
        step.error_msg = err
        self._step_rows[idx].refresh()
        self._disconnect_active_progress()
        print(f"❌ {step.display_label} failed: {err}", file=sys.stderr)
        QMessageBox.critical(
            self,
            "Queue Error",
            f"An error occurred during {step.display_label}:\n\n{err}\n\n"
            "Queue execution stopped.",
        )
        self._finalize_queue(success=False)

    def _finalize_queue(self, success: bool = True) -> None:
        """Reset UI lock, hide progress widgets, emit ``queue_finished``."""
        total_elapsed = time.time() - self._queue_total_start
        print(f"\n{'=' * 60}", file=sys.stderr)
        print(f"  🛒 Queue finished in {total_elapsed:.1f}s", file=sys.stderr)
        print(f"{'=' * 60}\n", file=sys.stderr)

        self._is_running = False
        self._active_step_idx = -1
        self._disconnect_active_progress()
        self.btn_run.setEnabled(True)
        self.btn_clear.setEnabled(True)
        self.preset_combo.setEnabled(True)
        for row in self._step_rows:
            row.btn_remove.setEnabled(row.step.status == StepStatus.PENDING)

        # Hide queue-level progress widgets (per-tab progress rows clear
        # themselves through the existing _bg.run lifecycle).
        self.toggle_progress.setVisible(False)
        self.toggle_progress.setRange(0, 0)
        self.body_progress_row.finish()
        self.body_progress_row.setVisible(False)

        self.queue_finished.emit()

    # ── Step dispatcher ────────────────────────────────────────────────

    def _dispatch_step(self, step: QueueStep, skip_existing: bool = False):
        """Dispatch ``step`` asynchronously; return the active widget.

        The returned widget exposes ``_bg.progress`` so the queue can
        drive the toggle-bar, body, and per-row progress bars in unison.
        Synchronous helpers (analysis tab) return ``None`` — the body
        bar stays in busy mode until on_done / on_failed fires.
        """
        cbs = {"on_done": self._step_done_cb, "on_failed": self._step_failed_cb}
        if step.step_type == StepType.TRAIN:
            return self._run_train(cbs)
        if step.step_type == StepType.SEGMENT:
            return self._run_segment(skip_existing, cbs)
        if step.step_type == StepType.CELLPOSE_SEGMENT:
            return self._run_cellpose_segment(step, skip_existing, cbs)
        if step.step_type == StepType.CELLPOSE_SAM_SEGMENT:
            return self._run_cellpose_sam_segment(step, skip_existing, cbs)
        if step.step_type == StepType.DEAD_MASK:
            return self._run_dead_mask(skip_existing, cbs)
        if step.step_type == StepType.TRACK:
            return self._run_track(skip_existing, cbs)
        if step.step_type == StepType.APOC_SEGMENT:
            return self._run_apoc_segment(step, skip_existing, cbs)
        if step.step_type == StepType.CONVPAINT_SEGMENT:
            return self._run_convpaint_segment(step, skip_existing, cbs)
        if step.step_type == StepType.FEATURE_EXTRACT:
            return self._run_feature_extract(skip_existing, cbs)
        if step.step_type == StepType.ACTIVE_KILLING:
            return self._run_active_killing(step, skip_existing, cbs)
        if step.step_type == StepType.FILTER:
            return self._run_filter(skip_existing, cbs)
        if step.step_type == StepType.DEATH_DYNAMICS:
            return self._run_death_dynamics(step, cbs)
        if step.step_type == StepType.MULTI_ORG_DEATH:
            return self._run_multi_org_death(step, cbs)
        if step.step_type == StepType.INTERACTION:
            return self._run_interaction(step, cbs)
        if step.step_type == StepType.MULTI_ORG_INTERACTION:
            return self._run_multi_org_interaction(step, cbs)
        if step.step_type == StepType.INVASIVENESS:
            return self._run_invasiveness(step, cbs)
        if step.step_type == StepType.SC_STATE_CLUSTER:
            return self._run_sc_state_cluster(step, cbs)
        if step.step_type == StepType.SC_TRAIN_STATE:
            return self._run_sc_train_state(step, cbs)
        if step.step_type == StepType.SC_APPLY_STATE:
            return self._run_sc_apply_state(step, cbs)
        if step.step_type == StepType.SC_TRACK_CLUSTER:
            return self._run_sc_track_cluster(step, cbs)
        # Unknown step — fail fast.
        QTimer.singleShot(
            0,
            lambda: self._step_failed_cb(f"Unknown step type: {step.step_type}"),
        )
        return None

    # ── Step runners (non-interactive) ─────────────────────────────────

    def _run_train(self, extra_callbacks):
        """Run classifier training (synchronous; queue advances on return)."""
        parent = self.parent()
        while parent and not hasattr(parent, 'tabs'):
            parent = parent.parent()
        if parent and hasattr(parent, 'tabs'):
            parent.tabs.setCurrentIndex(2)  # Segmentation tab

        pc_widget = self.segmentation_tab.pixel_classifier_page
        # ``run_train`` is intentionally synchronous; ``extra_callbacks``
        # fires before the call returns so we don't double-emit.
        pc_widget.run_train(interactive=False, extra_callbacks=extra_callbacks)
        return None  # synchronous — body bar stays busy

    def _run_segment(self, skip_existing, extra_callbacks):
        """Run batch segmentation asynchronously.

        Performs a pre-flight check: if any saved .joblib classifier has a
        feature count that does not match what make_features will produce for
        this dataset's pixel size, the step is aborted immediately with a
        clear diagnostic instead of failing deep inside a worker thread.
        """
        pc = self.segmentation_tab.pixel_classifier_page
        conflict_msg = self._preflight_pixel_classifier_features(pc)
        if conflict_msg:
            QMessageBox.critical(
                self,
                "Batch Segmentation Aborted",
                conflict_msg,
            )
            # Fire on_failed so the queue marks this step as errored and stops.
            from behav3d.napari._background_runner import fire_extra_callback
            fire_extra_callback(extra_callbacks, "on_failed", conflict_msg)
            return pc
        pc.run_batch_segmentation(
            interactive=False, skip_existing=skip_existing, block=False,
            extra_callbacks=extra_callbacks,
        )
        return pc

    def _preflight_pixel_classifier_features(self, pc_widget):
        """Return an error string if any .joblib classifier has a feature-count
        mismatch with the current make_features extractor, else return None."""
        try:
            import joblib
            import numpy as np
            from pathlib import Path
            from behav3d.preprocessing.segmentation.napari_pixelclassifier import make_features

            output_dir = self.metadata_loader.output_dir
            if not output_dir:
                return None
            pc_dir = Path(output_dir) / "images" / "PixelClassification"
            clf_files = list(pc_dir.glob("PixelClassifier_*.joblib")) if pc_dir.exists() else []
            if not clf_files:
                return None

            md = self.metadata_loader.metadata
            if md is None or md.empty:
                return None
            first_sample = md.iloc[0]
            pixel_xy = float(first_sample.get("pixel_distance_xy", 1.0))
            pixel_z  = float(first_sample.get("pixel_distance_z", pixel_xy))

            # Probe feature count: use a tiny dummy with the actual channel count.
            # We don't know how many channels were used at training time, so just
            # check the raw image channel count as an approximation.
            raw_shape = None
            try:
                from behav3d.io.images import load_image
                raw_path = Path(first_sample.get("raw_image_path", ""))
                zarr_path = Path(output_dir) / "images" / first_sample["sample_name"] / f"{first_sample['sample_name']}.zarr"
                probe_path = zarr_path if zarr_path.exists() else raw_path
                if probe_path.exists():
                    img = load_image(probe_path)
                    # shape: (T, C, Z, Y, X) or (T, C, Y, X)
                    n_ch = img.shape[1] if img.ndim >= 4 else 1
                    dummy = np.zeros((n_ch, 8, 8, 8), dtype=np.float32)
                    expected = make_features(pixel_xy, pixel_z, verbose=False)(dummy).shape[-1]
                    raw_shape = img.shape
            except Exception:
                return None  # can't probe → skip check

            if expected is None:
                return None

            mismatches = []
            for clf_path in clf_files:
                try:
                    clf = joblib.load(clf_path)
                    n_feat = getattr(clf, "n_features_in_", None)
                    if n_feat is not None and n_feat != expected:
                        mismatches.append((clf_path.name, int(n_feat), int(expected)))
                except Exception:
                    pass

            if not mismatches:
                return None

            lines = ["❌ Batch Segmentation aborted: classifier feature-count mismatch.\n"]
            for name, got, exp in mismatches:
                lines.append(
                    f"  • {name}: classifier expects {got} features "
                    f"but the current extractor produces {exp}."
                )
            lines.append(
                "\nThe classifier was likely trained with an older feature extractor.\n"
                "Please retrain it: go to the Pixel Classifier tab → "
                "\"Generate Training Data\" → \"Train Classifier\"."
            )
            return "\n".join(lines)
        except Exception:
            return None  # pre-flight errors must never block the user

    def _run_track(self, skip_existing, extra_callbacks):
        """Run batch tracking asynchronously."""
        self.tracking_tab.run_batch_tracking(
            interactive=False, skip_existing=skip_existing, block=False,
            extra_callbacks=extra_callbacks,
        )
        return self.tracking_tab

    def _run_feature_extract(self, skip_existing, extra_callbacks):
        """Run batch feature extraction asynchronously."""
        if self.feature_extraction_tab is None:
            QTimer.singleShot(
                0,
                lambda: self._step_failed_cb("Feature Extraction tab not wired to queue."),
            )
            return None
        self.feature_extraction_tab.run_batch_feature_extraction(
            interactive=False, skip_existing=skip_existing, block=False,
            extra_callbacks=extra_callbacks,
        )
        return self.feature_extraction_tab

    def _run_active_killing(self, step: QueueStep, skip_existing, extra_callbacks):
        """Run Active Killing analysis (synchronous)."""
        if self.feature_extraction_tab is None or not hasattr(
            self.feature_extraction_tab, "active_killing_panel"
        ):
            QTimer.singleShot(
                0,
                lambda: self._step_failed_cb("Feature Extraction tab not wired to queue."),
            )
            return None

        panel = self.feature_extraction_tab.active_killing_panel
        immune_type = step.params.get("immune_type", "")
        if immune_type:
            idx = panel.immune_combo.findText(str(immune_type))
            if idx >= 0:
                panel.immune_combo.setCurrentIndex(idx)

        self._pending_active_killing_params = dict(step.params)
        try:
            panel.run_analysis(interactive=False, extra_callbacks=extra_callbacks)
        finally:
            self._pending_active_killing_params = {}
        return None

    def _run_cellpose_segment(self, step: QueueStep, skip_existing, extra_callbacks):
        """Run cellpose segmentation asynchronously."""
        cp = self.segmentation_tab.cellpose_page
        cp.run_batch_cellpose(
            interactive=False, skip_existing=skip_existing, block=False,
            cell_type_override=step.params.get("cell_type"),
            model_path_override=step.params.get("model_path"),
            extra_callbacks=extra_callbacks,
        )
        return cp

    def _run_cellpose_sam_segment(self, step: QueueStep, skip_existing, extra_callbacks):
        """Run Cellpose-SAM segmentation from a queued params snapshot.

        ``step.params`` is the dict produced by
        :meth:`CellposeSAMWidget.get_queue_params`. As with APOC, the values are
        restored into the widget so the normal run path executes unchanged and
        the queue does not depend on whatever the user has since clicked.
        """
        cp = self.segmentation_tab.cellpose_sam_page
        p = step.params or {}

        model_name = p.get("model_name")
        if model_name and hasattr(cp, "combo_model"):
            idx = cp.combo_model.findText(str(model_name))
            if idx >= 0:
                cp.combo_model.setCurrentIndex(idx)

        # Torch device combo stores the device string as item data ("cuda:0"),
        # not as the display text ("GPU 0: NVIDIA ..."), so match on data.
        device = p.get("device")
        if device and hasattr(cp, "combo_gpu_device"):
            idx = cp.combo_gpu_device.findData(str(device))
            if idx >= 0:
                cp.combo_gpu_device.setCurrentIndex(idx)

        if hasattr(cp, "btn_force_cpu"):
            cp.btn_force_cpu.setChecked(bool(p.get("force_cpu", False)))

        # Restored before the run so a queued step keeps the power profile it was
        # queued with; a queue is exactly where an unattended overnight run lives,
        # which is where an abrupt shutdown costs the most.
        if hasattr(cp, "spin_n_threads"):
            cp.spin_n_threads.setValue(int(p.get("n_threads") or 0))
        if hasattr(cp, "spin_cooldown"):
            cp.spin_cooldown.setValue(float(p.get("cooldown_s") or 0.0))

        all_cell_types = bool(p.get("all_cell_types", False))
        if hasattr(cp, "check_all_cell_types"):
            cp.check_all_cell_types.setChecked(all_cell_types)

        # With all_cell_types restored above, leaving the override unset lets the
        # page expand to every detected cell type at run time.
        cell_type = None if all_cell_types else p.get("cell_type")
        cp.run_batch_cellpose_sam(
            interactive=False, skip_existing=skip_existing, block=False,
            cell_type_override=cell_type,
            extra_callbacks=extra_callbacks,
        )
        return cp

    def _run_dead_mask(self, skip_existing, extra_callbacks):
        """Run Otsu dead mask segmentation asynchronously."""
        cp = self.segmentation_tab.cellpose_page
        cp.run_otsu_threshold(
            interactive=False, block=False, extra_callbacks=extra_callbacks,
        )
        return cp

    def _run_filter(self, skip_existing, extra_callbacks):
        """Run batch filtering asynchronously."""
        if self.filtering_tab is None:
            QTimer.singleShot(
                0,
                lambda: self._step_failed_cb("Filtering tab not wired to queue."),
            )
            return None
        self.filtering_tab.run_batch_filtering(
            interactive=False, skip_existing=skip_existing, block=False,
            extra_callbacks=extra_callbacks,
        )
        return self.filtering_tab

    def _require_death_dynamics_tab(self):
        if self.analysis_tab is None or not hasattr(self.analysis_tab, "death_dynamics_tab"):
            raise RuntimeError("Analysis tab not wired to queue.")
        return self.analysis_tab.death_dynamics_tab

    def _run_sync_analysis(self, fn, extra_callbacks):
        """Run a synchronous analysis-tab call and fire ``extra_callbacks``.

        Wraps the existing ``death_dynamics_tab.run_*`` helpers so the
        state machine sees the same on_done / on_failed flow as the
        async batch methods.  These calls are synchronous on purpose —
        napari will briefly freeze for the duration of each step, but
        the queue still advances cleanly when ``fn`` returns.
        """
        try:
            result = fn()
            if extra_callbacks and extra_callbacks.get("on_done") is not None:
                extra_callbacks["on_done"](result)
        except Exception as e:
            traceback.print_exc()
            if extra_callbacks and extra_callbacks.get("on_failed") is not None:
                extra_callbacks["on_failed"](str(e))

    def _run_death_dynamics(self, step: QueueStep, extra_callbacks):
        dd = self._require_death_dynamics_tab()
        cell_types = step.params.get("cell_types") or []
        self._run_sync_analysis(
            lambda: dd.run_death_dynamics_for(cell_types, interactive=False),
            extra_callbacks,
        )
        return None

    def _run_multi_org_death(self, step: QueueStep, extra_callbacks):
        dd = self._require_death_dynamics_tab()
        cell_types = step.params.get("cell_types") or []
        self._run_sync_analysis(
            lambda: dd.run_multi_organoid_death_for(cell_types, interactive=False),
            extra_callbacks,
        )
        return None

    def _run_interaction(self, step: QueueStep, extra_callbacks):
        dd = self._require_death_dynamics_tab()
        cell_types = step.params.get("cell_types") or []
        interaction_cts = step.params.get("interaction_cell_types") or []
        group_by_line_condition = bool(
            step.params.get("group_by_line_condition", False)
        )
        self._run_sync_analysis(
            lambda: dd.run_interaction_for(
                cell_types,
                interaction_cts,
                group_by_line_condition=group_by_line_condition,
                interactive=False,
            ),
            extra_callbacks,
        )
        return None

    def _run_multi_org_interaction(self, step: QueueStep, extra_callbacks):
        dd = self._require_death_dynamics_tab()
        cell_types = step.params.get("cell_types") or []
        interaction_cts = step.params.get("interaction_cell_types") or []
        time_window_min = float(step.params.get("time_window_min", 60.0))
        group_by = step.params.get("group_by", "organoid_type")
        annotate_line_condition = bool(
            step.params.get("annotate_line_condition", False)
        )
        analysis_period_t = step.params.get("analysis_period_t")
        if analysis_period_t is None:
            analysis_period_t = step.params.get("analysis_period_min")
        self._run_sync_analysis(
            lambda: dd.run_multi_organoid_interaction_for(
                cell_types,
                interaction_cts,
                time_window_min=time_window_min,
                group_by=group_by,
                annotate_line_condition=annotate_line_condition,
                analysis_period_t=analysis_period_t,
                interactive=False,
            ),
            extra_callbacks,
        )
        return None

    def _run_invasiveness(self, step: QueueStep, extra_callbacks):
        dd = self._require_death_dynamics_tab()
        # Prefer the multi-immune list; fall back to the legacy single field
        # for older queued steps.
        immune = step.params.get("immune_cell_types")
        if not immune:
            single = step.params.get("immune_cell_type")
            immune = [single] if single else []
        targets = step.params.get("targets") or []
        summary_stat = step.params.get("summary_stat", "mean")
        group_by_line_condition = bool(
            step.params.get("group_by_line_condition", False)
        )
        analysis_period_t = step.params.get("analysis_period_t")
        if analysis_period_t is None:
            analysis_period_t = step.params.get("analysis_period_min")
        self._run_sync_analysis(
            lambda: dd.run_invasiveness_for(
                immune,
                targets,
                summary_stat=summary_stat,
                group_by_line_condition=group_by_line_condition,
                analysis_period_t=analysis_period_t,
                interactive=False,
            ),
            extra_callbacks,
        )
        return None

    def _run_apoc_segment(self, step: QueueStep, skip_existing, extra_callbacks):
        """Run APOC GPU batch segmentation from a queued params snapshot.

        ``step.params`` is the dict produced by :meth:`APOCWidget.get_queue_params`.
        The values are restored into the widget so the existing
        ``_on_run_segmentation`` code path can run unchanged; only attributes
        that actually exist on :class:`APOCWidget` are touched.
        Dispatches asynchronously and returns the widget so the queue
        can subscribe to its progress signal.
        """
        apoc_widget = self.segmentation_tab.apoc_page
        p = step.params or {}
        tw = getattr(apoc_widget, "_training_widget", None)

        # GPU device
        gpu_name = p.get("gpu_device_name")
        if gpu_name and hasattr(apoc_widget, "combo_gpu_device"):
            idx = apoc_widget.combo_gpu_device.findText(str(gpu_name))
            if idx >= 0:
                apoc_widget.combo_gpu_device.setCurrentIndex(idx)

        force_cpu = p.get("force_cpu", False)
        if hasattr(apoc_widget, "btn_force_cpu"):
            apoc_widget.btn_force_cpu.setChecked(force_cpu)

        # Strategy: snapshot stored ``strategy_name`` / ``strategy_index``.
        if tw is not None and tw.strategy_combo is not None:
            strat_name = p.get("strategy_name")
            if strat_name and tw.strategy_combo.findText(strat_name) >= 0:
                tw.strategy_combo.setCurrentText(strat_name)
            elif "strategy_index" in p:
                idx = int(p["strategy_index"])
                if 0 <= idx < tw.strategy_combo.count():
                    tw.strategy_combo.setCurrentIndex(idx)

        # Per-cell-type strategies in Advanced mode.
        per_ct_strategies = p.get("per_ct_strategies") or {}
        if tw is not None and per_ct_strategies:
            for ct, strat in per_ct_strategies.items():
                tab = tw.tabs.get(ct)
                if tab is None:
                    continue
                combo = getattr(tab, "_per_tab_strategy_combo", None)
                if combo is not None and combo.findText(str(strat)) >= 0:
                    combo.setCurrentText(str(strat))

        # Per-tab post-processing values.
        per_ct_params = p.get("per_ct_params") or {}
        if tw is not None and per_ct_params:
            for ct, cfg in per_ct_params.items():
                tab = tw.tabs.get(ct)
                if tab is None:
                    continue
                try:
                    tab.apply_config(cfg)
                except Exception:
                    pass

        # Worker / overwrite / timepoint settings.
        if "workers" in p and hasattr(apoc_widget, "spin_workers"):
            apoc_widget.spin_workers.setValue(int(p["workers"]))
        if hasattr(apoc_widget, "check_process_all"):
            if p.get("process_all", True):
                apoc_widget.check_process_all.setChecked(True)
            else:
                apoc_widget.check_process_all.setChecked(False)
                if hasattr(apoc_widget, "spin_t_start"):
                    apoc_widget.spin_t_start.setValue(int(p.get("t_start", 0)))
                if hasattr(apoc_widget, "spin_t_end"):
                    apoc_widget.spin_t_end.setValue(int(p.get("t_end", 100)))
        if hasattr(apoc_widget, "check_overwrite"):
            # skip_existing overrides any stored overwrite flag.
            apoc_widget.check_overwrite.setChecked(not skip_existing)

        apoc_widget._on_run_segmentation(
            interactive=False, block=False, extra_callbacks=extra_callbacks,
        )
        return apoc_widget

    def _run_convpaint_segment(self, step: QueueStep, skip_existing, extra_callbacks):
        """Run ConvPaint batch segmentation from a queued params snapshot.

        ``step.params`` is the dict produced by :meth:`ConvPaintWidget.get_queue_params`.
        Mirrors :meth:`_run_apoc_segment`: values are restored into the widget so
        the existing ``_on_run_segmentation`` code path runs unchanged.
        Dispatches asynchronously and returns the widget for progress subscription.
        """
        convpaint_widget = self.segmentation_tab.convpaint_page
        p = step.params or {}
        tw = getattr(convpaint_widget, "_training_widget", None)

        # Torch device
        device = p.get("device")
        if device and hasattr(convpaint_widget, "combo_gpu_device"):
            idx = convpaint_widget.combo_gpu_device.findData(str(device))
            if idx >= 0:
                convpaint_widget.combo_gpu_device.setCurrentIndex(idx)

        # Strategy
        if tw is not None and tw.strategy_combo is not None:
            strat_name = p.get("strategy_name")
            if strat_name and tw.strategy_combo.findText(strat_name) >= 0:
                tw.strategy_combo.setCurrentText(strat_name)
            elif "strategy_index" in p:
                idx = int(p["strategy_index"])
                if 0 <= idx < tw.strategy_combo.count():
                    tw.strategy_combo.setCurrentIndex(idx)

        # Per-cell-type strategies (Advanced mode)
        per_ct_strategies = p.get("per_ct_strategies") or {}
        if tw is not None and per_ct_strategies:
            for ct, strat in per_ct_strategies.items():
                tab = tw.tabs.get(ct)
                if tab is None:
                    continue
                combo = getattr(tab, "_per_tab_strategy_combo", None)
                if combo is not None and combo.findText(str(strat)) >= 0:
                    combo.setCurrentText(str(strat))

        # Per-tab post-processing values (apply via collect_params reverse path is
        # not available, so we apply individual widget values where present).
        per_ct_params = p.get("per_ct_params") or {}
        if tw is not None and per_ct_params:
            for ct, cfg in per_ct_params.items():
                tab = tw.tabs.get(ct)
                if tab is None or not isinstance(cfg, dict):
                    continue
                # ConvPaint per-tab keys are flat ``{ct}_<name>`` entries; the
                # tab itself exposes spinners/checkboxes named after the keys
                # it produced. Reverse-apply by looking for the corresponding
                # attribute on the tab.
                for key, value in cfg.items():
                    # Strip leading "{ct}_" prefix if present
                    short = key
                    prefix = f"{ct}_"
                    if short.startswith(prefix):
                        short = short[len(prefix):]
                    attr_map = {
                        "edt_threshold": "edt_threshold_spin",
                        "opening_nr_pixels": "opening_nr_pixels_spin",
                        "fill_holes": "fill_holes_cb",
                        "segment_size_min": "segment_size_min_spin",
                        "peak_min_distance": "peak_min_distance_spin",
                        "peak_min_ratio": "peak_min_ratio_spin",
                        "prob_mask_threshold": "prob_mask_threshold_spin",
                        "prob_seed_threshold": "prob_seed_threshold_spin",
                    }
                    attr_name = attr_map.get(short)
                    if attr_name is None:
                        continue
                    w = getattr(tab, attr_name, None)
                    if w is None:
                        continue
                    try:
                        if attr_name.endswith("_cb"):
                            w.setChecked(bool(value))
                        else:
                            w.setValue(float(value) if "threshold" in short or "ratio" in short or "distance" in short else int(value))
                    except Exception:
                        pass

        # Worker / overwrite / timepoint settings.
        if "workers" in p and hasattr(convpaint_widget, "spin_workers"):
            convpaint_widget.spin_workers.setValue(int(p["workers"]))
        if hasattr(convpaint_widget, "check_process_all"):
            if p.get("process_all", True):
                convpaint_widget.check_process_all.setChecked(True)
            else:
                convpaint_widget.check_process_all.setChecked(False)
                if hasattr(convpaint_widget, "spin_t_start"):
                    convpaint_widget.spin_t_start.setValue(int(p.get("t_start", 0)))
                if hasattr(convpaint_widget, "spin_t_end"):
                    convpaint_widget.spin_t_end.setValue(int(p.get("t_end", 100)))
        if hasattr(convpaint_widget, "check_overwrite"):
            convpaint_widget.check_overwrite.setChecked(not skip_existing)

        convpaint_widget._on_run_segmentation(
            interactive=False, block=False, extra_callbacks=extra_callbacks,
        )
        return convpaint_widget

    # ── Single Cell step runners ──────────────────────────────────────────

    def _require_single_cell_tab(self):
        """Return the SingleCellTab embedded in the AnalysisTab."""
        if self.analysis_tab is None:
            raise RuntimeError("Analysis tab not wired to queue.")
        sc = getattr(self.analysis_tab, "single_cell_tab", None)
        if sc is None:
            raise RuntimeError("Single Cell tab not found inside Analysis tab.")
        return sc

    def _run_sc_state_cluster(self, step: QueueStep, extra_callbacks):
        """Run HMM state clustering (asynchronous via BackgroundOperation)."""
        try:
            sc = self._require_single_cell_tab()
        except RuntimeError as e:
            QTimer.singleShot(0, lambda err=str(e): self._step_failed_cb(err))
            return None

        ct = step.params.get("cell_type") or ""
        if ct:
            combo = getattr(sc, "cell_type_combo", None)
            if combo is not None:
                idx = combo.findText(ct)
                if idx >= 0:
                    combo.setCurrentIndex(idx)

        sc.state_tab.run_state_classification(
            interactive=False, extra_callbacks=extra_callbacks
        )
        return sc.state_tab  # exposes ._bg.progress for queue progress bar

    def _run_sc_train_state(self, step: QueueStep, extra_callbacks):
        """Train state random-forest classifier (asynchronous)."""
        try:
            sc = self._require_single_cell_tab()
        except RuntimeError as e:
            QTimer.singleShot(0, lambda err=str(e): self._step_failed_cb(err))
            return None

        ct = step.params.get("cell_type") or ""
        if ct:
            combo = getattr(sc, "cell_type_combo", None)
            if combo is not None:
                idx = combo.findText(ct)
                if idx >= 0:
                    combo.setCurrentIndex(idx)

        sc.state_tab.run_train_state(
            interactive=False, extra_callbacks=extra_callbacks
        )
        return sc.state_tab

    def _run_sc_apply_state(self, step: QueueStep, extra_callbacks):
        """Apply state classifier to full dataset (asynchronous)."""
        try:
            sc = self._require_single_cell_tab()
        except RuntimeError as e:
            QTimer.singleShot(0, lambda err=str(e): self._step_failed_cb(err))
            return None

        ct = step.params.get("cell_type") or ""
        if ct:
            combo = getattr(sc, "cell_type_combo", None)
            if combo is not None:
                idx = combo.findText(ct)
                if idx >= 0:
                    combo.setCurrentIndex(idx)

        sc.state_tab.run_apply_state(
            interactive=False, extra_callbacks=extra_callbacks
        )
        return sc.state_tab

    def _run_sc_track_cluster(self, step: QueueStep, extra_callbacks):
        """Run DTW track clustering (asynchronous)."""
        try:
            sc = self._require_single_cell_tab()
        except RuntimeError as e:
            QTimer.singleShot(0, lambda err=str(e): self._step_failed_cb(err))
            return None

        ct = step.params.get("cell_type") or ""
        if ct:
            combo = getattr(sc, "cell_type_combo", None)
            if combo is not None:
                idx = combo.findText(ct)
                if idx >= 0:
                    combo.setCurrentIndex(idx)

        sc.track_tab.run_track_clustering(
            interactive=False, extra_callbacks=extra_callbacks
        )
        return sc.track_tab
