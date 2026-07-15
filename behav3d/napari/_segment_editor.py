"""Inline Qt widget for manual correction of tracked-segment volumes.

Embedded both in the BEHAV3D plugin's Visualization tab and in the
notebook-launched napari window.  Uses :class:`behav3d.editing.EditBuffer`
under the hood, so undo / redo / save / discard semantics are identical
across both surfaces.

Design notes
------------
* The editor binds to **one** Labels layer at a time and refuses any
  layer whose name does not end in ``" tracked segments"``.
* All operations stage in memory; the napari layer's ``.data`` array is
  rewritten in place after each commit so the user sees the result
  immediately.
* Save persists to the original ``*_tracked.zarr`` (only dirty frames)
  and regenerates the matching ``*_tracks.csv``.
* Tab-switch / window-close prompts are handled by ``request_exit()``,
  which the parent widget calls before tearing the editor down.
"""
from __future__ import annotations

import traceback
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import dask.array as da
import numpy as np
import pandas as pd
from qtpy.QtCore import Qt, QObject, QThread, QTimer, Signal
from qtpy.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QSpinBox,
    QStackedWidget,
    QTextEdit,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from behav3d.editing import (
    EditBuffer,
    create_label,
    delete_label,
    dilate_label,
    erode_label,
    lifetime_of,
    merge_labels,
    split_label,
)
from behav3d.napari._loaders import tracked_segments_paths_for
from behav3d.napari._tutorial_dialog import TutorialButton, TutorialStep

try:
    from napari.utils import progress as _NapariProgress
except Exception:
    _NapariProgress = None

try:
    from napari.utils.colormaps import DirectLabelColormap
except Exception:
    DirectLabelColormap = None


_SPLIT_SEEDS_LAYER = "__split_seeds__"
_HIGHLIGHT_LAYER = "__selected_label_highlight__"
_HIGHLIGHT_COLOR = "yellow"
_TRACKED_SUFFIX = " tracked segments"

_SEED_TUTORIAL_STEPS = [
    TutorialStep(
        "1. Select the label to edit",
        "seed_step1_select_label.png",
        "Click the cell you want to split/create from in the viewer "
        "(or type its ID above). It becomes the active label.",
    ),
    TutorialStep(
        "2. Click “Add seed points”",
        "seed_step2_add_seed_points.png",
        "This creates a temporary points layer and switches the view to "
        "2D so seed clicks land precisely.",
    ),
    TutorialStep(
        "3. Place seeds, one per sub-region",
        "seed_step3_place_seeds.png",
        "Click inside the label for each seed. Move the Z slider between "
        "clicks if you want seeds on different planes — seeds don't "
        "need to share a Z-slice, they just need to land inside the label.",
    ),
    TutorialStep(
        "4. Preview",
        "seed_step4_preview.png",
        "Once the seed counter shows enough seeds, click Preview to run "
        "the split/create on the current frame only, so you can check the "
        "result before committing to the full time range.",
    ),
    TutorialStep(
        "5. Apply",
        "seed_step5_apply.png",
        "Happy with the preview? Click Apply to commit it and propagate "
        "the result forward and backward across the track's lifetime.",
    ),
]


# ---------------------------------------------------------------------------
# Background materialisation worker
# ---------------------------------------------------------------------------
class _MaterialiseWorker(QObject):
    """Materialises a dask/zarr-backed labels array into RAM in a QThread.

    Emits ``done(np.ndarray)`` when finished, or ``failed(str)`` on error,
    followed in both cases by ``finished()`` (no arguments) so the owning
    QThread can be told to quit via a clean, zero-argument signal.
    The worker reads ``layer.data`` (whatever it is at start time) and emits
    the resulting numpy array — the caller is responsible for writing it back
    to the layer on the Qt thread.
    """

    done = Signal(object)   # np.ndarray
    failed = Signal(str)    # error message
    finished = Signal()     # always emitted last; used to quit the QThread

    def __init__(self, data) -> None:
        super().__init__()
        self._data = data

    def run(self) -> None:
        try:
            result = np.asarray(self._data)
            self.done.emit(result)
        except Exception as exc:
            self.failed.emit(str(exc))
        finally:
            self.finished.emit()


# ---------------------------------------------------------------------------
# Background edit operation worker
# ---------------------------------------------------------------------------
class _EditWorker(QObject):
    """Runs a single editing primitive in a background QThread.

    The callable ``fn`` receives ``*args, **kwargs`` and must return either an
    :class:`~behav3d.editing.OpResult` or a ``list[OpResult]`` (delete uses
    a list when multiple labels are removed in one shot).  All reads from the
    :class:`~behav3d.editing.EditBuffer` happen inside the thread; the caller
    is responsible for calling ``buf.apply()`` back on the Qt thread via the
    ``done`` signal so the buffer is never written from a background thread.
    """

    done = Signal(object)    # OpResult or list[OpResult]
    failed = Signal(str)     # error message
    finished = Signal()      # always emitted last; used to quit the QThread
    progress = Signal(int, int)  # (current_frame, total_frames)

    def __init__(self, fn, args, kwargs) -> None:
        super().__init__()
        self._fn = fn
        self._args = args
        self._kwargs = kwargs

    def run(self) -> None:
        try:
            result = self._fn(*self._args, **self._kwargs)
            self.done.emit(result)
        except Exception as exc:
            import traceback as _tb
            _tb.print_exc()
            self.failed.emit(str(exc))
        finally:
            self.finished.emit()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _is_tracked_layer_name(name: str) -> bool:
    return isinstance(name, str) and name.endswith(_TRACKED_SUFFIX)


def _cell_type_from_layer_name(name: str) -> Optional[str]:
    """Extract ``ct_name`` from ``"<sample> – <ct_name> tracked segments"``."""
    if not _is_tracked_layer_name(name):
        return None
    body = name[: -len(_TRACKED_SUFFIX)].rstrip()
    # body is "<sample> – <ct_name>" — split on the en-dash used elsewhere
    if " – " in body:
        return body.rsplit(" – ", 1)[1].strip()
    return None


def _sample_from_layer_name(name: str) -> Optional[str]:
    if not _is_tracked_layer_name(name):
        return None
    body = name[: -len(_TRACKED_SUFFIX)].rstrip()
    if " – " in body:
        return body.rsplit(" – ", 1)[0].strip()
    return None


# ---------------------------------------------------------------------------
# TrackedSegmentEditor
# ---------------------------------------------------------------------------
class TrackedSegmentEditor(QWidget):
    """Manual-correction editor bound to a single tracked-segments layer.

    Parameters
    ----------
    viewer:
        The napari viewer hosting the layers.
    metadata_loader:
        Object exposing ``.metadata`` (a ``pd.DataFrame``).  Used to
        resolve zarr/CSV paths and pixel sizes for the active sample.
    layer_name:
        Initial Labels layer to bind to.  Must end in
        ``" tracked segments"``.
    """

    def __init__(self, viewer, metadata_loader, layer_name: str, parent=None):
        super().__init__(parent)
        if not _is_tracked_layer_name(layer_name):
            raise ValueError(
                f"TrackedSegmentEditor only binds to '*{_TRACKED_SUFFIX}' "
                f"layers; got {layer_name!r}"
            )
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.layer_name = layer_name
        self.cell_type = _cell_type_from_layer_name(layer_name) or ""
        self.sample_name = _sample_from_layer_name(layer_name) or ""

        # Resolve metadata + paths.
        row = self._row_for_sample(self.sample_name)
        paths = tracked_segments_paths_for(row, self.cell_type) if row is not None else None
        if paths is None:
            raise ValueError(
                f"No tracked-segments path in metadata for sample={self.sample_name!r}, "
                f"cell_type={self.cell_type!r}"
            )
        self._zarr_path: Path = paths["zarr"]
        self._csv_path: Optional[Path] = paths["csv"]

        pix_xy = float(row.get("pixel_distance_xy") or 1.0)
        pix_z = float(row.get("pixel_distance_z") or 1.0)

        self.buffer = EditBuffer(
            zarr_path=self._zarr_path,
            csv_path=self._csv_path,
            pixel_size_xy=pix_xy,
            pixel_size_z=pix_z,
        )
        self.buffer.add_frames_changed_listener(self._on_frames_changed)

        # Selected labels (read from clicks on the layer).
        self._selected_labels: List[int] = []

        # Background materialisation state.
        self._materialise_thread: Optional[QThread] = None
        self._materialise_worker: Optional[_MaterialiseWorker] = None

        # Background edit operation state.
        self._edit_thread: Optional[QThread] = None
        self._edit_worker: Optional[_EditWorker] = None
        self._pending_on_success = None  # optional callback(op) after _commit_op
        self._pending_preview_tool: Optional[str] = None
        # Preview state — frames placed here are shown in the napari layer
        # without touching the EditBuffer or the undo stack.  Cleared before
        # every Apply / Undo / Redo / Discard / tool-switch.
        self._preview_dirty: Dict[int, np.ndarray] = {}
        self._preview_tool: Optional[str] = None
        # napari activity-dock progress bar (None when no operation is running).
        self._activity_progress = None
        # ndisplay saved while seed placement forces the view to 2D (None
        # when not currently overridden).
        self._ndisplay_before_seeds: Optional[int] = None

        self._build_ui()
        self._bind_to_layer()

    # ------------------------------------------------------------------
    # UI
    # ------------------------------------------------------------------
    def _build_ui(self) -> None:
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(4)

        group = QGroupBox("Manual edition")
        outer.addWidget(group)
        layout = QVBoxLayout(group)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        # ---- Header ---------------------------------------------------
        header = QHBoxLayout()
        title = QLabel(
            f"<b>{self.sample_name}</b> · {self.cell_type} <span style='color:#c62828'>● Editing</span>"
        )
        title.setTextFormat(Qt.RichText)
        header.addWidget(title)
        header.addStretch(1)
        layout.addLayout(header)

        # ---- Selected labels row -------------------------------------
        sel_row = QHBoxLayout()
        sel_row.addWidget(QLabel("Selected labels:"))
        self.lbl_selected = QLabel("(no labels selected)")
        self.lbl_selected.setStyleSheet("color:#1976D2; font-weight:bold;")
        sel_row.addWidget(self.lbl_selected, stretch=1)
        self.input_label_id = QLineEdit()
        self.input_label_id.setPlaceholderText("Add by ID")
        self.input_label_id.setFixedWidth(80)
        self.input_label_id.setToolTip(
            "Type a TrackID and press Enter to add it to the selection (use\n"
            "this when you already know the ID; clicking the layer also works)."
        )
        self.input_label_id.returnPressed.connect(self._add_label_by_id)
        sel_row.addWidget(self.input_label_id)
        self.btn_clear_sel = QPushButton("Clear")
        self.btn_clear_sel.clicked.connect(self._clear_selection)
        sel_row.addWidget(self.btn_clear_sel)
        layout.addLayout(sel_row)
        sel_hint = QLabel(
            "Tip: left-click adds a label, Shift-click adds another, right-click removes it."
        )
        sel_hint.setStyleSheet("color:#888; font-size:10px; font-style:italic;")
        layout.addWidget(sel_hint)

        # ---- Tool toolbar (drives the QStackedWidget) ----------------
        toolbar = QHBoxLayout()
        toolbar.setSpacing(2)
        self._tool_group = QButtonGroup(self)
        self._tool_group.setExclusive(True)
        self._tool_buttons: List[QToolButton] = []
        for i, (label, tip) in enumerate([
            ("Split", "Split a label into two using two or more seed points"),
            ("Merge", "Merge two or more labels into one"),
            ("Erode", "Erode the selected label by N pixels (3D)"),
            ("Dilate", "Expand the selected label by N pixels without crossing into other labels"),
            ("Delete", "Erase the selected label across the chosen time range"),
            ("Create", "Create new label(s) from seed points placed in unlabelled background pixels"),
        ]):
            btn = QToolButton()
            btn.setText(label)
            btn.setToolTip(tip)
            btn.setCheckable(True)
            btn.setStyleSheet(
                "QToolButton{padding:4px 8px;border:1px solid palette(mid);border-radius:3px;}"
                "QToolButton:checked{background:#1976D2;color:white;border-color:#0D47A1;}"
            )
            self._tool_group.addButton(btn, i)
            self._tool_buttons.append(btn)
            toolbar.addWidget(btn)
        toolbar.addStretch(1)
        layout.addLayout(toolbar)
        self._tool_group.idClicked.connect(self._on_tool_changed)

        # ---- Tool stack ---------------------------------------------
        self.tool_stack = QStackedWidget()
        layout.addWidget(self.tool_stack)
        self.tool_stack.addWidget(self._build_split_page())
        self.tool_stack.addWidget(self._build_merge_page())
        self.tool_stack.addWidget(self._build_erode_page())
        self.tool_stack.addWidget(self._build_dilate_page())
        self.tool_stack.addWidget(self._build_delete_page())
        self.tool_stack.addWidget(self._build_create_page())
        self._tool_buttons[0].setChecked(True)
        self.tool_stack.setCurrentIndex(0)

        # ---- Time range ---------------------------------------------
        rng_group = QGroupBox("Time range")
        rng_lay = QVBoxLayout(rng_group)
        rng_lay.setContentsMargins(6, 4, 6, 4)
        rng_lay.setSpacing(2)
        self.rb_full = QRadioButton("Full lifetime of the selected label(s) — recommended")
        self.rb_full.setChecked(True)
        self.rb_custom = QRadioButton("Custom range")
        rng_lay.addWidget(self.rb_full)
        rng_lay.addWidget(self.rb_custom)
        custom_row = QHBoxLayout()
        custom_row.addWidget(QLabel("Start T:"))
        self.spin_t0 = QSpinBox()
        self.spin_t0.setRange(0, max(0, int(self.buffer.shape[0]) - 1))
        self.spin_t0.setEnabled(False)
        custom_row.addWidget(self.spin_t0)
        custom_row.addSpacing(8)
        custom_row.addWidget(QLabel("End T:"))
        self.spin_t1 = QSpinBox()
        self.spin_t1.setRange(0, max(0, int(self.buffer.shape[0]) - 1))
        self.spin_t1.setValue(max(0, int(self.buffer.shape[0]) - 1))
        self.spin_t1.setEnabled(False)
        custom_row.addWidget(self.spin_t1)
        custom_row.addStretch(1)
        rng_lay.addLayout(custom_row)
        self.rb_custom.toggled.connect(
            lambda on: (self.spin_t0.setEnabled(on), self.spin_t1.setEnabled(on))
        )
        layout.addWidget(rng_group)

        # ---- History row ---------------------------------------------
        hist_row = QHBoxLayout()
        self.btn_undo = QPushButton("↶ Undo")
        self.btn_undo.setToolTip("Revert the latest layer change (in memory)")
        self.btn_undo.clicked.connect(self._on_undo)
        self.btn_redo = QPushButton("↷ Redo")
        self.btn_redo.clicked.connect(self._on_redo)
        hist_row.addWidget(self.btn_undo)
        hist_row.addWidget(self.btn_redo)
        self.lbl_history = QLabel("0 pending changes")
        self.lbl_history.setStyleSheet("color:#666; font-style:italic;")
        hist_row.addWidget(self.lbl_history, stretch=1)
        layout.addLayout(hist_row)

        # ---- Log -----------------------------------------------------
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(100)
        self.log.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(self.log)

        # ---- Progress bar (hidden until an operation is running) -----
        prog_row = QHBoxLayout()
        self.lbl_progress = QLabel("Running…")
        self.lbl_progress.setStyleSheet("color:#555; font-size:11px; min-width:120px;")
        self.lbl_progress.setVisible(False)
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 0)   # indeterminate by default
        self.progress_bar.setMaximumHeight(14)
        self.progress_bar.setTextVisible(True)
        self.progress_bar.setVisible(False)
        prog_row.addWidget(self.lbl_progress)
        prog_row.addWidget(self.progress_bar, stretch=1)
        layout.addLayout(prog_row)

        # ---- Footer --------------------------------------------------
        foot = QHBoxLayout()
        self.btn_save = QPushButton("💾 Save tracked segments")
        self.btn_save.setStyleSheet(
            "background-color:#28a745; color:white; font-weight:bold; padding:6px;"
        )
        self.btn_save.clicked.connect(self._on_save)
        self.btn_discard = QPushButton("✗ Discard changes")
        self.btn_discard.setStyleSheet(
            "background-color:#FB8C00; color:white; font-weight:bold; padding:6px;"
        )
        self.btn_discard.clicked.connect(self._on_discard)
        self.btn_stop = QPushButton("Stop editing")
        self.btn_stop.setToolTip("Close this editor (with a Save / Discard prompt if needed)")
        self.btn_stop.clicked.connect(self._on_stop_clicked)
        foot.addWidget(self.btn_save)
        foot.addWidget(self.btn_discard)
        foot.addWidget(self.btn_stop)
        layout.addLayout(foot)

        # Initial enable state.
        self._refresh_buttons()

    # ---- Tool pages ---------------------------------------------------
    def _build_split_page(self) -> QWidget:
        page = QWidget()
        lay = QVBoxLayout(page)
        lay.setContentsMargins(2, 2, 2, 2)
        info = QLabel(
            "Drop ≥ 2 seed points on the selected label, then Apply.\n"
            "Sub-labels propagate forward and backward across the track's lifetime."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color:#555; font-style:italic;")
        lay.addWidget(info)
        seed_row = QHBoxLayout()
        self.btn_seeds_layer = QPushButton("Add seed points")
        self.btn_seeds_layer.setToolTip(
            "Create / focus the temporary points layer used to seed the split."
        )
        self.btn_seeds_layer.clicked.connect(self._ensure_seeds_layer)
        seed_row.addWidget(self.btn_seeds_layer)
        self.lbl_seed_count = QLabel("0 seed(s) on current frame")
        self.lbl_seed_count.setStyleSheet("color:#1976D2;")
        seed_row.addWidget(self.lbl_seed_count, stretch=1)
        lay.addLayout(seed_row)
        seed_tutorial_row = QHBoxLayout()
        seed_tutorial_row.addWidget(
            TutorialButton("? How to place seeds", "Placing seed points", _SEED_TUTORIAL_STEPS)
        )
        seed_tutorial_row.addStretch(1)
        lay.addLayout(seed_tutorial_row)
        keep_row = QHBoxLayout()
        self.cb_keep_first = QCheckBox("First seed keeps the original TrackID")
        self.cb_keep_first.setChecked(True)
        keep_row.addWidget(self.cb_keep_first)
        keep_row.addStretch(1)
        lay.addLayout(keep_row)
        self.btn_split_preview = QPushButton("👁 Preview split (single timepoint)")
        self.btn_split_preview.setToolTip(
            "Run the split on the current frame only (no propagation) so you\n"
            "can verify the watershed result before committing to the full lifetime."
        )
        self.btn_split_preview.setStyleSheet(
            "QPushButton{color:#1565C0;border:1px solid #90CAF9;"
            "border-radius:3px;padding:4px;}"
            "QPushButton:hover{background:#E3F2FD;}"
            "QPushButton:disabled{color:#aaa;border-color:#ddd;}"
        )
        self.btn_split_preview.clicked.connect(self._preview_split)
        lay.addWidget(self.btn_split_preview)
        self.btn_split_apply = QPushButton("Apply split")
        self.btn_split_apply.clicked.connect(self._apply_split)
        lay.addWidget(self.btn_split_apply)
        lay.addStretch(1)
        return page

    def _build_merge_page(self) -> QWidget:
        page = QWidget()
        lay = QVBoxLayout(page)
        lay.setContentsMargins(2, 2, 2, 2)
        info = QLabel(
            "Pick 2 or more labels (left-click on the layer) and choose a target\n"
            "ID — every other selected label is rewritten to the target."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color:#555; font-style:italic;")
        lay.addWidget(info)
        row = QHBoxLayout()
        row.addWidget(QLabel("Target ID:"))
        self.combo_merge_target = QComboBox()
        self.combo_merge_target.setMinimumWidth(120)
        row.addWidget(self.combo_merge_target)
        row.addStretch(1)
        lay.addLayout(row)
        self.btn_merge_preview = QPushButton("👁 Preview merge (single timepoint)")
        self.btn_merge_preview.setToolTip(
            "Preview the merge on the current frame only so you can verify the\n"
            "result before committing to the full time range."
        )
        self.btn_merge_preview.setStyleSheet(
            "QPushButton{color:#1565C0;border:1px solid #90CAF9;"
            "border-radius:3px;padding:4px;}"
            "QPushButton:hover{background:#E3F2FD;}"
            "QPushButton:disabled{color:#aaa;border-color:#ddd;}"
        )
        self.btn_merge_preview.clicked.connect(self._preview_merge)
        lay.addWidget(self.btn_merge_preview)
        self.btn_merge_apply = QPushButton("Apply merge")
        self.btn_merge_apply.clicked.connect(self._apply_merge)
        lay.addWidget(self.btn_merge_apply)
        lay.addStretch(1)
        return page

    def _build_erode_page(self) -> QWidget:
        page = QWidget()
        lay = QVBoxLayout(page)
        lay.setContentsMargins(2, 2, 2, 2)
        info = QLabel("Erode the selected label with a 3D isotropic ball.")
        info.setWordWrap(True)
        info.setStyleSheet("color:#555; font-style:italic;")
        lay.addWidget(info)
        row = QHBoxLayout()
        row.addWidget(QLabel("Radius XY (px):"))
        self.spin_erode_xy = QSpinBox()
        self.spin_erode_xy.setRange(0, 50)
        self.spin_erode_xy.setValue(1)
        row.addWidget(self.spin_erode_xy)
        row.addSpacing(8)
        self.lbl_erode_info = QLabel()
        self.lbl_erode_info.setStyleSheet("color:#555;")
        row.addWidget(self.lbl_erode_info)
        row.addStretch(1)
        lay.addLayout(row)
        self.btn_erode_preview = QPushButton("👁 Preview erosion (single timepoint)")
        self.btn_erode_preview.setToolTip(
            "Preview the erosion on the current frame only so you can tune\n"
            "the radius before committing to the full lifetime."
        )
        self.btn_erode_preview.setStyleSheet(
            "QPushButton{color:#1565C0;border:1px solid #90CAF9;"
            "border-radius:3px;padding:4px;}"
            "QPushButton:hover{background:#E3F2FD;}"
            "QPushButton:disabled{color:#aaa;border-color:#ddd;}"
        )
        self.btn_erode_preview.clicked.connect(self._preview_erode)
        lay.addWidget(self.btn_erode_preview)
        self.btn_erode_apply = QPushButton("Apply erosion")
        self.btn_erode_apply.clicked.connect(self._apply_erode)
        lay.addWidget(self.btn_erode_apply)
        lay.addStretch(1)
        self.spin_erode_xy.valueChanged.connect(self._update_erode_info)
        self._update_erode_info()
        return page

    def _build_dilate_page(self) -> QWidget:
        page = QWidget()
        lay = QVBoxLayout(page)
        lay.setContentsMargins(2, 2, 2, 2)
        info = QLabel(
            "Dilate / expand the selected label with a 3D isotropic ball.\n"
            "Expansion never crosses into another label."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color:#555; font-style:italic;")
        lay.addWidget(info)
        row = QHBoxLayout()
        row.addWidget(QLabel("Radius XY (px):"))
        self.spin_dilate_xy = QSpinBox()
        self.spin_dilate_xy.setRange(0, 50)
        self.spin_dilate_xy.setValue(1)
        row.addWidget(self.spin_dilate_xy)
        row.addSpacing(8)
        self.lbl_dilate_info = QLabel()
        self.lbl_dilate_info.setStyleSheet("color:#555;")
        row.addWidget(self.lbl_dilate_info)
        row.addStretch(1)
        lay.addLayout(row)
        self.btn_dilate_preview = QPushButton("👁 Preview dilation (single timepoint)")
        self.btn_dilate_preview.setToolTip(
            "Preview the dilation on the current frame only so you can tune\n"
            "the radius before committing to the full lifetime."
        )
        self.btn_dilate_preview.setStyleSheet(
            "QPushButton{color:#1565C0;border:1px solid #90CAF9;"
            "border-radius:3px;padding:4px;}"
            "QPushButton:hover{background:#E3F2FD;}"
            "QPushButton:disabled{color:#aaa;border-color:#ddd;}"
        )
        self.btn_dilate_preview.clicked.connect(self._preview_dilate)
        lay.addWidget(self.btn_dilate_preview)
        self.btn_dilate_apply = QPushButton("Apply dilation")
        self.btn_dilate_apply.clicked.connect(self._apply_dilate)
        lay.addWidget(self.btn_dilate_apply)
        lay.addStretch(1)
        self.spin_dilate_xy.valueChanged.connect(self._update_dilate_info)
        self._update_dilate_info()
        return page

    def _build_delete_page(self) -> QWidget:
        page = QWidget()
        lay = QVBoxLayout(page)
        lay.setContentsMargins(2, 2, 2, 2)
        info = QLabel("Erase the selected label across the chosen time range.")
        info.setWordWrap(True)
        info.setStyleSheet("color:#555; font-style:italic;")
        lay.addWidget(info)
        self.cb_delete_confirm = QCheckBox("Yes, erase the selected label(s)")
        lay.addWidget(self.cb_delete_confirm)
        self.btn_delete_preview = QPushButton("👁 Preview delete (single timepoint)")
        self.btn_delete_preview.setToolTip(
            "Preview the deletion on the current frame only (confirmation checkbox\n"
            "not required) so you can verify before erasing the full time range."
        )
        self.btn_delete_preview.setStyleSheet(
            "QPushButton{color:#1565C0;border:1px solid #90CAF9;"
            "border-radius:3px;padding:4px;}"
            "QPushButton:hover{background:#E3F2FD;}"
            "QPushButton:disabled{color:#aaa;border-color:#ddd;}"
        )
        self.btn_delete_preview.clicked.connect(self._preview_delete)
        lay.addWidget(self.btn_delete_preview)
        self.btn_delete_apply = QPushButton("Apply delete")
        self.btn_delete_apply.clicked.connect(self._apply_delete)
        lay.addWidget(self.btn_delete_apply)
        lay.addStretch(1)
        return page

    def _build_create_page(self) -> QWidget:
        page = QWidget()
        lay = QVBoxLayout(page)
        lay.setContentsMargins(2, 2, 2, 2)
        info = QLabel(
            "Drop ≥ 1 seed point(s) on unlabelled (background) pixels, then Apply.\n"
            "Each seed creates one new label. Sub-labels propagate forward and\n"
            "backward across the chosen time range."
        )
        info.setWordWrap(True)
        info.setStyleSheet("color:#555; font-style:italic;")
        lay.addWidget(info)
        method_note = QLabel(
            "Note: segmentation uses Otsu thresholding + image-guided watershed "
            "on the selected channel (or Voronoi fallback when none is chosen). "
            "Tracking uses the same propagation-based approach as the rest of "
            "BEHAV3D — watershed from previous-frame markers, no appearance model."
        )
        method_note.setWordWrap(True)
        method_note.setStyleSheet("color:#888; font-size:10px; font-style:italic;")
        lay.addWidget(method_note)
        ch_row = QHBoxLayout()
        ch_row.addWidget(QLabel("Image channel:"))
        self.combo_create_channel = QComboBox()
        self.combo_create_channel.setToolTip(
            "Image channel used to guide the watershed and restrict the new label\n"
            "to bright (organoid-like) voxels via an automatic Otsu threshold.\n"
            "Select '(none)' to fall back to the Voronoi partition of all background."
        )
        ch_row.addWidget(self.combo_create_channel, stretch=1)
        lay.addLayout(ch_row)
        seed_row = QHBoxLayout()
        self.btn_create_seeds_layer = QPushButton("Add seed points")
        self.btn_create_seeds_layer.setToolTip(
            "Create / focus the temporary points layer used to seed the new label(s).\n"
            "Place seeds only on unlabelled voxels (bright ones when a channel is selected)."
        )
        self.btn_create_seeds_layer.clicked.connect(self._ensure_seeds_layer)
        seed_row.addWidget(self.btn_create_seeds_layer)
        self.lbl_create_seed_count = QLabel("0 seed(s) on current frame")
        self.lbl_create_seed_count.setStyleSheet("color:#1976D2;")
        seed_row.addWidget(self.lbl_create_seed_count, stretch=1)
        lay.addLayout(seed_row)
        create_tutorial_row = QHBoxLayout()
        create_tutorial_row.addWidget(
            TutorialButton("? How to place seeds", "Placing seed points", _SEED_TUTORIAL_STEPS)
        )
        create_tutorial_row.addStretch(1)
        lay.addLayout(create_tutorial_row)
        self.btn_create_preview = QPushButton("👁 Preview create (single timepoint)")
        self.btn_create_preview.setToolTip(
            "Run the watershed on the current frame only so you can verify\n"
            "the label boundaries before committing to the full time range."
        )
        self.btn_create_preview.setStyleSheet(
            "QPushButton{color:#1565C0;border:1px solid #90CAF9;"
            "border-radius:3px;padding:4px;}"
            "QPushButton:hover{background:#E3F2FD;}"
            "QPushButton:disabled{color:#aaa;border-color:#ddd;}"
        )
        self.btn_create_preview.clicked.connect(self._preview_create)
        lay.addWidget(self.btn_create_preview)
        self.btn_create_apply = QPushButton("Apply create")
        self.btn_create_apply.clicked.connect(self._apply_create)
        lay.addWidget(self.btn_create_apply)
        lay.addStretch(1)
        self._refresh_create_channels()
        return page

    def _refresh_create_channels(self) -> None:
        """Populate the channel combobox from the Image layers currently in the viewer."""
        if not hasattr(self, "combo_create_channel"):
            return
        try:
            import napari.layers as _nl
        except Exception:
            return
        self.combo_create_channel.clear()
        self.combo_create_channel.addItem("(none — Voronoi fallback)", userData=None)
        prefix = self.sample_name + " \u2013 "  # en-dash, same as layer name convention
        for layer in self.viewer.layers:
            if isinstance(layer, _nl.Image) and layer.name.startswith(prefix):
                self.combo_create_channel.addItem(layer.name, userData=layer.name)

    # ------------------------------------------------------------------
    # Background materialisation (no-op: kept for API compatibility)
    # ------------------------------------------------------------------
    def start_materialisation(self) -> None:
        """No-op: full RAM materialisation is no longer required.

        The layer is kept as a dask array throughout editing.  After each
        commit, undo, or redo, :meth:`_refresh_layer_data` rebuilds
        ``layer.data`` as a dask concatenation where only the dirty frames
        are served from in-memory numpy arrays; all other frames are read
        lazily from the on-disk zarr via :attr:`buffer._darr`.

        The method is preserved so that existing callers
        (e.g. ``_visualization._open_editor``) do not need to be changed.
        """

    def _on_operation_progress(self, current: int, total: int) -> None:
        """Update the progress bar on the Qt thread from the worker's signal."""
        if not hasattr(self, "progress_bar"):
            return
        if total > 0:
            self.progress_bar.setRange(0, total)
            self.progress_bar.setValue(current)
            if hasattr(self, "lbl_progress"):
                self.lbl_progress.setText(f"Frame {current}/{total}")
        self.progress_bar.setVisible(True)
        if hasattr(self, "lbl_progress"):
            self.lbl_progress.setVisible(True)

    # ------------------------------------------------------------------
    # napari activity-dock progress helpers
    # ------------------------------------------------------------------
    def _make_activity_progress(self, desc: str = "Editing…"):
        """Create an indeterminate napari activity-dock progress bar.

        ``total=0`` makes Qt render the bar as a spinning busy indicator until
        we know the real frame count.  The real total is set on the first
        progress signal via :meth:`_update_activity_progress`.
        """
        if _NapariProgress is None:
            return None
        try:
            pbr = _NapariProgress(total=0, desc=desc)
            # Try to expand the activity panel.  The attribute name for the
            # activity dock differs across napari versions; try each in turn.
            self._show_napari_activity_panel()
            return pbr
        except Exception:
            return None

    def _show_napari_activity_panel(self) -> None:
        """Best-effort attempt to make the napari activity panel visible."""
        try:
            win = self.viewer.window._qt_window
        except Exception:
            return
        for attr in (
            "_activity_dialog",
            "_qt_activity",
            "_activity_dock",
        ):
            obj = getattr(win, attr, None)
            if obj is not None and hasattr(obj, "show"):
                try:
                    obj.show()
                    obj.raise_()
                except Exception:
                    pass
                return
        try:
            from qtpy.QtWidgets import QDockWidget, QDialog
            for child in win.findChildren((QDockWidget, QDialog)):
                if "activity" in (child.objectName() or "").lower():
                    child.show()
                    child.raise_()
                    return
        except Exception:
            pass

    def _hide_napari_activity_panel(self) -> None:
        """Best-effort attempt to hide/collapse the napari activity panel."""
        try:
            win = self.viewer.window._qt_window
        except Exception:
            return
        for attr in (
            "_activity_dialog",
            "_qt_activity",
            "_activity_dock",
        ):
            obj = getattr(win, attr, None)
            if obj is not None and hasattr(obj, "hide"):
                try:
                    obj.hide()
                except Exception:
                    pass
                return
        try:
            from qtpy.QtWidgets import QDockWidget, QDialog
            for child in win.findChildren((QDockWidget, QDialog)):
                if "activity" in (child.objectName() or "").lower():
                    child.hide()
                    return
        except Exception:
            pass

    def _update_activity_progress(self, current: int, total: int) -> None:
        """Slot connected to the worker progress signal; updates the activity bar.

        On the first call the bar is re-initialised with the correct total via
        ``reset()`` so the Qt widget's maximum is set properly and the
        percentage is meaningful.  Subsequent calls use ``update(delta)`` —
        the canonical tqdm API — which calls ``display()`` internally and
        guarantees the Qt widget is refreshed.
        """
        pbr = self._activity_progress
        if pbr is None:
            return
        try:
            if total > 0 and pbr.total != total:
                # reset() sets n=0 and updates total in one shot; this also
                # updates the Qt progress bar's maximum correctly.
                try:
                    pbr.reset(total=total)
                except Exception:
                    pbr.total = total
                    pbr.n = 0
            delta = current - pbr.n
            if delta > 0:
                pbr.update(delta)
            elif delta < 0:
                # Shouldn't happen in normal use, but guard anyway.
                pbr.n = current
                pbr.refresh()
        except Exception:
            pass

    def _close_activity_progress(self) -> None:
        """Fill the bar to 100 %, wait 1 s, then close the activity entry."""
        pbr = self._activity_progress
        self._activity_progress = None
        if pbr is None:
            return
        try:
            # Drive to 100 % so the bar shows completion before the delay.
            if pbr.total and pbr.total > 0 and pbr.n < pbr.total:
                pbr.update(pbr.total - pbr.n)
        except Exception:
            pass

        def _do_close():
            try:
                pbr.close()
            except Exception:
                pass
            self._hide_napari_activity_panel()

        QTimer.singleShot(1000, _do_close)

    def _set_editing_enabled(self, enabled: bool) -> None:
        """Enable or disable all editing action buttons."""
        for btn in (
            getattr(self, "btn_split_preview", None),
            getattr(self, "btn_split_apply", None),
            getattr(self, "btn_merge_preview", None),
            getattr(self, "btn_merge_apply", None),
            getattr(self, "btn_erode_preview", None),
            getattr(self, "btn_erode_apply", None),
            getattr(self, "btn_dilate_preview", None),
            getattr(self, "btn_dilate_apply", None),
            getattr(self, "btn_delete_preview", None),
            getattr(self, "btn_delete_apply", None),
            getattr(self, "btn_create_preview", None),
            getattr(self, "btn_create_apply", None),
            getattr(self, "btn_undo", None),
            getattr(self, "btn_redo", None),
            getattr(self, "btn_save", None),
            getattr(self, "btn_seeds_layer", None),
            getattr(self, "btn_create_seeds_layer", None),
        ):
            if btn is not None:
                btn.setEnabled(enabled)
        combo = getattr(self, "combo_create_channel", None)
        if combo is not None:
            combo.setEnabled(enabled)
        # Show/hide progress bar.
        if not enabled:
            pb = getattr(self, "progress_bar", None)
            lbl = getattr(self, "lbl_progress", None)
            if pb is not None:
                pb.setRange(0, 0)   # indeterminate until first progress tick
                pb.setValue(0)
                pb.setVisible(True)
            if lbl is not None:
                lbl.setText("Running…")
                lbl.setVisible(True)
        else:
            pb = getattr(self, "progress_bar", None)
            lbl = getattr(self, "lbl_progress", None)
            if pb is not None:
                pb.setVisible(False)
            if lbl is not None:
                lbl.setVisible(False)
        # When re-enabling, restore the proper enable state (selection-aware).
        if enabled:
            self._refresh_buttons()

    # ------------------------------------------------------------------
    # Background edit operation management
    # ------------------------------------------------------------------
    def _run_operation_async(self, fn, args, kwargs=None, on_success=None,
                             desc: str = "Editing…") -> None:
        """Run ``fn(*args, **kwargs)`` in a background QThread.

        Disables all edit buttons while running and re-enables them when the
        thread finishes.  ``fn`` must return an
        :class:`~behav3d.editing.OpResult` or a list of them; the result is
        applied to :attr:`buffer` on the Qt thread via ``_on_edit_done``.

        Parameters
        ----------
        fn:
            Callable that performs the editing primitive.  Must not call
            ``buf.apply()`` itself — that is done by ``_on_edit_done``.
        args, kwargs:
            Positional and keyword arguments forwarded to ``fn``.
        on_success:
            Optional callable invoked with the result after every
            :meth:`_commit_op` succeeds (e.g. clear split seeds).
        desc:
            Short description shown in the napari activity-dock progress bar.
        """
        if kwargs is None:
            kwargs = {}
        self._pending_on_success = on_success
        self._activity_progress = self._make_activity_progress(desc)

        worker = _EditWorker(fn, args, kwargs)
        # Inject progress callback — primitives expose it as an optional kwarg.
        # We add it to the *same* dict that the worker already holds so the
        # worker will see it when it calls fn(**kwargs) on the thread.
        kwargs['progress_cb'] = worker.progress.emit
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.done.connect(self._on_edit_done)
        worker.failed.connect(self._on_edit_failed)
        worker.finished.connect(thread.quit)
        worker.progress.connect(self._on_operation_progress)
        worker.progress.connect(self._update_activity_progress)
        thread.finished.connect(self._on_edit_thread_finished)

        self._edit_thread = thread
        self._edit_worker = worker
        self._set_editing_enabled(False)
        thread.start()

    def _on_edit_done(self, result) -> None:
        """Called on the Qt thread when the operation has produced a result.

        ``result`` is either a single :class:`~behav3d.editing.OpResult` or a
        list of them (delete runs one per selected label).  Each is applied and
        committed here on the Qt thread so that ``EditBuffer`` is never written
        from a background thread.
        """
        self._close_activity_progress()
        ops = result if isinstance(result, list) else [result]
        any_committed = False
        for op in ops:
            try:
                self._commit_op(op)
                any_committed = True
            except MemoryError:
                n = len(op.new_frames)
                self._log_msg(
                    f"  ❌ Out of memory committing '{op.name}' "
                    f"({n} frame(s) × {tuple(self.buffer.shape[1:])}) — "
                    "use a smaller time range or reduce the radius."
                )
            except Exception as exc:
                self._log_msg(f"  ❌ Commit failed ({op.name}): {exc}")
        if self._pending_on_success is not None:
            try:
                self._pending_on_success(result)
            except Exception:
                pass
            self._pending_on_success = None
        self._set_editing_enabled(True)

    def _on_edit_failed(self, error: str) -> None:
        """Called on the Qt thread when the operation raised an exception."""
        self._close_activity_progress()
        self._log_msg(f"  ❌ Operation failed: {error}")
        self._set_editing_enabled(True)

    def _on_edit_thread_finished(self) -> None:
        """Called on the Qt thread when the worker OS thread has fully stopped.

        The only safe place to drop the Python references — clearing them
        earlier (inside ``_on_edit_done``/``_on_edit_failed``) would release
        the last reference to the QThread while the OS thread is still running.
        """
        self._edit_thread = None
        self._edit_worker = None

    # ------------------------------------------------------------------
    # Preview infrastructure
    # ------------------------------------------------------------------
    def _clear_preview(self) -> None:
        """Discard any active single-timepoint preview without touching the buffer."""
        if not self._preview_dirty:
            return
        frames = list(self._preview_dirty.keys())
        self._preview_dirty.clear()
        self._preview_tool = None
        self._refresh_layer_data(frames)

    def _run_preview_async(self, fn, args, kwargs=None, tool_name: str = "") -> None:
        """Run ``fn(*args, **kwargs)`` in a background QThread for preview only.

        The result is placed in ``_preview_dirty`` and shown in the napari layer
        without committing to the :class:`~behav3d.editing.EditBuffer` or the
        undo stack.  Buttons are disabled during execution (same as the real
        operations) to avoid race conditions.
        """
        if kwargs is None:
            kwargs = {}
        self._pending_preview_tool = tool_name
        self._activity_progress = self._make_activity_progress(
            f"Preview {tool_name}…" if tool_name else "Previewing…"
        )
        worker = _EditWorker(fn, args, kwargs)
        kwargs['progress_cb'] = worker.progress.emit
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.done.connect(self._on_preview_done)
        worker.failed.connect(self._on_preview_failed)
        worker.finished.connect(thread.quit)
        worker.progress.connect(self._on_operation_progress)
        worker.progress.connect(self._update_activity_progress)
        thread.finished.connect(self._on_edit_thread_finished)
        self._edit_thread = thread
        self._edit_worker = worker
        self._set_editing_enabled(False)
        thread.start()

    def _on_preview_done(self, result) -> None:
        self._close_activity_progress()
        tool_name = self._pending_preview_tool or "operation"
        self._pending_preview_tool = None
        self._preview_dirty.clear()
        ops = result if isinstance(result, list) else [result]
        frames_changed: List[int] = []
        for op in ops:
            for t, frame in op.new_frames.items():
                self._preview_dirty[int(t)] = np.asarray(frame).astype(
                    self.buffer.dtype, copy=False
                )
                frames_changed.append(int(t))
        self._preview_tool = tool_name
        if frames_changed:
            self._refresh_layer_data(frames_changed)
            ts = ", ".join(str(t) for t in sorted(frames_changed))
            self._log_msg(
                f"  👁 Preview {tool_name} at t={ts} — "
                f"click Apply to commit the full-lifetime operation, or switch "
                f"tools / undo to clear the preview."
            )
        else:
            self._log_msg(
                f"  ⚠️ Preview {tool_name}: no changes at this timepoint."
            )
        self._set_editing_enabled(True)

    def _on_preview_failed(self, error: str) -> None:
        self._close_activity_progress()
        self._pending_preview_tool = None
        self._log_msg(f"  ❌ Preview failed: {error}")
        self._set_editing_enabled(True)

    # ------------------------------------------------------------------
    # Layer binding
    # ------------------------------------------------------------------
    def _bind_to_layer(self) -> None:
        layer = self.viewer.layers[self.layer_name]
        layer.show_selected_label = False
        try:
            layer.n_edit_dimensions = 3  # never paint across time
        except Exception:
            pass
        try:
            # Default to pan/zoom so single clicks just select labels;
            # the user can still pick paint/erase from napari's controls.
            layer.mode = "pan_zoom"
        except Exception:
            pass
        layer.editable = True
        layer.mouse_drag_callbacks.append(self._on_layer_clicked)
        # The caller (visualization tab / notebook panel) is responsible
        # for loading the layer's data as a writable numpy array — we
        # deliberately do NOT replace ``layer.data`` here, because that
        # blanks the layer in some napari versions.  If the caller passed
        # a dask-backed layer we promote it on first edit (see
        # ``_refresh_layer_data``).
        try:
            self.viewer.dims.events.current_step.connect(
                lambda *_: self._refresh_seed_counter()
            )
            self.viewer.dims.events.current_step.connect(
                lambda *_: self._refresh_highlight_layer()
            )
        except Exception:
            pass
        try:
            self.viewer.layers.events.removed.connect(self._on_layer_removed)
        except Exception:
            pass
        try:
            self.viewer.layers.selection.active = self.viewer.layers[self.layer_name]
        except Exception:
            pass

    def _on_layer_removed(self, event) -> None:
        """Keep the tracked-segments layer active when the seed-points
        layer is removed.

        This covers deletion via the Apply/Cancel flow *and* manual
        deletion from napari's own layer list (e.g. the trash icon or
        the Delete key), which otherwise leaves napari to fall back to
        whatever layer happens to be first in the list.
        """
        removed = getattr(event, "value", None)
        if getattr(removed, "name", None) != _SPLIT_SEEDS_LAYER:
            return
        # The layer is already gone at this point, so the counter must be
        # refreshed explicitly — it won't get another data/current_step
        # event to react to.
        self._refresh_seed_counter()
        # Defer: napari's LayerList does its own "pick a new active layer"
        # bookkeeping in response to the same ``removed`` event, and
        # depending on connection order that can run *after* us — which
        # would leave both the tracked-segments layer and napari's default
        # pick selected. Running on the next event-loop tick guarantees we
        # go last and end up with exactly one layer selected.
        QTimer.singleShot(0, self._reselect_tracked_layer)

    def _reselect_tracked_layer(self) -> None:
        if self.layer_name in self.viewer.layers:
            layer = self.viewer.layers[self.layer_name]
            self.viewer.layers.selection.clear()
            self.viewer.layers.selection.add(layer)
            self.viewer.layers.selection.active = layer

    def _build_dask_with_dirty_frames(self) -> da.Array:
        """Return a dask array where dirty frames are served from memory.

        Clean frames delegate to the on-disk zarr via ``buffer._darr``;
        dirty frames are wrapped as single-frame dask chunks from the
        in-memory numpy arrays in ``buffer._dirty``.  Building the task
        graph is O(T) Python (no data is loaded), so it is safe to call
        after every commit, undo, or redo.
        """
        # Merge buffer dirty frames with any active preview frames; preview
        # takes precedence so the napari layer reflects the preview overlay.
        dirty = dict(self.buffer._dirty)
        dirty.update(self._preview_dirty)
        base = self.buffer._darr
        if not dirty:
            return base
        shape = self.buffer.shape
        slices: List[da.Array] = []
        for t in range(shape[0]):
            if t in dirty:
                slices.append(
                    da.from_array(
                        dirty[t][np.newaxis],
                        chunks=(1, *shape[1:]),
                    )
                )
            else:
                slices.append(base[t : t + 1])
        return da.concatenate(slices, axis=0)

    def _refresh_layer_data(self, frames: List[int]) -> None:
        if self.layer_name not in self.viewer.layers:
            return
        layer = self.viewer.layers[self.layer_name]
        try:
            layer.data = self._build_dask_with_dirty_frames()
            layer.refresh()
        except Exception:
            traceback.print_exc()

    # ------------------------------------------------------------------
    # Selection / interactions
    # ------------------------------------------------------------------
    def _on_layer_clicked(self, layer, event):
        # Only react to plain clicks, not drags.
        if event.type != "mouse_press":
            return
        try:
            value = layer.get_value(event.position, world=True)
        except Exception:
            return
        if value in (None, 0):
            return
        try:
            label_id = int(value)
        except (TypeError, ValueError):
            return
        if event.button == 1 and "Shift" not in (event.modifiers or []):
            # Replace selection with single label.
            self._set_selection([label_id])
        elif event.button == 1 and "Shift" in (event.modifiers or []):
            # Add to selection.
            cur = list(self._selected_labels)
            if label_id in cur:
                cur.remove(label_id)
            else:
                cur.append(label_id)
            self._set_selection(cur)
        elif event.button == 2:
            # Right-click removes.
            cur = [x for x in self._selected_labels if x != label_id]
            self._set_selection(cur)

    def _add_label_by_id(self) -> None:
        text = self.input_label_id.text().strip()
        if not text:
            return
        try:
            label_id = int(text)
        except ValueError:
            self._log_msg(f"Invalid label ID: {text!r}")
            return
        if label_id <= 0:
            self._log_msg("Label IDs must be positive integers.")
            return
        cur = list(self._selected_labels)
        if label_id not in cur:
            cur.append(label_id)
        self._set_selection(cur)
        self.input_label_id.clear()

    def _clear_selection(self) -> None:
        self._set_selection([])

    def _set_selection(self, labels: List[int]) -> None:
        self._selected_labels = [int(x) for x in labels]
        if self._selected_labels:
            self.lbl_selected.setText(", ".join(str(x) for x in self._selected_labels))
            try:
                layer = self.viewer.layers[self.layer_name]
                layer.selected_label = int(self._selected_labels[0])
            except Exception:
                pass
        else:
            self.lbl_selected.setText("(click on the tracked-segments layer)")
        # Refresh time-range defaults from the lifetime of the first label.
        if self._selected_labels:
            f, l = lifetime_of(self.buffer, self._selected_labels[0])
            if f is not None:
                self.spin_t0.setValue(int(f))
                self.spin_t1.setValue(int(l))
        # Update merge target combobox.
        self.combo_merge_target.clear()
        for lid in self._selected_labels:
            self.combo_merge_target.addItem(str(lid), userData=int(lid))
        self._refresh_buttons()
        self._refresh_seed_counter()
        self._refresh_highlight_layer()

    # ------------------------------------------------------------------
    # Tool dispatch
    # ------------------------------------------------------------------
    def _on_tool_changed(self, idx: int) -> None:
        self._clear_preview()
        self.tool_stack.setCurrentIndex(idx)
        # Seeds layer is shared by Split (0) and Create (5).
        is_seed_tool = idx in (0, 5)
        if not is_seed_tool and _SPLIT_SEEDS_LAYER in self.viewer.layers:
            self.viewer.layers[_SPLIT_SEEDS_LAYER].visible = False
        if not is_seed_tool:
            self._restore_ndisplay()
        if is_seed_tool and _SPLIT_SEEDS_LAYER in self.viewer.layers:
            self.viewer.layers[_SPLIT_SEEDS_LAYER].visible = True
        # Refresh channel list whenever the user switches to Create.
        if idx == 5:
            self._refresh_create_channels()

    def _selected_t_range(self, label_id: Optional[int] = None) -> Optional[Tuple[int, int]]:
        if self.rb_full.isChecked():
            if label_id is not None:
                f, l = lifetime_of(self.buffer, label_id)
                if f is None:
                    return None
                return (int(f), int(l))
            return None
        return (int(self.spin_t0.value()), int(self.spin_t1.value()))

    # ---- Split -------------------------------------------------------
    def _ensure_seeds_layer(self) -> None:
        # Split needs an existing label to seed the watershed from; Create
        # (index 5) deliberately has no such requirement.
        if self.tool_stack.currentIndex() == 0 and not self._selected_labels:
            QMessageBox.warning(
                self,
                "No label selected",
                "Select a label in the viewer (or type its ID) before adding seed points.",
            )
            return
        if self.viewer.dims.ndisplay == 3 and self._ndisplay_before_seeds is None:
            self._ndisplay_before_seeds = 3
            self.viewer.dims.ndisplay = 2
            self._log_msg(
                "  Switched to 2D view for precise seed placement — move "
                "the Z slider between clicks to place seeds at different depths."
            )
        if _SPLIT_SEEDS_LAYER not in self.viewer.layers:
            layer = self.viewer.add_points(
                np.empty((0, 4)),
                name=_SPLIT_SEEDS_LAYER,
                size=6,
                face_color="red",
                ndim=4,
            )
            # Prevent the Points layer's Text sub-visual from triggering
            # FreeType.  Napari has a bug where it never propagates
            # text.visible=False to the underlying vispy node, so we set
            # both the napari TextManager flag and the vispy node directly.
            layer.text.visible = False
            try:
                self.viewer.window._qt_viewer.layer_to_visual[layer].node.text.visible = False
            except Exception:
                pass
            layer.events.data.connect(lambda *_: self._refresh_seed_counter())
        layer = self.viewer.layers[_SPLIT_SEEDS_LAYER]
        layer.mode = "add"
        self.viewer.layers.selection.active = layer
        self._refresh_seed_counter()

    def _restore_ndisplay(self) -> None:
        if self._ndisplay_before_seeds is not None:
            try:
                self.viewer.dims.ndisplay = self._ndisplay_before_seeds
            except Exception:
                pass
            self._ndisplay_before_seeds = None

    def _refresh_seed_counter(self) -> None:
        n = 0
        if _SPLIT_SEEDS_LAYER in self.viewer.layers:
            try:
                seeds = np.asarray(self.viewer.layers[_SPLIT_SEEDS_LAYER].data)
                # Filter by current frame
                t_now = int(self.viewer.dims.current_step[0])
                if seeds.size:
                    n = int((seeds[:, 0].astype(int) == t_now).sum())
            except Exception:
                n = 0
        text = f"{n} seed(s) on current frame"
        if hasattr(self, "lbl_seed_count"):
            self.lbl_seed_count.setText(text)
        if hasattr(self, "lbl_create_seed_count"):
            self.lbl_create_seed_count.setText(text)

    def _refresh_highlight_layer(self) -> None:
        """Show the selected label(s) as a bright contour-only overlay.

        Uses a dedicated Labels layer rather than recoloring
        ``self.layer_name`` in place, so the rest of the tracked-segments
        layer keeps its normal per-object colors.  Whatever layer was
        active before this call (e.g. the seeds layer mid-placement, or
        the tracked-segments layer during label picking) is restored
        afterward so this never steals mouse interaction away from it.
        """
        try:
            prev_active = self.viewer.layers.selection.active
            if not self._selected_labels:
                if _HIGHLIGHT_LAYER in self.viewer.layers:
                    self.viewer.layers.remove(_HIGHLIGHT_LAYER)
                    if (
                        prev_active is not None
                        and getattr(prev_active, "name", None) != _HIGHLIGHT_LAYER
                    ):
                        self.viewer.layers.selection.active = prev_active
                return
            t_now = int(self.viewer.dims.current_step[0])
            mask = np.isin(self.buffer.peek(t_now), self._selected_labels).astype(np.uint8)
            if _HIGHLIGHT_LAYER in self.viewer.layers:
                layer = self.viewer.layers[_HIGHLIGHT_LAYER]
                layer.data = mask
                layer.visible = True
            else:
                layer = self.viewer.add_labels(mask, name=_HIGHLIGHT_LAYER, opacity=1.0)
                layer.editable = False
                layer.contour = 1
                if DirectLabelColormap is not None:
                    layer.colormap = DirectLabelColormap(
                        color_dict={1: _HIGHLIGHT_COLOR, None: "transparent"}
                    )
            if prev_active is not None and getattr(prev_active, "name", None) != _HIGHLIGHT_LAYER:
                self.viewer.layers.selection.active = prev_active
        except Exception:
            pass

    def _apply_split(self) -> None:
        self._clear_preview()
        if not self._selected_labels:
            self._log_msg("Pick a label first (click on the tracked-segments layer).")
            return
        label_id = self._selected_labels[0]
        if _SPLIT_SEEDS_LAYER not in self.viewer.layers:
            self._log_msg("Add seed points first (open the Split tool's seeds layer).")
            return
        seeds_all = np.asarray(self.viewer.layers[_SPLIT_SEEDS_LAYER].data)
        if seeds_all.size == 0:
            self._log_msg("No seeds added.")
            return
        t_ref = int(self.viewer.dims.current_step[0])
        seeds_t = seeds_all[seeds_all[:, 0].astype(int) == t_ref]
        if len(seeds_t) < 2:
            self._log_msg(
                f"Need ≥ 2 seeds on the current frame (t={t_ref}); got {len(seeds_t)}."
            )
            return
        seeds_zyx = [
            (int(round(s[1])), int(round(s[2])), int(round(s[3])))
            for s in seeds_t
        ]
        rng = self._selected_t_range(label_id)
        self._log_msg(f"  Running split of label {label_id} at t={t_ref}…")

        def _clear_seeds(_result):
            try:
                if _SPLIT_SEEDS_LAYER in self.viewer.layers:
                    self.viewer.layers.remove(_SPLIT_SEEDS_LAYER)
            except Exception:
                pass
            try:
                if self.layer_name in self.viewer.layers:
                    self.viewer.layers.selection.active = (
                        self.viewer.layers[self.layer_name]
                    )
            except Exception:
                pass
            self._restore_ndisplay()

        self._run_operation_async(
            split_label,
            (self.buffer,),
            dict(
                label_id=label_id,
                seeds_zyx=seeds_zyx,
                t_ref=t_ref,
                t_range=rng,
                keep_first_id=self.cb_keep_first.isChecked(),
            ),
            on_success=_clear_seeds,
            desc=f"Splitting label {label_id}…",
        )

    def _preview_split(self) -> None:
        if not self._selected_labels:
            self._log_msg("Pick a label first (click on the tracked-segments layer).")
            return
        label_id = self._selected_labels[0]
        if _SPLIT_SEEDS_LAYER not in self.viewer.layers:
            self._log_msg("Add seed points first (open the Split tool's seeds layer).")
            return
        seeds_all = np.asarray(self.viewer.layers[_SPLIT_SEEDS_LAYER].data)
        if seeds_all.size == 0:
            self._log_msg("No seeds added.")
            return
        t_ref = int(self.viewer.dims.current_step[0])
        seeds_t = seeds_all[seeds_all[:, 0].astype(int) == t_ref]
        if len(seeds_t) < 2:
            self._log_msg(
                f"Need ≥ 2 seeds on the current frame (t={t_ref}); got {len(seeds_t)}."
            )
            return
        seeds_zyx = [
            (int(round(s[1])), int(round(s[2])), int(round(s[3])))
            for s in seeds_t
        ]
        self._log_msg(f"  Previewing split of label {label_id} at t={t_ref}…")
        self._run_preview_async(
            split_label,
            (self.buffer,),
            dict(
                label_id=label_id,
                seeds_zyx=seeds_zyx,
                t_ref=t_ref,
                t_range=(t_ref, t_ref),
                keep_first_id=self.cb_keep_first.isChecked(),
            ),
            tool_name="split",
        )

    # ---- Merge -------------------------------------------------------
    def _apply_merge(self) -> None:
        self._clear_preview()
        if len(self._selected_labels) < 2:
            self._log_msg("Pick ≥ 2 labels (Shift-click to add).")
            return
        target = self.combo_merge_target.currentData()
        if target is None:
            target = min(self._selected_labels)
        # Pass t_range=None when "Full lifetime" is selected so that
        # merge_labels computes the union of all selected labels' lifetimes.
        # Using only the target's lifetime misses labels that exist exclusively
        # in timepoints where the target is absent (e.g. non-overlapping tracks).
        if self.rb_full.isChecked():
            rng = None
        else:
            rng = (int(self.spin_t0.value()), int(self.spin_t1.value()))
        self._log_msg(f"  Running merge of {self._selected_labels} → {target}…")
        self._run_operation_async(
            merge_labels,
            (self.buffer,),
            dict(
                label_ids=list(self._selected_labels),
                target_id=int(target),
                t_range=rng,
            ),
            desc=f"Merging labels → {target}…",
        )

    def _preview_merge(self) -> None:
        if len(self._selected_labels) < 2:
            self._log_msg("Pick ≥ 2 labels (Shift-click to add).")
            return
        target = self.combo_merge_target.currentData()
        if target is None:
            target = min(self._selected_labels)
        t_now = int(self.viewer.dims.current_step[0])
        self._log_msg(
            f"  Previewing merge of {self._selected_labels} → {target} at t={t_now}…"
        )
        self._run_preview_async(
            merge_labels,
            (self.buffer,),
            dict(
                label_ids=list(self._selected_labels),
                target_id=int(target),
                t_range=(t_now, t_now),
            ),
            tool_name="merge",
        )

    # ---- Erode / Dilate ---------------------------------------------
    def _iso_r_z(self, r_xy: int) -> int:
        """Z-slice radius that makes the erosion/dilation ball isotropic."""
        pz = self.buffer.pixel_size_z
        if pz <= 0:
            return r_xy
        return max(0, round(r_xy * self.buffer.pixel_size_xy / pz))

    def _update_erode_info(self) -> None:
        r_xy = int(self.spin_erode_xy.value())
        r_z = self._iso_r_z(r_xy)
        um = r_xy * self.buffer.pixel_size_xy
        self.lbl_erode_info.setText(f"= {um:.2f} µm  (Z: {r_z} slice(s))")

    def _update_dilate_info(self) -> None:
        r_xy = int(self.spin_dilate_xy.value())
        r_z = self._iso_r_z(r_xy)
        um = r_xy * self.buffer.pixel_size_xy
        self.lbl_dilate_info.setText(f"= {um:.2f} µm  (Z: {r_z} slice(s))")

    def _apply_erode(self) -> None:
        self._clear_preview()
        if not self._selected_labels:
            self._log_msg("Pick a label first.")
            return
        label_id = self._selected_labels[0]
        rng = self._selected_t_range(label_id)
        r_xy = int(self.spin_erode_xy.value())
        r_z = self._iso_r_z(r_xy)
        self._log_msg(f"  Running erode of label {label_id}…")
        self._run_operation_async(
            erode_label,
            (self.buffer,),
            dict(
                label_id=label_id,
                r_xy=r_xy,
                r_z=r_z,
                t_range=rng,
            ),
            desc=f"Eroding label {label_id}…",
        )

    def _preview_erode(self) -> None:
        if not self._selected_labels:
            self._log_msg("Pick a label first.")
            return
        label_id = self._selected_labels[0]
        t_now = int(self.viewer.dims.current_step[0])
        r_xy = int(self.spin_erode_xy.value())
        r_z = self._iso_r_z(r_xy)
        self._log_msg(f"  Previewing erosion of label {label_id} at t={t_now}…")
        self._run_preview_async(
            erode_label,
            (self.buffer,),
            dict(label_id=label_id, r_xy=r_xy, r_z=r_z, t_range=(t_now, t_now)),
            tool_name="erosion",
        )

    def _apply_dilate(self) -> None:
        self._clear_preview()
        if not self._selected_labels:
            self._log_msg("Pick a label first.")
            return
        label_id = self._selected_labels[0]
        rng = self._selected_t_range(label_id)
        r_xy = int(self.spin_dilate_xy.value())
        r_z = self._iso_r_z(r_xy)
        self._log_msg(f"  Running dilate of label {label_id}…")
        self._run_operation_async(
            dilate_label,
            (self.buffer,),
            dict(
                label_id=label_id,
                r_xy=r_xy,
                r_z=r_z,
                t_range=rng,
            ),
            desc=f"Dilating label {label_id}…",
        )

    def _preview_dilate(self) -> None:
        if not self._selected_labels:
            self._log_msg("Pick a label first.")
            return
        label_id = self._selected_labels[0]
        t_now = int(self.viewer.dims.current_step[0])
        r_xy = int(self.spin_dilate_xy.value())
        r_z = self._iso_r_z(r_xy)
        self._log_msg(f"  Previewing dilation of label {label_id} at t={t_now}…")
        self._run_preview_async(
            dilate_label,
            (self.buffer,),
            dict(label_id=label_id, r_xy=r_xy, r_z=r_z, t_range=(t_now, t_now)),
            tool_name="dilation",
        )

    # ---- Delete ------------------------------------------------------
    def _apply_delete(self) -> None:
        self._clear_preview()
        if not self._selected_labels:
            self._log_msg("Pick a label first.")
            return
        if not self.cb_delete_confirm.isChecked():
            self._log_msg("Tick the confirmation box to delete.")
            return
        label_ids = list(self._selected_labels)
        # Resolve time ranges on the Qt thread (accesses Qt widgets).
        ranges = {lid: self._selected_t_range(lid) for lid in label_ids}
        self._log_msg(f"  Running delete of label(s) {label_ids}…")

        def _do_delete_all(_buf, progress_cb=None):
            n = len(label_ids)
            results = []
            for i, lid in enumerate(label_ids):
                # For a single label forward progress_cb directly so the bar
                # tracks per-frame; for multiple labels track per-label.
                if n == 1:
                    results.append(delete_label(_buf, label_id=lid, t_range=ranges[lid],
                                                progress_cb=progress_cb))
                else:
                    if progress_cb:
                        progress_cb(i, n)
                    results.append(delete_label(_buf, label_id=lid, t_range=ranges[lid]))
                    if progress_cb:
                        progress_cb(i + 1, n)
            return results

        def _uncheck(_result):
            self.cb_delete_confirm.setChecked(False)

        self._run_operation_async(
            _do_delete_all,
            (self.buffer,),
            on_success=_uncheck,
            desc=f"Deleting label(s) {label_ids}…",
        )

    def _preview_delete(self) -> None:
        if not self._selected_labels:
            self._log_msg("Pick a label first.")
            return
        label_ids = list(self._selected_labels)
        t_now = int(self.viewer.dims.current_step[0])
        self._log_msg(f"  Previewing delete of label(s) {label_ids} at t={t_now}…")

        def _do_delete_preview(_buf, progress_cb=None):
            return [
                delete_label(_buf, label_id=lid, t_range=(t_now, t_now),
                             progress_cb=progress_cb)
                for lid in label_ids
            ]

        self._run_preview_async(
            _do_delete_preview,
            (self.buffer,),
            tool_name="delete",
        )

    # ---- Create ------------------------------------------------------
    def _get_create_image_ref(self, t_ref: int) -> Optional[np.ndarray]:
        """Return the ZYX frame at ``t_ref`` from the selected channel, or None."""
        layer_name = (
            self.combo_create_channel.currentData()
            if hasattr(self, "combo_create_channel")
            else None
        )
        if not layer_name:
            return None
        try:
            return np.asarray(self.viewer.layers[layer_name].data[t_ref])
        except Exception:
            return None

    def _get_create_image_stack(self):
        """Return the full T×Z×Y×X channel array for the selected channel, or None.

        Used by :meth:`_apply_create` to pass per-frame intensity guidance to
        :func:`create_label` so that propagation frames use the same
        Otsu + negated-image watershed as the reference frame.
        """
        layer_name = (
            self.combo_create_channel.currentData()
            if hasattr(self, "combo_create_channel")
            else None
        )
        if not layer_name:
            return None
        try:
            return self.viewer.layers[layer_name].data
        except Exception:
            return None

    def _apply_create(self) -> None:
        self._clear_preview()
        if _SPLIT_SEEDS_LAYER not in self.viewer.layers:
            self._log_msg("Add seed points first (open the Create tool's seeds layer).")
            return
        seeds_all = np.asarray(self.viewer.layers[_SPLIT_SEEDS_LAYER].data)
        if seeds_all.size == 0:
            self._log_msg("No seeds added.")
            return
        t_ref = int(self.viewer.dims.current_step[0])
        seeds_t = seeds_all[seeds_all[:, 0].astype(int) == t_ref]
        if len(seeds_t) < 1:
            self._log_msg(
                f"No seeds on the current frame (t={t_ref}). "
                "Navigate to the frame where you placed seeds."
            )
            return
        seeds_zyx = [
            (int(round(s[1])), int(round(s[2])), int(round(s[3])))
            for s in seeds_t
        ]
        rng = self._selected_t_range(label_id=None)
        image_ref = self._get_create_image_ref(t_ref)
        ch_info = (
            f" (channel: {self.combo_create_channel.currentText()})"
            if image_ref is not None and hasattr(self, "combo_create_channel")
            else " (Voronoi fallback)"
        )
        self._log_msg(
            f"  Creating {len(seeds_zyx)} new label(s) from seeds at t={t_ref}{ch_info}…"
        )

        def _clear_seeds(_result):
            try:
                if _SPLIT_SEEDS_LAYER in self.viewer.layers:
                    self.viewer.layers.remove(_SPLIT_SEEDS_LAYER)
            except Exception:
                pass
            try:
                if self.layer_name in self.viewer.layers:
                    self.viewer.layers.selection.active = (
                        self.viewer.layers[self.layer_name]
                    )
            except Exception:
                pass
            self._restore_ndisplay()

        self._run_operation_async(
            create_label,
            (self.buffer,),
            dict(
                seeds_zyx=seeds_zyx,
                t_ref=t_ref,
                t_range=rng,
                image_ref=image_ref,
                image_stack=self._get_create_image_stack(),
            ),
            on_success=_clear_seeds,
            desc=f"Creating {len(seeds_zyx)} label(s)…",
        )

    def _preview_create(self) -> None:
        if _SPLIT_SEEDS_LAYER not in self.viewer.layers:
            self._log_msg("Add seed points first (open the Create tool's seeds layer).")
            return
        seeds_all = np.asarray(self.viewer.layers[_SPLIT_SEEDS_LAYER].data)
        if seeds_all.size == 0:
            self._log_msg("No seeds added.")
            return
        t_ref = int(self.viewer.dims.current_step[0])
        seeds_t = seeds_all[seeds_all[:, 0].astype(int) == t_ref]
        if len(seeds_t) < 1:
            self._log_msg(
                f"No seeds on the current frame (t={t_ref})."
            )
            return
        seeds_zyx = [
            (int(round(s[1])), int(round(s[2])), int(round(s[3])))
            for s in seeds_t
        ]
        image_ref = self._get_create_image_ref(t_ref)
        ch_info = (
            f" (channel: {self.combo_create_channel.currentText()})"
            if image_ref is not None and hasattr(self, "combo_create_channel")
            else " (Voronoi fallback)"
        )
        self._log_msg(
            f"  Previewing create of {len(seeds_zyx)} label(s) at t={t_ref}{ch_info}…"
        )
        self._run_preview_async(
            create_label,
            (self.buffer,),
            dict(
                seeds_zyx=seeds_zyx,
                t_ref=t_ref,
                t_range=(t_ref, t_ref),
                image_ref=image_ref,
            ),
            tool_name="create",
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _unique_darr(darr: da.Array) -> da.Array:
        """Return *darr* re-tagged with a unique dask task-graph token.

        Two ``da.from_zarr`` calls on the same zarr path produce arrays with
        the same auto-generated task-graph keys.  If napari or dask has any
        per-key slice cache, the second array would reuse results from the
        first.  Wrapping with ``map_blocks(identity)`` and a UUID name gives
        the array an unambiguous identity so all caches treat it as fresh.
        """
        import uuid
        try:
            return darr.map_blocks(
                lambda b: b,
                dtype=darr.dtype,
                name=f"zarr_saved_{uuid.uuid4().hex}",
            )
        except Exception:
            return darr

    # ------------------------------------------------------------------
    # Commit / undo / redo / save / discard
    # ------------------------------------------------------------------
    def _commit_op(self, op) -> None:
        if not op.new_frames:
            self._log_msg(f"  ⚠️ {op.summary or op.name}: no frames changed.")
            self._refresh_buttons()
            return
        large = len(op.new_frames) > self.buffer.max_undoable_frames
        self.buffer.apply(op)
        if large:
            self._log_msg(
                f"  ✅ {op.summary}  "
                f"⚠️ Large operation ({len(op.new_frames)} frames) — undo not available."
            )
        else:
            self._log_msg(f"  ✅ {op.summary}")
        self._refresh_buttons()

    def _on_undo(self) -> None:
        self._clear_preview()
        msg = self.buffer.undo()
        self._log_msg(f"  ↶ {msg or 'nothing to undo'}")
        self._refresh_buttons()

    def _on_redo(self) -> None:
        self._clear_preview()
        msg = self.buffer.redo()
        self._log_msg(f"  ↷ {msg or 'nothing to redo'}")
        self._refresh_buttons()

    def _on_save(self) -> None:
        try:
            n, csv_path = self.buffer.save()
        except Exception as exc:
            traceback.print_exc()
            self._log_msg(f"❌ Save failed: {exc}")
            return
        if n == 0:
            self._log_msg("Nothing to save.")
            return
        self._log_msg(
            f"💾 Saved {n} frame(s) → {self._zarr_path.name}"
            + (f"; rewrote CSV {csv_path.name}" if csv_path else "")
        )
        # After saving, _dirty is cleared and buffer._darr is reloaded from the
        # freshly-written zarr.  We need to:
        #   1. Give the new dask array a unique token so napari / dask cannot
        #      mistake it for the pre-save array (same zarr path → same default
        #      task-graph keys → stale slice cache).
        #   2. Assign the array to the layer AND force napari's slicer to
        #      re-fire by doing a brief step-bounce; layer.refresh() alone only
        #      repaints the current OpenGL frame without re-slicing.
        if self.layer_name in self.viewer.layers:
            try:
                layer = self.viewer.layers[self.layer_name]
                fresh_darr = self._unique_darr(self.buffer._darr)
                layer.data = fresh_darr
                # Step-bounce: move one frame away and back so napari's async
                # slicer requests a fresh chunk from the new dask graph for
                # every subsequent frame navigation.
                try:
                    axis = 0
                    T = int(self.buffer.shape[0])
                    t_now = int(self.viewer.dims.current_step[axis])
                    t_other = (t_now + 1) % T
                    self.viewer.dims.set_current_step(axis, t_other)
                    self.viewer.dims.set_current_step(axis, t_now)
                except Exception:
                    pass
                # Explicit cache-clearing calls (best-effort, API varies).
                for method in ("_reset_cache", "_reset_labels_cache",
                               "set_view_slice", "_set_view_slice"):
                    fn = getattr(layer, method, None)
                    if callable(fn):
                        try:
                            fn()
                        except Exception:
                            pass
                        break
                layer.refresh()
            except Exception:
                traceback.print_exc()
        self._refresh_buttons()

    def _on_discard(self) -> None:
        self._clear_preview()
        if not self.buffer.is_dirty():
            self._log_msg("Nothing to discard.")
            return
        ret = QMessageBox.question(
            self,
            "Discard all edits?",
            "All in-memory edits will be lost.  Continue?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No,
        )
        if ret != QMessageBox.Yes:
            return
        self.buffer.discard()
        # buffer.discard() calls _emit() which triggers _refresh_layer_data()
        # and resets layer.data to buffer._darr — no extra refresh needed here.
        self._log_msg("Discarded all in-memory edits.")
        self._refresh_buttons()

    def _on_stop_clicked(self) -> None:
        if self.request_exit():
            # The Visualization tab owns the lifecycle, so hide ourselves
            # only after it has unhooked us.  Emit a custom Qt signal-like
            # callback if attached.
            if hasattr(self, "stop_callback") and callable(self.stop_callback):
                self.stop_callback()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def request_exit(self) -> bool:
        """Returns True if it is safe to tear the editor down.

        If there are unsaved edits, prompts Save / Discard / Cancel and
        acts on the user's choice.
        """
        if self._edit_thread is not None and self._edit_thread.isRunning():
            QMessageBox.information(
                self,
                "Operation in progress",
                "The current editing operation has not finished yet.\n\n"
                "Please wait for it to complete, then close.",
            )
            return False
        if not self.buffer.is_dirty():
            self._cleanup()
            return True
        box = QMessageBox(self)
        box.setWindowTitle("Unsaved tracked-segment edits")
        box.setText(
            f"You have {self.buffer.dirty_count()} unsaved frame edit(s) for "
            f"{self.cell_type}.  Save them or discard?"
        )
        b_save = box.addButton("Save and close", QMessageBox.AcceptRole)
        b_disc = box.addButton("Discard and close", QMessageBox.DestructiveRole)
        b_canc = box.addButton("Cancel", QMessageBox.RejectRole)
        box.setDefaultButton(b_canc)
        box.exec_()
        clicked = box.clickedButton()
        if clicked is b_canc:
            return False
        if clicked is b_save:
            try:
                self.buffer.save()
            except Exception as exc:
                traceback.print_exc()
                QMessageBox.critical(self, "Save failed", str(exc))
                return False
        else:
            self.buffer.discard()
        self._cleanup()
        return True

    def _cleanup(self) -> None:
        # Stop any in-progress background edit operation.
        if self._edit_thread is not None:
            try:
                self._edit_thread.quit()
                # Allow up to 60 s for long operations (e.g. 500-frame split).
                # If the thread is still alive after that, terminate it forcefully
                # to avoid "QThread: Destroyed while thread is still running" and
                # the resulting fatal crash (Windows 0xC0000409).
                if not self._edit_thread.wait(60_000):
                    self._edit_thread.terminate()
                    self._edit_thread.wait(3_000)
            except Exception:
                pass
            self._edit_thread = None
            self._edit_worker = None
        # Stop any in-progress background materialisation (legacy path).
        if self._materialise_thread is not None:
            try:
                self._materialise_thread.quit()
                self._materialise_thread.wait(2000)  # ms
            except Exception:
                pass
            self._materialise_thread = None
            self._materialise_worker = None
        # Remove the seeds layer if present.
        try:
            if _SPLIT_SEEDS_LAYER in self.viewer.layers:
                self.viewer.layers.remove(_SPLIT_SEEDS_LAYER)
        except Exception:
            pass
        # Remove the selection-highlight overlay if present.
        try:
            if _HIGHLIGHT_LAYER in self.viewer.layers:
                self.viewer.layers.remove(_HIGHLIGHT_LAYER)
        except Exception:
            pass
        # Remove our mouse callback from the bound layer (best-effort).
        try:
            layer = self.viewer.layers[self.layer_name]
            if self._on_layer_clicked in layer.mouse_drag_callbacks:
                layer.mouse_drag_callbacks.remove(self._on_layer_clicked)
        except Exception:
            pass

    # ------------------------------------------------------------------
    # Misc
    # ------------------------------------------------------------------
    def _on_frames_changed(self, frames: List[int]) -> None:
        self._refresh_layer_data(frames)
        self._refresh_highlight_layer()

    def _refresh_buttons(self) -> None:
        has_sel = bool(self._selected_labels)
        for btn in (
            getattr(self, "btn_split_preview", None),
            getattr(self, "btn_split_apply", None),
            getattr(self, "btn_erode_preview", None),
            getattr(self, "btn_erode_apply", None),
            getattr(self, "btn_dilate_preview", None),
            getattr(self, "btn_dilate_apply", None),
            getattr(self, "btn_delete_preview", None),
            getattr(self, "btn_delete_apply", None),
        ):
            if btn is not None:
                btn.setEnabled(has_sel)
        needs_two = len(self._selected_labels) >= 2
        if hasattr(self, "btn_merge_preview"):
            self.btn_merge_preview.setEnabled(needs_two)
        if hasattr(self, "btn_merge_apply"):
            self.btn_merge_apply.setEnabled(needs_two)
        # Create does not need a label selection — always enabled.
        for btn in (
            getattr(self, "btn_create_preview", None),
            getattr(self, "btn_create_apply", None),
        ):
            if btn is not None:
                btn.setEnabled(True)
        if hasattr(self, "btn_undo"):
            self.btn_undo.setEnabled(self.buffer.can_undo())
        if hasattr(self, "btn_redo"):
            self.btn_redo.setEnabled(self.buffer.can_redo())
        n = self.buffer.dirty_count()
        if hasattr(self, "lbl_history"):
            self.lbl_history.setText(
                f"{n} pending frame change(s); {len(self.buffer.history())} step(s) in history"
            )

    def _log_msg(self, msg: str) -> None:
        self.log.append(msg)
        sb = self.log.verticalScrollBar()
        sb.setValue(sb.maximum())
        print(msg, flush=True)

    def _row_for_sample(self, sample_name: str) -> Optional[pd.Series]:
        md = getattr(self.metadata_loader, "metadata", None)
        if md is None or md.empty:
            return None
        rows = md[md["sample_name"].astype(str) == str(sample_name)]
        if rows.empty:
            return None
        return rows.iloc[0]