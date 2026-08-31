"""
BEHAV3D napari plugin – Single Cell Analysis Tab.

Provides two inner sub-tabs:

- **🔬 State Classification** — HMM-based behavioral state clustering,
  intrinsic + full cluster renaming, reports (Step 3), and inline
  backprojection (Step 4).
- **🛤️ Track Classification** — DTW-based trajectory clustering, rename,
  RF train/apply (Step 3), create plots (Step 4), and inline
  backprojection (Step 5).

Both sub-tabs show a warning banner + disable their steps when the
prerequisite file for that cell type is not found on disk. The "Apply"
sections at the top remain enabled so existing outputs can always be loaded.
"""
from __future__ import annotations

import datetime
import traceback
from pathlib import Path
from typing import Optional

import numpy as np
import yaml
from qtpy.QtCore import Qt, QTimer
from qtpy.QtGui import QDesktopServices
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDialog,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMenu,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QStackedWidget,
    QTabWidget,
    QTextEdit,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from behav3d.napari._analysis import CollapsibleSection, DualListGroupSelector
from behav3d.napari._background_runner import (
    BackgroundOperation,
    ProgressBarRow,
    ThreadSafeLogger,
)
from behav3d.napari._pdf_view import open_pdf_in_napari
from behav3d.napari._rename_dialog import RenameClusterDialog
from behav3d.napari._preview_dims import (
    clear_viewer_layers,
    disconnect_all_preview_dims_listeners,
    register_preview_dims_listener,
    stop_dim_playback,
    unregister_preview_dims_listener,
)
from behav3d.core.qt_help import HelpButton, make_help_row, reset_scroll_on_page_change
from behav3d.core.utils import rmtree_ignore_missing, format_timepoints_as_time
from behav3d.core.metadata import resolve_metadata_csv_path


# ═══════════════════════════════════════════════════════════════════════════
# Shared helpers
# ═══════════════════════════════════════════════════════════════════════════

def _detect_sc_cell_types(metadata_loader) -> list[str]:
    """Return immune + other cell types (non-multicolor), plus any
    ``cell_type_groups`` (yml-based, post-filtering) groups categorized as
    immune/other."""
    if metadata_loader is None or getattr(metadata_loader, "metadata", None) is None:
        return []
    try:
        from behav3d.core.metadata import (
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            detect_merged_cell_types_from_metadata,
            filter_multicolor_inputs,
        )
        from behav3d.widgets.utils import detect_cell_type_category
        from behav3d.analysis.grouping import list_cell_type_groups, group_category

        md = metadata_loader.metadata
        imm = list(filter_multicolor_inputs(detect_immune_cell_types_from_metadata(md)))
        oth = list(filter_multicolor_inputs(detect_other_cell_types_from_metadata(md)))

        for ct in detect_merged_cell_types_from_metadata(md):
            try:
                cat = detect_cell_type_category(ct, md)
            except Exception:
                cat = "immune"
            if cat == "immune" and ct not in imm:
                imm.append(ct)
            elif cat == "other" and ct not in oth:
                oth.append(ct)

        # Every group is offered here regardless of its best-effort display
        # category — groups exist specifically for Death Dynamics/Single
        # Cell use (see behav3d/analysis/grouping.py), so an "organoid"
        # category (picked from the group's first member, for icon/bucket
        # purposes only) must never cause a group to be silently dropped.
        params = getattr(metadata_loader, "behav3d_parameters", {}) or {}
        for group_id in list_cell_type_groups(params):
            cat = group_category(params, md, group_id)
            if cat == "other":
                if group_id not in oth:
                    oth.append(group_id)
            elif group_id not in imm:
                imm.append(group_id)

        return imm + oth
    except Exception:
        traceback.print_exc()
        return []


def _group_tracked_status_for(metadata_loader, cell_type):
    """Return ``{"built", "total"}`` tracked-segment status if ``cell_type``
    is a group, else ``None`` (not a group, or status can't be determined).

    Training (HMM/DTW) only needs the group's merged track-features CSV
    (already produced at group-creation time), so this is only relevant for
    gating backprojection/exemplar features, which need per-sample tracked
    zarrs built via ``create_group_tracked_segments``.
    """
    if metadata_loader is None or not cell_type:
        return None
    try:
        from behav3d.analysis.grouping import (
            group_tracked_segments_status,
            list_cell_type_groups,
        )

        params = getattr(metadata_loader, "behav3d_parameters", {}) or {}
        if cell_type not in list_cell_type_groups(params):
            return None
        output_dir = getattr(metadata_loader, "output_dir", None)
        metadata = getattr(metadata_loader, "metadata", None)
        if not output_dir or metadata is None:
            return None
        return group_tracked_segments_status(output_dir, cell_type, metadata)
    except Exception:
        return None


def _apply_group_tracked_gate(warning_label, grp_bp, metadata_loader, cell_type):
    """Disable ``grp_bp`` (Backprojection) and warn if ``cell_type`` is a
    group with no tracked segments built yet.

    Only touches the backprojection group box — training group boxes are
    left alone since training doesn't need tracked segments.
    """
    status = _group_tracked_status_for(metadata_loader, cell_type)
    if status is None or status["built"] > 0:
        return
    grp_bp.setEnabled(False)
    extra = (
        f"⚠ Tracked segments not built for group '{cell_type}' — "
        "Backprojection/Exemplar export unavailable.\n"
        "Build them from the Group Builder dialog (Analysis tab) → "
        "'Build Tracked Segments'."
    )
    # Only stack onto an already-visible warning set earlier in this same
    # prerequisite check; a hidden label may hold stale text from a
    # previous cell-type selection.
    if warning_label.isVisible() and warning_label.text():
        warning_label.setText(f"{warning_label.text()}\n\n{extra}")
    else:
        warning_label.setText(extra)
    warning_label.show()


def _log_backprojection_manifest(log_fn, manifest):
    """Log a backprojection export manifest's success/skip/error counts.

    ``manifest`` is the dict returned by
    ``export_behavioral_state_backprojection_zarrs`` (and, transitively,
    ``export_track_cluster_backprojection``): keys ``output_paths``,
    ``skipped_samples``, ``errors``. Per-sample skip/error reasons
    previously only reached ``warnings.warn`` (invisible in the GUI) — this
    surfaces them in the log instead, which matters most for a group whose
    tracked segments haven't been built for every sample yet.
    """
    if not isinstance(manifest, dict):
        log_fn("✅ Backprojection export done.")
        return
    n_ok = len(manifest.get("output_paths", {}) or {})
    skipped = manifest.get("skipped_samples", []) or []
    errors = manifest.get("errors", []) or []
    log_fn(f"✅ Backprojection export done: {n_ok} sample(s) written.")
    if skipped:
        preview = "; ".join(
            f"{s.get('sample_name', '?')} ({s.get('reason', 'unknown')})"
            for s in skipped[:5]
        )
        more = f" (+{len(skipped) - 5} more)" if len(skipped) > 5 else ""
        log_fn(f"⚠ {len(skipped)} sample(s) skipped: {preview}{more}")
    if errors:
        preview = "; ".join(
            f"{e.get('sample_name', '?')} ({e.get('error', 'unknown')})"
            for e in errors[:5]
        )
        more = f" (+{len(errors) - 5} more)" if len(errors) > 5 else ""
        log_fn(f"❌ {len(errors)} sample(s) errored: {preview}{more}")


def _style_primary(btn: QPushButton):
    btn.setStyleSheet(
        "QPushButton { background: #28a745; color: white; font-weight: bold; "
        "border-radius: 4px; padding: 5px 10px; }"
        "QPushButton:hover { background: #218838; }"
        "QPushButton:disabled { background: #3a3a3a; color: #666; }"
    )


def _style_secondary(btn: QPushButton):
    btn.setStyleSheet(
        "QPushButton { background: #2d5a8e; color: #ddd; font-weight: bold; "
        "border-radius: 4px; padding: 5px 10px; }"
        "QPushButton:hover { background: #3a72b3; }"
        "QPushButton:disabled { background: #3a3a3a; color: #666; }"
    )


def _style_rename(btn: QPushButton):
    btn.setStyleSheet(
        "QPushButton { background: #5a3e7a; color: #ddd; font-weight: bold; "
        "border-radius: 4px; padding: 5px 10px; }"
        "QPushButton:hover { background: #7a5ea0; }"
        "QPushButton:disabled { background: #3a3a3a; color: #666; }"
    )


def _make_queue_btn() -> QPushButton:
    btn = QPushButton("+🛒")
    btn.setFixedSize(36, 28)
    btn.setToolTip("Add to Processing Queue")
    btn.setStyleSheet(
        "QPushButton { background: #1a1a2e; color: #ffc107; "
        "border: 1px solid #ffc107; border-radius: 4px; "
        "font-size: 11px; font-weight: bold; }"
        "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
    )
    return btn


def _make_view_btn() -> QPushButton:
    btn = QPushButton("👁")
    btn.setFixedSize(36, 28)
    btn.setToolTip("View result (not yet available)")
    btn.setEnabled(False)
    btn.setStyleSheet(
        "QPushButton { background: #1a2e1a; color: #6dd56d; "
        "border: 1px solid #6dd56d; border-radius: 4px; "
        "font-size: 12px; font-weight: bold; }"
        "QPushButton:hover { background: #6dd56d; color: #1a2e1a; }"
        "QPushButton:disabled { background: #1a1a1a; color: #555; "
        "border: 1px solid #333; }"
    )
    return btn


def _wire_view_btn(btn: QPushButton, handler, kind: str) -> None:
    """Connect a 👁 button to *handler*, passing the button along.

    The handler needs the button so a multi-result chooser menu can be popped
    up underneath it; a bare ``QMenu.exec_()`` is equivalent to ``exec_(pos())``
    and a never-shown menu's ``pos()`` is (0, 0), i.e. the screen corner.
    """
    btn.clicked.connect(lambda: handler(kind, btn))


def _make_timepoint_time_label() -> QLabel:
    """Return a small muted QLabel for showing a timepoint-window's real-time equivalent."""
    lbl = QLabel("")
    lbl.setStyleSheet("color: #999; font-style: italic; font-size: 10px;")
    return lbl


def _make_info_label(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    lbl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
    lbl.setAlignment(Qt.AlignLeft | Qt.AlignTop)
    return lbl


def _fit_list_widget_height(list_widget, min_rows: int = 2, max_rows: int = 5) -> None:
    """Size a QListWidget's height to its current item count instead of a fixed box,
    so a short candidate-column list doesn't leave a lot of empty space."""
    n = max(min_rows, min(list_widget.count(), max_rows))
    row_h = list_widget.sizeHintForRow(0) if list_widget.count() else 18
    list_widget.setFixedHeight(n * row_h + 2 * list_widget.frameWidth() + 4)


def _condition_levels_for_column(col, *, metadata_loader=None, adata=None, h5ad_path=None) -> list[str]:
    """Cheap, synchronous discovery of the distinct values of a metadata/obs column, for
    populating a "group conditions" level picker without triggering a full h5ad load.

    Tries, in order: the already-loaded metadata.csv DataFrame (covers exp_nr, well, and
    every *_line_condition column), an already-loaded AnnData's .obs (covers columns like
    origin_cell_type once a model is loaded in memory), then a targeted h5py read of just
    that one obs column on disk (same trick as `_read_h5ad_obs_column`/
    `_max_available_track_length` use for the fast synchronous phase of `_reload()`).
    """
    if not col:
        return []
    md = getattr(metadata_loader, "metadata", None) if metadata_loader else None
    if md is not None and col in md.columns:
        return sorted(md[col].dropna().astype(str).unique().tolist())
    if adata is not None and col in adata.obs.columns:
        return sorted(adata.obs[col].dropna().astype(str).unique().tolist())
    if h5ad_path is not None:
        try:
            import h5py
            with h5py.File(str(h5ad_path), "r") as f:
                if "obs" in f and col in f["obs"]:
                    values = TrackClassificationSubTab._read_h5ad_obs_column(f, col)
                    return sorted({str(v) for v in values if v is not None})
        except Exception:
            pass
    return []


def _add_group_axis_row(layout: QVBoxLayout, label_text: str, combo: QComboBox):
    """Adds a "Label: [combo]" row to `layout`, followed by a "Group conditions" checkbox
    and an initially-hidden DualListGroupSelector beneath it - for merging that axis's
    metadata-column levels into two custom groups. Returns (checkbox, selector); the
    caller is responsible for connecting the checkbox's `toggled` signal and populating
    the selector's levels."""
    row = QHBoxLayout()
    row.addWidget(QLabel(label_text))
    row.addWidget(combo, stretch=1)
    layout.addLayout(row)
    checkbox = QCheckBox("Group conditions")
    checkbox.setToolTip(
        "Pool this axis's levels into two custom groups (left/right), instead of one "
        "panel per level."
    )
    layout.addWidget(checkbox)
    selector = DualListGroupSelector()
    selector.setVisible(False)
    layout.addWidget(selector)
    return checkbox, selector


def _make_browse_row(
    le: QLineEdit,
    btn: QPushButton,
    help_title: str = "",
    help_desc: str = "",
) -> QHBoxLayout:
    """Return HBoxLayout: [line_edit | browse_btn | ? help]."""
    row = QHBoxLayout()
    row.setContentsMargins(0, 0, 0, 0)
    row.addWidget(le, stretch=1)
    row.addWidget(btn)
    if help_title:
        row.addWidget(HelpButton(help_title, help_desc))
    return row


def _make_chk_help_row(chk: QCheckBox, title: str, desc: str) -> QHBoxLayout:
    """Return HBoxLayout: [checkbox | ? help]."""
    row = QHBoxLayout()
    row.setContentsMargins(0, 0, 0, 0)
    row.addWidget(chk)
    row.addWidget(HelpButton(title, desc))
    row.addStretch()
    return row


# ── Default selections (used when no prior config exists) ───────────────────
_DEFAULT_TIMEPOINT_FEATURES = {"speed", "sphericity", "elongation", "extent", "solidity"}
_DEFAULT_WINDOW_FEATURES    = {"net_displacement"}
_DEFAULT_LOG_SCALE_FEATURES = {"speed"}

_TECHNICAL_OBS_COLS = {
    "position_t", "sample_name", "TrackID",
    "intrinsic_behavioral_cluster", "full_behavioral_cluster",
    "hmm_intrinsic_behavioral_state", "behavioral_state",
}

_CHANNEL_COLORS = ["cyan", "yellow", "green", "red", "blue", "magenta"]


def _bp_add_raw_channels(viewer, raw_img, sample_name, saved_channels, dim_order="TCZYX"):
    """Add raw_img to viewer as one Image layer per channel. Returns list of layer names."""
    import dask.array as da
    if not isinstance(raw_img, da.Array):
        raw_img = da.from_array(raw_img)
    if dim_order != "TCZYX" and len(dim_order) == 5:
        try:
            axes = [dim_order.index(d) for d in "TCZYX"]
            raw_img = da.transpose(raw_img, axes)
        except ValueError:
            pass
    n_channels = raw_img.shape[1] if raw_img.ndim >= 5 else 1
    added = []
    for c in range(n_channels):
        ch_data = raw_img[:, c, :, :, :] if raw_img.ndim >= 5 else raw_img
        saved = saved_channels.get(c) or saved_channels.get(str(c))
        color = saved["colormap"] if saved and "colormap" in saved else _CHANNEL_COLORS[c % len(_CHANNEL_COLORS)]
        layer_name = f"{sample_name} – Ch{c}"
        add_kwargs = dict(name=layer_name, colormap=color, blending="additive", visible=True)
        if saved and "contrast_limits" in saved:
            add_kwargs["contrast_limits"] = tuple(saved["contrast_limits"])
        viewer.add_image(ch_data, **add_kwargs)
        added.append(layer_name)
    return added


def _bp_save_channel_display(viewer, metadata_loader, out_dir_fn):
    """Snapshot current channel display settings and persist to YAML."""
    params = getattr(metadata_loader, "behav3d_parameters", None)
    if not isinstance(params, dict):
        return
    channels_cfg = params.setdefault("viewer_display", {}).setdefault("channels", {})
    for layer in viewer.layers:
        if " – Ch" not in layer.name:
            continue
        try:
            ch_idx = int(layer.name.rsplit("Ch", 1)[1])
        except (ValueError, IndexError):
            continue
        if hasattr(layer, "contrast_limits") and hasattr(layer, "colormap"):
            channels_cfg[ch_idx] = {
                "colormap": str(layer.colormap.name),
                "contrast_limits": [float(layer.contrast_limits[0]), float(layer.contrast_limits[1])],
            }
    _save_behav3d_params(metadata_loader, out_dir_fn)


def _save_behav3d_params(metadata_loader, out_dir_fn):
    params = getattr(metadata_loader, "behav3d_parameters", None)
    if not isinstance(params, dict):
        return
    cfg_path = getattr(metadata_loader, "behav3d_parameters_path", None)
    if cfg_path is None:
        out = out_dir_fn()
        if out:
            cfg_path = out / "behav3d_parameters.yml"
    if cfg_path is not None:
        try:
            with open(cfg_path, "w", encoding="utf-8") as f:
                yaml.safe_dump(params, f, sort_keys=False)
        except Exception:
            pass


# ═══════════════════════════════════════════════════════════════════════════
# State Classification sub-tab
# ═══════════════════════════════════════════════════════════════════════════

class StateClassificationSubTab(QWidget):
    """Steps 1-4 for HMM-based behavioral state classification + backprojection."""

    def __init__(self, viewer=None, metadata_loader=None, cell_type_getter=None,
                 parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._get_cell_type = cell_type_getter
        self._model_adata = None
        self._hmm_model = None
        self._hmm_model_cell_type = None
        self._bg = BackgroundOperation(self)
        self._preload_bg = BackgroundOperation(self)
        self._last_features_key: tuple = ()
        # Current-timepoint-only state backprojection preview (see
        # `_refresh_state_bp_layer`): recomputed on every dims scrub instead
        # of writing a full per-sample zarr up front.
        self._state_bp_preview: Optional[dict] = None
        self._state_bp_dims_callback = None
        # Populated by _populate_dynamic_features; reused by the log-scaling
        # "Preview feature distributions" histogram button.
        self._logscale_candidate_cols: list = []
        self._logscale_csv_path: Optional[Path] = None

        self._display_save_timer = QTimer(self)
        self._display_save_timer.setSingleShot(True)
        self._display_save_timer.setInterval(1000)
        self._display_save_timer.timeout.connect(self._persist_bp_viewer_display)

        # Queue buttons (wired by _widget.py)
        self.btn_queue_state_cluster = _make_queue_btn()
        self.btn_queue_train_state = _make_queue_btn()
        self.btn_queue_apply_state = _make_queue_btn()

        self._init_ui()

    # ── UI ──────────────────────────────────────────────────────────────

    def _init_ui(self):
        from behav3d.napari._guided import make_back_header

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # Pinned above the scroll area (not inside it) so it stays visible
        # regardless of how far the plotting page's content scrolls. Only
        # shown while the plotting page (outer page 1) is active.
        self._plotting_header, self._pipeline_title = make_back_header(
            on_back=self._plotting_back, label="‹ Back"
        )
        self._plotting_header.hide()
        outer.addWidget(self._plotting_header)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        content = QWidget()
        content_lay = QVBoxLayout(content)
        content_lay.setContentsMargins(0, 0, 0, 0)
        content_lay.setSpacing(0)
        scroll.setWidget(content)

        # Outer per-subtab stack: page 0 = clustering settings (Steps 1/2/4 +
        # the "Generate analysis and plots" trigger), page 1 = the fully
        # isolated plotting page, reached only via that trigger.
        self._subtab_stack = QStackedWidget()
        content_lay.addWidget(self._subtab_stack)

        settings_page = QWidget()
        lay = QVBoxLayout(settings_page)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        self._subtab_stack.addWidget(settings_page)  # outer page 0

        # ── Warning label (prerequisite check) ──────────────────────────
        self.warning_label = QLabel("")
        self.warning_label.setWordWrap(True)
        self.warning_label.setStyleSheet(
            "QLabel { background: #3d2200; color: #ffaa44; border-radius: 4px; "
            "padding: 6px 8px; font-size: 11px; }"
        )
        self.warning_label.hide()
        lay.addWidget(self.warning_label)

        # ── Apply Existing Checkbox ──────────────────────────────────────
        self.chk_apply_existing = QCheckBox("Apply existing behavioral state classification")
        self.chk_apply_existing.setChecked(False)
        lay.addWidget(self.chk_apply_existing)

        # ── Group: Apply Existing ────────────────────────────────────────
        self.grp_apply = QGroupBox("Apply saved HMM artifact")
        apply_lay = QVBoxLayout(self.grp_apply)
        apply_lay.setSpacing(4)

        apply_lay.addWidget(_make_info_label("Load an existing HMM deployment artifact (.pkl) to apply "
                                  "the pre-trained classification to the full dataset."))

        self.le_hmm_artifact = QLineEdit()
        self.le_hmm_artifact.setPlaceholderText("Path to saved HMM deployment artifact .pkl")
        self.btn_browse_hmm = QPushButton("Browse…")
        apply_lay.addLayout(_make_browse_row(
            self.le_hmm_artifact, self.btn_browse_hmm,
            "HMM artifact path",
            "Path to the .pkl HMM deployment artifact generated by a previous "
            "run of State Classification."
        ))

        self.btn_apply_hmm = QPushButton("▶ Apply saved HMM artifact")
        _style_primary(self.btn_apply_hmm)
        apply_lay.addWidget(self.btn_apply_hmm)
        self.grp_apply.hide()
        lay.addWidget(self.grp_apply)

        # ── Group: Training (Step 1) ─────────────────────────────────────
        self.grp_train = QGroupBox("Step 1 — State Clustering")
        train_lay = QVBoxLayout(self.grp_train)
        train_lay.setSpacing(4)

        # Feature Selection
        feat_sel_sec = CollapsibleSection("Feature Selection", expanded=True)
        self.feat_sel_scroll = QScrollArea()
        self.feat_sel_scroll.setWidgetResizable(True)
        self.feat_sel_scroll.setMinimumHeight(150)
        feat_sel_content = QWidget()
        self.feat_sel_lay = QVBoxLayout(feat_sel_content)
        self.feat_sel_lay.setSpacing(2)

        self.feat_sel_lay.addWidget(QLabel("<b>Window features</b>"))
        win_feat_form = QFormLayout()
        self.spin_window_size = QSpinBox()
        self.spin_window_size.setRange(1, 500)
        self.spin_window_size.setValue(5)
        win_feat_form.addRow("Window size (timepoints):", make_help_row(
            self.spin_window_size, "Window size",
            "Size of the rolling window (in timepoints/frames) used to compute window "
            "features (net_displacement, straightness, mean_square_displacement). The line "
            "below shows the equivalent real time, computed from this sample's "
            "time_interval/time_unit metadata."
        ))
        self.lbl_window_size_time = _make_timepoint_time_label()
        win_feat_form.addRow("", self.lbl_window_size_time)
        self.spin_window_size.valueChanged.connect(self._update_timepoint_time_labels)
        self.feat_sel_lay.addLayout(win_feat_form)

        self.chk_net_disp = QCheckBox("net_displacement")
        self.chk_straight = QCheckBox("straightness")
        self.chk_msd = QCheckBox("mean_square_displacement")
        self.feat_sel_lay.addWidget(self.chk_net_disp)
        self.feat_sel_lay.addWidget(self.chk_straight)
        self.feat_sel_lay.addWidget(self.chk_msd)

        self.feat_sel_lay.addWidget(QLabel("<b>Timepoint features</b>"))
        self.timepoint_features_lay = QVBoxLayout()
        self.feat_sel_lay.addLayout(self.timepoint_features_lay)

        self.feat_sel_scroll.setWidget(feat_sel_content)
        feat_sel_sec.addWidget(self.feat_sel_scroll)
        train_lay.addWidget(feat_sel_sec)

        # Feature Processing
        feat_proc_sec = CollapsibleSection("Feature Processing", expanded=False)
        self.feat_proc_scroll = QScrollArea()
        self.feat_proc_scroll.setWidgetResizable(True)
        self.feat_proc_scroll.setMinimumHeight(150)
        feat_proc_content = QWidget()
        self.feat_proc_lay = QVBoxLayout(feat_proc_content)

        # Log Scaling inside Feature Processing
        log_scale_sec = CollapsibleSection("Log scaling", expanded=False)
        # Distribution preview: sits just below the "Log scaling" header and
        # above the per-feature checkboxes, so the user can inspect which
        # features are skewed (and would benefit from log scaling) first.
        self.btn_preview_distributions = QPushButton("\U0001F4CA  Preview feature distributions")
        self.btn_preview_distributions.setToolTip(
            "Show histograms of the continuous features for the selected cell "
            "type, to judge which ones benefit from log scaling before applying it."
        )
        self.btn_preview_distributions.clicked.connect(self._on_preview_feature_distributions)
        log_scale_sec.addWidget(self.btn_preview_distributions)
        self.log_scale_lay = QVBoxLayout()
        log_scale_container = QWidget()
        log_scale_container.setLayout(self.log_scale_lay)
        log_scale_sec.addWidget(log_scale_container)
        self.feat_proc_lay.addWidget(log_scale_sec)

        proc_form = QFormLayout()
        self.spin_hmm_feature_smoothing_window = QSpinBox()
        self.spin_hmm_feature_smoothing_window.setRange(1, 100)
        self.spin_hmm_feature_smoothing_window.setValue(5)
        proc_form.addRow("Smooth window (timepoints):", make_help_row(
            self.spin_hmm_feature_smoothing_window, "Feature smoothing window",
            "Size of the rolling smoothing window (in timepoints/frames) applied to features "
            "before HMM fitting. 1 = no smoothing. The line below shows the equivalent real "
            "time, computed from this sample's time_interval/time_unit metadata."
        ))
        self.lbl_hmm_smoothing_time = _make_timepoint_time_label()
        proc_form.addRow("", self.lbl_hmm_smoothing_time)
        self.spin_hmm_feature_smoothing_window.valueChanged.connect(self._update_timepoint_time_labels)

        self.spin_quant_lo = QDoubleSpinBox()
        self.spin_quant_lo.setRange(0.0, 0.5)
        self.spin_quant_lo.setSingleStep(0.01)
        self.spin_quant_lo.setValue(0.0)
        proc_form.addRow("Low percentile cap:", make_help_row(
            self.spin_quant_lo, "Low percentile cap",
            "Quantile at which to clip the lower end of feature values (0 = no clipping)."
        ))

        self.spin_quant_hi = QDoubleSpinBox()
        self.spin_quant_hi.setRange(0.5, 1.0)
        self.spin_quant_hi.setSingleStep(0.01)
        self.spin_quant_hi.setValue(1.0)
        proc_form.addRow("High percentile cap:", make_help_row(
            self.spin_quant_hi, "High percentile cap",
            "Quantile at which to clip the upper end of feature values (0.99 = clip top 1%)."
        ))
        self.feat_proc_lay.addLayout(proc_form)

        self.feat_proc_scroll.setWidget(feat_proc_content)
        feat_proc_sec.addWidget(self.feat_proc_scroll)
        train_lay.addWidget(feat_proc_sec)

        # Binary Group Selection
        bin_grp_sec = CollapsibleSection("Binary Group Selection", expanded=False)
        self.bin_grp_scroll = QScrollArea()
        self.bin_grp_scroll.setWidgetResizable(True)
        self.bin_grp_scroll.setMinimumHeight(100)
        bin_grp_content = QWidget()
        self.bin_grp_lay = QVBoxLayout(bin_grp_content)
        self.bin_grp_scroll.setWidget(bin_grp_content)
        bin_grp_sec.addWidget(self.bin_grp_scroll)
        train_lay.addWidget(bin_grp_sec)

        # Advanced Settings
        adv_sec = CollapsibleSection("⚙ Advanced Configuration", expanded=False)
        adv_form = QFormLayout()

        self.combo_hmm_n_states_mode = QComboBox()
        self.combo_hmm_n_states_mode.addItems(["fixed", "auto"])
        adv_form.addRow("State selection mode:", make_help_row(
            self.combo_hmm_n_states_mode, "State selection mode",
            "'fixed' uses n_states directly. 'auto' fits a separate HMM for every k "
            "between k_min and k_max and selects the model with the lowest BIC "
            "(Bayesian Information Criterion)."
        ))

        self.spin_hmm_k_min = QSpinBox()
        self.spin_hmm_k_min.setRange(2, 50)
        self.spin_hmm_k_min.setValue(2)
        self.row_k_min = QLabel("k_min (Auto mode):")
        k_min_row = make_help_row(
            self.spin_hmm_k_min, "k_min",
            "Minimum number of states to test in auto mode."
        )
        self.help_k_min = k_min_row.itemAt(1).widget()
        adv_form.addRow(self.row_k_min, k_min_row)

        self.spin_hmm_k_max = QSpinBox()
        self.spin_hmm_k_max.setRange(2, 50)
        self.spin_hmm_k_max.setValue(8)
        self.row_k_max = QLabel("k_max (Auto mode):")
        k_max_row = make_help_row(
            self.spin_hmm_k_max, "k_max",
            "Maximum number of states to test in auto mode."
        )
        self.help_k_max = k_max_row.itemAt(1).widget()
        adv_form.addRow(self.row_k_max, k_max_row)

        self.spin_hmm_start_offset = QSpinBox()
        self.spin_hmm_start_offset.setRange(0, 100000)
        self.spin_hmm_start_offset.setValue(1)
        adv_form.addRow("Start offset:", make_help_row(
            self.spin_hmm_start_offset, "Start offset",
            "Number of initial timepoints to skip per track (e.g. dead frames). "
            "0 = use all timepoints."
        ))

        self.combo_hmm_start_offset_fill_mode = QComboBox()
        self.combo_hmm_start_offset_fill_mode.addItems(["backfill", "leave_unassigned"])
        adv_form.addRow("Skipped frames:", make_help_row(
            self.combo_hmm_start_offset_fill_mode, "Skipped frames fill mode",
            "'backfill' assigns the first predicted state to skipped frames. "
            "'leave_unassigned' leaves them as NaN."
        ))

        self.combo_hmm_covariance_type = QComboBox()
        self.combo_hmm_covariance_type.addItems(["full", "diag", "spherical", "tied"])
        adv_form.addRow("Covariance type:", make_help_row(
            self.combo_hmm_covariance_type, "Covariance type",
            "Covariance matrix structure for the HMM Gaussian emissions. "
            "'full' is most flexible; 'diag' is faster and avoids overfitting on small data."
        ))

        self.spin_hmm_n_iter = QSpinBox()
        self.spin_hmm_n_iter.setRange(10, 10000)
        self.spin_hmm_n_iter.setValue(200)
        adv_form.addRow("n_iter:", make_help_row(
            self.spin_hmm_n_iter, "Max EM iterations",
            "Maximum number of Expectation-Maximisation iterations for HMM fitting."
        ))

        self.spin_hmm_tol = QDoubleSpinBox()
        self.spin_hmm_tol.setRange(1e-6, 1.0)
        self.spin_hmm_tol.setSingleStep(1e-3)
        self.spin_hmm_tol.setDecimals(6)
        self.spin_hmm_tol.setValue(1e-3)
        adv_form.addRow("tol:", make_help_row(
            self.spin_hmm_tol, "Convergence tolerance",
            "Convergence threshold for the EM algorithm. Smaller = more precise but slower."
        ))

        self.spin_hmm_min_covar = QDoubleSpinBox()
        self.spin_hmm_min_covar.setRange(1e-6, 1.0)
        self.spin_hmm_min_covar.setSingleStep(1e-3)
        self.spin_hmm_min_covar.setDecimals(6)
        self.spin_hmm_min_covar.setValue(1e-3)
        adv_form.addRow("min_covar:", make_help_row(
            self.spin_hmm_min_covar, "Minimum covariance",
            "Floor added to covariance diagonal to prevent numerical instability."
        ))

        self.chk_hmm_sticky = QCheckBox("Sticky HMM")
        adv_form.addRow("", _make_chk_help_row(
            self.chk_hmm_sticky, "Sticky HMM",
            "Adds a self-transition bias (kappa) to make the model prefer staying "
            "in the same state, producing longer, more stable behavioral bouts."
        ))

        self.spin_hmm_stickiness_kappa = QDoubleSpinBox()
        self.spin_hmm_stickiness_kappa.setRange(0.0, 100.0)
        self.spin_hmm_stickiness_kappa.setValue(8.0)
        self.row_kappa = QLabel("kappa (Sticky):")
        kappa_row = make_help_row(
            self.spin_hmm_stickiness_kappa, "Stickiness kappa",
            "Self-transition bias strength. Higher = longer bouts before switching states."
        )
        self.help_kappa = kappa_row.itemAt(1).widget()
        adv_form.addRow(self.row_kappa, kappa_row)

        self.spin_hmm_transmat_alpha = QDoubleSpinBox()
        self.spin_hmm_transmat_alpha.setRange(0.0, 100.0)
        self.spin_hmm_transmat_alpha.setValue(1.0)
        self.row_alpha = QLabel("alpha (Sticky):")
        alpha_row = make_help_row(
            self.spin_hmm_transmat_alpha, "Transition matrix alpha",
            "Dirichlet prior concentration for the transition matrix in sticky HMM mode."
        )
        self.help_alpha = alpha_row.itemAt(1).widget()
        adv_form.addRow(self.row_alpha, alpha_row)

        self.spin_seed = QSpinBox()
        self.spin_seed.setRange(0, 99999)
        self.spin_seed.setValue(123)
        self.spin_seed.setMaximumWidth(90)
        adv_form.addRow("Random seed:", make_help_row(
            self.spin_seed, "Random seed",
            "Seed for reproducibility of HMM initialisation."
        ))

        adv_sec.addLayout(adv_form)
        train_lay.addWidget(adv_sec)

        # n_states
        nstates_form = QFormLayout()
        self.spin_hmm_n_states = QSpinBox()
        self.spin_hmm_n_states.setRange(2, 50)
        self.spin_hmm_n_states.setValue(4)
        nstates_form.addRow("n_states:", make_help_row(
            self.spin_hmm_n_states, "Number of states",
            "Number of HMM behavioral states. Ignored in 'auto' mode (use k_min/k_max). "
            "Typical values: 3–8."
        ))
        train_lay.addLayout(nstates_form)

        # Configuration summary
        self.config_summary_label = QLabel("")
        self.config_summary_label.setWordWrap(True)
        self.config_summary_label.setStyleSheet(
            "QLabel { background: #1a2e2e; color: #66ddcc; border-radius: 4px; "
            "padding: 8px 10px; font-size: 11px; font-family: monospace; }"
        )
        self.config_summary_label.hide()
        train_lay.addWidget(self.config_summary_label)

        # Run button
        row_run = QHBoxLayout()
        self.btn_run_state = QPushButton("▶ Run State Classification")
        _style_primary(self.btn_run_state)
        row_run.addWidget(self.btn_run_state, stretch=1)
        row_run.addWidget(self.btn_queue_state_cluster)
        self.btn_view_state = _make_view_btn()
        row_run.addWidget(self.btn_view_state)
        train_lay.addLayout(row_run)

        lay.addWidget(self.grp_train)

        # ── Step 2: Rename Clusters ──────────────────────────────────────
        self.grp2 = QGroupBox("Step 2 — Rename Clusters")
        g2 = QVBoxLayout(self.grp2)
        g2.setSpacing(4)

        self.rename_status_lbl = QLabel("ℹ Run state classification first to enable renaming.")
        self.rename_status_lbl.setStyleSheet("color: #999; font-size: 11px;")
        g2.addWidget(self.rename_status_lbl)

        self.btn_rename_intrinsic = QPushButton("✏  Rename Primary Dynamic State Clusters")
        _style_rename(self.btn_rename_intrinsic)
        self.btn_rename_intrinsic.setEnabled(False)
        g2.addWidget(self.btn_rename_intrinsic)

        self.btn_rename_full = QPushButton("✏  Rename Full Behavioral Clusters (Binary Groups)")
        _style_rename(self.btn_rename_full)
        self.btn_rename_full.setEnabled(False)
        g2.addWidget(self.btn_rename_full)
        lay.addWidget(self.grp2)

        # ── Step 3: Reports (built directly into per-pipeline group boxes) ──
        self.grp_state_diagnostics = QGroupBox("Diagnostics")
        g_state_diagnostics = QVBoxLayout(self.grp_state_diagnostics)
        g_state_diagnostics.setSpacing(4)
        g_state_diagnostics.addWidget(_make_info_label(
            "Runs quality-control diagnostics on the HMM state clustering "
            "(state means, transitions, feature distributions)."
        ))
        state_diag_row = QHBoxLayout()
        self.btn_state_diagnostics = QPushButton("▶ Create Diagnostics")
        _style_secondary(self.btn_state_diagnostics)
        state_diag_row.addWidget(self.btn_state_diagnostics, stretch=1)
        self.btn_view_state_diagnostics = _make_view_btn()
        state_diag_row.addWidget(self.btn_view_state_diagnostics)
        g_state_diagnostics.addLayout(state_diag_row)

        self.grp_state_composition = QGroupBox("State Composition Report")
        g_state_composition = QVBoxLayout(self.grp_state_composition)
        g_state_composition.setSpacing(4)
        self.combo_composition_group_x = QComboBox()
        self.combo_composition_group_x.setMinimumWidth(160)
        self.combo_composition_group_x.currentTextChanged.connect(
            self._on_composition_group_x_changed
        )
        self.chk_composition_group_x_conditions, self.group_selector_composition_x = _add_group_axis_row(
            g_state_composition, "Group in X:", self.combo_composition_group_x
        )
        self.chk_composition_group_x_conditions.toggled.connect(
            self._on_composition_group_x_conditions_toggled
        )
        self.combo_composition_group_y = QComboBox()
        self.combo_composition_group_y.setMinimumWidth(160)
        self.combo_composition_group_y.currentTextChanged.connect(
            self._on_composition_group_y_changed
        )
        self.chk_composition_group_y_conditions, self.group_selector_composition_y = _add_group_axis_row(
            g_state_composition, "Group in Y:", self.combo_composition_group_y
        )
        self.chk_composition_group_y_conditions.toggled.connect(
            self._on_composition_group_y_conditions_toggled
        )
        g_state_composition.addWidget(QLabel("Group per page (Ctrl/Cmd click for multiple):"))
        self.list_composition_group_cols = QListWidget()
        self.list_composition_group_cols.setSelectionMode(QAbstractItemView.ExtendedSelection)
        g_state_composition.addWidget(self.list_composition_group_cols)
        g_state_composition.addWidget(_make_info_label(
            "Plots the proportion of each behavioral state per sample. Tick \"Group conditions\" "
            "under Group in X/Y to pool that axis's levels into two custom groups."
        ))
        comp_row = QHBoxLayout()
        self.btn_state_composition = QPushButton("▶ State Composition Report")
        _style_secondary(self.btn_state_composition)
        comp_row.addWidget(self.btn_state_composition, stretch=1)
        self.btn_view_composition = _make_view_btn()
        comp_row.addWidget(self.btn_view_composition)
        g_state_composition.addLayout(comp_row)

        self.grp_state_transition = QGroupBox("State Transition Report")
        g_state_transition = QVBoxLayout(self.grp_state_transition)
        g_state_transition.setSpacing(4)
        g_state_transition.addWidget(_make_info_label(
            "Plots state-to-state transition probabilities."
        ))
        trans_row = QHBoxLayout()
        self.btn_state_transition = QPushButton("▶ State Transition Report")
        _style_secondary(self.btn_state_transition)
        trans_row.addWidget(self.btn_state_transition, stretch=1)
        self.btn_view_transition = _make_view_btn()
        trans_row.addWidget(self.btn_view_transition)
        g_state_transition.addLayout(trans_row)

        self.grp_state_comparison = QGroupBox("Condition Comparison Report")
        g_state_comparison = QVBoxLayout(self.grp_state_comparison)
        g_state_comparison.setSpacing(4)
        g_state_comparison.addWidget(_make_info_label(
            "Condition comparison (overall proportions, Welch's t-test): each row is one pairwise "
            "comparison of \"Compare condition\"'s levels; \"Group in X\" splits it into side-by-side "
            "columns from another condition. Tick \"Group conditions\" to instead pool levels into "
            "two custom groups and compare just those."
        ))
        comparison_form = QFormLayout()
        comparison_form.setSpacing(3)
        self.combo_comparison_condition_col = QComboBox()
        self.combo_comparison_condition_col.setMinimumWidth(200)
        self.combo_comparison_condition_col.currentTextChanged.connect(
            self._on_comparison_condition_col_changed
        )
        comparison_form.addRow("Compare condition:", self.combo_comparison_condition_col)
        g_state_comparison.addLayout(comparison_form)

        self.combo_comparison_group_x = QComboBox()
        self.combo_comparison_group_x.setMinimumWidth(160)
        self.combo_comparison_group_x.currentTextChanged.connect(
            self._on_comparison_group_x_changed
        )
        self.chk_comparison_group_x_conditions, self.group_selector_comparison_x = _add_group_axis_row(
            g_state_comparison, "Group in X:", self.combo_comparison_group_x
        )
        self.chk_comparison_group_x_conditions.toggled.connect(
            self._on_comparison_group_x_conditions_toggled
        )

        self.line_comparison_group_y = QLineEdit()
        self.line_comparison_group_y.setMinimumWidth(160)
        self.line_comparison_group_y.setReadOnly(True)
        y_row = QHBoxLayout()
        y_row.addWidget(QLabel("Group in Y:"))
        y_row.addWidget(self.line_comparison_group_y, stretch=1)
        g_state_comparison.addLayout(y_row)
        self.chk_comparison_group_conditions = QCheckBox("Group conditions")
        self.chk_comparison_group_conditions.setToolTip(
            "Pool the compared column's levels into two custom groups (left/right) and run a "
            "single comparison between them, instead of every pairwise level combination."
        )
        self.chk_comparison_group_conditions.toggled.connect(
            self._on_comparison_group_conditions_toggled
        )
        g_state_comparison.addWidget(self.chk_comparison_group_conditions)
        self.group_selector_comparison = DualListGroupSelector()
        self.group_selector_comparison.setVisible(False)
        g_state_comparison.addWidget(self.group_selector_comparison)

        g_state_comparison.addWidget(QLabel("Group per page (Ctrl/Cmd click for multiple):"))
        self.list_comparison_group_cols = QListWidget()
        self.list_comparison_group_cols.setSelectionMode(QAbstractItemView.ExtendedSelection)
        g_state_comparison.addWidget(self.list_comparison_group_cols)
        comparison_row = QHBoxLayout()
        self.btn_condition_comparison = QPushButton("▶ Condition Comparison Report")
        _style_secondary(self.btn_condition_comparison)
        comparison_row.addWidget(self.btn_condition_comparison, stretch=1)
        self.btn_view_condition_comparison = _make_view_btn()
        comparison_row.addWidget(self.btn_view_condition_comparison)
        g_state_comparison.addLayout(comparison_row)

        from behav3d.napari._guided import GuidedPanel
        from behav3d.napari.analysis_guided_copy import (
            STATE_REPORTS_ENTRY, STATE_REPORT_PIPELINES,
        )

        # Entry trigger stays inline on the clustering-settings page; Start
        # navigates to the fully isolated plotting page (outer page 1) below.
        lay.addWidget(
            GuidedPanel([STATE_REPORTS_ENTRY],
                        start_cb=lambda _id: self._open_plotting_page())
        )

        plotting_page = QWidget()
        plotting_lay = QVBoxLayout(plotting_page)
        plotting_lay.setContentsMargins(0, 0, 0, 0)
        plotting_lay.setSpacing(0)

        # Inner stack: page 0 = report list, page 1 = a report's settings.
        self._plots_stack = QStackedWidget()

        pipeline_page = QWidget()
        pipeline_lay = QVBoxLayout(pipeline_page)
        pipeline_lay.setContentsMargins(0, 0, 0, 0)
        pipeline_lay.setSpacing(0)
        pipeline_lay.addWidget(
            GuidedPanel(STATE_REPORT_PIPELINES, start_cb=self._on_pipeline_start)
        )
        self._plots_stack.addWidget(pipeline_page)  # inner page 0

        report_settings_page = QWidget()
        report_settings_lay = QVBoxLayout(report_settings_page)
        report_settings_lay.setContentsMargins(0, 0, 0, 0)
        report_settings_lay.setSpacing(0)
        self._pipeline_scroll = QScrollArea()
        self._pipeline_scroll.setWidgetResizable(True)
        self._pipeline_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        report_settings_lay.addWidget(self._pipeline_scroll)
        self._plots_stack.addWidget(report_settings_page)  # inner page 1

        pipeline_content = QWidget()
        pipeline_content_lay = QVBoxLayout(pipeline_content)
        pipeline_content_lay.setContentsMargins(6, 6, 6, 6)
        pipeline_content_lay.setSpacing(6)
        pipeline_content_lay.addWidget(self.grp_state_diagnostics)
        pipeline_content_lay.addWidget(self.grp_state_composition)
        pipeline_content_lay.addWidget(self.grp_state_transition)
        pipeline_content_lay.addWidget(self.grp_state_comparison)
        pipeline_content_lay.addStretch()
        self._pipeline_scroll.setWidget(pipeline_content)
        reset_scroll_on_page_change(self._plots_stack)

        self._plots_stack.setCurrentIndex(0)
        plotting_lay.addWidget(self._plots_stack)
        self._subtab_stack.addWidget(plotting_page)  # outer page 1
        reset_scroll_on_page_change(self._subtab_stack)

        self._subtab_stack.setCurrentIndex(0)

        # ── Step 4: Backprojection ───────────────────────────────────────
        self.grp_bp = QGroupBox("Step 4 — Backprojection")
        g_bp = QVBoxLayout(self.grp_bp)
        g_bp.setSpacing(6)

        # Live napari view
        grp_view = QGroupBox("Live Napari Layer Backprojection")
        g_view = QVBoxLayout(grp_view)
        g_view.setSpacing(4)

        info_bp = QLabel(
            "Loads the state adata and overlays colored state labels onto "
            "the selected sample image in napari."
        )
        info_bp.setStyleSheet("color: #999; font-size: 10px;")
        info_bp.setWordWrap(True)
        g_view.addWidget(info_bp)

        bp_form = QFormLayout()
        bp_form.setSpacing(3)

        self.sample_combo = QComboBox()
        self.sample_combo.addItem("— All samples —")
        self.sample_combo.setMinimumWidth(120)
        self.sample_combo.setToolTip(
            "Select a sample for backprojection visualisation. "
            "'— All samples —' uses the first available sample."
        )
        bp_form.addRow("Sample:", self.sample_combo)

        self.combo_state_color_by = QComboBox()
        self.combo_state_color_by.addItems(
            ["full_behavioral_cluster", "intrinsic_behavioral_cluster", "raw_hmm_state"]
        )
        self.combo_state_color_by.setMinimumWidth(220)
        bp_form.addRow("Color by:", make_help_row(
            self.combo_state_color_by, "Color by",
            "Which clustering label to use for coloring the backprojection overlay. "
            "'full_behavioral_cluster' includes binary grouping; "
            "'intrinsic_behavioral_cluster' shows only the HMM states; "
            "'raw_hmm_state' shows the HMM output of Step 1 before any renaming."
        ))

        self.spin_state_opacity = QSpinBox()
        self.spin_state_opacity.setRange(10, 100)
        self.spin_state_opacity.setValue(80)
        self.spin_state_opacity.setSuffix(" %")
        self.spin_state_opacity.setMaximumWidth(90)
        bp_form.addRow("Opacity:", make_help_row(
            self.spin_state_opacity, "Opacity",
            "Opacity of the colored state overlay layer in napari (10–100%)."
        ))

        self.chk_show_state_trajectories = QCheckBox("Show trajectories")
        self.chk_show_state_trajectories.setChecked(True)
        bp_form.addRow("Trajectories:", _make_chk_help_row(
            self.chk_show_state_trajectories, "Show trajectories",
            "Overlay each track's full path as a line whose color changes over "
            "time to match its state at each timepoint. On by default — adds "
            "one napari Tracks layer per state."
        ))
        g_view.addLayout(bp_form)

        view_row = QHBoxLayout()
        self.btn_show_state_bp = QPushButton("▶ Show State Backprojection in Napari")
        _style_primary(self.btn_show_state_bp)
        view_row.addWidget(self.btn_show_state_bp, stretch=1)
        g_view.addLayout(view_row)

        export_run_row = QHBoxLayout()
        self.btn_export_state_bp = QPushButton("▶ Export State Backprojection")
        _style_secondary(self.btn_export_state_bp)
        export_run_row.addWidget(self.btn_export_state_bp, stretch=1)
        g_view.addLayout(export_run_row)
        g_bp.addWidget(grp_view)

        lay.addWidget(self.grp_bp)

        # ── Progress + Log ───────────────────────────────────────────────
        self.progress_row = ProgressBarRow()
        lay.addWidget(self.progress_row)

        lay.addWidget(QLabel("Log"))
        self.log_edit = QTextEdit()
        self.log_edit.setReadOnly(True)
        self.log_edit.setMaximumHeight(120)
        self.log_edit.setStyleSheet(
            "background: #111; color: #ccc; font-family: monospace; font-size: 10px;"
        )
        lay.addWidget(self.log_edit)
        lay.addStretch(1)

        self._setup_signals()

    def _setup_signals(self):
        self.chk_apply_existing.toggled.connect(self._toggle_apply_mode)
        self.combo_hmm_n_states_mode.currentTextChanged.connect(self._toggle_n_states_mode)
        self.chk_hmm_sticky.toggled.connect(self._toggle_sticky_hmm)
        self.btn_run_state.clicked.connect(self._on_run_state)
        _wire_view_btn(self.btn_view_state, self._on_view, "state_qc")
        self.btn_rename_intrinsic.clicked.connect(self._on_rename_intrinsic)
        self.btn_rename_full.clicked.connect(self._on_rename_full)
        self.btn_state_diagnostics.clicked.connect(self._on_state_diagnostics)
        _wire_view_btn(self.btn_view_state_diagnostics, self._on_view, "state_diagnostics")
        self.btn_state_composition.clicked.connect(self._on_state_composition)
        _wire_view_btn(self.btn_view_composition, self._on_view, "state_composition")
        self.btn_state_transition.clicked.connect(self._on_state_transition)
        _wire_view_btn(self.btn_view_transition, self._on_view, "state_transition")
        self.btn_condition_comparison.clicked.connect(self._on_condition_comparison)
        _wire_view_btn(self.btn_view_condition_comparison, self._on_view, "state_condition_comparison")
        self.btn_browse_hmm.clicked.connect(self._browse_hmm_artifact)
        self.btn_apply_hmm.clicked.connect(self._on_apply_hmm)
        self.btn_show_state_bp.clicked.connect(self._on_show_state_bp)
        self.btn_export_state_bp.clicked.connect(self._on_export_state_bp)

        self.spin_window_size.valueChanged.connect(self._update_config_summary)
        self.chk_net_disp.toggled.connect(self._update_config_summary)
        self.chk_straight.toggled.connect(self._update_config_summary)
        self.chk_msd.toggled.connect(self._update_config_summary)
        self.spin_hmm_feature_smoothing_window.valueChanged.connect(self._update_config_summary)
        self.spin_quant_lo.valueChanged.connect(self._update_config_summary)
        self.spin_quant_hi.valueChanged.connect(self._update_config_summary)

        self._toggle_n_states_mode(self.combo_hmm_n_states_mode.currentText())
        self._toggle_sticky_hmm(self.chk_hmm_sticky.isChecked())

    # ── Prerequisite check ───────────────────────────────────────────────

    def _check_prerequisites(self) -> bool:
        """Check if track-features CSV exists; show warning and disable steps if not."""
        ct = self._cell_type()
        out = self._out_dir()
        if not ct or not out:
            self.warning_label.hide()
            return True

        csv_filtered = out / "analysis" / ct / "track_features" / f"BEHAV3D_{ct}_combined_track_features_filtered.csv"
        csv_unfiltered = out / "analysis" / ct / "track_features" / f"BEHAV3D_{ct}_combined_track_features.csv"

        if csv_filtered.exists():
            self.warning_label.hide()
            for grp in [self.grp_train, self.grp2, self._subtab_stack, self.grp_bp]:
                grp.setEnabled(True)
            _apply_group_tracked_gate(self.warning_label, self.grp_bp, self.metadata_loader, ct)
            return True
        elif csv_unfiltered.exists():
            self.warning_label.setText(
                f"⚠ Filtered data not found for cell type '{ct}'.\n"
                "Run Filtering first for full functionality.\n"
                "Using unfiltered data as a fallback."
            )
            self.warning_label.show()
            for grp in [self.grp_train, self.grp2, self._subtab_stack, self.grp_bp]:
                grp.setEnabled(True)
            _apply_group_tracked_gate(self.warning_label, self.grp_bp, self.metadata_loader, ct)
            return True
        else:
            self.warning_label.setText(
                f"⚠ Feature data not found for cell type '{ct}'.\n"
                "Run Filtering first, then return here."
            )
            self.warning_label.show()
            for grp in [self.grp_train, self.grp2, self._subtab_stack, self.grp_bp]:
                grp.setEnabled(False)
            return False

    # ── Guided pipeline dispatch (Step 3 reports) ────────────────────────
    _STATE_PIPELINE_RUN_BUTTONS = {
        "state_diagnostics": "btn_state_diagnostics",
        "state_transition": "btn_state_transition",
    }

    def _on_pipeline_start(self, pipeline_id: str):
        """Start from a pipeline card: run directly (no params) or open its settings."""
        from behav3d.napari.analysis_guided_copy import STATE_REPORT_PIPELINES
        spec = next((s for s in STATE_REPORT_PIPELINES if s["id"] == pipeline_id), None)
        if spec is not None and not spec.get("has_params", True):
            btn_name = self._STATE_PIPELINE_RUN_BUTTONS.get(pipeline_id)
            if btn_name is not None:
                getattr(self, btn_name).click()
            return
        self._focus_pipeline(pipeline_id)
        self._plots_stack.setCurrentIndex(1)
        self._pipeline_scroll.verticalScrollBar().setValue(0)

    def _focus_pipeline(self, pipeline_id: str):
        """Show only the group box relevant to ``pipeline_id``, hide the rest."""
        from behav3d.napari.analysis_guided_copy import STATE_REPORT_PIPELINES
        visible = {
            "state_diagnostics": {self.grp_state_diagnostics},
            "state_composition": {self.grp_state_composition},
            "state_transition": {self.grp_state_transition},
            "state_comparison": {self.grp_state_comparison},
        }.get(pipeline_id, set())
        for group in (
            self.grp_state_diagnostics, self.grp_state_composition,
            self.grp_state_transition, self.grp_state_comparison,
        ):
            group.setVisible(group in visible)
        title = next(
            (s["title"] for s in STATE_REPORT_PIPELINES if s["id"] == pipeline_id), ""
        )
        self._pipeline_title.setText(title)

    def _open_plotting_page(self):
        """Entry card Start: jump to the isolated plotting page, reset to the
        report list."""
        self._plots_stack.setCurrentIndex(0)
        self._pipeline_title.setText("")
        self._subtab_stack.setCurrentIndex(1)
        self._plotting_header.show()

    def _plotting_back(self):
        """Contextual Back inside the plotting page: report settings → report
        list; report list → clustering settings."""
        if self._plots_stack.currentIndex() == 1:
            self._plots_stack.setCurrentIndex(0)
            self._pipeline_title.setText("")
        else:
            self._subtab_stack.setCurrentIndex(0)
            self._plotting_header.hide()

    # ── Toggle helpers ───────────────────────────────────────────────────

    def _toggle_apply_mode(self, checked):
        """Hide Steps 1+2 when apply-existing mode is active."""
        self.grp_train.setVisible(not checked)
        self.grp2.setVisible(not checked)
        self.grp_apply.setVisible(checked)

    def _toggle_n_states_mode(self, mode):
        auto = mode == "auto"
        self.spin_hmm_n_states.setEnabled(not auto)
        self.spin_hmm_k_min.setVisible(auto)
        self.row_k_min.setVisible(auto)
        self.help_k_min.setVisible(auto)
        self.spin_hmm_k_max.setVisible(auto)
        self.row_k_max.setVisible(auto)
        self.help_k_max.setVisible(auto)

    def _toggle_sticky_hmm(self, checked):
        self.spin_hmm_stickiness_kappa.setVisible(checked)
        self.row_kappa.setVisible(checked)
        self.help_kappa.setVisible(checked)
        self.spin_hmm_transmat_alpha.setVisible(checked)
        self.row_alpha.setVisible(checked)
        self.help_alpha.setVisible(checked)

    def _update_config_summary(self):
        tp_feats = sorted(
            f for f, cb in getattr(self, "_timepoint_checkboxes", {}).items()
            if cb.isChecked()
        )
        log_feats = sorted(
            f for f, cb in getattr(self, "_logscale_checkboxes", {}).items()
            if cb.isChecked()
        )
        bin_groups = sorted(
            f for f, cb in getattr(self, "_bingrp_checkboxes", {}).items()
            if cb.isChecked()
        )
        win_feats = [
            n for n, chk in [
                ("net_displacement", self.chk_net_disp),
                ("straightness", self.chk_straight),
                ("mean_square_displacement", self.chk_msd),
            ] if chk.isChecked()
        ]
        lo = self.spin_quant_lo.value()
        hi = self.spin_quant_hi.value()
        capping = f"lower={lo:.2f}, upper={hi:.2f}" if (lo > 0 or hi < 1.0) else "none"
        lines = [
            f"Timepoint features ({len(tp_feats)}): {', '.join(tp_feats) or '—'}",
            f"Smoothing window: {self.spin_hmm_feature_smoothing_window.value()}",
            f"Window features: {', '.join(win_feats) or '—'}  |  Window size: {self.spin_window_size.value()}",
            f"Log-scaled: {', '.join(log_feats) or '—'}",
            f"Binary groups: {', '.join(bin_groups) or '—'}",
            f"Quantile capping: {capping}",
        ]
        self.config_summary_label.setText("\n".join(lines))
        self.config_summary_label.show()

    def _browse_hmm_artifact(self):
        fpath, _ = QFileDialog.getOpenFileName(
            self, "Select HMM Deployment Artifact", "", "Pickle Files (*.pkl)"
        )
        if fpath:
            self.le_hmm_artifact.setText(fpath)

    def _on_apply_hmm(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return

        pkl_path = self.le_hmm_artifact.text().strip()
        if not pkl_path or not Path(pkl_path).exists():
            QMessageBox.warning(self, "Error", "Please select a valid .pkl file.")
            return

        out = self._out_dir()
        self._log(f"▶ Applying existing classification for '{ct}'...")

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import apply_hmm_deployment_artifact_to_full_dataset
            apply_hmm_deployment_artifact_to_full_dataset(
                output_dir=str(out),
                cell_type=ct,
                hmm_deployment_artifact=Path(pkl_path),
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Apply classification ({ct})...",
            progress_row=self.progress_row,
            buttons=[self.btn_apply_hmm],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Application done for '{ct}'."),
                self._reload()
            ),
            on_failed=lambda e: self._log(f"❌ Application failed: {e}"),
        )

    # ── Metadata / reload ────────────────────────────────────────────────

    def on_metadata_updated(self):
        self._reload()

    def _update_timepoint_time_labels(self, *_):
        metadata = getattr(self.metadata_loader, "metadata", None)
        self.lbl_window_size_time.setText(
            format_timepoints_as_time(self.spin_window_size.value(), metadata)
        )
        self.lbl_hmm_smoothing_time.setText(
            format_timepoints_as_time(self.spin_hmm_feature_smoothing_window.value(), metadata)
        )

    def _reload(self):
        try:
            self._update_timepoint_time_labels()
            ct = self._cell_type()
            if not ct:
                self._model_adata = None
                self._refresh_buttons()
                return

            self._populate_dynamic_features(ct)

            self._check_prerequisites()
            self._update_view_buttons()
            self._update_bp_buttons()

            path = self._model_adata_path(ct)
            if not path or not path.exists():
                self._model_adata = None
                self._refresh_buttons()
                return

            if self._preload_bg.is_running():
                return

            self._model_adata = None
            self._refresh_buttons()

            def _load():
                import anndata as ad
                return ad.read_h5ad(str(path))

            def _on_done(result):
                self._model_adata = result
                self._refresh_buttons()

            def _on_failed(err):
                print(f"[BEHAV3D] State model adata failed to load: {err}")
                self._model_adata = None
                self._refresh_buttons()

            self._preload_bg.run(
                fn=_load,
                inject_progress=False,
                on_done=_on_done,
                on_failed=_on_failed,
            )
        except Exception:
            traceback.print_exc()

    def _rebuild_log_scale_features(self, state=None):
        while self.log_scale_lay.count():
            child = self.log_scale_lay.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
                
        selected_feats = [f for f, cb in getattr(self, "_timepoint_checkboxes", {}).items() if cb.isChecked()]
        ct = self._cell_type()
        cfg = getattr(self.metadata_loader, "behav3d_parameters", {}).get("state_classification", {}).get(ct, {}) if ct else {}
        saved_log_raw = cfg.get("log_scale_features", None)
        saved_log = set(saved_log_raw) if saved_log_raw is not None else (
            _DEFAULT_LOG_SCALE_FEATURES & set(selected_feats)
        )
        prev_checked = {f for f, cb in getattr(self, "_logscale_checkboxes", {}).items() if cb.isChecked()}
        
        self._logscale_checkboxes = {}
        if not selected_feats:
            self.log_scale_lay.addWidget(QLabel("<i>No features selected.</i>"))
            self._update_config_summary()
            return

        grid_log = QGridLayout()
        for i, f in enumerate(selected_feats):
            cb = QCheckBox(f)
            if f in prev_checked or f in saved_log:
                cb.setChecked(True)
            cb.stateChanged.connect(self._update_config_summary)
            self._logscale_checkboxes[f] = cb
            grid_log.addWidget(cb, i // 3, i % 3)
        self.log_scale_lay.addLayout(grid_log)
        self._update_config_summary()

    def _on_preview_feature_distributions(self, *_):
        """Show histograms of the continuous features for the selected cell
        type so the user can judge which benefit from log scaling."""
        from behav3d.napari._pdf_view import show_matplotlib_figure
        from behav3d.widgets.utils import build_feature_distribution_figure

        ct = self._cell_type()
        csv_path = getattr(self, "_logscale_csv_path", None)
        candidates = list(getattr(self, "_logscale_candidate_cols", []) or [])
        if not csv_path or not Path(csv_path).exists() or not candidates:
            QMessageBox.information(
                self,
                "Feature distributions",
                "No feature table is loaded yet for this cell type. "
                "Select a cell type with extracted track features first.",
            )
            return

        # Prefer the features the user has actually selected (those shown in the
        # log-scaling list); fall back to all continuous candidates.
        selected = [f for f, cb in getattr(self, "_timepoint_checkboxes", {}).items()
                    if cb.isChecked() and f in candidates]
        feats = selected or candidates
        title = f"Feature distributions — {ct}"
        fig, truncated = build_feature_distribution_figure(csv_path, feats, title=title)
        if fig is None:
            QMessageBox.warning(
                self, "Feature distributions",
                f"Could not read features from:\n{csv_path}",
            )
            return
        if truncated:
            title += "  (showing first 36)"
        show_matplotlib_figure(fig, title=title, parent=self)

    def _populate_dynamic_features(self, ct):
        from behav3d.widgets.utils import behav3d_calculated_features, excluded_non_behavior_columns
        from behav3d.core.utils import expand_column_patterns
        import pandas as pd
        out = self._out_dir()
        if not out:
            return

        base = out / "analysis" / ct / "track_features"
        csv_path = base / f"BEHAV3D_{ct}_combined_track_features_filtered.csv"
        if not csv_path.exists():
            csv_path = base / f"BEHAV3D_{ct}_combined_track_features.csv"

        try:
            mtime = csv_path.stat().st_mtime if csv_path.exists() else None
        except OSError:
            mtime = None
        features_key = (ct, str(csv_path), mtime)
        if features_key == self._last_features_key:
            return

        self._timepoint_checkboxes = {}
        self._logscale_checkboxes = {}
        self._bingrp_checkboxes = {}

        self.setUpdatesEnabled(False)
        try:
            for lay in [self.timepoint_features_lay, self.log_scale_lay, self.bin_grp_lay]:
                while lay.count():
                    child = lay.takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()

            if not csv_path.exists():
                self.timepoint_features_lay.addWidget(_make_info_label(
                    "<i>No track-features CSV found for this cell type.\n Run feature extraction first.</i>"
                ))
                self.log_scale_lay.addWidget(_make_info_label(
                    "<i>No track-features CSV found.</i>"
                ))
                self.bin_grp_lay.addWidget(_make_info_label(
                    "<i>No binary columns detected yet.</i>"
                ))
                for w in [self.spin_window_size, self.chk_net_disp, self.chk_straight, self.chk_msd,
                          self.spin_hmm_feature_smoothing_window, self.spin_quant_lo, self.spin_quant_hi]:
                    w.setEnabled(False)
                self._last_features_key = features_key
                return

            # Read at least 5 rows so pandas can infer column dtypes properly.
            # If nrows=0, many columns default to dtype 'object'.
            for w in [self.spin_window_size, self.chk_net_disp, self.chk_straight, self.chk_msd,
                      self.spin_hmm_feature_smoothing_window, self.spin_quant_lo, self.spin_quant_hi]:
                w.setEnabled(True)

            ct = self._cell_type()
            cfg = getattr(self.metadata_loader, "behav3d_parameters", {}).get("state_classification", {}).get(ct, {})

            # ── Features ────────────────────────────────────────────────────────
            saved_features_raw = cfg.get("selected_features", cfg.get("features", None))
            saved_features = set(saved_features_raw) if saved_features_raw is not None else _DEFAULT_TIMEPOINT_FEATURES
            saved_bin = set(cfg.get("binary_features_to_group", []))

            # ── Window features ──────────────────────────────────────────────────
            saved_window_raw = cfg.get("additional_window_features", None)
            saved_window = set(saved_window_raw) if saved_window_raw is not None else _DEFAULT_WINDOW_FEATURES
            self.chk_net_disp.setChecked("net_displacement" in saved_window)
            self.chk_straight.setChecked("straightness" in saved_window)
            self.chk_msd.setChecked("mean_square_displacement" in saved_window)

            # ── Feature Processing ───────────────────────────────────────────────
            self.spin_window_size.setValue(int(cfg.get("window_features_window", self.spin_window_size.value())))
            self.spin_hmm_feature_smoothing_window.setValue(int(cfg.get("feature_smoothing_window", self.spin_hmm_feature_smoothing_window.value())))
            lo = cfg.get("lower_quantile_cap", None)
            self.spin_quant_lo.setValue(float(lo) if lo is not None else 0.0)
            hi = cfg.get("upper_quantile_cap", None)
            self.spin_quant_hi.setValue(float(hi) if hi is not None else 1.0)

            # ── HMM Parameters ───────────────────────────────────────────────────
            n_states_cfg = cfg.get("n_states", None)
            if isinstance(n_states_cfg, int):
                self.combo_hmm_n_states_mode.setCurrentText("fixed")
                self.spin_hmm_n_states.setValue(n_states_cfg)
            elif n_states_cfg == "auto":
                self.combo_hmm_n_states_mode.setCurrentText("auto")
            if "k_min" in cfg:
                self.spin_hmm_k_min.setValue(int(cfg["k_min"]))
            if "k_max" in cfg:
                self.spin_hmm_k_max.setValue(int(cfg["k_max"]))
            if "start_offset" in cfg:
                self.spin_hmm_start_offset.setValue(int(cfg["start_offset"]))
            if "start_offset_fill_mode" in cfg:
                self.combo_hmm_start_offset_fill_mode.setCurrentText(cfg["start_offset_fill_mode"])
            if "covariance_type" in cfg:
                self.combo_hmm_covariance_type.setCurrentText(cfg["covariance_type"])
            if "n_iter" in cfg:
                self.spin_hmm_n_iter.setValue(int(cfg["n_iter"]))
            if "tol" in cfg:
                self.spin_hmm_tol.setValue(float(cfg["tol"]))
            if "min_covar" in cfg:
                self.spin_hmm_min_covar.setValue(float(cfg["min_covar"]))
            if "sticky" in cfg:
                self.chk_hmm_sticky.setChecked(bool(cfg["sticky"]))
            if "stickiness_kappa" in cfg:
                self.spin_hmm_stickiness_kappa.setValue(float(cfg["stickiness_kappa"]))
            if "transmat_alpha" in cfg:
                self.spin_hmm_transmat_alpha.setValue(float(cfg["transmat_alpha"]))
            if "random_state" in cfg:
                self.spin_seed.setValue(int(cfg["random_state"]))

            cols = list(pd.read_csv(csv_path, nrows=0).columns)
            md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
            excluded = excluded_non_behavior_columns(cols, metadata=md)
            usable_cols = [c for c in cols if c not in excluded]
            # Value-based binary detection over the full CSV (see
            # behav3d.widgets.base_state_classification.detect_binary_columns_from_csv).
            # The previous 5-row dtype heuristic mis-classified numeric feature
            # columns as "binary" whenever the sampled rows were NaN/blank, so
            # switching cell types could dump every feature into the binary list.
            from behav3d.widgets.base_state_classification import (
                detect_binary_columns_from_csv,
                detect_non_numeric_columns_from_csv,
            )
            bin_cols = detect_binary_columns_from_csv(Path(csv_path), usable_cols)
            bin_set = set(bin_cols)
            # Columns that aren't binary and can't be parsed as continuous numbers
            # (e.g. "touching_27ts" holding comma-separated contact-ID lists, or
            # unit/label columns) must not be offered as selectable HMM features --
            # picking one silently breaks the .h5ad write later on.
            non_feature_candidates = [c for c in usable_cols if c not in bin_set]
            non_numeric_cols = set(
                detect_non_numeric_columns_from_csv(Path(csv_path), non_feature_candidates)
            )
            feat_cols = [c for c in non_feature_candidates if c not in non_numeric_cols]
            # Continuous features eligible for log scaling; reused by the
            # "Preview feature distributions" histogram button.
            self._logscale_candidate_cols = list(feat_cols)
            self._logscale_csv_path = Path(csv_path)

            from copy import deepcopy
            base_groups = deepcopy(behav3d_calculated_features)
            matched = set()
            for gname, patterns in base_groups.items():
                vals = []
                for pat in patterns:
                    vals.extend(expand_column_patterns(pat, feat_cols))
                clean_vals = sorted({x for x in vals if x in feat_cols})
                if clean_vals:
                    group_sec = CollapsibleSection(gname, expanded=False)
                    group_content = QWidget()
                    grid = QGridLayout(group_content)
                    for i, f in enumerate(clean_vals):
                        cb = QCheckBox(f)
                        if f in saved_features:
                            cb.setChecked(True)
                        cb.stateChanged.connect(self._rebuild_log_scale_features)
                        self._timepoint_checkboxes[f] = cb
                        grid.addWidget(cb, i // 3, i % 3)
                    group_sec.addWidget(group_content)
                    self.timepoint_features_lay.addWidget(group_sec)
                    matched.update(clean_vals)

            other = sorted([c for c in feat_cols if c not in matched])
            if other:
                group_sec = CollapsibleSection("other", expanded=False)
                group_content = QWidget()
                grid = QGridLayout(group_content)
                for i, f in enumerate(other):
                    cb = QCheckBox(f)
                    if f in saved_features:
                        cb.setChecked(True)
                    cb.stateChanged.connect(self._rebuild_log_scale_features)
                    self._timepoint_checkboxes[f] = cb
                    grid.addWidget(cb, i // 3, i % 3)
                group_sec.addWidget(group_content)
                self.timepoint_features_lay.addWidget(group_sec)

            self._rebuild_log_scale_features()

            if not bin_cols:
                self.bin_grp_lay.addWidget(QLabel("<i>No binary columns detected yet.</i>"))
            else:
                for b in bin_cols:
                    cb = QCheckBox(b)
                    if b in saved_bin:
                        cb.setChecked(True)
                    cb.stateChanged.connect(self._update_config_summary)
                    self._bingrp_checkboxes[b] = cb
                    self.bin_grp_lay.addWidget(cb)
                self._update_config_summary()

            self._last_features_key = features_key
        finally:
            self.setUpdatesEnabled(True)

    def _model_adata_path(self, cell_type: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return (
            out / "analysis" / cell_type / "behavioral_states" / "processing"
            / f"BEHAV3D_{cell_type}_behavioral_states_modeldata.h5ad"
        )

    def _full_adata_path(self, cell_type: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return (
            out / "analysis" / cell_type / "behavioral_states"
            / f"BEHAV3D_{cell_type}_behavioral_states.h5ad"
        )

    def _classifier_path(self, cell_type: str, kind="full") -> Optional[Path]:
        """Deprecated: kept for view-button compatibility."""
        out = self._out_dir()
        if not out:
            return None
        return (
            out / "analysis" / cell_type / "behavioral_states" / "processing"
            / "full_behavioral_classification"
            / f"state_classification_random_forest_{kind}_{cell_type}.pkl"
        )

    def _report_path(self, cell_type: str, report: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        state_base = out / "analysis" / cell_type / "behavioral_states"
        _REPORT_PATHS = {
            "state_composition_report": state_base / "state_composition" / "state_composition_report.pdf",
            "state_transition_report": state_base / "state_transitions" / "transition_matrix" / "transition_matrix_heatmap.pdf",
        }
        return _REPORT_PATHS.get(report)

    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

    def _state_run_qc_pdfs(self, ct: str) -> list:
        """The two QC PDFs written by a State Classification run itself.

        `run_hmm_state_clustering` puts these under the ``raw/`` sub-folder of
        the QC directory (see `_resolve_hmm_quality_control_outdir`), which is
        what separates them from the re-generated PDFs that the Step 3
        "Create Diagnostics" button writes to the QC folder root.
        """
        qc_dir = self._state_diagnostics_qc_dir(ct)
        if not qc_dir:
            return []
        raw = qc_dir / "raw"
        return [
            ("Clustering diagnostics", raw / "behavioral_clustering_diagnostics.pdf"),
            ("Feature distributions", raw / "behavioral_clustering_feature_distributions.pdf"),
        ]

    def _state_diagnostics_qc_dir(self, ct: str) -> Optional[Path]:
        """Path to the HMM diagnostics QC folder, without creating it."""
        out = self._out_dir()
        if not out or not ct:
            return None
        return (
            out / "analysis" / ct / "behavioral_states" / "processing"
            / "hmm_behavioral_classification" / "quality_control"
        )

    def _resolve_hmm_model_for_state_diagnostics(self, ct: str):
        """Return the cached HMM model for `ct`, or load it from the saved deployment artifact."""
        from behav3d.analysis.behavior.state.classification import (
            load_hmm_deployment_artifact,
            _resolve_hmm_deployment_artifact_path,
        )
        hmm_model = (
            getattr(self, "_hmm_model", None)
            if getattr(self, "_hmm_model_cell_type", None) == ct
            else None
        )
        if hmm_model is None:
            out = self._out_dir()
            artifact_path = _resolve_hmm_deployment_artifact_path(
                output_dir=str(out) if out else "", cell_type=ct
            )
            if artifact_path.exists():
                try:
                    stored = load_hmm_deployment_artifact(str(artifact_path))
                    hmm_model = stored.get("model")
                    self._hmm_model = hmm_model
                    self._hmm_model_cell_type = ct
                except Exception:
                    pass
        return hmm_model

    def _regenerate_state_diagnostics(self, ct: str, *, model_adata=None, verbose: bool = True):
        """(Re)generate the HMM state-classification diagnostics PDF/CSVs for `ct`."""
        from behav3d.analysis.behavior.state.classification import (
            INTRINSIC_STATE_COL,
            save_hmm_quality_control_outputs,
            _resolve_hmm_quality_control_outdir,
        )
        adata = model_adata if model_adata is not None else self._model_adata
        if adata is None:
            raise ValueError("No model adata loaded.")
        out = self._out_dir()
        preprocessing_meta = adata.uns.get("preprocessing", {})
        if not isinstance(preprocessing_meta, dict):
            preprocessing_meta = {}
        feature_cols = preprocessing_meta.get(
            "continuous_feature_cols",
            preprocessing_meta.get("kept_features", list(adata.var_names)),
        )
        feature_cols = [] if feature_cols is None else list(feature_cols)
        feature_cols = [str(c) for c in feature_cols if str(c) in adata.var_names]
        scaler_meta = preprocessing_meta.get("scaler", {}) if isinstance(preprocessing_meta, dict) else {}
        qc_dir = _resolve_hmm_quality_control_outdir(output_dir=str(out) if out else "", cell_type=ct)
        hmm_model = self._resolve_hmm_model_for_state_diagnostics(ct)
        return save_hmm_quality_control_outputs(
            adata,
            feature_cols=feature_cols,
            output_dir=qc_dir,
            model=hmm_model,
            selection_df=None,
            cluster_col=INTRINSIC_STATE_COL,
            scaler_mean=scaler_meta.get("mean", None) if isinstance(scaler_meta, dict) else None,
            scaler_scale=scaler_meta.get("scale", None) if isinstance(scaler_meta, dict) else None,
            title=f"all_data | hmm | curated {INTRINSIC_STATE_COL}",
            preprocessing_params=preprocessing_meta,
            verbose=verbose,
        )

    def _load_model_adata(self, path: Path):
        try:
            import anndata as ad
            self._model_adata = ad.read_h5ad(str(path))
        except Exception:
            traceback.print_exc()
            self._model_adata = None

    # ── Button state management ──────────────────────────────────────────

    def _refresh_buttons(self):
        has_model = self._model_adata is not None
        has_intrinsic = has_model and "intrinsic_behavioral_cluster" in self._model_adata.obs.columns
        has_full = has_model and "full_behavioral_cluster" in self._model_adata.obs.columns

        self.btn_rename_intrinsic.setEnabled(has_intrinsic)
        self.btn_rename_full.setEnabled(has_full)
        self.btn_state_diagnostics.setEnabled(has_intrinsic)

        if has_intrinsic:
            n_intr = self._model_adata.obs["intrinsic_behavioral_cluster"].astype(str).nunique()
            n_full = (
                self._model_adata.obs["full_behavioral_cluster"].astype(str).nunique()
                if has_full else 0
            )
            self.rename_status_lbl.setText(
                f"✅ Model loaded: {n_intr} intrinsic clusters"
                + (f", {n_full} full clusters." if has_full else ".")
            )
        else:
            self.rename_status_lbl.setText("ℹ Run state classification first to enable renaming.")

        self._refresh_composition_group_cols()

    def _composition_candidate_columns(self):
        cols = []
        added = set()

        # Prefer columns from the metadata CSV, but only surface the known
        # experimental-design grouping columns (exp_nr, well, *_line_condition).
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        if md is not None and hasattr(md, "columns"):
            for col in md.columns:
                if col in ("exp_nr", "well") or col.endswith("_line_condition"):
                    cols.append(col)
                    added.add(col)

        # Metadata CSV never carries per-cell columns like `origin_cell_type`
        # (stamped during multi-population track merging), so also check the
        # adata/h5ad obs columns for those, whether or not metadata was loaded.
        if self._model_adata is not None:
            obs_cols = list(self._model_adata.obs.columns)
        else:
            ct = self._cell_type()
            full_path = self._full_adata_path(ct) if ct else None
            if not full_path or not full_path.exists():
                return cols
            try:
                import h5py
                with h5py.File(str(full_path), "r") as f:
                    obs_cols = list(f.get("obs", {}).keys())
            except Exception:
                return cols
        for col in obs_cols:
            if col in added or col.startswith("_") or col in _TECHNICAL_OBS_COLS:
                continue
            if col in ("exp_nr", "well", "origin_cell_type") or col.endswith("_line_condition"):
                cols.append(col)
                added.add(col)
        return cols

    def _refresh_composition_group_cols(self):
        candidate_cols = self._composition_candidate_columns()
        self.list_composition_group_cols.clear()
        for col in candidate_cols:
            self.list_composition_group_cols.addItem(col)
        _fit_list_widget_height(self.list_composition_group_cols)

        for combo in (self.combo_composition_group_x, self.combo_composition_group_y):
            prev = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItem("(none)")
            combo.addItems(candidate_cols)
            combo.setCurrentText(prev if prev in (["(none)"] + candidate_cols) else "(none)")
            combo.blockSignals(False)
        self._refresh_composition_group_x_levels()
        self._refresh_composition_group_y_levels()

        prev_cond = self.combo_comparison_condition_col.currentText()
        self.combo_comparison_condition_col.blockSignals(True)
        self.combo_comparison_condition_col.clear()
        self.combo_comparison_condition_col.addItems(candidate_cols)
        if prev_cond in candidate_cols:
            self.combo_comparison_condition_col.setCurrentText(prev_cond)
        self.combo_comparison_condition_col.blockSignals(False)

        self.list_comparison_group_cols.clear()
        for col in candidate_cols:
            self.list_comparison_group_cols.addItem(col)
        _fit_list_widget_height(self.list_comparison_group_cols)

        prev = self.combo_comparison_group_x.currentText()
        self.combo_comparison_group_x.blockSignals(True)
        self.combo_comparison_group_x.clear()
        self.combo_comparison_group_x.addItem("(none)")
        self.combo_comparison_group_x.addItems(candidate_cols)
        self.combo_comparison_group_x.setCurrentText(
            prev if prev in (["(none)"] + candidate_cols) else "(none)"
        )
        self.combo_comparison_group_x.blockSignals(False)

        self._sync_comparison_group_y_text()
        self._refresh_comparison_group_levels()
        self._refresh_comparison_group_x_levels()

    def _on_composition_group_x_changed(self, _text):
        self._refresh_composition_group_x_levels()

    def _on_composition_group_y_changed(self, _text):
        self._refresh_composition_group_y_levels()

    def _refresh_composition_group_x_levels(self):
        col = self.combo_composition_group_x.currentText()
        col = None if col in ("", "(none)") else col
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._model_adata,
            h5ad_path=self._model_adata_path(ct) if ct else None,
        ) if col else []
        self.group_selector_composition_x.set_items(levels)

    def _refresh_composition_group_y_levels(self):
        col = self.combo_composition_group_y.currentText()
        col = None if col in ("", "(none)") else col
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._model_adata,
            h5ad_path=self._model_adata_path(ct) if ct else None,
        ) if col else []
        self.group_selector_composition_y.set_items(levels)

    def _on_composition_group_x_conditions_toggled(self, checked: bool):
        self.group_selector_composition_x.setVisible(checked)
        if checked:
            self._refresh_composition_group_x_levels()

    def _on_composition_group_y_conditions_toggled(self, checked: bool):
        self.group_selector_composition_y.setVisible(checked)
        if checked:
            self._refresh_composition_group_y_levels()

    def _sync_comparison_group_y_text(self):
        self.line_comparison_group_y.setText(self.combo_comparison_condition_col.currentText())

    def _on_comparison_condition_col_changed(self, _text):
        self._sync_comparison_group_y_text()
        self._refresh_comparison_group_levels()

    def _refresh_comparison_group_levels(self):
        col = self.combo_comparison_condition_col.currentText()
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._model_adata,
            h5ad_path=self._model_adata_path(ct) if ct else None,
        )
        self.group_selector_comparison.set_items(levels)

    def _on_comparison_group_conditions_toggled(self, checked: bool):
        self.group_selector_comparison.setVisible(checked)
        if checked:
            self._refresh_comparison_group_levels()

    def _on_comparison_group_x_changed(self, _text):
        self._refresh_comparison_group_x_levels()

    def _refresh_comparison_group_x_levels(self):
        col = self.combo_comparison_group_x.currentText()
        col = None if col in ("", "(none)") else col
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._model_adata,
            h5ad_path=self._model_adata_path(ct) if ct else None,
        ) if col else []
        self.group_selector_comparison_x.set_items(levels)

    def _on_comparison_group_x_conditions_toggled(self, checked: bool):
        self.group_selector_comparison_x.setVisible(checked)
        if checked:
            self._refresh_comparison_group_x_levels()

    def _update_view_buttons(self):
        ct = self._cell_type()
        if not ct:
            for btn in (
                self.btn_view_state, self.btn_view_state_diagnostics, self.btn_view_composition,
                self.btn_view_transition, self.btn_view_condition_comparison,
            ):
                btn.setEnabled(False)
            return
        self.btn_view_state.setEnabled(
            any(p.exists() for _lbl, p in self._state_run_qc_pdfs(ct))
        )
        qc_dir = self._state_diagnostics_qc_dir(ct)
        self.btn_view_state_diagnostics.setEnabled(
            bool(qc_dir and qc_dir.exists() and any(qc_dir.glob("*.pdf")))
        )
        comp = self._report_path(ct, "state_composition_report")
        self.btn_view_composition.setEnabled(bool(comp and comp.exists()))
        trans = self._report_path(ct, "state_transition_report")
        self.btn_view_transition.setEnabled(bool(trans and trans.exists()))
        self.btn_view_condition_comparison.setEnabled(
            len(self._get_view_candidates("state_condition_comparison", ct)) > 0
        )

    def _update_bp_buttons(self):
        ct = self._cell_type()
        full_path = self._full_adata_path(ct) if ct else None
        model_path = self._model_adata_path(ct) if ct else None
        has_full = bool(full_path and full_path.exists())
        has_model = bool(model_path and model_path.exists())
        self.btn_show_state_bp.setEnabled(has_full or has_model)
        self.btn_export_state_bp.setEnabled(has_full)

    # ── Click handlers ───────────────────────────────────────────────────

    def _cell_type(self) -> str:
        return self._get_cell_type() if self._get_cell_type else ""

    def _sample(self) -> str:
        txt = self.sample_combo.currentText()
        return "" if txt.startswith("—") else txt

    def _log(self, msg: str):
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        formatted = f"[{ts}] {msg}"
        self.log_edit.append(formatted)
        print(formatted)

    # ── Open-folder pop-up ─────────────────────────────────────────────
    def _offer_open_results_folder(self, folder: Path, what: str = "State Classification"):
        """Pop a dialog offering to open the results folder in the OS file manager.

        Silently no-ops if ``folder`` does not exist (e.g. analysis failed).
        """
        try:
            if not folder or not Path(folder).exists():
                return
            box = QMessageBox(self)
            box.setIcon(QMessageBox.Information)
            box.setWindowTitle(f"{what} complete")
            box.setText(f"{what} analysis is complete.")
            box.setInformativeText(
                f"Results have been saved to:\n{folder}\n\n"
                "Would you like to open the results folder?"
            )
            btn_open = box.addButton("Open folder", QMessageBox.AcceptRole)
            box.addButton("Close", QMessageBox.RejectRole)
            box.setDefaultButton(btn_open)
            box.exec_()
            if box.clickedButton() is btn_open:
                from qtpy.QtCore import QUrl
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))
        except Exception as e:
            self._log(f"Could not show results-folder dialog: {e}")

    def _on_run_state(self):
        ct = self._cell_type()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        self._log(f"▶ Running state classification for '{ct}'…")
        params = self._collect_state_params(ct)

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import run_hmm_state_clustering
            res = run_hmm_state_clustering(**params, verbose=True, return_details=True)
            self._hmm_model = res["hmm_model"]
            self._hmm_model_cell_type = ct
            return res["model_adata"]

        self._bg.run(
            fn=_run,
            desc=f"State classification ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_run_state],
            viewer=self.viewer,
            inject_progress=False,
            on_done=self._on_state_done,
            on_failed=lambda e: self._log(f"❌ State classification failed: {e}"),
        )

    def _collect_state_params(self, ct: str) -> dict:
        feats = [f for f, cb in getattr(self, "_timepoint_checkboxes", {}).items() if cb.isChecked()]
        bin_grps = [b for b, cb in getattr(self, "_bingrp_checkboxes", {}).items() if cb.isChecked()]
        log_feats = [f for f, cb in getattr(self, "_logscale_checkboxes", {}).items() if cb.isChecked()]

        # Collect window-feature checkboxes and pass under the correct backend key.
        win_feats = []
        if self.chk_net_disp.isChecked():
            win_feats.append("net_displacement")
        if self.chk_straight.isChecked():
            win_feats.append("straightness")
        if self.chk_msd.isChecked():
            win_feats.append("mean_square_displacement")

        out_dir = self._out_dir()
        lower_cap = self.spin_quant_lo.value()
        upper_cap = self.spin_quant_hi.value()
        return {
            "features": feats,
            "additional_window_features": win_feats if win_feats else None,
            "binary_features_to_group": bin_grps,
            "output_dir": str(out_dir) if out_dir else "",
            "cell_type": ct,
            "log_scale_features": log_feats,
            "window_features_window": self.spin_window_size.value(),
            "feature_smoothing_window": self.spin_hmm_feature_smoothing_window.value(),
            "lower_quantile_cap": lower_cap if lower_cap > 0 else None,
            "upper_quantile_cap": upper_cap if upper_cap < 1.0 else None,
            "n_states": (
                self.spin_hmm_n_states.value()
                if self.combo_hmm_n_states_mode.currentText() == "fixed"
                else "auto"
            ),
            "k_min": self.spin_hmm_k_min.value(),
            "k_max": self.spin_hmm_k_max.value(),
            "start_offset": self.spin_hmm_start_offset.value(),
            "start_offset_fill_mode": self.combo_hmm_start_offset_fill_mode.currentText(),
            "covariance_type": self.combo_hmm_covariance_type.currentText(),
            "n_iter": self.spin_hmm_n_iter.value(),
            "tol": self.spin_hmm_tol.value(),
            "min_covar": self.spin_hmm_min_covar.value(),
            "sticky": self.chk_hmm_sticky.isChecked(),
            "stickiness_kappa": self.spin_hmm_stickiness_kappa.value(),
            "transmat_alpha": self.spin_hmm_transmat_alpha.value(),
            "random_state": self.spin_seed.value(),
        }

    def _persist_state_cfg(self, ct: str):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return
        collected = self._collect_state_params(ct)
        collected["selected_features"] = collected.pop("features", collected.get("selected_features", []))
        params.setdefault("state_classification", {})[ct] = collected
        _save_behav3d_params(self.metadata_loader, self._out_dir)

    def _on_state_done(self, result, interactive: bool = True):
        ct = self._cell_type()
        self._persist_state_cfg(ct)
        self._log(f"✅ State classification complete for '{ct}'.")
        path = self._model_adata_path(ct)
        if path and path.exists():
            self._load_model_adata(path)
        if self._hmm_model is not None and ct:
            out = self._out_dir()
            if out:
                from behav3d.analysis.behavior.state.classification import (
                    save_hmm_deployment_artifact,
                    _resolve_hmm_deployment_artifact_path,
                )
                _art_path = _resolve_hmm_deployment_artifact_path(output_dir=str(out), cell_type=ct)
                try:
                    save_hmm_deployment_artifact(
                        output_path=_art_path,
                        model_adata=self._model_adata,
                        hmm_model=self._hmm_model,
                        verbose=False,
                    )
                except Exception:
                    pass
        self._refresh_buttons()
        self._update_view_buttons()
        self._update_bp_buttons()
        self._notify_results()
        QTimer.singleShot(0, lambda _ct=ct: self._apply_to_full_dataset_after_rename(_ct))
        if interactive:
            diagnostics_pdf = None
            try:
                diagnostics_pdf = result.uns.get("clustering", {}).get("diagnostics_pdf")
            except Exception:
                diagnostics_pdf = None
            if diagnostics_pdf:
                self._offer_open_results_folder(
                    Path(diagnostics_pdf).parent, what="State Classification"
                )

    def run_state_classification(self, interactive=True, extra_callbacks=None):
        """Called from queue runner (_queue.py)."""
        ct = self._cell_type()
        if not ct:
            err = "No cell type selected for state classification."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        if self._bg.is_running():
            err = "BackgroundOperation already running."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return

        params = self._collect_state_params(ct)
        on_done_cb = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_cb = extra_callbacks.get("on_failed") if extra_callbacks else None

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import run_hmm_state_clustering
            res = run_hmm_state_clustering(**params, verbose=True, return_details=True)
            self._hmm_model = res["hmm_model"]
            self._hmm_model_cell_type = ct
            return res["model_adata"]

        def _done(result):
            self._on_state_done(result, interactive=interactive)
            if on_done_cb:
                on_done_cb(result)

        def _fail(err):
            self._log(f"❌ {err}")
            if on_fail_cb:
                on_fail_cb(err)

        self._bg.run(
            fn=_run,
            desc=f"State classification ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_run_state],
            viewer=self.viewer,
            inject_progress=False,
            on_done=_done,
            on_failed=_fail,
        )

    def _on_rename_intrinsic(self):
        ct = self._cell_type()
        if self._model_adata is None:
            QMessageBox.warning(self, "No model", "Run state classification first.")
            return
        dlg = RenameClusterDialog(
            mode="intrinsic",
            model_adata=self._model_adata,
            adata_path=self._model_adata_path(ct),
            parent=self,
        )
        dlg.clusters_renamed.connect(
            lambda mapping: self._log(f"✅ Intrinsic clusters renamed: {mapping}")
        )
        if dlg.exec_() == QDialog.Accepted:
            self._refresh_buttons()
            self._update_view_buttons()
            self._apply_to_full_dataset_after_rename(ct)

    def _on_rename_full(self):
        ct = self._cell_type()
        if self._model_adata is None:
            QMessageBox.warning(self, "No model", "Run state classification first.")
            return
        if "full_behavioral_cluster" not in self._model_adata.obs.columns:
            QMessageBox.warning(
                self, "No full clusters",
                "Run intrinsic cluster renaming first (full clusters are created when "
                "intrinsic states are combined with binary groups)."
            )
            return
        dlg = RenameClusterDialog(
            mode="full",
            model_adata=self._model_adata,
            adata_path=self._model_adata_path(ct),
            parent=self,
        )
        dlg.clusters_renamed.connect(
            lambda mapping: self._log(f"✅ Full behavioral clusters renamed: {mapping}")
        )
        if dlg.exec_() == QDialog.Accepted:
            self._refresh_buttons()
            self._update_view_buttons()
            self._apply_to_full_dataset_after_rename(ct)

    def _apply_to_full_dataset_after_rename(self, ct: str):
        """Rebuild the pkl artifact with renamed labels and apply it to the full dataset."""
        out = self._out_dir()
        if not out or self._model_adata is None or not ct:
            return
        if self._bg.is_running():
            self._log("⚠ Background job running — skipping auto-apply to full dataset.")
            return

        from behav3d.analysis.behavior.state.classification import (
            INTRINSIC_STATE_COL,
            FULL_STATE_COL,
            save_hmm_deployment_artifact,
            apply_hmm_deployment_artifact_to_full_dataset,
            _resolve_hmm_deployment_artifact_path,
        )
        from behav3d.analysis.behavior.state.utils import (
            _get_classification_state_colors,
            _get_classification_state_order,
            _set_classification_state_colors,
            _set_classification_state_order,
        )

        artifact_path = _resolve_hmm_deployment_artifact_path(output_dir=str(out), cell_type=ct)

        hmm_model = self._resolve_hmm_model_for_state_diagnostics(ct)

        if hmm_model is None:
            self._log(
                "⚠ Cannot apply states to full dataset: HMM model not in memory and no saved artifact found. "
                "Re-run Step 1 to regenerate."
            )
            return

        # Sync columns so the artifact builder uses the renamed labels.
        if "intrinsic_behavioral_cluster" in self._model_adata.obs.columns:
            self._model_adata.obs[INTRINSIC_STATE_COL] = (
                self._model_adata.obs["intrinsic_behavioral_cluster"].copy()
            )
        if "full_behavioral_cluster" in self._model_adata.obs.columns:
            self._model_adata.obs[FULL_STATE_COL] = (
                self._model_adata.obs["full_behavioral_cluster"].copy()
            )

        # Sync saved colors/order too — the rename dialog saves them under the raw
        # dialog obs-column names, not the classification-constant column names the
        # deployment artifact/full-dataset pipeline looks them up under.
        raw_intrinsic_colors = _get_classification_state_colors(self._model_adata, "intrinsic_behavioral_cluster")
        raw_intrinsic_order = _get_classification_state_order(self._model_adata, "intrinsic_behavioral_cluster")
        if raw_intrinsic_colors:
            _set_classification_state_colors(self._model_adata, INTRINSIC_STATE_COL, raw_intrinsic_colors)
        if raw_intrinsic_order:
            _set_classification_state_order(self._model_adata, INTRINSIC_STATE_COL, raw_intrinsic_order)

        raw_full_colors = _get_classification_state_colors(self._model_adata, "full_behavioral_cluster")
        raw_full_order = _get_classification_state_order(self._model_adata, "full_behavioral_cluster")
        if raw_full_colors:
            _set_classification_state_colors(self._model_adata, FULL_STATE_COL, raw_full_colors)
        if raw_full_order:
            _set_classification_state_order(self._model_adata, FULL_STATE_COL, raw_full_order)

        try:
            save_hmm_deployment_artifact(
                output_path=artifact_path,
                model_adata=self._model_adata,
                hmm_model=hmm_model,
                cell_type=ct,
                output_dir=str(out),
                verbose=False,
            )
        except Exception as exc:
            self._log(f"⚠ Could not save updated HMM artifact: {exc}")
            return

        self._log(f"▶ Applying state labels to full dataset for '{ct}'…")
        _artifact_path = artifact_path
        _out = out
        model_adata = self._model_adata

        def _run(**kw):
            apply_hmm_deployment_artifact_to_full_dataset(
                output_dir=str(_out),
                cell_type=ct,
                hmm_deployment_artifact=_artifact_path,
                verbose=True,
            )
            diagnostics_warning = None
            try:
                self._regenerate_state_diagnostics(ct, model_adata=model_adata, verbose=True)
            except Exception as exc:
                diagnostics_warning = str(exc)
            return {"diagnostics_warning": diagnostics_warning}

        self._bg.run(
            fn=_run,
            desc=f"Apply states to full dataset ({ct})…",
            progress_row=self.progress_row,
            buttons=[],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ State labels applied to full dataset for '{ct}'."),
                self._log(f"⚠ Could not refresh state diagnostics: {r.get('diagnostics_warning')}")
                if isinstance(r, dict) and r.get("diagnostics_warning")
                else self._log("✅ State diagnostics refreshed."),
                self._update_bp_buttons(),
                self._update_view_buttons(),
            ),
            on_failed=lambda e: self._log(f"❌ Apply to full dataset failed: {e}"),
        )

    def _on_state_diagnostics(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._model_adata is None:
            QMessageBox.warning(self, "No model", "Run state classification first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        self._log(f"▶ Creating state diagnostics for '{ct}'…")
        model_adata = self._model_adata

        def _run(**kw):
            return self._regenerate_state_diagnostics(ct, model_adata=model_adata, verbose=True)

        self._bg.run(
            fn=_run,
            desc=f"State diagnostics ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_state_diagnostics],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Diagnostics done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Diagnostics failed: {e}"),
        )

    def _on_state_composition(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        full_path = self._full_adata_path(ct)
        if not full_path or not full_path.exists():
            QMessageBox.warning(self, "No data", "Run state classification first.")
            return
        out = self._out_dir()
        self._log(f"▶ Generating state composition report for '{ct}'…")
        selected_cols = [item.text() for item in self.list_composition_group_cols.selectedItems()]
        group_x = self.combo_composition_group_x.currentText()
        group_x = None if group_x in ("", "(none)") else group_x
        group_y = self.combo_composition_group_y.currentText()
        group_y = None if group_y in ("", "(none)") else group_y
        group_x_levels_map = None
        if self.chk_composition_group_x_conditions.isChecked():
            if not group_x or not self.group_selector_composition_x.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one Group X level into each group, or untick "
                    "\"Group conditions\" under Group in X.",
                )
                return
            left = self.group_selector_composition_x.left_items()
            right = self.group_selector_composition_x.right_items()
            group_x_levels_map = {lvl: "+".join(left) for lvl in left}
            group_x_levels_map.update({lvl: "+".join(right) for lvl in right})
        group_y_levels_map = None
        if self.chk_composition_group_y_conditions.isChecked():
            if not group_y or not self.group_selector_composition_y.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one Group Y level into each group, or untick "
                    "\"Group conditions\" under Group in Y.",
                )
                return
            left = self.group_selector_composition_y.left_items()
            right = self.group_selector_composition_y.right_items()
            group_y_levels_map = {lvl: "+".join(left) for lvl in left}
            group_y_levels_map.update({lvl: "+".join(right) for lvl in right})
        _raw_md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        md_snapshot = _raw_md.copy() if _raw_md is not None else None

        def _run(**kw):
            import anndata as ad
            from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
            from behav3d.analysis.behavior.state.utils import _resolve_state_paths
            from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
                save_state_composition_report,
            )
            adata = ad.read_h5ad(str(full_path))

            # Inject known metadata grouping columns from metadata CSV into obs (for
            # existing h5ad files that predate Fix 1 in classification.py).
            if md_snapshot is not None and "sample_name" in md_snapshot.columns:
                _known_meta = ["exp_nr", "well"] + [
                    c for c in md_snapshot.columns if c.endswith("_line_condition")
                ]
                cols_to_inject = [
                    c for c in _known_meta
                    if c in md_snapshot.columns and c not in adata.obs.columns
                ]
                if cols_to_inject:
                    meta_map = (
                        md_snapshot[["sample_name"] + cols_to_inject]
                        .drop_duplicates("sample_name")
                        .set_index("sample_name")
                    )
                    for col in cols_to_inject:
                        adata.obs[col] = (
                            adata.obs["sample_name"]
                            .map(meta_map[col])
                            .astype(str)
                            .fillna("(unknown)")
                        )

            from behav3d.analysis.behavior.state.utils import (
                _get_classification_state_colors,
                _get_classification_state_order,
            )
            composition_dir = _resolve_state_paths(out, ct).state_composition_outdir
            composition_dir.mkdir(parents=True, exist_ok=True)
            return save_state_composition_report(
                adata=adata,
                output_pdf_path=composition_dir / "state_composition_report.pdf",
                output_csv_path=composition_dir / "state_composition_report.csv",
                time_col="position_t",
                state_col=FULL_STATE_COL,
                sample_col="sample_name",
                include_pooled_summary=True,
                group_cols=selected_cols,
                group_x=group_x,
                group_y=group_y,
                group_x_levels_map=group_x_levels_map,
                group_y_levels_map=group_y_levels_map,
                state_colors=_get_classification_state_colors(adata, FULL_STATE_COL),
                state_order=_get_classification_state_order(adata, FULL_STATE_COL),
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"State composition report ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_state_composition],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Composition report done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Composition report failed: {e}"),
        )

    def _on_state_transition(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        full_path = self._full_adata_path(ct)
        if not full_path or not full_path.exists():
            QMessageBox.warning(self, "No data", "Run state classification first.")
            return
        out = self._out_dir()
        self._log(f"▶ Generating state transition report for '{ct}'…")

        def _run(**kw):
            import anndata as ad
            from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
            from behav3d.analysis.behavior.state.utils import _resolve_state_paths
            from behav3d.analysis.behavior.state.visualization.plots.state_transitions import (
                save_state_transition_report,
            )
            from behav3d.analysis.behavior.state.utils import (
                _get_classification_state_colors,
                _get_classification_state_order,
            )
            adata = ad.read_h5ad(str(full_path))
            transition_dir = _resolve_state_paths(out, ct).state_transitions_outdir
            transition_dir.mkdir(parents=True, exist_ok=True)
            return save_state_transition_report(
                adata=adata,
                output_dir=transition_dir,
                state_col=FULL_STATE_COL,
                time_col="position_t",
                state_colors=_get_classification_state_colors(adata, FULL_STATE_COL),
                state_order=_get_classification_state_order(adata, FULL_STATE_COL),
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"State transition report ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_state_transition],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Transition report done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Transition report failed: {e}"),
        )

    def _on_condition_comparison(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        full_path = self._full_adata_path(ct)
        if not full_path or not full_path.exists():
            QMessageBox.warning(self, "No data", "Run state classification first.")
            return
        condition_col = self.combo_comparison_condition_col.currentText()
        group_cols = [item.text() for item in self.list_comparison_group_cols.selectedItems()]
        group_x = self.combo_comparison_group_x.currentText()
        group_x = None if group_x in ("", "(none)") else group_x
        if not condition_col:
            QMessageBox.warning(self, "Missing selection", "Select a condition column to compare.")
            return
        condition_groups = None
        if self.chk_comparison_group_conditions.isChecked():
            if not self.group_selector_comparison.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one condition level into each group, or untick "
                    "\"Group conditions\".",
                )
                return
            left = self.group_selector_comparison.left_items()
            right = self.group_selector_comparison.right_items()
            condition_groups = {lvl: "+".join(left) for lvl in left}
            condition_groups.update({lvl: "+".join(right) for lvl in right})
        group_x_levels_map = None
        if self.chk_comparison_group_x_conditions.isChecked():
            if not group_x or not self.group_selector_comparison_x.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one Group X level into each group, or untick "
                    "\"Group conditions\" under Group in X.",
                )
                return
            left = self.group_selector_comparison_x.left_items()
            right = self.group_selector_comparison_x.right_items()
            group_x_levels_map = {lvl: "+".join(left) for lvl in left}
            group_x_levels_map.update({lvl: "+".join(right) for lvl in right})
        out = self._out_dir()
        self._log(f"▶ Generating condition comparison report for '{ct}'…")

        def _run(**kw):
            import anndata as ad
            from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
            from behav3d.analysis.behavior.state.utils import (
                _resolve_state_paths,
                _get_classification_state_colors,
                _get_classification_state_order,
            )
            from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
                save_state_condition_comparison_report,
            )
            from behav3d.analysis.behavior.utils import _sanitize_filename_token
            adata = ad.read_h5ad(str(full_path))
            state_paths = _resolve_state_paths(out, ct)
            comp_dir = state_paths.state_composition_outdir / "behavior_proportions"
            comp_dir.mkdir(parents=True, exist_ok=True)
            cond_token = _sanitize_filename_token(condition_col, fallback="condition")
            out_pdf = comp_dir / f"condition_comparison_{cond_token}.pdf"
            return save_state_condition_comparison_report(
                adata=adata,
                output_pdf_path=out_pdf,
                output_csv_path=out_pdf.with_suffix(".csv"),
                state_col=FULL_STATE_COL,
                sample_col="sample_name",
                condition_col=condition_col,
                group_cols=group_cols or None,
                group_x=group_x,
                group_x_levels_map=group_x_levels_map,
                condition_groups=condition_groups,
                state_colors=_get_classification_state_colors(adata, FULL_STATE_COL),
                state_order=_get_classification_state_order(adata, FULL_STATE_COL),
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Condition comparison report ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_condition_comparison],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Condition comparison report done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Condition comparison report failed: {e}"),
        )

    # ── Backprojection ───────────────────────────────────────────────────

    def _on_bp_layer_display_changed(self, event=None):
        self._display_save_timer.start()

    def _persist_bp_viewer_display(self):
        _bp_save_channel_display(self.viewer, self.metadata_loader, self._out_dir)

    def _current_viewer_frame(self) -> int:
        if self.viewer is None:
            return 0
        try:
            return int(self.viewer.dims.current_step[0])
        except Exception:
            return 0

    def _teardown_state_bp_preview(self):
        if self._state_bp_dims_callback is not None and self.viewer is not None:
            try:
                self.viewer.dims.events.current_step.disconnect(self._state_bp_dims_callback)
            except Exception:
                pass
        unregister_preview_dims_listener(self.viewer, self)
        self._state_bp_dims_callback = None
        self._state_bp_preview = None

    def _connect_state_bp_dims_listener(self):
        if self._state_bp_dims_callback is not None and self.viewer is not None:
            try:
                self.viewer.dims.events.current_step.disconnect(self._state_bp_dims_callback)
            except Exception:
                pass
            unregister_preview_dims_listener(self.viewer, self)
            self._state_bp_dims_callback = None
        if self.viewer is None:
            return

        def _on_step(*_):
            self._refresh_state_bp_layer()

        try:
            self.viewer.dims.events.current_step.connect(_on_step)
            self._state_bp_dims_callback = _on_step
            register_preview_dims_listener(self.viewer, self, _on_step)
        except Exception:
            self._state_bp_dims_callback = None

    def _refresh_state_bp_layer(self):
        """Recompute the state-class overlay for whichever timepoint napari's
        time slider is currently on. Mirrors the Feature Backprojection tab's
        ``_refresh_feature_layer`` (see ``behav3d.napari._feature_backprojection``)
        so opening the preview never requires writing a full per-sample zarr."""
        from behav3d.analysis.behavior.state.visualization.backprojection import (
            backproject_state_at_timepoint,
        )
        from behav3d.io.images import load_image_timepoint

        preview = self._state_bp_preview
        if preview is None or self.viewer is None:
            return

        t = self._current_viewer_frame()
        try:
            labels_frame = np.asarray(load_image_timepoint(preview["tracked_path"], t))
        except Exception as exc:
            self._log(f"❌ Could not read frame {t}: {exc}")
            return
        if labels_frame.ndim == 4:
            labels_frame = labels_frame[0]

        mapped, _ids_with_value = backproject_state_at_timepoint(
            labels_frame=labels_frame,
            state_code_lookup=preview["code_lookup"],
            time_index=t,
        )

        name = preview["layer_name"]
        try:
            layer = self.viewer.layers[name]
            layer.data = mapped
            layer.refresh()
        except (KeyError, ValueError):
            layer = self.viewer.add_labels(mapped, name=name, opacity=preview["opacity"])
            from behav3d.analysis.behavior.state.visualization.backprojection import (
                _apply_state_code_colors_to_layer,
            )
            _apply_state_code_colors_to_layer(layer, preview.get("code_colors", {}))

    def _on_show_state_bp(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        color_by = self.combo_state_color_by.currentText()
        if color_by == "raw_hmm_state":
            state_path = self._model_adata_path(ct)
            err_hint = "Model adata not found. Run State Classification (Step 1) first."
        else:
            state_path = self._full_adata_path(ct)
            err_hint = f"State adata not found:\n{state_path}\n\nRun State Classification first."
        if not state_path or not state_path.exists():
            QMessageBox.warning(self, "No state adata", err_hint)
            return
        opacity = self.spin_state_opacity.value() / 100.0
        self._log(f"▶ Loading state backprojection for '{ct}' / sample '{sample}'…")
        try:
            import scanpy as sc
            from behav3d.analysis.behavior.state.visualization.backprojection import (
                _resolve_raw_image_path,
                _resolve_tracked_image_path,
                _build_code_map,
                _build_state_code_color_map,
                _add_mapping_dock_widget,
                _build_state_mapping_text,
                _align_labels_to_raw_shape_for_view,
                _validate_required_obs_columns,
                prepare_state_code_lookup,
            )
            from behav3d.analysis.behavior.state.utils import (
                _get_classification_state_colors,
                _get_classification_state_order,
                _normalize_label_color_map,
            )
            from behav3d.analysis.backprojection import filter_track_image_to_ids
            from behav3d.io.images import load_image
            out_dir = self._out_dir()
            if not out_dir:
                raise ValueError("No output directory set.")
            adata = sc.read_h5ad(str(state_path))
            if self.chk_show_state_trajectories.isChecked():
                try:
                    from behav3d.analysis.behavior.track.visualization.plots.exemplar_coordinate_utils import (
                        ensure_exemplar_coordinate_columns,
                    )
                    ensure_exemplar_coordinate_columns(
                        adata, output_dir=out_dir, cell_type=ct, require_pixel_for_video=True,
                    )
                except Exception as exc:
                    self._log(f"⚠️ Could not prepare trajectory positions: {exc}")
            resolved_col = "hmm_intrinsic_behavioral_state_raw" if color_by == "raw_hmm_state" else color_by
            state_col = resolved_col if (resolved_col and resolved_col in adata.obs.columns) else "full_behavioral_cluster"
            obs_samples = adata.obs["sample_name"].astype(str)
            sample_name = sample if sample else obs_samples.iloc[0]
            sample_adata = adata[obs_samples == str(sample_name)]
            if sample_adata.n_obs == 0:
                raise ValueError(f"No rows found for sample '{sample_name}' in state adata.")
            _validate_required_obs_columns(sample_adata.obs, required_cols=["TrackID", "position_t", state_col])

            state_order = _get_classification_state_order(sample_adata, state_col)
            code_map = _build_code_map(sample_adata.obs, state_col=state_col, state_order=state_order)
            if len(code_map) == 0:
                raise ValueError(f"'{state_col}' has no non-empty labels for sample '{sample_name}'.")
            state_colors = _normalize_label_color_map(
                code_map.keys(), colors=_get_classification_state_colors(sample_adata, state_col)
            )
            code_colors = _build_state_code_color_map(code_map, state_colors=state_colors)
            label_map = {str(code): str(label) for label, code in code_map.items()}
            state_code_lookup = prepare_state_code_lookup(
                sample_adata.obs, state_col=state_col, code_map=code_map,
            )

            gui_metadata_csv_path = resolve_metadata_csv_path(self.metadata_loader)
            raw_path = _resolve_raw_image_path(
                out_dir, sample_name, verbose=False, metadata_csv_path=gui_metadata_csv_path
            )
            if raw_path is None or not Path(raw_path).exists():
                raise FileNotFoundError(f"Raw image not found for sample '{sample_name}'.")
            tracked_path = _resolve_tracked_image_path(
                out_dir, sample_name, ct, verbose=False, metadata_csv_path=gui_metadata_csv_path
            )
            if tracked_path is None or not Path(tracked_path).exists():
                raise FileNotFoundError(
                    f"Tracked image not found for sample '{sample_name}', cell_type '{ct}'."
                )
            raw_img = load_image(raw_path)
            tracked_img = load_image(tracked_path)
            tracked_view = _align_labels_to_raw_shape_for_view(tracked_img, raw_img, "TrackID", verbose=False)
            keep_ids = state_code_lookup["TrackID"].unique()
            tracked_view = filter_track_image_to_ids(tracked_view, keep_ids)

            # Full viewer reset — not just the specific layer names this
            # preview itself uses. Whatever the Visualization tab (raw
            # channels, Segments, Tracked Segments, Tracks) or another
            # backprojection preview (Track Classification, Feature
            # Backprojection) had loaded should not linger under this
            # preview's overlay, and any of their still-connected dims
            # listeners must be dropped before we clear layers so none of
            # them can fire reentrantly mid-clear.
            self._teardown_state_bp_preview()
            stop_dim_playback(self.viewer)
            disconnect_all_preview_dims_listeners(self.viewer)
            clear_viewer_layers(self.viewer)
            saved_channels = (
                getattr(self.metadata_loader, "behav3d_parameters", {})
                .get("viewer_display", {})
                .get("channels", {})
            )
            try:
                md = self.metadata_loader.metadata
                row = md[md["sample_name"] == sample_name].iloc[0]
                dim_order = str(row.get("dimension_order", "TCZYX")).strip() or "TCZYX"
            except Exception:
                dim_order = "TCZYX"
            ch_names = _bp_add_raw_channels(self.viewer, raw_img, sample_name, saved_channels, dim_order)
            for lname in ch_names:
                try:
                    layer = self.viewer.layers[lname]
                    layer.events.contrast_limits.connect(self._on_bp_layer_display_changed)
                    layer.events.colormap.connect(self._on_bp_layer_display_changed)
                except (KeyError, IndexError):
                    pass
            self.viewer.add_labels(tracked_view, name="filtered TrackID", visible=False, opacity=opacity)

            self._state_bp_preview = {
                "tracked_path": Path(tracked_path),
                "code_lookup": state_code_lookup,
                "layer_name": color_by,
                "opacity": opacity,
                "code_colors": code_colors,
            }
            self._refresh_state_bp_layer()
            self._connect_state_bp_dims_listener()

            if self.chk_show_state_trajectories.isChecked():
                try:
                    from behav3d.analysis.behavior.track.visualization.backprojection import (
                        add_track_cluster_trajectory_layers,
                    )
                    from behav3d.analysis.behavior.state.visualization.backprojection import (
                        prepare_state_trajectory_data,
                    )
                    trajectory_data = prepare_state_trajectory_data(
                        sample_adata.obs, state_col=state_col,
                    )
                    add_track_cluster_trajectory_layers(
                        self.viewer,
                        trajectory_data=trajectory_data,
                        code_colors=code_colors,
                        label_map=label_map,
                        output_col=state_col,
                        tail_length=int(tracked_img.shape[0]),
                    )
                except Exception as exc:
                    self._log(f"⚠️ Could not add state trajectory layers: {exc}")

            mapping_text = _build_state_mapping_text(label_map, code_colors)
            _existing_dock = getattr(self, "_state_mapping_dock", None)
            if _existing_dock is not None:
                try:
                    self.viewer.window.remove_dock_widget(_existing_dock)
                except Exception:
                    pass
            self._state_mapping_dock = _add_mapping_dock_widget(
                self.viewer,
                mapping_text=mapping_text,
                label_map=label_map,
                code_colors=code_colors,
                title="State Class Mapping",
            )
            self._log("✅ State backprojection loaded (current timepoint; updates as you scrub).")
        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Backprojection failed: {e}")

    def _on_export_state_bp(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        state_path = self._full_adata_path(ct)
        if not state_path or not state_path.exists():
            QMessageBox.warning(self, "No state adata", "Run State Classification first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        color_by = self.combo_state_color_by.currentText()
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Exporting state backprojection for '{ct}'…")

        def _run(**kw):
            import scanpy as sc
            from behav3d.analysis.behavior.state.visualization.backprojection import (
                export_behavioral_state_backprojection_zarrs,
            )
            adata = sc.read_h5ad(str(state_path))
            state_col = color_by if color_by else "full_behavioral_cluster"
            sample_name = sample if sample else None
            export_adata = (
                adata[adata.obs["sample_name"].astype(str) == str(sample_name)]
                if sample_name is not None else adata
            )
            logger("▶ Writing state backprojection zarrs…")
            return export_behavioral_state_backprojection_zarrs(
                adata=export_adata,
                output_dir=out,
                cell_type=ct,
                state_col=state_col,
                enforce_time_coverage=False,
                n_workers=1,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Exporting state backprojection ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_export_state_bp, self.btn_show_state_bp],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: _log_backprojection_manifest(self._log, r),
            on_failed=lambda e: self._log(f"❌ Export failed: {e}"),
        )

    # ── View helpers ─────────────────────────────────────────────────────

    def _on_view(self, kind: str, btn: Optional[QPushButton] = None):
        ct = self._cell_type()
        if not ct:
            return
        candidates = self._get_view_candidates(kind, ct)
        existing = [(lbl, p) for lbl, p in candidates if p.exists()]
        if not existing:
            return
        if len(existing) == 1:
            self._open_pdf(existing[0][1])
        else:
            menu = QMenu(self)
            for lbl, path in existing:
                act = menu.addAction(f"👁  {lbl}")
                act.triggered.connect(lambda _=False, _p=path: self._open_pdf(_p))
            if btn is not None:
                menu.exec_(btn.mapToGlobal(btn.rect().bottomLeft()))
            else:
                menu.exec_()

    def _get_view_candidates(self, kind: str, ct: str):
        out = self._out_dir()
        if not out:
            return []
        if kind == "state_qc":
            return self._state_run_qc_pdfs(ct)
        if kind == "state_diagnostics":
            qc_dir = self._state_diagnostics_qc_dir(ct)
            if not qc_dir or not qc_dir.exists():
                return []
            return [(f.stem, f) for f in sorted(qc_dir.glob("*.pdf"))]
        if kind == "state_composition":
            p = self._report_path(ct, "state_composition_report")
            return [(f"Composition report ({ct})", p)] if p else []
        if kind == "state_transition":
            p = self._report_path(ct, "state_transition_report")
            return [(f"Transition report ({ct})", p)] if p else []
        if kind == "state_condition_comparison":
            comp_dir = out / "analysis" / ct / "behavioral_states" / "state_composition" / "behavior_proportions"
            if not comp_dir.exists():
                return []
            return [
                (p.stem.replace("condition_comparison_", ""), p)
                for p in sorted(comp_dir.glob("condition_comparison_*.pdf"))
            ]
        return []

    def _open_pdf(self, path: Path):
        try:
            panel = self._results_panel()
            dpi = int(panel.spin_dpi.value()) if panel else 150
            reuse = getattr(panel, "chk_reuse", None)
            reuse = reuse.isChecked() if reuse else True
            open_pdf_in_napari(path, dpi=dpi, reuse=reuse)
        except Exception as e:
            traceback.print_exc()
            self._log(f"Could not open in napari: {e}")

    def _results_panel(self):
        node = self.parent()
        while node and not hasattr(node, "results_panel"):
            node = node.parent()
        return getattr(node, "results_panel", None) if node else None

    def _notify_results(self):
        panel = self._results_panel()
        if panel:
            try:
                panel.refresh()
            except Exception:
                pass


# ═══════════════════════════════════════════════════════════════════════════
# Track Classification sub-tab
# ═══════════════════════════════════════════════════════════════════════════

class TrackClassificationSubTab(QWidget):
    """Steps 1-5 for DTW-based track trajectory classification + backprojection."""

    def __init__(self, viewer=None, metadata_loader=None, cell_type_getter=None,
                 parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._get_cell_type = cell_type_getter
        self._track_adata = None
        self._track_adata_load_error: Optional[str] = None
        self._bg = BackgroundOperation(self)
        self._preload_bg = BackgroundOperation(self)
        # Current-timepoint-only track-cluster backprojection preview (see
        # `_refresh_track_bp_layer`): recomputed on every dims scrub instead
        # of writing a full per-sample zarr up front.
        self._track_bp_preview: Optional[dict] = None
        self._track_bp_dims_callback = None

        self._display_save_timer = QTimer(self)
        self._display_save_timer.setSingleShot(True)
        self._display_save_timer.setInterval(1000)
        self._display_save_timer.timeout.connect(self._persist_bp_viewer_display)

        # Queue buttons (wired by _widget.py)
        self.btn_queue_track_cluster = _make_queue_btn()
        self.btn_queue_train_track = _make_queue_btn()
        self.btn_queue_apply_track = _make_queue_btn()

        self._init_ui()

    # ── UI ───────────────────────────────────────────────────────────────

    def _init_ui(self):
        from behav3d.napari._guided import make_back_header

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # Pinned above the scroll area (not inside it) so it stays visible
        # regardless of how far the plotting page's content scrolls. Only
        # shown while the plotting page (outer page 1) is active.
        self._plotting_header, self._pipeline_title = make_back_header(
            on_back=self._plotting_back, label="‹ Back"
        )
        self._plotting_header.hide()
        outer.addWidget(self._plotting_header)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        content = QWidget()
        content_lay = QVBoxLayout(content)
        content_lay.setContentsMargins(0, 0, 0, 0)
        content_lay.setSpacing(0)
        scroll.setWidget(content)

        # Outer per-subtab stack: page 0 = clustering settings (Steps 1/2/4 +
        # the "Generate analysis and plots" trigger), page 1 = the fully
        # isolated plotting page, reached only via that trigger.
        self._subtab_stack = QStackedWidget()
        content_lay.addWidget(self._subtab_stack)

        settings_page = QWidget()
        lay = QVBoxLayout(settings_page)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        self._subtab_stack.addWidget(settings_page)  # outer page 0

        # ── Warning label (prerequisite check) ──────────────────────────
        self.warning_label = QLabel("")
        self.warning_label.setWordWrap(True)
        self.warning_label.setStyleSheet(
            "QLabel { background: #3d2200; color: #ffaa44; border-radius: 4px; "
            "padding: 6px 8px; font-size: 11px; }"
        )
        self.warning_label.hide()
        lay.addWidget(self.warning_label)

        # ── Apply Pretrained Checkbox ────────────────────────────────────
        self.chk_apply_pretrained = QCheckBox("Apply pretrained trajectory classifier")
        self.chk_apply_pretrained.setChecked(False)
        lay.addWidget(self.chk_apply_pretrained)

        # ── Group: Apply Pretrained ──────────────────────────────────────
        self.grp_apply_pretrained = QGroupBox("Apply pretrained trajectory classifier")
        ap_lay = QVBoxLayout(self.grp_apply_pretrained)
        ap_lay.setSpacing(5)

        ap_lay.addWidget(_make_info_label("Load an existing trained classifier (.pkl) and behavioral states file (.h5ad)"
                                           " to classify trajectories without running Steps 1–3."))

        self.le_pretrained_clf_path = QLineEdit()
        self.le_pretrained_clf_path.setPlaceholderText("Path to trained classifier .pkl")
        self.btn_browse_pretrained_clf = QPushButton("Browse…")
        ap_lay.addLayout(_make_browse_row(
            self.le_pretrained_clf_path, self.btn_browse_pretrained_clf,
            "Classifier path",
            "Path to the trained Random Forest classifier .pkl file."
        ))

        self.le_pretrained_states_path = QLineEdit()
        self.le_pretrained_states_path.setPlaceholderText("Path to behavioral states .h5ad")
        self.btn_browse_pretrained_states = QPushButton("Browse…")
        ap_lay.addLayout(_make_browse_row(
            self.le_pretrained_states_path, self.btn_browse_pretrained_states,
            "States h5ad path",
            "Path to the behavioral states .h5ad file (output of State Classification). "
            "Auto-fills if State Classification output exists for the selected cell type."
        ))

        self.btn_run_apply_pretrained = QPushButton("▶ Apply Pretrained Classifier")
        _style_primary(self.btn_run_apply_pretrained)
        ap_lay.addWidget(self.btn_run_apply_pretrained)
        self.grp_apply_pretrained.hide()
        lay.addWidget(self.grp_apply_pretrained)

        # ── Step 1: Track Clustering ─────────────────────────────────────
        self.grp1 = QGroupBox("Step 1 — Track Clustering")
        g1 = QVBoxLayout(self.grp1)
        g1.setSpacing(4)

        # "Use original" at top (shown only when original mode is active)
        self.chk_use_original_top = QCheckBox(
            "Use original feature-based BEHAV3D DTW Clustering"
        )
        self.chk_use_original_top.setChecked(False)
        self.chk_use_original_top.hide()
        g1.addWidget(self.chk_use_original_top)

        # Shared: Trajectory size + N clusters (always visible)
        basic_form = QFormLayout()
        basic_form.setSpacing(3)
        self.spin_traj_size = QSpinBox()
        self.spin_traj_size.setRange(1, 9999)
        self.spin_traj_size.setValue(100)
        self.spin_traj_size.setMaximumWidth(120)
        basic_form.addRow("Trajectory size (timepoints):", make_help_row(
            self.spin_traj_size, "Trajectory size",
            "Number of timepoints to include per track. With the default 'last' trim mode "
            "this keeps each track's final N timepoints (dropping its earlier ones); "
            "tracks shorter than this are discarded. Match the value used in your notebook. "
            "The line below shows the equivalent real time, computed from this sample's "
            "time_interval/time_unit metadata."
        ))
        self.lbl_traj_size_time = _make_timepoint_time_label()
        basic_form.addRow("", self.lbl_traj_size_time)
        self.spin_traj_size.valueChanged.connect(self._update_timepoint_time_labels)
        self.spin_n_clusters = QSpinBox()
        self.spin_n_clusters.setRange(2, 200)
        self.spin_n_clusters.setValue(6)
        self.spin_n_clusters.setMaximumWidth(90)
        basic_form.addRow("N clusters:", make_help_row(
            self.spin_n_clusters, "N clusters",
            "Number of trajectory clusters to create via hierarchical clustering of "
            "DTW distances."
        ))
        g1.addLayout(basic_form)

        # Trajectory basis + clustering method — the two main choices, kept always
        # visible (not buried in Advanced Configuration) since they each gate a
        # different set of downstream parameter fields.
        self._basis_method_frame = QFrame()
        basis_method_form = QFormLayout(self._basis_method_frame)
        basis_method_form.setSpacing(3)

        self.combo_trajectory_basis = QComboBox()
        self.combo_trajectory_basis.addItems(["dtw", "bouts"])
        self.combo_trajectory_basis.setMaximumWidth(130)
        basis_method_form.addRow("Basis:", make_help_row(
            self.combo_trajectory_basis, "Trajectory basis",
            "'dtw' clusters tracks by dynamic time warping distance over their raw "
            "per-timepoint state sequences. 'bouts' instead describes each track with "
            "bout/proportion features (fraction of time per state, bout counts/lengths, "
            "state-transition probabilities, n-grams) and clusters those feature vectors."
        ))

        self.combo_clustering_method = QComboBox()
        self.combo_clustering_method.addItems(["agglomerative", "leiden"])
        self.combo_clustering_method.setMaximumWidth(130)
        basis_method_form.addRow("Method:", make_help_row(
            self.combo_clustering_method, "Clustering method",
            "Algorithm used to group tracks from the chosen basis above. "
            "'agglomerative' uses a fixed cluster count (N clusters, Linkage below). "
            "'leiden' finds density-based communities instead, so the number of "
            "clusters is emergent (set via Leiden resolution below)."
        ))
        g1.addWidget(self._basis_method_frame)

        # UMAP parameters (only in original mode)
        self._umap_frame = QFrame()
        self._umap_frame.setVisible(False)
        umap_form = QFormLayout(self._umap_frame)
        umap_form.setSpacing(3)
        self.spin_umap_neighbors = QSpinBox()
        self.spin_umap_neighbors.setRange(2, 200)
        self.spin_umap_neighbors.setValue(15)
        self.spin_umap_neighbors.setMaximumWidth(90)
        umap_form.addRow("UMAP n_neighbors:", make_help_row(
            self.spin_umap_neighbors, "UMAP n_neighbors",
            "Number of neighbours used to build the UMAP graph. Higher values capture "
            "more global structure."
        ))
        self.spin_umap_min_dist = QDoubleSpinBox()
        self.spin_umap_min_dist.setRange(0.001, 1.0)
        self.spin_umap_min_dist.setSingleStep(0.05)
        self.spin_umap_min_dist.setDecimals(3)
        self.spin_umap_min_dist.setValue(0.1)
        self.spin_umap_min_dist.setMaximumWidth(90)
        umap_form.addRow("UMAP min_dist:", make_help_row(
            self.spin_umap_min_dist, "UMAP min_dist",
            "Minimum distance between points in the UMAP embedding. Smaller values "
            "produce tighter clusters."
        ))
        g1.addWidget(self._umap_frame)

        # Advanced Configuration (hidden in original mode, contains "Use original" checkbox)
        self.adv1 = CollapsibleSection("⚙ Advanced Configuration", expanded=False)

        # Trim mode + divide-long-tracks apply to both bases (DTW and Bouts/proportions
        # both truncate/split tracks the same way before clustering), so they stay in a
        # plain shared form rather than a basis-toggled frame.
        trim_form = QFormLayout()
        trim_form.setSpacing(3)
        self.combo_trim = QComboBox()
        self.combo_trim.addItems(["last", "first"])
        self.combo_trim.setMaximumWidth(130)
        trim_form.addRow("Trim mode:", make_help_row(
            self.combo_trim, "Trim mode",
            "How to trim each track to Trajectory size: "
            "'last' keeps each track's final N timepoints (removes leading/early ones); "
            "'first' keeps each track's first N timepoints (removes trailing/late ones)."
        ))
        self.chk_split_long_tracks = QCheckBox("Divide long tracks")
        self.chk_split_long_tracks.setChecked(False)
        trim_form.addRow("", _make_chk_help_row(
            self.chk_split_long_tracks, "Divide long tracks",
            "Split tracks longer than Trajectory size into non-overlapping full-length "
            "analysis windows. Leftover timepoints are discarded. Original TrackID values "
            "are preserved for backprojection."
        ))
        self.adv1.addLayout(trim_form)

        # Linkage (DTW) and Bouts linkage each get their own self-contained frame —
        # not a shared QFormLayout row — so hiding one via setVisible() also hides its
        # label and help button together instead of leaving them stranded.
        self._dtw_linkage_frame = QFrame()
        dtw_linkage_form = QFormLayout(self._dtw_linkage_frame)
        dtw_linkage_form.setSpacing(3)
        self.combo_linkage = QComboBox()
        self.combo_linkage.addItems(["average", "complete", "single"])
        self.combo_linkage.setMaximumWidth(130)
        dtw_linkage_form.addRow("Linkage:", make_help_row(
            self.combo_linkage, "Linkage",
            "Agglomerative clustering linkage method. 'average' is the most commonly used; "
            "'complete' produces more compact clusters. Only used with the DTW basis and "
            "'agglomerative' Method."
        ))
        self.adv1.addWidget(self._dtw_linkage_frame)

        self._bouts_linkage_frame = QFrame()
        bouts_linkage_form = QFormLayout(self._bouts_linkage_frame)
        bouts_linkage_form.setSpacing(3)
        self.combo_bouts_linkage = QComboBox()
        self.combo_bouts_linkage.addItems(["ward", "complete", "average", "single"])
        self.combo_bouts_linkage.setMaximumWidth(130)
        bouts_linkage_form.addRow("Bouts linkage:", make_help_row(
            self.combo_bouts_linkage, "Bouts linkage",
            "Agglomerative clustering linkage method for the Bouts/proportions basis. "
            "'ward' (default) minimizes within-cluster variance and is the standard choice "
            "for feature vectors. Only used with the Bouts basis and 'agglomerative' Method."
        ))
        self.adv1.addWidget(self._bouts_linkage_frame)

        self._leiden_frame = QFrame()
        leiden_form = QFormLayout(self._leiden_frame)
        leiden_form.setSpacing(3)
        self.spin_leiden_neighbors = QSpinBox()
        self.spin_leiden_neighbors.setRange(2, 200)
        self.spin_leiden_neighbors.setValue(15)
        self.spin_leiden_neighbors.setMaximumWidth(90)
        leiden_form.addRow("Leiden neighbors:", make_help_row(
            self.spin_leiden_neighbors, "Leiden neighbors",
            "Number of nearest neighbours used to build the graph Leiden clusters on. "
            "Higher values smooth over noise but can merge distinct groups."
        ))
        self.spin_leiden_resolution = QDoubleSpinBox()
        self.spin_leiden_resolution.setRange(0.01, 20.0)
        self.spin_leiden_resolution.setSingleStep(0.1)
        self.spin_leiden_resolution.setDecimals(2)
        self.spin_leiden_resolution.setValue(1.0)
        self.spin_leiden_resolution.setMaximumWidth(90)
        leiden_form.addRow("Leiden resolution:", make_help_row(
            self.spin_leiden_resolution, "Leiden resolution",
            "Controls how many Leiden clusters are found: higher values produce more, "
            "smaller clusters; lower values produce fewer, larger ones."
        ))
        self.adv1.addWidget(self._leiden_frame)

        # Bouts/proportions-only feature toggles (which track-describing features feed
        # the clustering) + PCA + optional exemplar-PDF generation during Run, since the
        # bouts pipeline has no separable "diagnostics only" step like the DTW one does.
        self._bouts_frame = QFrame()
        bouts_form = QFormLayout(self._bouts_frame)
        bouts_form.setSpacing(3)

        self.chk_bouts_use_fractions = QCheckBox("State fractions")
        self.chk_bouts_use_fractions.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_use_fractions, "State fractions",
            "Include each state's fraction of total track time as a clustering feature."
        ))
        self.chk_bouts_use_bout_stats = QCheckBox("Bout stats")
        self.chk_bouts_use_bout_stats.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_use_bout_stats, "Bout stats",
            "Include per-state bout count, mean length, and max length as clustering features."
        ))
        self.chk_bouts_use_transitions = QCheckBox("Transitions")
        self.chk_bouts_use_transitions.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_use_transitions, "Transitions",
            "Include state-to-state transition probabilities as clustering features."
        ))
        self.chk_bouts_use_ngrams = QCheckBox("Bigrams/trigrams")
        self.chk_bouts_use_ngrams.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_use_ngrams, "Bigrams/trigrams",
            "Include state bigram and trigram counts as clustering features."
        ))
        self.chk_bouts_do_pca = QCheckBox("PCA before clustering")
        self.chk_bouts_do_pca.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_do_pca, "PCA before clustering",
            "Reduce the feature matrix with PCA (95% variance retained) before clustering."
        ))
        self.chk_bouts_use_clr = QCheckBox("Log-ratio transform proportions (CLR)")
        self.chk_bouts_use_clr.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_use_clr, "Log-ratio transform proportions (CLR)",
            "State fractions, bout-count shares, and each transition row sum to 1 by "
            "construction (compositional data). CLR opens them up before PCA/clustering so "
            "Euclidean distance doesn't manufacture spurious negative correlations between "
            "states purely from that sum-to-1 constraint."
        ))
        self.chk_bouts_log_bout_length = QCheckBox("Log-transform bout lengths")
        self.chk_bouts_log_bout_length.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_log_bout_length, "Log-transform bout lengths",
            "Bout/pause durations are typically heavy-tailed; log1p compresses that tail so "
            "a few very long bouts don't dominate Euclidean distance."
        ))
        self.chk_bouts_block_scaling = QCheckBox("Balance feature blocks (MFA)")
        self.chk_bouts_block_scaling.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_block_scaling, "Balance feature blocks (MFA)",
            "Fractions/bout-stats/transitions/n-grams have very different column counts. "
            "Without balancing, PCA gives more weight to whichever block simply has more "
            "columns. This divides each block by the leading singular value of its own PCA "
            "(Multiple Factor Analysis) so every block contributes comparable spread."
        ))
        self.chk_bouts_drop_redundant = QCheckBox("Drop redundant/constant features")
        self.chk_bouts_drop_redundant.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_drop_redundant, "Drop redundant/constant features",
            "Drops near-constant features and one of any pair of highly correlated features "
            "(threshold 0.95) before clustering -- catches exact collinearities such as K "
            "fractions that must sum to 1."
        ))
        self.chk_bouts_plot_exemplars = QCheckBox("Also generate exemplar PDFs")
        self.chk_bouts_plot_exemplars.setChecked(True)
        bouts_form.addRow("", _make_chk_help_row(
            self.chk_bouts_plot_exemplars, "Also generate exemplar PDFs",
            "Generate exemplar-track PDFs during Run, matching the DTW basis's automatic "
            "exemplar overview. Bouts diagnostics/exemplar PDFs are produced inline here; "
            "the separate Diagnostics/Exemplar PDF buttons below are DTW-only."
        ))
        self.adv1.addWidget(self._bouts_frame)

        # Parallel computation / save-distance-matrix are DTW-distance-matrix-specific
        # (bouts clusters feature vectors directly, no pairwise distance matrix exists),
        # so group them in their own frame for clean hide/show.
        self._dtw_technical_frame = QFrame()
        dtw_technical_lay = QVBoxLayout(self._dtw_technical_frame)
        dtw_technical_lay.setContentsMargins(0, 0, 0, 0)
        dtw_technical_lay.setSpacing(3)

        self.chk_parallel = QCheckBox("Parallel computation")
        self.chk_parallel.setChecked(True)
        dtw_technical_lay.addLayout(_make_chk_help_row(
            self.chk_parallel, "Parallel computation",
            "Use parallel computing for DTW distance matrix calculation. "
            "Faster on multi-core machines."
        ))

        self.chk_save_dist = QCheckBox("Save distance matrix CSV")
        self.chk_save_dist.setChecked(False)
        dtw_technical_lay.addLayout(_make_chk_help_row(
            self.chk_save_dist, "Save distance matrix",
            "Save the full DTW pairwise distance matrix to a CSV file. "
            "Can be large for many tracks."
        ))
        self.adv1.addWidget(self._dtw_technical_frame)

        self.spin_seed = QSpinBox()
        self.spin_seed.setRange(0, 99999)
        self.spin_seed.setValue(123)
        self.spin_seed.setMaximumWidth(90)
        seed_form = QFormLayout()
        seed_form.setSpacing(3)
        seed_form.addRow("Random seed:", make_help_row(
            self.spin_seed, "Random seed",
            "Seed for reproducibility of clustering and UMAP."
        ))
        self.adv1.addLayout(seed_form)

        # Separator + "Use original" checkbox at the bottom of Advanced Config
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet("color: #444;")
        self.adv1.addWidget(sep)

        self.chk_use_original = QCheckBox(
            "Use original feature-based BEHAV3D DTW Clustering"
        )
        self.chk_use_original.setChecked(False)
        self.adv1.addLayout(_make_chk_help_row(
            self.chk_use_original, "Use original BEHAV3D DTW",
            "Switch to the original feature-based BEHAV3D DTW pipeline (run_tcell_analysis). "
            "When checked, the checkbox moves to the top of Step 1, Advanced Configuration "
            "is hidden, and UMAP parameters appear."
        ))

        g1.addWidget(self.adv1)

        # Run button (text changes with mode)
        run_row = QHBoxLayout()
        self.btn_run_track = QPushButton("▶ Run Track Clustering")
        _style_primary(self.btn_run_track)
        run_row.addWidget(self.btn_run_track, stretch=1)
        run_row.addWidget(self.btn_queue_track_cluster)
        g1.addLayout(run_row)

        lay.addWidget(self.grp1)

        # ── Step 2: Rename ───────────────────────────────────────────────
        self.grp2 = QGroupBox("Step 2 — Rename Track Clusters")
        g2 = QVBoxLayout(self.grp2)
        g2.setSpacing(4)

        self.rename_track_status = QLabel("ℹ Run clustering first to enable renaming.")
        self.rename_track_status.setStyleSheet("color: #999; font-size: 11px;")
        self.rename_track_status.setWordWrap(True)
        g2.addWidget(self.rename_track_status)

        self.btn_rename_track = QPushButton("✏  Rename Track Clusters")
        _style_rename(self.btn_rename_track)
        self.btn_rename_track.setEnabled(False)
        g2.addWidget(self.btn_rename_track)
        lay.addWidget(self.grp2)

        # ── Train Track Classifier (last step) ───────────────────────────
        self.grp3 = CollapsibleSection("Train track classifier", expanded=False)
        self.grp3._toggle.setStyleSheet(
            "QToolButton { border: 2px solid #28a745; border-radius: 5px; font-weight: bold; "
            "color: #28a745; font-size: 12px; background-color: rgba(40, 167, 69, 0.06); "
            "padding: 4px 8px; text-align: left; }"
            "QToolButton:hover { background-color: rgba(40, 167, 69, 0.12); }"
            "QToolButton:checked { border-bottom-left-radius: 0px; border-bottom-right-radius: 0px; }"
        )
        self.grp3._content.setStyleSheet(
            "QFrame { border: 2px solid #28a745; border-top: none; "
            "border-bottom-left-radius: 5px; border-bottom-right-radius: 5px; "
            "background-color: rgba(40, 167, 69, 0.06); padding: 4px; }"
        )
        g3 = self.grp3.contentLayout()
        g3.setSpacing(6)

        # Sub-section A: Train RF Classifier
        grp_train_rf = QGroupBox("Train RF Classifier on Named Clusters")
        g_train_rf = QVBoxLayout(grp_train_rf)
        g_train_rf.setSpacing(4)

        rf_form = QFormLayout()
        rf_form.setSpacing(3)

        self.spin_track_n_est = QSpinBox()
        self.spin_track_n_est.setRange(10, 2000)
        self.spin_track_n_est.setValue(300)
        self.spin_track_n_est.setMaximumWidth(90)
        rf_form.addRow("Number of trees:", make_help_row(
            self.spin_track_n_est, "Number of trees",
            "Number of decision trees in the Random Forest classifier. "
            "More trees = more stable predictions but slower training. Default: 300."
        ))

        self.spin_track_max_depth = QSpinBox()
        self.spin_track_max_depth.setRange(0, 100)
        self.spin_track_max_depth.setValue(0)
        self.spin_track_max_depth.setSpecialValueText("Unlimited (0)")
        self.spin_track_max_depth.setMaximumWidth(130)
        rf_form.addRow("Max depth:", make_help_row(
            self.spin_track_max_depth, "Max tree depth",
            "Maximum depth of each decision tree (0 = unlimited). Deeper trees fit more "
            "complex patterns but may overfit on small datasets."
        ))

        self.spin_track_test_pct = QDoubleSpinBox()
        self.spin_track_test_pct.setRange(0.0, 50.0)
        self.spin_track_test_pct.setSingleStep(5.0)
        self.spin_track_test_pct.setDecimals(1)
        self.spin_track_test_pct.setValue(20.0)
        self.spin_track_test_pct.setSuffix(" %")
        self.spin_track_test_pct.setMaximumWidth(110)
        rf_form.addRow("Test holdout %:", make_help_row(
            self.spin_track_test_pct, "Test holdout %",
            "Percentage of data held out for classifier validation (0 = OOB-only mode, "
            "no holdout). Typical value: 20%."
        ))
        g_train_rf.addLayout(rf_form)

        # Advanced RF parameters
        adv3 = CollapsibleSection("⚙ Advanced Parameters", expanded=False)
        adv3_form = QFormLayout()
        adv3_form.setSpacing(3)

        self.spin_track_min_leaf = QSpinBox()
        self.spin_track_min_leaf.setRange(1, 100)
        self.spin_track_min_leaf.setValue(2)
        self.spin_track_min_leaf.setMaximumWidth(90)
        adv3_form.addRow("Min samples leaf:", make_help_row(
            self.spin_track_min_leaf, "Min samples leaf",
            "Minimum number of samples required to be at a leaf node. Higher values "
            "regularise the tree."
        ))

        self.spin_track_min_split = QSpinBox()
        self.spin_track_min_split.setRange(2, 100)
        self.spin_track_min_split.setValue(2)
        self.spin_track_min_split.setMaximumWidth(90)
        adv3_form.addRow("Min samples split:", make_help_row(
            self.spin_track_min_split, "Min samples split",
            "Minimum number of samples required to split an internal node."
        ))

        self.combo_track_max_features = QComboBox()
        self.combo_track_max_features.addItems(["sqrt", "log2"])
        self.combo_track_max_features.setMaximumWidth(100)
        adv3_form.addRow("Max features:", make_help_row(
            self.combo_track_max_features, "Max features",
            "Number of features considered when looking for the best split: "
            "'sqrt' (recommended) or 'log2'."
        ))

        self.spin_track_n_jobs = QSpinBox()
        self.spin_track_n_jobs.setRange(-1, 64)
        self.spin_track_n_jobs.setValue(-1)
        self.spin_track_n_jobs.setSpecialValueText("All cores (-1)")
        self.spin_track_n_jobs.setMaximumWidth(130)
        adv3_form.addRow("n_jobs:", make_help_row(
            self.spin_track_n_jobs, "n_jobs",
            "Number of parallel jobs for fitting trees (-1 = use all available CPU cores)."
        ))
        adv3.addLayout(adv3_form)
        g_train_rf.addWidget(adv3)

        train_rf_run_row = QHBoxLayout()
        self.btn_train_track = QPushButton("▶ Train RF Classifier")
        _style_primary(self.btn_train_track)
        self.btn_train_track.setEnabled(False)
        train_rf_run_row.addWidget(self.btn_train_track, stretch=1)
        train_rf_run_row.addWidget(self.btn_queue_train_track)
        g_train_rf.addLayout(train_rf_run_row)
        g3.addWidget(grp_train_rf)

        # Sub-section B: Apply Classifier to New Data
        grp_apply_clf = QGroupBox("Apply Classifier to New Data")
        g_apply_clf = QVBoxLayout(grp_apply_clf)
        g_apply_clf.setSpacing(4)

        self.le_apply_clf_path = QLineEdit()
        self.le_apply_clf_path.setPlaceholderText("Path to trained classifier .pkl (auto-fills after training)")
        self.btn_browse_apply_clf = QPushButton("Browse…")
        g_apply_clf.addLayout(_make_browse_row(
            self.le_apply_clf_path, self.btn_browse_apply_clf,
            "Classifier path",
            "Path to a trained Random Forest track classifier .pkl. "
            "Auto-fills after training or when the expected output file exists on disk."
        ))

        self.le_apply_states_path = QLineEdit()
        self.le_apply_states_path.setPlaceholderText("Path to behavioral states .h5ad (auto-fills from State Classification)")
        self.btn_browse_apply_states = QPushButton("Browse…")
        g_apply_clf.addLayout(_make_browse_row(
            self.le_apply_states_path, self.btn_browse_apply_states,
            "States h5ad path",
            "Path to the behavioral states .h5ad file produced by State Classification. "
            "Auto-fills when the expected file exists for the selected cell type."
        ))

        apply_clf_run_row = QHBoxLayout()
        self.btn_apply_track = QPushButton("▶ Apply Classifier")
        _style_secondary(self.btn_apply_track)
        self.btn_apply_track.setEnabled(False)
        apply_clf_run_row.addWidget(self.btn_apply_track, stretch=1)
        apply_clf_run_row.addWidget(self.btn_queue_apply_track)
        g_apply_clf.addLayout(apply_clf_run_row)
        g3.addWidget(grp_apply_clf)

        # ── Step 3: Create Plots (built directly into per-pipeline group boxes) ──
        self.grp_diag = QGroupBox("Diagnostic")
        g_diag = QVBoxLayout(self.grp_diag)
        g_diag.setSpacing(4)
        g_diag.addWidget(_make_info_label(
            "Runs quality-control diagnostics on the trajectory clustering "
            "(silhouette/distance plots)."
        ))
        diag_row = QHBoxLayout()
        self.btn_diagnostics = QPushButton("▶ Create Diagnostics")
        _style_secondary(self.btn_diagnostics)
        diag_row.addWidget(self.btn_diagnostics, stretch=1)
        self.btn_view_diagnostics = _make_view_btn()
        diag_row.addWidget(self.btn_view_diagnostics)
        g_diag.addLayout(diag_row)

        self.grp_track_proportions = QGroupBox("Track Proportions")
        g_prop = QVBoxLayout(self.grp_track_proportions)
        g_prop.setSpacing(4)
        g_prop.addWidget(_make_info_label(
            "Plots how track cluster proportions vary across samples. Tick \"Group conditions\" "
            "under Group in X/Y to pool that axis's levels into two custom groups."
        ))
        self.combo_track_proportion_group_x = QComboBox()
        self.combo_track_proportion_group_x.setMinimumWidth(160)
        self.combo_track_proportion_group_x.currentTextChanged.connect(
            self._on_track_proportion_group_x_changed
        )
        self.chk_track_proportion_group_x_conditions, self.group_selector_track_proportion_x = _add_group_axis_row(
            g_prop, "Group in X:", self.combo_track_proportion_group_x
        )
        self.chk_track_proportion_group_x_conditions.toggled.connect(
            self._on_track_proportion_group_x_conditions_toggled
        )
        self.combo_track_proportion_group_y = QComboBox()
        self.combo_track_proportion_group_y.setMinimumWidth(160)
        self.combo_track_proportion_group_y.currentTextChanged.connect(
            self._on_track_proportion_group_y_changed
        )
        self.chk_track_proportion_group_y_conditions, self.group_selector_track_proportion_y = _add_group_axis_row(
            g_prop, "Group in Y:", self.combo_track_proportion_group_y
        )
        self.chk_track_proportion_group_y_conditions.toggled.connect(
            self._on_track_proportion_group_y_conditions_toggled
        )
        g_prop.addWidget(QLabel("Group per page (Ctrl/Cmd click for multiple):"))
        self.list_track_proportion_group_cols = QListWidget()
        self.list_track_proportion_group_cols.setSelectionMode(QAbstractItemView.ExtendedSelection)
        g_prop.addWidget(self.list_track_proportion_group_cols)
        prop_row = QHBoxLayout()
        self.btn_track_proportions = QPushButton("▶ Create Track Proportion Plots")
        _style_secondary(self.btn_track_proportions)
        prop_row.addWidget(self.btn_track_proportions, stretch=1)
        self.btn_view_track_proportions = _make_view_btn()
        prop_row.addWidget(self.btn_view_track_proportions)
        g_prop.addLayout(prop_row)

        self.grp_window_transitions = QGroupBox("Window Transitions")
        g_wintrans = QVBoxLayout(self.grp_window_transitions)
        g_wintrans.setSpacing(4)
        g_wintrans.addWidget(_make_info_label(
            "Sankey diagram of how each track's trajectory cluster changes from one "
            "window to the next - reconnects the sub-tracks 'Divide long tracks' split "
            "a track into, back through their shared parent TrackID. One diagram per "
            "sample, plus a pooled page. Needs 'Divide long tracks' to have been used."
        ))
        wintrans_row = QHBoxLayout()
        self.btn_window_transitions = QPushButton("▶ Create Window Transition Sankey")
        _style_secondary(self.btn_window_transitions)
        wintrans_row.addWidget(self.btn_window_transitions, stretch=1)
        self.btn_view_window_transitions = _make_view_btn()
        wintrans_row.addWidget(self.btn_view_window_transitions)
        g_wintrans.addLayout(wintrans_row)

        self.grp_track_comparison = QGroupBox("Condition Comparison Report")
        g_track_comparison = QVBoxLayout(self.grp_track_comparison)
        g_track_comparison.setSpacing(4)
        g_track_comparison.addWidget(_make_info_label(
            "Condition comparison (overall proportions, Welch's t-test): each row is one pairwise "
            "comparison of \"Compare condition\"'s levels; \"Group in X\" splits it into side-by-side "
            "columns from another condition. Tick \"Group conditions\" to instead pool levels into "
            "two custom groups and compare just those."
        ))
        track_comparison_form = QFormLayout()
        track_comparison_form.setSpacing(3)
        self.combo_track_comparison_condition_col = QComboBox()
        self.combo_track_comparison_condition_col.setMinimumWidth(200)
        self.combo_track_comparison_condition_col.currentTextChanged.connect(
            self._on_track_comparison_condition_col_changed
        )
        track_comparison_form.addRow("Compare condition:", self.combo_track_comparison_condition_col)
        g_track_comparison.addLayout(track_comparison_form)

        self.combo_track_comparison_group_x = QComboBox()
        self.combo_track_comparison_group_x.setMinimumWidth(160)
        self.combo_track_comparison_group_x.currentTextChanged.connect(
            self._on_track_comparison_group_x_changed
        )
        self.chk_track_comparison_group_x_conditions, self.group_selector_track_comparison_x = _add_group_axis_row(
            g_track_comparison, "Group in X:", self.combo_track_comparison_group_x
        )
        self.chk_track_comparison_group_x_conditions.toggled.connect(
            self._on_track_comparison_group_x_conditions_toggled
        )

        self.line_track_comparison_group_y = QLineEdit()
        self.line_track_comparison_group_y.setMinimumWidth(160)
        self.line_track_comparison_group_y.setReadOnly(True)
        track_y_row = QHBoxLayout()
        track_y_row.addWidget(QLabel("Group in Y:"))
        track_y_row.addWidget(self.line_track_comparison_group_y, stretch=1)
        g_track_comparison.addLayout(track_y_row)
        self.chk_track_comparison_group_conditions = QCheckBox("Group conditions")
        self.chk_track_comparison_group_conditions.setToolTip(
            "Pool the compared column's levels into two custom groups (left/right) and run a "
            "single comparison between them, instead of every pairwise level combination."
        )
        self.chk_track_comparison_group_conditions.toggled.connect(
            self._on_track_comparison_group_conditions_toggled
        )
        g_track_comparison.addWidget(self.chk_track_comparison_group_conditions)
        self.group_selector_track_comparison = DualListGroupSelector()
        self.group_selector_track_comparison.setVisible(False)
        g_track_comparison.addWidget(self.group_selector_track_comparison)

        g_track_comparison.addWidget(QLabel("Group per page (Ctrl/Cmd click for multiple):"))
        self.list_track_comparison_group_cols = QListWidget()
        self.list_track_comparison_group_cols.setSelectionMode(QAbstractItemView.ExtendedSelection)
        g_track_comparison.addWidget(self.list_track_comparison_group_cols)
        track_comparison_row = QHBoxLayout()
        self.btn_track_condition_comparison = QPushButton("▶ Condition Comparison Report")
        _style_secondary(self.btn_track_condition_comparison)
        track_comparison_row.addWidget(self.btn_track_condition_comparison, stretch=1)
        self.btn_view_track_condition_comparison = _make_view_btn()
        track_comparison_row.addWidget(self.btn_view_track_condition_comparison)
        g_track_comparison.addLayout(track_comparison_row)

        self.grp_contact_analysis = QGroupBox("Contact analysis")
        g_contact = QVBoxLayout(self.grp_contact_analysis)
        g_contact.setSpacing(4)
        g_contact.addWidget(QLabel(
            "Contact-based grouping (tracks with vs. without a sufficiently long contact bout):"
        ))
        g_contact.addWidget(_make_info_label(
            "Tick \"Group conditions\" under Group in X/Y to pool that axis's levels into two "
            "custom groups."
        ))
        contact_form = QFormLayout()
        contact_form.setSpacing(3)
        self.combo_contact_col = QComboBox()
        self.combo_contact_col.setMinimumWidth(200)
        contact_form.addRow("Contact column:", self.combo_contact_col)
        self.spin_contact_min_bout = QSpinBox()
        self.spin_contact_min_bout.setRange(1, 100000)
        self.spin_contact_min_bout.setValue(5)
        self.spin_contact_min_bout.setMaximumWidth(100)
        contact_form.addRow("Min. contiguous bout:", make_help_row(
            self.spin_contact_min_bout, "Min. contiguous contact bout (timepoints)",
            "A track counts as 'contact' if it has at least one unbroken run of this many "
            "consecutive contact timepoints; otherwise it is 'no_contact'."
        ))
        self.chk_use_target_class = QCheckBox("Use contact cell classification")
        self.chk_use_target_class.setChecked(False)
        contact_form.addRow("", self.chk_use_target_class)
        self.combo_target_class_source = QComboBox()
        self.combo_target_class_source.addItem("State classification", "state")
        self.combo_target_class_source.addItem("Track classification", "track")
        contact_form.addRow("Target classification:", self.combo_target_class_source)
        self.combo_target_state_col = QComboBox()
        self.combo_target_state_col.addItems(
            ["full_behavioral_cluster", "intrinsic_behavioral_cluster", "raw_hmm_state"]
        )
        contact_form.addRow("Target state column:", self.combo_target_state_col)
        self.label_target_class_warning = QLabel("")
        self.label_target_class_warning.setWordWrap(True)
        self.label_target_class_warning.setStyleSheet(
            "QLabel { background: #3d2200; color: #ffaa44; border-radius: 4px; "
            "padding: 6px 8px; font-size: 11px; }"
        )
        self.label_target_class_warning.hide()
        contact_form.addRow("", self.label_target_class_warning)
        g_contact.addLayout(contact_form)
        self.combo_contact_group_x = QComboBox()
        self.combo_contact_group_x.setMinimumWidth(160)
        self.combo_contact_group_x.currentTextChanged.connect(
            self._on_contact_group_x_changed
        )
        self.chk_contact_group_x_conditions, self.group_selector_contact_x = _add_group_axis_row(
            g_contact, "Group in X:", self.combo_contact_group_x
        )
        self.chk_contact_group_x_conditions.toggled.connect(
            self._on_contact_group_x_conditions_toggled
        )
        self.combo_contact_group_y = QComboBox()
        self.combo_contact_group_y.setMinimumWidth(160)
        self.combo_contact_group_y.currentTextChanged.connect(
            self._on_contact_group_y_changed
        )
        self.chk_contact_group_y_conditions, self.group_selector_contact_y = _add_group_axis_row(
            g_contact, "Group in Y:", self.combo_contact_group_y
        )
        self.chk_contact_group_y_conditions.toggled.connect(
            self._on_contact_group_y_conditions_toggled
        )
        g_contact.addWidget(QLabel("Group per page (Ctrl/Cmd click for multiple):"))
        self.list_contact_group_cols = QListWidget()
        self.list_contact_group_cols.setSelectionMode(QAbstractItemView.ExtendedSelection)
        g_contact.addWidget(self.list_contact_group_cols)
        contact_row = QHBoxLayout()
        self.btn_contact_analysis = QPushButton("▶ Run Contact-vs-No-Contact Analysis")
        _style_secondary(self.btn_contact_analysis)
        contact_row.addWidget(self.btn_contact_analysis, stretch=1)
        self.btn_view_contact_analysis = _make_view_btn()
        contact_row.addWidget(self.btn_view_contact_analysis)
        g_contact.addLayout(contact_row)

        g_contact.addWidget(QLabel(
            "Contact duration comparison (how long tracks stay in contact with each touched class "
            "of 'Use contact cell classification', pairwise and each class vs. the rest):"
        ))
        duration_form = QFormLayout()
        duration_form.setSpacing(3)
        self.combo_duration_test_mode = QComboBox()
        self.combo_duration_test_mode.addItem("Welch's t-test (unpaired)", "welch")
        self.combo_duration_test_mode.addItem("Paired t-test", "paired")
        self.combo_duration_test_mode.currentIndexChanged.connect(self._on_duration_test_mode_changed)
        duration_form.addRow("Test:", make_help_row(
            self.combo_duration_test_mode, "Comparison test",
            "'Welch's t-test': compares the raw per-touch durations between the two groups without "
            "assuming equal variances. 'Paired t-test': averages durations within each 'Pairing "
            "column' value first, then pairs the two groups' averages for the same value (e.g. the "
            "same sample) — pairing units missing either side are dropped."
        ))
        self.combo_duration_pairing_col = QComboBox()
        duration_form.addRow("Pairing column:", self.combo_duration_pairing_col)
        self.spin_duration_comparisons_per_page = QSpinBox()
        self.spin_duration_comparisons_per_page.setRange(1, 100)
        self.spin_duration_comparisons_per_page.setValue(12)
        self.spin_duration_comparisons_per_page.setMaximumWidth(100)
        duration_form.addRow("Comparisons per page:", self.spin_duration_comparisons_per_page)
        g_contact.addLayout(duration_form)
        duration_row = QHBoxLayout()
        self.btn_duration_comparison = QPushButton("▶ Create Contact Duration Comparison")
        _style_secondary(self.btn_duration_comparison)
        duration_row.addWidget(self.btn_duration_comparison, stretch=1)
        self.btn_view_duration_comparison = _make_view_btn()
        duration_row.addWidget(self.btn_view_duration_comparison)
        g_contact.addLayout(duration_row)
        self._on_duration_test_mode_changed()

        g_contact.addWidget(QLabel("State-shift analysis (behavioral state before → after contact):"))
        shift_form = QFormLayout()
        shift_form.setSpacing(3)
        self.combo_contact_shift_window_mode = QComboBox()
        self.combo_contact_shift_window_mode.addItems(["Fixed", "Full"])
        shift_form.addRow("Window mode:", make_help_row(
            self.combo_contact_shift_window_mode, "Before/after window mode",
            "'Fixed': compare a fixed number of timepoints immediately before the contact bout "
            "starts vs. immediately after it ends. 'Full': compare everything before the bout to "
            "everything after it, within the track's classified window."
        ))
        self.spin_contact_shift_window_length = QSpinBox()
        self.spin_contact_shift_window_length.setRange(1, 100000)
        self.spin_contact_shift_window_length.setValue(10)
        self.spin_contact_shift_window_length.setMaximumWidth(100)
        shift_form.addRow("Fixed window length:", make_help_row(
            self.spin_contact_shift_window_length, "Fixed window length (timepoints)",
            "Only used in 'Fixed' window mode: number of timepoints immediately before the bout "
            "start / after the bout end to compare."
        ))
        self.combo_contact_shift_state_col = QComboBox()
        self.combo_contact_shift_state_col.addItems(
            ["full_behavioral_cluster", "intrinsic_behavioral_cluster", "raw_hmm_state"]
        )
        shift_form.addRow("State column:", self.combo_contact_shift_state_col)
        g_contact.addLayout(shift_form)
        shift_row = QHBoxLayout()
        self.btn_contact_state_shift = QPushButton("▶ Run State-Shift Analysis")
        _style_secondary(self.btn_contact_state_shift)
        shift_row.addWidget(self.btn_contact_state_shift, stretch=1)
        self.btn_view_contact_state_shift = _make_view_btn()
        shift_row.addWidget(self.btn_view_contact_state_shift)
        g_contact.addLayout(shift_row)

        g_contact.addWidget(QLabel(
            "Track contact overview (full state trajectory + contact-bout markers, one page per sample):"
        ))
        overview_form = QFormLayout()
        overview_form.setSpacing(3)
        self.spin_track_overview_rows_per_page = QSpinBox()
        self.spin_track_overview_rows_per_page.setRange(1, 100)
        self.spin_track_overview_rows_per_page.setValue(6)
        self.spin_track_overview_rows_per_page.setMaximumWidth(100)
        overview_form.addRow("Tracks per page:", self.spin_track_overview_rows_per_page)
        g_contact.addLayout(overview_form)
        overview_row = QHBoxLayout()
        self.btn_track_contact_overview = QPushButton("▶ Create Track Contact Overview")
        _style_secondary(self.btn_track_contact_overview)
        overview_row.addWidget(self.btn_track_contact_overview, stretch=1)
        self.btn_view_track_contact_overview = _make_view_btn()
        overview_row.addWidget(self.btn_view_track_contact_overview)
        g_contact.addLayout(overview_row)

        self.grp_exemplar = QGroupBox("Exemplar Tracks")
        g_exemplar = QVBoxLayout(self.grp_exemplar)
        g_exemplar.setSpacing(4)
        g_exemplar.addWidget(_make_info_label(
            "Selects and renders representative example tracks per cluster."
        ))
        ex_form = QFormLayout()
        ex_form.setSpacing(3)
        self.spin_n_per_cluster = QSpinBox()
        self.spin_n_per_cluster.setRange(1, 100)
        self.spin_n_per_cluster.setValue(10)
        self.spin_n_per_cluster.setMaximumWidth(80)
        ex_form.addRow("Exemplars / cluster:", make_help_row(
            self.spin_n_per_cluster, "Exemplars per cluster",
            "Number of exemplar tracks to select per trajectory cluster."
        ))
        g_exemplar.addLayout(ex_form)

        opts_row = QHBoxLayout()
        self.chk_overview_statebars = QCheckBox("Overview statebars")
        self.chk_overview_statebars.setChecked(True)
        opts_row.addWidget(self.chk_overview_statebars)
        opts_row.addWidget(HelpButton("Overview statebars",
            "Generate a PDF with stacked state bar plots for all exemplar tracks."))

        self.chk_backproj_pdf = QCheckBox("PDFs")
        self.chk_backproj_pdf.setChecked(False)
        opts_row.addWidget(self.chk_backproj_pdf)
        opts_row.addWidget(HelpButton("PDFs",
            "Generate backprojection PDFs showing exemplar tracks overlaid on raw images."))

        self.chk_backproj_mp4 = QCheckBox("MP4")
        self.chk_backproj_mp4.setChecked(False)
        opts_row.addWidget(self.chk_backproj_mp4)
        opts_row.addWidget(HelpButton("MP4",
            "Generate backprojection video (MP4) for exemplar tracks."))
        opts_row.addStretch()
        g_exemplar.addLayout(opts_row)

        exemplar_row = QHBoxLayout()
        self.btn_exemplars = QPushButton("▶ Create Exemplar PDFs")
        _style_secondary(self.btn_exemplars)
        exemplar_row.addWidget(self.btn_exemplars, stretch=1)
        self.btn_view_exemplars = _make_view_btn()
        exemplar_row.addWidget(self.btn_view_exemplars)
        g_exemplar.addLayout(exemplar_row)

        from behav3d.napari._guided import GuidedPanel
        from behav3d.napari.analysis_guided_copy import (
            TRACK_PLOTS_ENTRY, TRACK_PLOT_PIPELINES,
        )

        # Entry trigger stays inline on the clustering-settings page; Start
        # navigates to the fully isolated plotting page (outer page 1) below.
        lay.addWidget(
            GuidedPanel([TRACK_PLOTS_ENTRY],
                        start_cb=lambda _id: self._open_plotting_page())
        )

        plotting_page = QWidget()
        plotting_lay = QVBoxLayout(plotting_page)
        plotting_lay.setContentsMargins(0, 0, 0, 0)
        plotting_lay.setSpacing(0)

        # Inner stack: page 0 = report list, page 1 = a report's settings.
        self._plots_stack = QStackedWidget()

        pipeline_page = QWidget()
        pipeline_lay = QVBoxLayout(pipeline_page)
        pipeline_lay.setContentsMargins(0, 0, 0, 0)
        pipeline_lay.setSpacing(0)
        pipeline_lay.addWidget(
            GuidedPanel(TRACK_PLOT_PIPELINES, start_cb=self._on_pipeline_start)
        )
        self._plots_stack.addWidget(pipeline_page)  # inner page 0

        report_settings_page = QWidget()
        report_settings_lay = QVBoxLayout(report_settings_page)
        report_settings_lay.setContentsMargins(0, 0, 0, 0)
        report_settings_lay.setSpacing(0)
        self._pipeline_scroll = QScrollArea()
        self._pipeline_scroll.setWidgetResizable(True)
        self._pipeline_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        report_settings_lay.addWidget(self._pipeline_scroll)
        self._plots_stack.addWidget(report_settings_page)  # inner page 1

        pipeline_content = QWidget()
        pipeline_content_lay = QVBoxLayout(pipeline_content)
        pipeline_content_lay.setContentsMargins(6, 6, 6, 6)
        pipeline_content_lay.setSpacing(6)
        pipeline_content_lay.addWidget(self.grp_diag)
        pipeline_content_lay.addWidget(self.grp_track_proportions)
        pipeline_content_lay.addWidget(self.grp_window_transitions)
        pipeline_content_lay.addWidget(self.grp_track_comparison)
        pipeline_content_lay.addWidget(self.grp_contact_analysis)
        pipeline_content_lay.addWidget(self.grp_exemplar)
        pipeline_content_lay.addStretch()
        self._pipeline_scroll.setWidget(pipeline_content)
        reset_scroll_on_page_change(self._plots_stack)

        self._plots_stack.setCurrentIndex(0)
        plotting_lay.addWidget(self._plots_stack)
        # Kept as an attribute so `_check_prerequisites` can disable *only* this
        # plotting page (Steps 3/4 create-plots) when state adata is missing,
        # without disabling the whole `_subtab_stack` — which would also grey out
        # page 0's Step 1 (grp1) and its "Run Original BEHAV3D DTW" button.
        self._plotting_page = plotting_page
        self._subtab_stack.addWidget(plotting_page)  # outer page 1
        reset_scroll_on_page_change(self._subtab_stack)

        self._subtab_stack.setCurrentIndex(0)

        # ── Step 5: Backprojection ───────────────────────────────────────
        self.grp_bp = QGroupBox("Step 4 — Backprojection")
        g_bp5 = QVBoxLayout(self.grp_bp)
        g_bp5.setSpacing(6)

        # Live napari view
        grp_track_view = QGroupBox("Live Napari Layer Backprojection (Tracks)")
        g_track_view = QVBoxLayout(grp_track_view)
        g_track_view.setSpacing(4)

        info_tbp = QLabel(
            "Loads the track cluster adata and overlays colored trajectory cluster "
            "labels onto the selected sample image in napari."
        )
        info_tbp.setStyleSheet("color: #999; font-size: 10px;")
        info_tbp.setWordWrap(True)
        g_track_view.addWidget(info_tbp)

        tbp_form = QFormLayout()
        tbp_form.setSpacing(3)

        self.sample_combo = QComboBox()
        self.sample_combo.addItem("— All samples —")
        self.sample_combo.setMinimumWidth(120)
        self.sample_combo.setToolTip(
            "Select a sample for backprojection visualisation. "
            "'— All samples —' uses the first available sample."
        )
        tbp_form.addRow("Sample:", self.sample_combo)

        self.combo_track_color_by = QComboBox()
        # populated dynamically from actual obs columns when track adata is loaded
        self.combo_track_color_by.setMinimumWidth(220)
        tbp_form.addRow("Color by:", make_help_row(
            self.combo_track_color_by, "Color by",
            "Which track cluster label to use for coloring the backprojection overlay."
        ))

        self.spin_track_opacity = QSpinBox()
        self.spin_track_opacity.setRange(10, 100)
        self.spin_track_opacity.setValue(80)
        self.spin_track_opacity.setSuffix(" %")
        self.spin_track_opacity.setMaximumWidth(90)
        tbp_form.addRow("Opacity:", make_help_row(
            self.spin_track_opacity, "Opacity",
            "Opacity of the colored track overlay layer in napari (10–100%)."
        ))

        self.chk_show_track_trajectories = QCheckBox("Show trajectories")
        self.chk_show_track_trajectories.setChecked(True)
        tbp_form.addRow("Trajectories:", _make_chk_help_row(
            self.chk_show_track_trajectories, "Show trajectories",
            "Overlay each track's full path as a colored line, matching its class "
            "color below. On by default — adds one napari Tracks layer per class."
        ))
        g_track_view.addLayout(tbp_form)

        track_view_row = QHBoxLayout()
        self.btn_show_track_bp = QPushButton("▶ Show Track Backprojection in Napari")
        _style_primary(self.btn_show_track_bp)
        track_view_row.addWidget(self.btn_show_track_bp, stretch=1)
        g_track_view.addLayout(track_view_row)

        track_export_run_row = QHBoxLayout()
        self.btn_export_track_bp = QPushButton("▶ Export Track Backprojection")
        _style_secondary(self.btn_export_track_bp)
        track_export_run_row.addWidget(self.btn_export_track_bp, stretch=1)
        g_track_view.addLayout(track_export_run_row)
        g_bp5.addWidget(grp_track_view)
        lay.addWidget(self.grp_bp)
        lay.addWidget(self.grp3)

        # ── Progress + Log ───────────────────────────────────────────────
        self.progress_row = ProgressBarRow()
        lay.addWidget(self.progress_row)

        lay.addWidget(QLabel("Log"))
        self.log_edit = QTextEdit()
        self.log_edit.setReadOnly(True)
        self.log_edit.setMaximumHeight(120)
        self.log_edit.setStyleSheet(
            "background: #111; color: #ccc; font-family: monospace; font-size: 10px;"
        )
        lay.addWidget(self.log_edit)
        lay.addStretch(1)

        self._setup_signals()
        self._apply_clustering_controls_mode()

    def _setup_signals(self):
        self.chk_apply_pretrained.toggled.connect(self._toggle_pretrained_mode)
        self.combo_trajectory_basis.currentTextChanged.connect(self._apply_clustering_controls_mode)
        self.combo_clustering_method.currentTextChanged.connect(self._apply_clustering_controls_mode)
        self.chk_use_original.toggled.connect(self._on_original_toggled_from_adv)
        self.chk_use_original_top.toggled.connect(self._on_original_toggled_from_top)
        self.btn_run_track.clicked.connect(self._on_run_cluster)
        self.btn_rename_track.clicked.connect(self._on_rename_track)
        self.btn_train_track.clicked.connect(self._on_train_track)
        self.btn_apply_track.clicked.connect(self._on_apply_track)
        self.btn_exemplars.clicked.connect(self._on_exemplars)
        _wire_view_btn(self.btn_view_exemplars, self._on_view, "track_exemplars")
        self.btn_diagnostics.clicked.connect(self._on_diagnostics)
        _wire_view_btn(self.btn_view_diagnostics, self._on_view, "track_diagnostics")
        self.btn_track_proportions.clicked.connect(self._on_track_proportions)
        _wire_view_btn(self.btn_view_track_proportions, self._on_view, "track_proportions")
        self.btn_window_transitions.clicked.connect(self._on_window_transitions)
        _wire_view_btn(self.btn_view_window_transitions, self._on_view, "window_transitions")
        self.btn_track_condition_comparison.clicked.connect(self._on_track_condition_comparison)
        _wire_view_btn(self.btn_view_track_condition_comparison, self._on_view, "track_condition_comparison")
        self.btn_contact_analysis.clicked.connect(self._on_contact_analysis)
        _wire_view_btn(self.btn_view_contact_analysis, self._on_view, "contact_analysis")
        self.combo_contact_col.currentTextChanged.connect(self._sync_contact_target_class_controls)
        self.chk_use_target_class.stateChanged.connect(self._sync_contact_target_class_controls)
        self.combo_target_class_source.currentTextChanged.connect(self._sync_contact_target_class_controls)
        self.btn_duration_comparison.clicked.connect(self._on_duration_comparison)
        _wire_view_btn(self.btn_view_duration_comparison, self._on_view, "contact_duration")
        self.btn_contact_state_shift.clicked.connect(self._on_contact_state_shift_analysis)
        _wire_view_btn(self.btn_view_contact_state_shift, self._on_view, "contact_state_shift")
        self.btn_track_contact_overview.clicked.connect(self._on_track_contact_overview)
        _wire_view_btn(self.btn_view_track_contact_overview, self._on_view, "track_contact_overview")
        self.btn_browse_pretrained_clf.clicked.connect(self._browse_pretrained_clf)
        self.btn_browse_pretrained_states.clicked.connect(self._browse_pretrained_states)
        self.btn_run_apply_pretrained.clicked.connect(self._on_apply_pretrained)
        self.btn_browse_apply_clf.clicked.connect(self._browse_apply_clf)
        self.btn_browse_apply_states.clicked.connect(self._browse_apply_states)
        self.btn_show_track_bp.clicked.connect(self._on_show_track_bp)
        self.btn_export_track_bp.clicked.connect(self._on_export_track_bp)

    # ── Prerequisite check ───────────────────────────────────────────────

    def _check_prerequisites(self) -> bool:
        """Check if behavioral states h5ad exists; conditionally enable/disable steps.

        - If state adata is absent: Step 1 (grp1) stays enabled but is locked into
          'original BEHAV3D DTW' mode (chk_use_original forced True + disabled).
          Steps 2-5 are disabled since they all require state adata.
        - If state adata is present: all steps enabled; if the checkbox was previously
          force-locked, it is automatically unchecked and re-enabled (standard mode).
        """
        ct = self._cell_type()
        out = self._out_dir()
        if not ct or not out:
            self.warning_label.hide()
            return True

        states_path = self._state_adata_path(ct)
        if not states_path or not states_path.exists():
            self.warning_label.setText(
                f"⚠ Behavioral states not found for cell type '{ct}'.\n"
                "Run State Classification first to unlock all steps.\n"
                "You can still run the original BEHAV3D DTW clustering (Step 1) below."
            )
            self.warning_label.show()

            # Step 1 stays available, but lock into 'original DTW' mode.
            # NB: disable only the plotting page, not the whole `_subtab_stack`
            # — grp1 lives on page 0 of that stack, so disabling the stack would
            # also grey out Step 1's "Run Original BEHAV3D DTW" button despite
            # the setEnabled(True) above.
            self.grp1.setEnabled(True)
            for grp in [self.grp2, self.grp3, self._plotting_page, self.grp_bp]:
                grp.setEnabled(False)

            # Force 'use original' on and prevent the user from unchecking it.
            for chk in (self.chk_use_original, self.chk_use_original_top):
                chk.blockSignals(True)
                chk.setChecked(True)
                chk.setEnabled(False)
                chk.blockSignals(False)
            self._apply_original_mode(True)

            return False
        else:
            self.warning_label.hide()
            for grp in [self.grp1, self.grp2, self.grp3, self._plotting_page, self.grp_bp]:
                grp.setEnabled(True)

            # Only revert to standard mode if the checkbox was previously force-locked
            # (i.e. disabled because state adata was missing). Otherwise preserve the
            # user's/config's current choice.
            was_force_locked = not self.chk_use_original.isEnabled()
            for chk in (self.chk_use_original, self.chk_use_original_top):
                chk.blockSignals(True)
                if was_force_locked:
                    chk.setChecked(False)
                chk.setEnabled(True)
                chk.blockSignals(False)
            self._apply_original_mode(self.chk_use_original.isChecked())

            _apply_group_tracked_gate(self.warning_label, self.grp_bp, self.metadata_loader, ct)
            return True

    @staticmethod
    def _read_h5ad_obs_column(f, col):
        """Read one obs column from an open h5ad file, handling both plain
        arrays and AnnData's categorical (categories/codes) group encoding."""
        import h5py
        node = f["obs"][col]
        if isinstance(node, h5py.Group):
            categories = node["categories"][:]
            categories = [c.decode() if isinstance(c, bytes) else c for c in categories]
            codes = node["codes"][:]
            return [categories[c] if c >= 0 else None for c in codes]
        data = node[:]
        return [v.decode() if isinstance(v, bytes) else v for v in data]

    def _max_available_track_length(self, ct: str):
        """Longest (sample_name, TrackID) track in this cell type's behavioral-states
        h5ad — the same file/keys the dtaidistance clustering step reads — so the
        Trajectory size spinbox can never be pushed past what any track can supply."""
        states_path = self._state_adata_path(ct) if ct else None
        if not states_path or not states_path.exists():
            return None
        try:
            import h5py
            import pandas as pd
            with h5py.File(str(states_path), "r") as f:
                obs = f.get("obs", {})
                if "sample_name" not in obs or "TrackID" not in obs:
                    return None
                sample_names = self._read_h5ad_obs_column(f, "sample_name")
                track_ids = self._read_h5ad_obs_column(f, "TrackID")
            if not sample_names or len(sample_names) != len(track_ids):
                return None
            counts = pd.DataFrame({"sample_name": sample_names, "TrackID": track_ids}).groupby(
                ["sample_name", "TrackID"]
            ).size()
            max_len = int(counts.max()) if len(counts) else 0
            return max_len if max_len > 0 else None
        except Exception:
            return None

    def _apply_max_trajectory_size(self, ct: str):
        """Cap the Trajectory size spinbox at the longest available track,
        so a value guaranteed to filter out every track can't be selected."""
        max_len = self._max_available_track_length(ct)
        self.spin_traj_size.setMaximum(max_len if max_len else 9999)

    # ── Guided pipeline dispatch (Step 3 create plots) ───────────────────
    _TRACK_PIPELINE_RUN_BUTTONS = {
        "track_diagnostics": "btn_diagnostics",
        "track_window_transitions": "btn_window_transitions",
    }

    def _on_pipeline_start(self, pipeline_id: str):
        """Start from a pipeline card: run directly (no params) or open its settings."""
        from behav3d.napari.analysis_guided_copy import TRACK_PLOT_PIPELINES
        spec = next((s for s in TRACK_PLOT_PIPELINES if s["id"] == pipeline_id), None)
        if spec is not None and not spec.get("has_params", True):
            btn_name = self._TRACK_PIPELINE_RUN_BUTTONS.get(pipeline_id)
            if btn_name is not None:
                getattr(self, btn_name).click()
            return
        self._focus_pipeline(pipeline_id)
        self._plots_stack.setCurrentIndex(1)
        self._pipeline_scroll.verticalScrollBar().setValue(0)

    def _focus_pipeline(self, pipeline_id: str):
        """Show only the group box relevant to ``pipeline_id``, hide the rest."""
        from behav3d.napari.analysis_guided_copy import TRACK_PLOT_PIPELINES
        visible = {
            "track_diagnostics": {self.grp_diag},
            "track_proportions": {self.grp_track_proportions},
            "track_window_transitions": {self.grp_window_transitions},
            "track_comparison": {self.grp_track_comparison},
            "track_contact": {self.grp_contact_analysis},
            "track_exemplars": {self.grp_exemplar},
        }.get(pipeline_id, set())
        for group in (self.grp_diag, self.grp_track_proportions, self.grp_window_transitions,
                      self.grp_track_comparison, self.grp_contact_analysis, self.grp_exemplar):
            group.setVisible(group in visible)
        title = next(
            (s["title"] for s in TRACK_PLOT_PIPELINES if s["id"] == pipeline_id), ""
        )
        self._pipeline_title.setText(title)

    def _open_plotting_page(self):
        """Entry card Start: jump to the isolated plotting page, reset to the
        report list."""
        self._plots_stack.setCurrentIndex(0)
        self._pipeline_title.setText("")
        self._subtab_stack.setCurrentIndex(1)
        self._plotting_header.show()

    def _plotting_back(self):
        """Contextual Back inside the plotting page: report settings → report
        list; report list → clustering settings."""
        if self._plots_stack.currentIndex() == 1:
            self._plots_stack.setCurrentIndex(0)
            self._pipeline_title.setText("")
        else:
            self._subtab_stack.setCurrentIndex(0)
            self._plotting_header.hide()

    # ── Toggle helpers ───────────────────────────────────────────────────

    def _toggle_pretrained_mode(self, checked: bool):
        """Hide Steps 1–2 and the train classifier box when apply-pretrained mode is active."""
        self.grp1.setVisible(not checked)
        self.grp2.setVisible(not checked)
        self.grp3.setVisible(not checked)
        self.grp_apply_pretrained.setVisible(checked)

    def _on_original_toggled_from_adv(self, checked: bool):
        """Checkbox in Advanced Config toggled → sync top checkbox + apply mode."""
        self.chk_use_original_top.blockSignals(True)
        self.chk_use_original_top.setChecked(checked)
        self.chk_use_original_top.blockSignals(False)
        self._apply_original_mode(checked)
        if checked:
            self._show_original_dtw_disclaimer()

    def _on_original_toggled_from_top(self, checked: bool):
        """Top checkbox toggled (user unchecks from top) → sync adv checkbox + apply mode."""
        self.chk_use_original.blockSignals(True)
        self.chk_use_original.setChecked(checked)
        self.chk_use_original.blockSignals(False)
        self._apply_original_mode(checked)
        if checked:
            self._show_original_dtw_disclaimer()

    def _apply_clustering_controls_mode(self, *_args):
        """Update visibility of basis/method-specific controls (Linkage vs Bouts
        linkage vs Leiden neighbors/resolution; DTW-only technical controls vs
        bouts feature toggles)."""
        is_bouts = self.combo_trajectory_basis.currentText() == "bouts"
        is_leiden = self.combo_clustering_method.currentText() == "leiden"

        self._leiden_frame.setVisible(is_leiden)
        self._dtw_linkage_frame.setVisible(not is_leiden and not is_bouts)
        self._bouts_linkage_frame.setVisible(not is_leiden and is_bouts)

        # Trim mode / divide-long-tracks apply to both bases, so they stay visible
        # regardless of is_bouts. Parallel/save-distance-matrix are DTW-distance-matrix
        # specific and bouts' own feature toggles are bouts-only.
        self._dtw_technical_frame.setVisible(not is_bouts)
        self._bouts_frame.setVisible(is_bouts)

        if not self.chk_use_original.isChecked():
            self.btn_run_track.setText(
                "▶ Run Bout/Proportion Clustering" if is_bouts else "▶ Run Track Clustering"
            )

    def _apply_original_mode(self, checked: bool):
        """Update visibility of UI sections for original vs dtaidistance/bouts mode."""
        self.chk_use_original_top.setVisible(checked)
        self.adv1.setVisible(not checked)
        self._basis_method_frame.setVisible(not checked)
        self._umap_frame.setVisible(checked)
        if checked:
            self.btn_run_track.setText("▶ Run Original BEHAV3D DTW")
        else:
            self._apply_clustering_controls_mode()

    def _show_original_dtw_disclaimer(self):
        """Warn the user that the original BEHAV3D DTW pipeline requires equal-length tracks."""
        if getattr(self, '_hide_original_dtw_disclaimer', False):
            return

        msg_box = QMessageBox(self)
        msg_box.setIcon(QMessageBox.Warning)
        msg_box.setWindowTitle("Track length requirement")
        msg_box.setText(
            "The original feature-based BEHAV3D DTW pipeline requires that all tracks "
            "be of the same length for a reliable result.\n\n"
            "If your tracks are not all the same length, please go back to the "
            "Filtering tab and filter/trim them to a uniform length before running "
            "this analysis."
        )
        
        cb = QCheckBox("Don't show again")
        msg_box.setCheckBox(cb)
        msg_box.exec_()
        
        if cb.isChecked():
            self._hide_original_dtw_disclaimer = True

    # ── Browse helpers ───────────────────────────────────────────────────

    def _browse_pretrained_clf(self):
        fpath, _ = QFileDialog.getOpenFileName(
            self, "Select Classifier", "", "Pickle Files (*.pkl)"
        )
        if fpath:
            self.le_pretrained_clf_path.setText(fpath)

    def _browse_pretrained_states(self):
        fpath, _ = QFileDialog.getOpenFileName(
            self, "Select Behavioral States", "", "H5AD Files (*.h5ad)"
        )
        if fpath:
            self.le_pretrained_states_path.setText(fpath)

    def _browse_apply_clf(self):
        fpath, _ = QFileDialog.getOpenFileName(
            self, "Select Classifier", "", "Pickle Files (*.pkl)"
        )
        if fpath:
            self.le_apply_clf_path.setText(fpath)

    def _browse_apply_states(self):
        fpath, _ = QFileDialog.getOpenFileName(
            self, "Select Behavioral States", "", "H5AD Files (*.h5ad)"
        )
        if fpath:
            self.le_apply_states_path.setText(fpath)

    # ── Config persistence / restoration ────────────────────────────────────

    def _collect_track_params(self, ct: str) -> dict:
        return {
            "behavioral_trajectory_size": int(self.spin_traj_size.value()),
            "n_clusters":                 int(self.spin_n_clusters.value()),
            "trajectory_basis":           self.combo_trajectory_basis.currentText(),
            "linkage":                    self.combo_linkage.currentText(),
            "clustering_method":          self.combo_clustering_method.currentText(),
            "leiden_n_neighbors":         int(self.spin_leiden_neighbors.value()),
            "leiden_resolution":          float(self.spin_leiden_resolution.value()),
            "bouts_linkage":              self.combo_bouts_linkage.currentText(),
            "bouts_use_fractions":        self.chk_bouts_use_fractions.isChecked(),
            "bouts_use_bout_stats":       self.chk_bouts_use_bout_stats.isChecked(),
            "bouts_use_transitions":      self.chk_bouts_use_transitions.isChecked(),
            "bouts_use_ngrams":           self.chk_bouts_use_ngrams.isChecked(),
            "bouts_do_pca":               self.chk_bouts_do_pca.isChecked(),
            "bouts_use_clr":              self.chk_bouts_use_clr.isChecked(),
            "bouts_log_bout_length":      self.chk_bouts_log_bout_length.isChecked(),
            "bouts_block_scaling":        self.chk_bouts_block_scaling.isChecked(),
            "bouts_drop_redundant":       self.chk_bouts_drop_redundant.isChecked(),
            "bouts_plot_exemplars":       self.chk_bouts_plot_exemplars.isChecked(),
            "trajectory_trim_mode":       self.combo_trim.currentText(),
            "split_long_tracks":          self.chk_split_long_tracks.isChecked(),
            "parallel":                   self.chk_parallel.isChecked(),
            "save_distance_matrix":       self.chk_save_dist.isChecked(),
            "random_state":               int(self.spin_seed.value()),
            "use_original":               self.chk_use_original.isChecked(),
            "umap_n_neighbors":           int(self.spin_umap_neighbors.value()),
            "umap_min_dist":              float(self.spin_umap_min_dist.value()),
            "rf_n_estimators":            int(self.spin_track_n_est.value()),
            "rf_max_depth":               int(self.spin_track_max_depth.value()),
            "rf_test_size_pct":           float(self.spin_track_test_pct.value()),
            "rf_min_samples_leaf":        int(self.spin_track_min_leaf.value()),
            "rf_min_samples_split":       int(self.spin_track_min_split.value()),
            "rf_max_features":            self.combo_track_max_features.currentText(),
            "rf_n_jobs":                  int(self.spin_track_n_jobs.value()),
        }

    def _persist_track_cfg(self, ct: str):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if not isinstance(params, dict):
            return
        params.setdefault("track_classification", {})[ct] = self._collect_track_params(ct)
        _save_behav3d_params(self.metadata_loader, self._out_dir)

    def _populate_track_settings(self, ct: str):
        cfg = getattr(self.metadata_loader, "behav3d_parameters", {}).get("track_classification", {}).get(ct, {})
        if not cfg:
            return
        if "behavioral_trajectory_size" in cfg:
            self.spin_traj_size.setValue(int(cfg["behavioral_trajectory_size"]))
        if "n_clusters" in cfg:
            self.spin_n_clusters.setValue(int(cfg["n_clusters"]))
        if "trajectory_basis" in cfg:
            self.combo_trajectory_basis.setCurrentText(cfg["trajectory_basis"])
        if "linkage" in cfg:
            self.combo_linkage.setCurrentText(cfg["linkage"])
        if "clustering_method" in cfg:
            self.combo_clustering_method.setCurrentText(cfg["clustering_method"])
        if "leiden_n_neighbors" in cfg:
            self.spin_leiden_neighbors.setValue(int(cfg["leiden_n_neighbors"]))
        if "leiden_resolution" in cfg:
            self.spin_leiden_resolution.setValue(float(cfg["leiden_resolution"]))
        if "bouts_linkage" in cfg:
            self.combo_bouts_linkage.setCurrentText(cfg["bouts_linkage"])
        if "bouts_use_fractions" in cfg:
            self.chk_bouts_use_fractions.setChecked(bool(cfg["bouts_use_fractions"]))
        if "bouts_use_bout_stats" in cfg:
            self.chk_bouts_use_bout_stats.setChecked(bool(cfg["bouts_use_bout_stats"]))
        if "bouts_use_transitions" in cfg:
            self.chk_bouts_use_transitions.setChecked(bool(cfg["bouts_use_transitions"]))
        if "bouts_use_ngrams" in cfg:
            self.chk_bouts_use_ngrams.setChecked(bool(cfg["bouts_use_ngrams"]))
        if "bouts_do_pca" in cfg:
            self.chk_bouts_do_pca.setChecked(bool(cfg["bouts_do_pca"]))
        if "bouts_use_clr" in cfg:
            self.chk_bouts_use_clr.setChecked(bool(cfg["bouts_use_clr"]))
        if "bouts_log_bout_length" in cfg:
            self.chk_bouts_log_bout_length.setChecked(bool(cfg["bouts_log_bout_length"]))
        if "bouts_block_scaling" in cfg:
            self.chk_bouts_block_scaling.setChecked(bool(cfg["bouts_block_scaling"]))
        if "bouts_drop_redundant" in cfg:
            self.chk_bouts_drop_redundant.setChecked(bool(cfg["bouts_drop_redundant"]))
        if "bouts_plot_exemplars" in cfg:
            self.chk_bouts_plot_exemplars.setChecked(bool(cfg["bouts_plot_exemplars"]))
        self._apply_clustering_controls_mode()
        if "trajectory_trim_mode" in cfg:
            self.combo_trim.setCurrentText(cfg["trajectory_trim_mode"])
        if "split_long_tracks" in cfg:
            self.chk_split_long_tracks.setChecked(bool(cfg["split_long_tracks"]))
        if "parallel" in cfg:
            self.chk_parallel.setChecked(bool(cfg["parallel"]))
        if "save_distance_matrix" in cfg:
            self.chk_save_dist.setChecked(bool(cfg["save_distance_matrix"]))
        if "random_state" in cfg:
            self.spin_seed.setValue(int(cfg["random_state"]))
        if "umap_n_neighbors" in cfg:
            self.spin_umap_neighbors.setValue(int(cfg["umap_n_neighbors"]))
        if "umap_min_dist" in cfg:
            self.spin_umap_min_dist.setValue(float(cfg["umap_min_dist"]))
        if "rf_n_estimators" in cfg:
            self.spin_track_n_est.setValue(int(cfg["rf_n_estimators"]))
        if "rf_max_depth" in cfg:
            self.spin_track_max_depth.setValue(int(cfg["rf_max_depth"]))
        if "rf_test_size_pct" in cfg:
            self.spin_track_test_pct.setValue(float(cfg["rf_test_size_pct"]))
        if "rf_min_samples_leaf" in cfg:
            self.spin_track_min_leaf.setValue(int(cfg["rf_min_samples_leaf"]))
        if "rf_min_samples_split" in cfg:
            self.spin_track_min_split.setValue(int(cfg["rf_min_samples_split"]))
        if "rf_max_features" in cfg:
            self.combo_track_max_features.setCurrentText(cfg["rf_max_features"])
        if "rf_n_jobs" in cfg:
            self.spin_track_n_jobs.setValue(int(cfg["rf_n_jobs"]))
        if "use_original" in cfg:
            val = bool(cfg["use_original"])
            for chk in (self.chk_use_original, self.chk_use_original_top):
                chk.blockSignals(True)
                chk.setChecked(val)
                chk.blockSignals(False)
            self._apply_original_mode(val)

    # ── Metadata / reload ────────────────────────────────────────────────

    def on_metadata_updated(self):
        self._reload()

    def _update_timepoint_time_labels(self, *_):
        metadata = getattr(self.metadata_loader, "metadata", None)
        self.lbl_traj_size_time.setText(
            format_timepoints_as_time(self.spin_traj_size.value(), metadata)
        )

    def _reload(self):
        self._update_timepoint_time_labels()
        ct = self._cell_type()
        if not ct:
            self._track_adata = None
            self._track_adata_load_error = None
            self._refresh_buttons()
            return

        # ── Fast synchronous phase (main thread) ──────────────────────────
        path = self._track_adata_path(ct)

        # Auto-fill classifier path if .pkl exists
        clf_path = self._track_classifier_path(ct)
        if clf_path and clf_path.exists():
            self.le_apply_clf_path.setText(str(clf_path))
            self.le_pretrained_clf_path.setText(str(clf_path))

        # Auto-fill states path from State Classification output
        states_path = self._state_adata_path(ct)
        if states_path and states_path.exists():
            self.le_apply_states_path.setText(str(states_path))
            self.le_pretrained_states_path.setText(str(states_path))

        self._populate_track_settings(ct)
        self._apply_max_trajectory_size(ct)
        self._check_prerequisites()
        self._update_view_buttons()
        self._update_bp_buttons()

        # Reset adata so buttons reflect "not yet loaded" state immediately.
        self._track_adata = None
        self._track_adata_load_error = None
        self._refresh_buttons()

        if not path or not path.exists():
            return

        # ── Slow async phase (background thread) ──────────────────────────
        # Skip if a previous preload is still in flight (e.g. rapid tab switches).
        if self._preload_bg.is_running():
            return

        def _load():
            import h5py
            import anndata as ad
            if not h5py.is_hdf5(str(path)):
                raise ValueError(f"Not a valid HDF5 file: {path.name}")
            return ad.read_h5ad(str(path))

        def _on_done(result):
            self._track_adata = result
            self._refresh_buttons()
            self._sync_track_cluster_combo(result)

        def _on_failed(err):
            print(f"[BEHAV3D] Track adata failed to load: {err}")
            self._track_adata_load_error = err
            self._refresh_buttons()

        self._preload_bg.run(
            fn=_load,
            inject_progress=False,
            on_done=_on_done,
            on_failed=_on_failed,
        )

    def _autofill_paths(self, ct: str):
        """Lightweight path auto-fill called on tab switch.

        Only checks whether the states h5ad and classifier pkl exist and fills
        the corresponding line-edit fields.  Avoids the full _reload() cost
        (widget rebuild, prerequisite check, async h5ad load) that is
        unnecessary when switching tabs after a previous _reload already ran.
        """
        if not ct:
            return
        clf_path = self._track_classifier_path(ct)
        if clf_path and clf_path.exists():
            self.le_apply_clf_path.setText(str(clf_path))
            self.le_pretrained_clf_path.setText(str(clf_path))
        states_path = self._state_adata_path(ct)
        if states_path and states_path.exists():
            self.le_apply_states_path.setText(str(states_path))
            self.le_pretrained_states_path.setText(str(states_path))

    # ── Path helpers ─────────────────────────────────────────────────────

    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

    def _track_adata_path(self, ct: str) -> Optional[Path]:
        from behav3d.analysis.behavior.track.utils import get_dtaidistance_track_trajectories_filename
        out = self._out_dir()
        if not out:
            return None
        d = out / "analysis" / ct / "behavorial_trajectories"
        canonical = d / get_dtaidistance_track_trajectories_filename(ct)
        if canonical.exists():
            return canonical
        if d.exists():
            for f in sorted(d.glob("*.h5ad")):
                return f
        return None

    def _track_diagnostics_qc_dirs(self, traj_dir: Optional[Path]) -> list:
        """Both possible diagnostics output folders: dtaidistance's shared
        quality_control/, and the original-BEHAV3D method's own
        original_behav3d/ folder (root — Create Diagnostics writes there)."""
        if not traj_dir:
            return []
        return [
            traj_dir / "quality_control",
            traj_dir / "original_behav3d",
        ]

    def _sync_track_cluster_combo(self, adata_tracks) -> None:
        """Repopulate the 'Color by' combo from actual obs columns in adata_tracks."""
        from behav3d.napari._rename_dialog import _track_cluster_col
        primary = _track_cluster_col(adata_tracks)
        if primary is None:
            return
        options = [primary]
        original_col = f"{primary}_original"
        if original_col in adata_tracks.obs.columns:
            options.append(original_col)
        current = self.combo_track_color_by.currentText()
        self.combo_track_color_by.blockSignals(True)
        self.combo_track_color_by.clear()
        self.combo_track_color_by.addItems(options)
        if current in options:
            self.combo_track_color_by.setCurrentText(current)
        self.combo_track_color_by.blockSignals(False)

    def _track_classifier_path(self, ct: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return out / "analysis" / ct / "behavorial_trajectories" / f"classifier_{ct}.pkl"

    def _state_adata_path(self, ct: str) -> Optional[Path]:
        """Path to the behavioral states h5ad produced by State Classification."""
        out = self._out_dir()
        if not out:
            return None
        return (
            out / "analysis" / ct / "behavioral_states"
            / f"BEHAV3D_{ct}_behavioral_states.h5ad"
        )

    def _cell_type(self) -> str:
        return self._get_cell_type() if self._get_cell_type else ""

    def _sample(self) -> str:
        txt = self.sample_combo.currentText()
        return "" if txt.startswith("—") else txt

    def _log(self, msg: str):
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        formatted = f"[{ts}] {msg}"
        self.log_edit.append(formatted)
        print(formatted)

    # ── Button state management ──────────────────────────────────────────

    def _refresh_buttons(self):
        has_adata = self._track_adata is not None
        self.btn_rename_track.setEnabled(has_adata)
        self.btn_train_track.setEnabled(has_adata)
        self.btn_apply_track.setEnabled(has_adata)
        self.btn_track_proportions.setEnabled(has_adata)
        if has_adata:
            self.rename_track_status.setText(
                f"✅ Track adata loaded: {self._track_adata.n_obs} rows."
            )
        elif self._track_adata_load_error:
            self.rename_track_status.setText(
                f"⚠ Track adata failed to load — re-run clustering. ({self._track_adata_load_error})"
            )
        else:
            self.rename_track_status.setText("ℹ Run clustering first to enable renaming.")
        self._refresh_track_proportion_group_cols()
        self._refresh_contact_columns()

    def _track_proportion_candidate_columns(self):
        cols = []
        added = set()
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        if md is not None and hasattr(md, "columns"):
            for col in md.columns:
                if col in ("exp_nr", "well") or col.endswith("_line_condition"):
                    cols.append(col)
                    added.add(col)
        if self._track_adata is not None:
            for col in self._track_adata.obs.columns:
                if col in added or col.startswith("_") or col in _TECHNICAL_OBS_COLS:
                    continue
                if col in ("exp_nr", "well", "origin_cell_type") or col.endswith("_line_condition"):
                    cols.append(col)
                    added.add(col)
        return cols

    def _refresh_track_proportion_group_cols(self):
        candidate_cols = self._track_proportion_candidate_columns()
        self.list_track_proportion_group_cols.clear()
        for col in candidate_cols:
            self.list_track_proportion_group_cols.addItem(col)
        _fit_list_widget_height(self.list_track_proportion_group_cols)

        self.list_contact_group_cols.clear()
        for col in candidate_cols:
            self.list_contact_group_cols.addItem(col)
        _fit_list_widget_height(self.list_contact_group_cols)

        self.list_track_comparison_group_cols.clear()
        for col in candidate_cols:
            self.list_track_comparison_group_cols.addItem(col)
        _fit_list_widget_height(self.list_track_comparison_group_cols)

        for combo in (
            self.combo_track_proportion_group_x, self.combo_track_proportion_group_y,
            self.combo_contact_group_x, self.combo_contact_group_y,
            self.combo_track_comparison_group_x,
        ):
            prev = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItem("(none)")
            combo.addItems(candidate_cols)
            combo.setCurrentText(prev if prev in (["(none)"] + candidate_cols) else "(none)")
            combo.blockSignals(False)
        self._refresh_track_proportion_group_x_levels()
        self._refresh_track_proportion_group_y_levels()
        self._refresh_contact_group_x_levels()
        self._refresh_contact_group_y_levels()
        self._refresh_track_comparison_group_x_levels()

        prev_cond = self.combo_track_comparison_condition_col.currentText()
        self.combo_track_comparison_condition_col.blockSignals(True)
        self.combo_track_comparison_condition_col.clear()
        self.combo_track_comparison_condition_col.addItems(candidate_cols)
        if prev_cond in candidate_cols:
            self.combo_track_comparison_condition_col.setCurrentText(prev_cond)
        self.combo_track_comparison_condition_col.blockSignals(False)
        self._sync_track_comparison_group_y_text()
        self._refresh_track_comparison_group_levels()

        # Pairing column for the paired contact-duration test: "sample_name" is always offered
        # first (the common case, matching a typical R paired-t-test workflow) even though it's
        # not itself one of the condition-like candidate_cols above.
        pairing_options = ["sample_name"] + [c for c in candidate_cols if c != "sample_name"]
        prev_pairing = self.combo_duration_pairing_col.currentText()
        self.combo_duration_pairing_col.blockSignals(True)
        self.combo_duration_pairing_col.clear()
        self.combo_duration_pairing_col.addItems(pairing_options)
        self.combo_duration_pairing_col.setCurrentText(
            prev_pairing if prev_pairing in pairing_options else "sample_name"
        )
        self.combo_duration_pairing_col.blockSignals(False)

    def _on_track_proportion_group_x_changed(self, _text):
        self._refresh_track_proportion_group_x_levels()

    def _on_track_proportion_group_y_changed(self, _text):
        self._refresh_track_proportion_group_y_levels()

    def _refresh_track_proportion_group_x_levels(self):
        col = self.combo_track_proportion_group_x.currentText()
        col = None if col in ("", "(none)") else col
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._track_adata,
            h5ad_path=self._track_adata_path(ct) if ct else None,
        ) if col else []
        self.group_selector_track_proportion_x.set_items(levels)

    def _refresh_track_proportion_group_y_levels(self):
        col = self.combo_track_proportion_group_y.currentText()
        col = None if col in ("", "(none)") else col
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._track_adata,
            h5ad_path=self._track_adata_path(ct) if ct else None,
        ) if col else []
        self.group_selector_track_proportion_y.set_items(levels)

    def _on_track_proportion_group_x_conditions_toggled(self, checked: bool):
        self.group_selector_track_proportion_x.setVisible(checked)
        if checked:
            self._refresh_track_proportion_group_x_levels()

    def _on_track_proportion_group_y_conditions_toggled(self, checked: bool):
        self.group_selector_track_proportion_y.setVisible(checked)
        if checked:
            self._refresh_track_proportion_group_y_levels()

    def _on_contact_group_x_changed(self, _text):
        self._refresh_contact_group_x_levels()

    def _on_contact_group_y_changed(self, _text):
        self._refresh_contact_group_y_levels()

    def _refresh_contact_group_x_levels(self):
        col = self.combo_contact_group_x.currentText()
        col = None if col in ("", "(none)") else col
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._track_adata,
            h5ad_path=self._track_adata_path(ct) if ct else None,
        ) if col else []
        self.group_selector_contact_x.set_items(levels)

    def _refresh_contact_group_y_levels(self):
        col = self.combo_contact_group_y.currentText()
        col = None if col in ("", "(none)") else col
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._track_adata,
            h5ad_path=self._track_adata_path(ct) if ct else None,
        ) if col else []
        self.group_selector_contact_y.set_items(levels)

    def _on_contact_group_x_conditions_toggled(self, checked: bool):
        self.group_selector_contact_x.setVisible(checked)
        if checked:
            self._refresh_contact_group_x_levels()

    def _on_contact_group_y_conditions_toggled(self, checked: bool):
        self.group_selector_contact_y.setVisible(checked)
        if checked:
            self._refresh_contact_group_y_levels()

    def _sync_track_comparison_group_y_text(self):
        self.line_track_comparison_group_y.setText(self.combo_track_comparison_condition_col.currentText())

    def _on_track_comparison_condition_col_changed(self, _text):
        self._sync_track_comparison_group_y_text()
        self._refresh_track_comparison_group_levels()

    def _refresh_track_comparison_group_levels(self):
        col = self.combo_track_comparison_condition_col.currentText()
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._track_adata,
            h5ad_path=self._track_adata_path(ct) if ct else None,
        )
        self.group_selector_track_comparison.set_items(levels)

    def _on_track_comparison_group_conditions_toggled(self, checked: bool):
        self.group_selector_track_comparison.setVisible(checked)
        if checked:
            self._refresh_track_comparison_group_levels()

    def _on_track_comparison_group_x_changed(self, _text):
        self._refresh_track_comparison_group_x_levels()

    def _refresh_track_comparison_group_x_levels(self):
        col = self.combo_track_comparison_group_x.currentText()
        col = None if col in ("", "(none)") else col
        ct = self._cell_type()
        levels = _condition_levels_for_column(
            col,
            metadata_loader=self.metadata_loader,
            adata=self._track_adata,
            h5ad_path=self._track_adata_path(ct) if ct else None,
        ) if col else []
        self.group_selector_track_comparison_x.set_items(levels)

    def _on_track_comparison_group_x_conditions_toggled(self, checked: bool):
        self.group_selector_track_comparison_x.setVisible(checked)
        if checked:
            self._refresh_track_comparison_group_x_levels()

    def _track_features_csv_path(self, ct: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        base = out / "analysis" / ct / "track_features"
        filtered = base / f"BEHAV3D_{ct}_combined_track_features_filtered.csv"
        if filtered.exists():
            return filtered
        return base / f"BEHAV3D_{ct}_combined_track_features.csv"

    def _contact_target_cell_type(self) -> Optional[str]:
        from behav3d.analysis.behavior.track.contact_grouping import contact_col_target_cell_type
        contact_col = self.combo_contact_col.currentText()
        if not contact_col:
            return None
        try:
            return contact_col_target_cell_type(contact_col)
        except ValueError:
            return None

    def _target_class_availability(self):
        """(available, warning_text) for the currently selected target-classification option.
        ``available`` is always True when the checkbox is unchecked (nothing to block)."""
        import pandas as pd
        from behav3d.analysis.behavior.track.contact_grouping import touching_column_name
        if not self.chk_use_target_class.isChecked():
            return True, ""

        target_ct = self._contact_target_cell_type()
        if not target_ct:
            return False, "Select a contact column first."

        ct = self._cell_type()
        csv_path = self._track_features_csv_path(ct) if ct else None
        touching_col = touching_column_name(target_ct)
        has_touching = False
        if csv_path and csv_path.exists():
            try:
                has_touching = touching_col in pd.read_csv(csv_path, nrows=0).columns
            except Exception:
                has_touching = False
        if not has_touching:
            return False, (
                f"Per-cell contact identity isn't available for '{self.combo_contact_col.currentText()}' "
                f"('{touching_col}' column missing). Recompute contact features with the pixel/mask-based "
                f"method to enable this."
            )

        if self.combo_target_class_source.currentData() == "track":
            if self._track_adata_path(target_ct) is None:
                return False, (
                    f"No track classification found for '{target_ct}'. Run Track Classification for "
                    f"'{target_ct}' first."
                )
        else:
            state_path = self._state_adata_path(target_ct)
            if not state_path or not state_path.exists():
                return False, (
                    f"No behavioral state classification found for '{target_ct}'. Run State "
                    f"Classification for '{target_ct}' first."
                )
        return True, ""

    def _sync_contact_target_class_controls(self):
        use_target = self.chk_use_target_class.isChecked()
        self.combo_target_class_source.setVisible(use_target)
        self.combo_target_state_col.setVisible(
            use_target and self.combo_target_class_source.currentData() == "state"
        )
        available, warning = self._target_class_availability()
        if use_target and warning:
            self.label_target_class_warning.setText(warning)
            self.label_target_class_warning.show()
        else:
            self.label_target_class_warning.hide()
        contact_cols_present = self.combo_contact_col.count() > 0
        self.btn_contact_analysis.setEnabled(
            contact_cols_present
            and self._track_adata is not None
            and not (use_target and not available)
        )
        # Duration-by-class comparison is meaningless without a target classification to split by.
        self.btn_duration_comparison.setEnabled(
            contact_cols_present
            and self._track_adata is not None
            and use_target
            and available
        )

    def _on_duration_test_mode_changed(self, *_args):
        is_paired = self.combo_duration_test_mode.currentData() == "paired"
        self.combo_duration_pairing_col.setVisible(is_paired)

    def _refresh_contact_columns(self):
        import pandas as pd
        from behav3d.analysis.behavior.track.contact_grouping import list_available_contact_columns
        ct = self._cell_type()
        contact_cols = []
        if ct:
            csv_path = self._track_features_csv_path(ct)
            if csv_path and csv_path.exists():
                try:
                    contact_cols = list_available_contact_columns(pd.read_csv(csv_path, nrows=0))
                except Exception:
                    contact_cols = []
        prev = self.combo_contact_col.currentText()
        self.combo_contact_col.blockSignals(True)
        self.combo_contact_col.clear()
        self.combo_contact_col.addItems(contact_cols)
        if prev in contact_cols:
            self.combo_contact_col.setCurrentText(prev)
        self.combo_contact_col.blockSignals(False)
        self.combo_contact_col.setEnabled(len(contact_cols) > 0)
        self._sync_contact_target_class_controls()
        state_adata_path = self._state_adata_path(ct) if ct else None
        has_states = bool(state_adata_path and state_adata_path.exists())
        self.btn_contact_state_shift.setEnabled(
            len(contact_cols) > 0 and self._track_adata is not None and has_states
        )
        self.btn_track_contact_overview.setEnabled(
            len(contact_cols) > 0 and self._track_adata is not None and has_states
        )

    def _update_view_buttons(self):
        ct = self._cell_type()
        if not ct:
            for btn in (
                self.btn_view_exemplars,
                self.btn_view_diagnostics, self.btn_view_track_proportions,
                self.btn_view_window_transitions,
                self.btn_view_track_condition_comparison, self.btn_view_contact_analysis,
                self.btn_view_contact_state_shift, self.btn_view_track_contact_overview,
                self.btn_view_duration_comparison,
            ):
                btn.setEnabled(False)
            return
        out = self._out_dir()
        traj_dir = (out / "analysis" / ct / "behavorial_trajectories") if out else None
        self.btn_view_exemplars.setEnabled(
            bool(traj_dir and any(traj_dir.glob("exemplar_tracks*.pdf")))
        )
        qc_dirs = self._track_diagnostics_qc_dirs(traj_dir)
        self.btn_view_diagnostics.setEnabled(
            any(any(qc_dir.glob("*.pdf")) for qc_dir in qc_dirs)
        )
        proportions_dir = traj_dir / "behavior_proportions" if traj_dir else None
        self.btn_view_track_proportions.setEnabled(
            bool(
                proportions_dir
                and proportions_dir.exists()
                and any(proportions_dir.glob("*.pdf"))
            )
        )
        window_transitions_dir = traj_dir / "window_transitions" if traj_dir else None
        self.btn_view_window_transitions.setEnabled(
            bool(
                window_transitions_dir
                and window_transitions_dir.exists()
                and any(window_transitions_dir.glob("*.pdf"))
            )
        )
        comparisons_dir = traj_dir / "behavior_comparisons" if traj_dir else None
        self.btn_view_track_condition_comparison.setEnabled(
            bool(
                comparisons_dir
                and comparisons_dir.exists()
                and any(comparisons_dir.glob("condition_comparison_*.pdf"))
            )
        )
        contact_dir = traj_dir / "contact_analysis" if traj_dir else None
        self.btn_view_contact_analysis.setEnabled(
            bool(
                contact_dir
                and contact_dir.exists()
                and any(contact_dir.glob("*/*.pdf"))
            )
        )
        self.btn_view_contact_state_shift.setEnabled(
            bool(
                contact_dir
                and contact_dir.exists()
                and any(contact_dir.glob("*/contact_state_shift.pdf"))
            )
        )
        self.btn_view_track_contact_overview.setEnabled(
            bool(
                contact_dir
                and contact_dir.exists()
                and any(contact_dir.glob("*/track_contact_overview.pdf"))
            )
        )
        self.btn_view_duration_comparison.setEnabled(
            bool(
                contact_dir
                and contact_dir.exists()
                and any(contact_dir.glob("*/contact_duration_comparison.pdf"))
            )
        )

    def _update_bp_buttons(self):
        ct = self._cell_type()
        track_path = self._track_adata_path(ct) if ct else None
        state_path = self._state_adata_path(ct) if ct else None
        has_track = bool(track_path and track_path.exists())
        has_state = bool(state_path and state_path.exists())
        self.btn_show_track_bp.setEnabled(has_track and has_state)
        self.btn_export_track_bp.setEnabled(has_track and has_state)

    # ── Click handlers ───────────────────────────────────────────────────

    def _on_run_cluster(self):
        """Dispatch to original or dtaidistance based on current mode."""
        if self.chk_use_original.isChecked():
            self._on_run_original()
        else:
            self._on_run_track()

    def _on_run_track(self):
        ct = self._cell_type()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        self._log(f"▶ Running track clustering (dtaidistance) for '{ct}'…")
        self._dispatch_track_cluster(ct)

    def _dispatch_track_cluster(self, ct: str, extra_callbacks=None):
        out = self._out_dir()
        is_bouts = self.combo_trajectory_basis.currentText() == "bouts"

        if is_bouts:
            params = {
                "output_dir": str(out) if out else "",
                "cell_type": ct,
                "behavioral_trajectory_size": int(self.spin_traj_size.value()),
                "trajectory_trim_mode": self.combo_trim.currentText(),
                "split_long_tracks": self.chk_split_long_tracks.isChecked(),
                "use_fractions": self.chk_bouts_use_fractions.isChecked(),
                "use_bout_stats": self.chk_bouts_use_bout_stats.isChecked(),
                "use_transitions": self.chk_bouts_use_transitions.isChecked(),
                "use_bigrams": self.chk_bouts_use_ngrams.isChecked(),
                "use_trigrams": self.chk_bouts_use_ngrams.isChecked(),
                "do_pca": self.chk_bouts_do_pca.isChecked(),
                "use_clr_transform": self.chk_bouts_use_clr.isChecked(),
                "log_transform_bout_lengths": self.chk_bouts_log_bout_length.isChecked(),
                "do_block_scaling": self.chk_bouts_block_scaling.isChecked(),
                "drop_highly_correlated": self.chk_bouts_drop_redundant.isChecked(),
                "drop_low_variance": self.chk_bouts_drop_redundant.isChecked(),
                "clustering_method": self.combo_clustering_method.currentText(),
                "n_clusters": int(self.spin_n_clusters.value()),
                "agglomerative_linkage": self.combo_bouts_linkage.currentText(),
                "n_neighbors": int(self.spin_leiden_neighbors.value()),
                "leiden_resolution": float(self.spin_leiden_resolution.value()),
                "cluster_key": "ClusterID",
                "plot_results": True,
                "plot_exemplars": self.chk_bouts_plot_exemplars.isChecked(),
                "random_state": int(self.spin_seed.value()),
                # Share the DTW basis's output folder so both bases' diagnostics/exemplar
                # PDFs and the canonical model .h5ad land in the same place on disk.
                "output_subdir_name": "behavorial_trajectories",
            }
        else:
            params = {
                "output_dir": str(out) if out else "",
                "cell_type": ct,
                "behavioral_trajectory_size": int(self.spin_traj_size.value()),
                "n_clusters": int(self.spin_n_clusters.value()),
                "random_state": int(self.spin_seed.value()),
                "linkage": self.combo_linkage.currentText(),
                "clustering_method": self.combo_clustering_method.currentText(),
                "leiden_n_neighbors": int(self.spin_leiden_neighbors.value()),
                "leiden_resolution": float(self.spin_leiden_resolution.value()),
                "trajectory_trim_mode": self.combo_trim.currentText(),
                "split_long_tracks": self.chk_split_long_tracks.isChecked(),
                "parallel": self.chk_parallel.isChecked(),
                "save_distance_matrix": self.chk_save_dist.isChecked(),
                "plot_results": True,
            }

        on_done_ext = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_ext = extra_callbacks.get("on_failed") if extra_callbacks else None

        def _run(**kw):
            _traj_dir = out / "analysis" / ct / "behavorial_trajectories"
            if _traj_dir.exists():
                rmtree_ignore_missing(_traj_dir)
            if is_bouts:
                from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
                from behav3d.analysis.behavior.track.bouts import run_state_based_analysis
                from behav3d.analysis.behavior.track.utils import (
                    get_dtaidistance_track_trajectories_filename,
                )
                result = run_state_based_analysis(state_col=FULL_STATE_COL, verbose=True, **params)
                # Write to the same canonical path the dtaidistance branch uses so rename,
                # backprojection, exemplar training, contact analysis, etc. keep working
                # unchanged regardless of which basis produced the model.
                model_path = _traj_dir / get_dtaidistance_track_trajectories_filename(ct)
                model_path.parent.mkdir(parents=True, exist_ok=True)
                result.write(model_path, compression="gzip")
                return result
            from behav3d.analysis.behavior.track.state_dtw import (
                run_categorical_dtaidistance_trajectory_clustering,
            )
            return run_categorical_dtaidistance_trajectory_clustering(**params, verbose=True)

        def _done(r):
            if is_bouts:
                self._log(f"✅ Bout/proportion clustering done for '{ct}'.")
            else:
                umap_error = (r.uns.get("visualization", {}) or {}).get("umap_error")
                if umap_error:
                    self._log(f"⚠ Track clustering done for '{ct}', but UMAP was skipped: {umap_error}")
                else:
                    self._log(f"✅ Track clustering done for '{ct}'.")
            self._persist_track_cfg(ct)
            self._reload()
            self._notify_results()

            # The auto exemplar-overview step below relies on save_dtaidistance_exemplar_overview /
            # save_dtaidistance_medoid_overview, which need a DTW pairwise distance matrix — not
            # applicable to a bouts-derived model. The bouts pipeline already writes its own
            # diagnostics/exemplar PDFs inline (via plot_results/plot_exemplars above).
            if is_bouts:
                if on_done_ext:
                    on_done_ext(r)
                return

            _track_adata = r
            _n_per = int(self.spin_n_per_cluster.value())
            _seed = int(self.spin_seed.value())
            _log = self._log

            def _run_overview(**kw):
                from behav3d.analysis.behavior.track.state_dtw import (
                    save_dtaidistance_exemplar_overview,
                    save_dtaidistance_medoid_overview,
                )
                from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
                _raw_dir = _resolve_dtaidistance_paths(str(out), ct)["quality_control_outfolder"] / "raw"
                _raw_dir.mkdir(parents=True, exist_ok=True)
                _result = save_dtaidistance_exemplar_overview(
                    _track_adata,
                    output_dir=str(out),
                    cell_type=ct,
                    n_per_cluster=_n_per,
                    random_state=_seed,
                    outfolder=_raw_dir,
                    verbose=True,
                )
                try:
                    save_dtaidistance_medoid_overview(
                        _track_adata,
                        output_dir=str(out),
                        cell_type=ct,
                        outfolder=_raw_dir,
                        verbose=True,
                    )
                except Exception as _exc:
                    _log(f"⚠ Could not generate medoid overview: {_exc}")
                return _result

            def _overview_done(_):
                _log("✅ Exemplar overview done.")
                QTimer.singleShot(0, self._update_view_buttons)

            def _start_overview():
                self._bg.run(
                    fn=_run_overview,
                    desc=f"Exemplar overview ({ct})…",
                    progress_row=self.progress_row,
                    buttons=[],
                    viewer=self.viewer,
                    inject_progress=False,
                    on_done=_overview_done,
                    on_failed=lambda e: self._log(f"⚠ Exemplar overview failed: {e}"),
                )
            QTimer.singleShot(0, _start_overview)
            if on_done_ext:
                on_done_ext(r)

        def _fail(e):
            self._log(f"❌ Track clustering failed: {e}")
            if on_fail_ext:
                on_fail_ext(e)

        self._bg.run(
            fn=_run,
            desc=f"Track clustering ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_run_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=_done,
            on_failed=_fail,
        )

    def run_track_clustering(self, interactive=True, extra_callbacks=None):
        """Called from queue runner."""
        ct = self._cell_type()
        if not ct:
            err = "No cell type selected."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        if self._bg.is_running():
            err = "BackgroundOperation already running."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        self._dispatch_track_cluster(ct, extra_callbacks=extra_callbacks)

    def _on_run_original(self):
        ct = self._cell_type()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        out = self._out_dir()
        self._log(f"▶ Running original BEHAV3D feature-based DTW for '{ct}'…")
        n_clust = int(self.spin_n_clusters.value())
        traj_size = int(self.spin_traj_size.value())
        n_neigh = int(self.spin_umap_neighbors.value())
        min_dist = float(self.spin_umap_min_dist.value())

        def _run(**kw):
            import shutil
            _traj_dir = out / "analysis" / ct / "behavorial_trajectories"
            if _traj_dir.exists():
                rmtree_ignore_missing(_traj_dir)
            from behav3d.analysis.behavior.track.feature_dtw import run_tcell_analysis
            return run_tcell_analysis(
                output_dir=str(out) if out else "",
                cell_type=ct,
                nr_of_clusters=n_clust,
                min_track_length=traj_size,
                max_track_length=traj_size,
                umap_n_neighbors=n_neigh,
                umap_minimal_distance=min_dist,
                feature_scaling_preset="original_behav3d",
                output_subdir_name="behavorial_trajectories/original_behav3d/raw",
            )

        def _done(_):
            from behav3d.analysis.behavior.track.feature_dtw import _create_original_behav3d_adata
            try:
                _create_original_behav3d_adata(str(out), ct)
            except Exception as e:
                self._log(f"⚠ Could not create h5ad from original BEHAV3D results: {e}")
            self._log(f"✅ Original BEHAV3D clustering done for '{ct}'.")
            self._persist_track_cfg(ct)
            self._reload()
            self._notify_results()

        self._bg.run(
            fn=_run,
            desc=f"Original BEHAV3D DTW ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_run_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=_done,
            on_failed=lambda e: self._log(f"❌ Original BEHAV3D failed: {e}"),
        )

    def _on_rename_track(self):
        ct = self._cell_type()
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        path = self._track_adata_path(ct)
        dlg = RenameClusterDialog(
            mode="track",
            model_adata=self._track_adata,
            adata_path=path,
            parent=self,
        )
        dlg.clusters_renamed.connect(
            lambda mapping: self._log(f"✅ Track clusters renamed: {mapping}")
        )
        if dlg.exec_() == QDialog.Accepted:
            self._refresh_buttons()
            self._update_view_buttons()
            if self._track_adata is not None:
                self._sync_track_cluster_combo(self._track_adata)
                self._regenerate_track_reports_after_rename(ct)

    def _regenerate_track_reports_after_rename(self, ct: str):
        """Re-run diagnostics + proportion plots so plots reflect the just-renamed labels/colors."""
        if self._track_adata is None or self._bg.is_running():
            return
        out = self._out_dir()
        track_adata = self._track_adata
        method = (track_adata.uns.get("dtai_trajectory_clustering", {}) or {}).get("method")
        if method == "original_behav3d_feature_dtw":
            # Known limitation: renaming feature-DTW clusters via this generic dialog updates
            # adata.uns["classification"], but diagnostics for this method read colors/order from
            # a separate YAML store (see `_save_feature_dtw_quality_control`) that isn't refreshed
            # here, so its own reports can go stale after a rename. Regenerate manually via
            # "Create diagnostics" if needed.
            return

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import save_dtaidistance_diagnostics
            from behav3d.analysis.behavior.track.visualization.plots.reports import (
                generate_track_clustering_report_pdfs,
                save_track_class_proportions_by_sample_plot,
            )
            from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
            from behav3d.napari._rename_dialog import _track_cluster_col
            cluster_col = _track_cluster_col(track_adata) or "ClusterID"
            paths = _resolve_dtaidistance_paths(str(out) if out else "", ct)
            if method == "bouts_feature_clustering":
                # adata.X is a per-track feature matrix here, not a pairwise DTW
                # distance matrix - regenerate via the same generator bouts.py
                # itself uses at clustering time instead of the DTW-only path.
                diag = generate_track_clustering_report_pdfs(
                    adata_tracks=track_adata,
                    outfolder=paths["clustering_outfolder"],
                    cluster_key=cluster_col,
                    verbose=True,
                )
            else:
                diag = save_dtaidistance_diagnostics(
                    adata_tracks=track_adata,
                    output_dir=str(out) if out else "",
                    cell_type=ct,
                    verbose=True,
                )
            prop = save_track_class_proportions_by_sample_plot(
                track_adata,
                paths["behavior_proportions_outfolder"],
                sample_col="sample_name",
                class_col=cluster_col,
                verbose=True,
            )
            result = {"diagnostics": diag, "proportions": prop}
            if "trajectory_window_id" in track_adata.obs.columns:
                # Only meaningful when 'Divide long tracks' was used - most runs
                # don't have this column, so skip silently rather than erroring
                # the whole refresh for the common case. A real failure here is
                # caught locally so it doesn't take down the reports above,
                # which already succeeded by this point.
                try:
                    from behav3d.analysis.behavior.track.visualization.plots.window_transitions import (
                        save_window_transition_report,
                    )
                    result["window_transitions"] = save_window_transition_report(
                        track_adata,
                        output_dir=str(out) if out else "",
                        cell_type=ct,
                        cluster_key=cluster_col,
                        verbose=True,
                    )
                except Exception as exc:
                    result["window_transitions_error"] = str(exc)
            return result

        self._bg.run(
            fn=_run,
            desc=f"Refreshing track reports after rename ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_rename_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log("✅ Track reports refreshed after rename."),
                self._update_view_buttons(),
            ),
            on_failed=lambda e: self._log(f"⚠ Could not refresh track reports after rename: {e}"),
        )

    def _on_train_track(self):
        ct = self._cell_type()
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        out = self._out_dir()
        n_est = int(self.spin_track_n_est.value())
        max_depth_val = int(self.spin_track_max_depth.value())
        max_depth = None if max_depth_val == 0 else max_depth_val
        test_size = float(self.spin_track_test_pct.value()) / 100.0
        min_leaf = int(self.spin_track_min_leaf.value())
        min_split = int(self.spin_track_min_split.value())
        max_features = self.combo_track_max_features.currentText()
        n_jobs = int(self.spin_track_n_jobs.value())
        random_state = int(self.spin_seed.value())
        self._log(f"▶ Training track classifier for '{ct}'…")

        # Capture adata reference to avoid closure over self
        track_adata = self._track_adata

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import (
                train_dtaidistance_trajectory_classifier,
            )
            return train_dtaidistance_trajectory_classifier(
                output_dir=str(out) if out else "",
                cell_type=ct,
                model_adata=track_adata,
                classifier_n_estimators=n_est,
                classifier_max_depth=max_depth,
                validation_test_size=test_size,
                classifier_min_samples_leaf=min_leaf,
                classifier_min_samples_split=min_split,
                classifier_max_features=max_features,
                classifier_n_jobs=n_jobs,
                random_state=random_state,
                verbose=True,
            )

        def _done(r):
            self._log(f"✅ Track classifier trained for '{ct}'.")
            self._persist_track_cfg(ct)
            # Auto-fill classifier path
            clf_path = None
            if isinstance(r, dict):
                clf_path = r.get("classifier_path")
            if not clf_path:
                p = self._track_classifier_path(ct)
                clf_path = str(p) if p and p.exists() else None
            if clf_path:
                self.le_apply_clf_path.setText(clf_path)
                self.le_pretrained_clf_path.setText(clf_path)
            self._update_view_buttons()
            self._notify_results()

        self._bg.run(
            fn=_run,
            desc=f"Training track classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_train_track, self.btn_apply_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=_done,
            on_failed=lambda e: self._log(f"❌ Train track failed: {e}"),
        )

    def run_train_track(self, interactive=True, extra_callbacks=None):
        """Called from queue runner."""
        ct = self._cell_type()
        if not ct or self._track_adata is None:
            err = "No track adata available — run track clustering first."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        if self._bg.is_running():
            err = "BackgroundOperation already running."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        on_done_cb = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_cb = extra_callbacks.get("on_failed") if extra_callbacks else None
        out = self._out_dir()
        n_est = int(self.spin_track_n_est.value())
        max_depth_val = int(self.spin_track_max_depth.value())
        max_depth = None if max_depth_val == 0 else max_depth_val
        test_size = float(self.spin_track_test_pct.value()) / 100.0
        min_leaf = int(self.spin_track_min_leaf.value())
        min_split = int(self.spin_track_min_split.value())
        max_features = self.combo_track_max_features.currentText()
        n_jobs = int(self.spin_track_n_jobs.value())
        random_state = int(self.spin_seed.value())
        track_adata = self._track_adata

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import (
                train_dtaidistance_trajectory_classifier,
            )
            return train_dtaidistance_trajectory_classifier(
                output_dir=str(out) if out else "",
                cell_type=ct,
                model_adata=track_adata,
                classifier_n_estimators=n_est,
                classifier_max_depth=max_depth,
                validation_test_size=test_size,
                classifier_min_samples_leaf=min_leaf,
                classifier_min_samples_split=min_split,
                classifier_max_features=max_features,
                classifier_n_jobs=n_jobs,
                random_state=random_state,
                verbose=True,
            )

        def _done(r):
            self._log(f"✅ Track classifier trained for '{ct}'.")
            self._persist_track_cfg(ct)
            self._update_view_buttons()
            self._notify_results()
            if on_done_cb:
                on_done_cb(r)

        def _fail(e):
            self._log(f"❌ Train track failed: {e}")
            if on_fail_cb:
                on_fail_cb(e)

        self._bg.run(
            fn=_run,
            desc=f"Training track classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_train_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=_done,
            on_failed=_fail,
        )

    def _on_apply_track(self):
        ct = self._cell_type()
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        clf_path_str = self.le_apply_clf_path.text().strip()
        states_path_str = self.le_apply_states_path.text().strip()
        if not clf_path_str or not Path(clf_path_str).exists():
            QMessageBox.warning(self, "No classifier",
                                "Select a valid classifier .pkl file first.")
            return
        if not states_path_str or not Path(states_path_str).exists():
            QMessageBox.warning(self, "No states file",
                                "Select a valid behavioral states .h5ad file first.")
            return
        out = self._out_dir()
        self._log(f"▶ Applying track classifier for '{ct}'…")
        clf_path = Path(clf_path_str)
        states_path = Path(states_path_str)

        def _run(**kw):
            from behav3d.analysis.behavior.track.bouts import apply_track_classifier_to_subtracks
            return apply_track_classifier_to_subtracks(
                output_dir=str(out) if out else "",
                cell_type=ct,
                classifier_artifact_or_path=str(clf_path),
                adata_full_path=str(states_path),
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Applying track classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_apply_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Track classifier applied for '{ct}'."),
                self._reload(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Apply track failed: {e}"),
        )

    def run_apply_track(self, interactive=True, extra_callbacks=None):
        """Called from queue runner."""
        ct = self._cell_type()
        if not ct:
            err = "No cell type selected."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        if self._bg.is_running():
            err = "BackgroundOperation already running."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        clf_path_str = self.le_apply_clf_path.text().strip()
        states_path_str = self.le_apply_states_path.text().strip()
        if not clf_path_str or not Path(clf_path_str).exists():
            err = "Track classifier .pkl not found — browse or train first."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        if not states_path_str or not Path(states_path_str).exists():
            err = "Behavioral states .h5ad not found — check the path."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return

        on_done_cb = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_cb = extra_callbacks.get("on_failed") if extra_callbacks else None
        out = self._out_dir()
        clf_path = Path(clf_path_str)
        states_path = Path(states_path_str)

        def _run(**kw):
            from behav3d.analysis.behavior.track.bouts import apply_track_classifier_to_subtracks
            return apply_track_classifier_to_subtracks(
                output_dir=str(out) if out else "",
                cell_type=ct,
                classifier_artifact_or_path=str(clf_path),
                adata_full_path=str(states_path),
                verbose=True,
            )

        def _done(r):
            self._log(f"✅ Track classifier applied for '{ct}'.")
            self._reload()
            self._notify_results()
            if on_done_cb:
                on_done_cb(r)

        def _fail(e):
            self._log(f"❌ Apply track failed: {e}")
            if on_fail_cb:
                on_fail_cb(e)

        self._bg.run(
            fn=_run,
            desc=f"Applying track classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_apply_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=_done,
            on_failed=_fail,
        )

    def _on_apply_pretrained(self):
        """Apply a pre-trained classifier + pre-existing states file."""
        ct = self._cell_type()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        clf_path_str = self.le_pretrained_clf_path.text().strip()
        states_path_str = self.le_pretrained_states_path.text().strip()
        if not clf_path_str or not Path(clf_path_str).exists():
            QMessageBox.warning(self, "No classifier",
                                "Select a valid classifier .pkl file.")
            return
        if not states_path_str or not Path(states_path_str).exists():
            QMessageBox.warning(self, "No states file",
                                "Select a valid behavioral states .h5ad file.")
            return
        out = self._out_dir()
        clf_path = Path(clf_path_str)
        states_path = Path(states_path_str)
        self._log(f"▶ Applying pretrained classifier for '{ct}'…")

        def _run(**kw):
            from behav3d.analysis.behavior.track.bouts import apply_track_classifier_to_subtracks
            return apply_track_classifier_to_subtracks(
                output_dir=str(out) if out else "",
                cell_type=ct,
                classifier_artifact_or_path=str(clf_path),
                adata_full_path=str(states_path),
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Applying pretrained classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_run_apply_pretrained],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Pretrained classifier applied for '{ct}'."),
                self._reload(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Apply pretrained failed: {e}"),
        )

    def _on_exemplars(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        out = self._out_dir()
        n_per = int(self.spin_n_per_cluster.value())
        make_overview = self.chk_overview_statebars.isChecked()
        make_bp_pdf = self.chk_backproj_pdf.isChecked()
        make_bp_mp4 = self.chk_backproj_mp4.isChecked()
        self._log(f"▶ Creating exemplar PDFs for '{ct}'…")
        track_adata = self._track_adata
        state_adata_path = self._state_adata_path(ct)

        def _run(**kw):
            from pathlib import Path as _Path
            from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
                plot_exemplar_tracks_by_cluster,
                save_exemplar_statebar_track_pdf_per_cluster,
                save_exemplar_statebar_backprojection_pdf,
                save_exemplar_statebar_backprojection_video_per_cluster,
            )
            from matplotlib.backends.backend_pdf import PdfPages
            import matplotlib.pyplot as plt
            import anndata as _ad
            from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
            from behav3d.napari._rename_dialog import _track_cluster_col

            cluster_col = _track_cluster_col(track_adata) or "ClusterID"

            if state_adata_path is None or not state_adata_path.exists():
                raise FileNotFoundError(
                    f"Behavioral states h5ad not found at {state_adata_path}. "
                    "Run State Classification first."
                )
            full_adata = _ad.read_h5ad(str(state_adata_path))

            out_path = _Path(out) if out else _Path(".")
            exemplar_root = (
                out_path / "analysis" / ct / "behavorial_trajectories" / "example_tracks"
            )
            exemplar_root.mkdir(parents=True, exist_ok=True)
            results = {}

            if make_overview:
                fig, _, _ = plot_exemplar_tracks_by_cluster(
                    full_adata,
                    track_adata,
                    n_per_cluster=n_per,
                    sample_key="sample_name",
                    track_key="TrackID",
                    time_key="position_t",
                    state_key=FULL_STATE_COL,
                    cluster_key=cluster_col,
                    tmin_key="position_t_min",
                    tmax_key="position_t_max",
                )
                overview_pdf = exemplar_root / "example_tracks_overview.pdf"
                with PdfPages(overview_pdf) as pdf:
                    pdf.savefig(fig, bbox_inches="tight", dpi=300)
                plt.close(fig)
                results["overview_pdf"] = str(overview_pdf)

                statebar_out = save_exemplar_statebar_track_pdf_per_cluster(
                    adata_full=full_adata,
                    out_dir=exemplar_root,
                    adata_tracks=track_adata,
                    n_per_cluster=n_per,
                    sample_key="sample_name",
                    track_key="TrackID",
                    time_key="position_t",
                    state_key=FULL_STATE_COL,
                    cluster_key=cluster_col,
                    tmin_key="position_t_min",
                    tmax_key="position_t_max",
                )
                results["statebar_pdfs"] = statebar_out

            if make_bp_pdf:
                bp_pdf_out = save_exemplar_statebar_backprojection_pdf(
                    full_adata,
                    output_dir=str(out_path),
                    cell_type=ct,
                    out_dir=exemplar_root / "backprojection",
                    adata_tracks=track_adata,
                    n_per_cluster=n_per,
                    sample_key="sample_name",
                    track_key="TrackID",
                    time_key="position_t",
                    state_key=FULL_STATE_COL,
                    cluster_key=cluster_col,
                    tmin_key="position_t_min",
                    tmax_key="position_t_max",
                    verbose=True,
                )
                results["backprojection_pdf"] = bp_pdf_out

            if make_bp_mp4:
                bp_mp4_out = save_exemplar_statebar_backprojection_video_per_cluster(
                    full_adata,
                    output_dir=str(out_path),
                    cell_type=ct,
                    out_dir=exemplar_root / "backprojection",
                    adata_tracks=track_adata,
                    n_per_cluster=n_per,
                    sample_key="sample_name",
                    track_key="TrackID",
                    time_key="position_t",
                    state_key=FULL_STATE_COL,
                    cluster_key=cluster_col,
                    tmin_key="position_t_min",
                    tmax_key="position_t_max",
                    verbose=True,
                )
                results["backprojection_mp4"] = bp_mp4_out

            return results

        self._bg.run(
            fn=_run,
            desc=f"Exemplar PDFs ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_exemplars],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Exemplar PDFs done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Exemplar PDFs failed: {e}"),
        )

    def _on_diagnostics(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        out = self._out_dir()
        self._log(f"▶ Creating diagnostics for '{ct}'…")
        track_adata = self._track_adata

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import (
                save_dtaidistance_diagnostics,
                save_dtaidistance_exemplar_overview,
                save_dtaidistance_medoid_overview,
            )
            method = (track_adata.uns.get("dtai_trajectory_clustering", {}) or {}).get("method")
            if method == "original_behav3d_feature_dtw":
                from behav3d.analysis.behavior.track.feature_dtw import (
                    _save_feature_dtw_quality_control,
                )
                from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
                result = _save_feature_dtw_quality_control(
                    output_dir=str(out) if out else "",
                    cell_type=ct,
                    proportions_outfolder=_resolve_dtaidistance_paths(
                        str(out) if out else "", ct
                    )["behavior_proportions_outfolder"],
                )
            elif method == "bouts_feature_clustering":
                # adata.X is a per-track feature matrix here, not a pairwise DTW
                # distance matrix - regenerate via the same generator bouts.py
                # itself uses at clustering time instead of the DTW-only path.
                from behav3d.analysis.behavior.track.visualization.plots.reports import (
                    generate_track_clustering_report_pdfs,
                )
                from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
                from behav3d.napari._rename_dialog import _track_cluster_col
                paths = _resolve_dtaidistance_paths(str(out) if out else "", ct)
                cluster_col = _track_cluster_col(track_adata) or "ClusterID"
                result = generate_track_clustering_report_pdfs(
                    adata_tracks=track_adata,
                    outfolder=paths["clustering_outfolder"],
                    cluster_key=cluster_col,
                    verbose=True,
                )
            else:
                result = save_dtaidistance_diagnostics(
                    adata_tracks=track_adata,
                    output_dir=str(out) if out else "",
                    cell_type=ct,
                    verbose=True,
                )
            try:
                save_dtaidistance_exemplar_overview(
                    track_adata,
                    output_dir=str(out) if out else "",
                    cell_type=ct,
                    verbose=True,
                )
            except Exception as _exc:
                print(f"[BEHAV3D] Could not generate exemplar overview: {_exc}")
            try:
                save_dtaidistance_medoid_overview(
                    track_adata,
                    output_dir=str(out) if out else "",
                    cell_type=ct,
                    verbose=True,
                )
            except Exception as _exc:
                print(f"[BEHAV3D] Could not generate medoid overview: {_exc}")
            return result

        self._bg.run(
            fn=_run,
            desc=f"Track diagnostics ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_diagnostics],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: self._on_diagnostics_done(r, ct),
            on_failed=lambda e: self._log(f"❌ Diagnostics failed: {e}"),
        )

    def _on_diagnostics_done(self, result, ct):
        umap_error = (result or {}).get("umap_error")
        if umap_error:
            self._log(f"⚠ Diagnostics done for '{ct}', but UMAP was skipped: {umap_error}")
        else:
            self._log(f"✅ Diagnostics done for '{ct}'.")
        self._update_view_buttons()
        self._notify_results()

    def _on_track_proportions(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        out = self._out_dir()
        self._log(f"▶ Creating track proportion plots for '{ct}'…")
        track_adata = self._track_adata
        selected_cols = [item.text() for item in self.list_track_proportion_group_cols.selectedItems()]
        group_x = self.combo_track_proportion_group_x.currentText()
        group_x = None if group_x in ("", "(none)") else group_x
        group_y = self.combo_track_proportion_group_y.currentText()
        group_y = None if group_y in ("", "(none)") else group_y
        group_x_levels_map = None
        if self.chk_track_proportion_group_x_conditions.isChecked():
            if not group_x or not self.group_selector_track_proportion_x.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one Group X level into each group, or untick "
                    "\"Group conditions\" under Group in X.",
                )
                return
            left = self.group_selector_track_proportion_x.left_items()
            right = self.group_selector_track_proportion_x.right_items()
            group_x_levels_map = {lvl: "+".join(left) for lvl in left}
            group_x_levels_map.update({lvl: "+".join(right) for lvl in right})
        group_y_levels_map = None
        if self.chk_track_proportion_group_y_conditions.isChecked():
            if not group_y or not self.group_selector_track_proportion_y.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one Group Y level into each group, or untick "
                    "\"Group conditions\" under Group in Y.",
                )
                return
            left = self.group_selector_track_proportion_y.left_items()
            right = self.group_selector_track_proportion_y.right_items()
            group_y_levels_map = {lvl: "+".join(left) for lvl in left}
            group_y_levels_map.update({lvl: "+".join(right) for lvl in right})
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None

        def _run(**kw):
            from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
            from behav3d.analysis.behavior.track.visualization.plots.reports import (
                save_track_class_proportions_by_sample_plot,
            )
            from behav3d.core.metadata import merge_condition_columns_into_obs
            from behav3d.napari._rename_dialog import _track_cluster_col
            cluster_col = _track_cluster_col(track_adata) or "ClusterID"
            prop_dir = _resolve_dtaidistance_paths(str(out) if out else "", ct)["behavior_proportions_outfolder"]
            all_cols = selected_cols + [c for c in (group_x, group_y) if c]
            cols_to_merge = [c for c in all_cols if c not in track_adata.obs.columns]
            if cols_to_merge and md is not None:
                merge_condition_columns_into_obs(track_adata, md, cols_to_merge)
            return save_track_class_proportions_by_sample_plot(
                track_adata,
                prop_dir,
                sample_col="sample_name",
                class_col=cluster_col,
                group_cols=selected_cols or None,
                group_x=group_x,
                group_y=group_y,
                group_x_levels_map=group_x_levels_map,
                group_y_levels_map=group_y_levels_map,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Track proportion plots ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_track_proportions],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Track proportion plots done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Track proportion plots failed: {e}"),
        )

    def _on_window_transitions(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if "trajectory_window_id" not in self._track_adata.obs.columns:
            QMessageBox.warning(
                self, "No sub-track windows",
                "This report needs tracks split into windows. Enable 'Divide long "
                "tracks' in Step 1's Advanced Configuration and re-run track clustering.",
            )
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        out = self._out_dir()
        self._log(f"▶ Creating window transition Sankey for '{ct}'…")
        track_adata = self._track_adata

        def _run(**kw):
            from behav3d.analysis.behavior.track.visualization.plots.window_transitions import (
                save_window_transition_report,
            )
            from behav3d.napari._rename_dialog import _track_cluster_col
            cluster_col = _track_cluster_col(track_adata) or "ClusterID"
            return save_window_transition_report(
                track_adata,
                output_dir=str(out) if out else "",
                cell_type=ct,
                cluster_key=cluster_col,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Window transition Sankey ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_window_transitions],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Window transition Sankey done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Window transition Sankey failed: {e}"),
        )

    def _on_track_condition_comparison(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        condition_col = self.combo_track_comparison_condition_col.currentText()
        group_cols = [item.text() for item in self.list_track_comparison_group_cols.selectedItems()]
        group_x = self.combo_track_comparison_group_x.currentText()
        group_x = None if group_x in ("", "(none)") else group_x
        if not condition_col:
            QMessageBox.warning(self, "Missing selection", "Select a condition column to compare.")
            return
        condition_groups = None
        if self.chk_track_comparison_group_conditions.isChecked():
            if not self.group_selector_track_comparison.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one condition level into each group, or untick "
                    "\"Group conditions\".",
                )
                return
            left = self.group_selector_track_comparison.left_items()
            right = self.group_selector_track_comparison.right_items()
            condition_groups = {lvl: "+".join(left) for lvl in left}
            condition_groups.update({lvl: "+".join(right) for lvl in right})
        group_x_levels_map = None
        if self.chk_track_comparison_group_x_conditions.isChecked():
            if not group_x or not self.group_selector_track_comparison_x.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one Group X level into each group, or untick "
                    "\"Group conditions\" under Group in X.",
                )
                return
            left = self.group_selector_track_comparison_x.left_items()
            right = self.group_selector_track_comparison_x.right_items()
            group_x_levels_map = {lvl: "+".join(left) for lvl in left}
            group_x_levels_map.update({lvl: "+".join(right) for lvl in right})
        out = self._out_dir()
        self._log(f"▶ Creating track condition comparison plot for '{ct}'…")
        track_adata = self._track_adata
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None

        def _run(**kw):
            from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
            from behav3d.analysis.behavior.track.visualization.plots.reports import (
                save_track_condition_comparison_report,
            )
            from behav3d.core.metadata import merge_condition_columns_into_obs
            from behav3d.napari._rename_dialog import _track_cluster_col
            cluster_col = _track_cluster_col(track_adata) or "ClusterID"
            comparison_dir = _resolve_dtaidistance_paths(str(out) if out else "", ct)["behavior_comparisons_outfolder"]
            all_cols = [condition_col] + group_cols + [c for c in (group_x,) if c]
            cols_to_merge = [c for c in all_cols if c not in track_adata.obs.columns]
            if cols_to_merge and md is not None:
                merge_condition_columns_into_obs(track_adata, md, cols_to_merge)
            return save_track_condition_comparison_report(
                track_adata,
                comparison_dir,
                sample_col="sample_name",
                class_col=cluster_col,
                condition_col=condition_col,
                group_cols=group_cols or None,
                group_x=group_x,
                group_x_levels_map=group_x_levels_map,
                condition_groups=condition_groups,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Track condition comparison ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_track_condition_comparison],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Track condition comparison plot done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Track condition comparison failed: {e}"),
        )

    def _on_contact_analysis(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        contact_col = self.combo_contact_col.currentText()
        if not contact_col:
            QMessageBox.warning(self, "Missing selection", "Select a contact column to group tracks by.")
            return
        min_bout_length = int(self.spin_contact_min_bout.value())
        csv_path = self._track_features_csv_path(ct)
        if not csv_path or not csv_path.exists():
            QMessageBox.warning(self, "No data", "Track-features CSV not found. Run feature extraction first.")
            return
        use_target_class = self.chk_use_target_class.isChecked()
        target_available, target_warning = self._target_class_availability()
        if use_target_class and not target_available:
            QMessageBox.warning(self, "Target classification unavailable", target_warning)
            return
        target_ct = self._contact_target_cell_type() if use_target_class else None
        target_source = self.combo_target_class_source.currentData() if use_target_class else None
        target_state_choice = self.combo_target_state_col.currentText()
        out = self._out_dir()
        self._log(f"▶ Running contact-vs-no-contact analysis for '{ct}'…")
        track_adata = self._track_adata
        selected_extra_cols = [item.text() for item in self.list_contact_group_cols.selectedItems()]
        group_x = self.combo_contact_group_x.currentText()
        group_x = None if group_x in ("", "(none)") else group_x
        group_y = self.combo_contact_group_y.currentText()
        group_y = None if group_y in ("", "(none)") else group_y
        group_x_levels_map = None
        if self.chk_contact_group_x_conditions.isChecked():
            if not group_x or not self.group_selector_contact_x.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one Group X level into each group, or untick "
                    "\"Group conditions\" under Group in X.",
                )
                return
            left = self.group_selector_contact_x.left_items()
            right = self.group_selector_contact_x.right_items()
            group_x_levels_map = {lvl: "+".join(left) for lvl in left}
            group_x_levels_map.update({lvl: "+".join(right) for lvl in right})
        group_y_levels_map = None
        if self.chk_contact_group_y_conditions.isChecked():
            if not group_y or not self.group_selector_contact_y.is_configured():
                QMessageBox.warning(
                    self, "Missing selection",
                    "Move at least one Group Y level into each group, or untick "
                    "\"Group conditions\" under Group in Y.",
                )
                return
            left = self.group_selector_contact_y.left_items()
            right = self.group_selector_contact_y.right_items()
            group_y_levels_map = {lvl: "+".join(left) for lvl in left}
            group_y_levels_map.update({lvl: "+".join(right) for lvl in right})
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None

        def _run(**kw):
            import pandas as pd
            import anndata as ad
            from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
            from behav3d.analysis.behavior.track.contact_grouping import (
                touching_column_name,
                build_target_class_lookup_from_state_adata,
                build_target_class_lookup_from_track_adata,
            )
            from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
            from behav3d.analysis.behavior.track.visualization.plots.reports import (
                save_track_contact_group_analysis,
            )
            from behav3d.core.metadata import merge_condition_columns_into_obs
            from behav3d.napari._rename_dialog import _track_cluster_col
            cluster_col = _track_cluster_col(track_adata) or "ClusterID"
            df_timepoints = pd.read_csv(csv_path)
            contact_dir = _resolve_dtaidistance_paths(str(out) if out else "", ct)["outfolder"]
            all_extra_cols = selected_extra_cols + [c for c in (group_x, group_y) if c]
            cols_to_merge = [c for c in all_extra_cols if c not in track_adata.obs.columns]
            if cols_to_merge and md is not None:
                merge_condition_columns_into_obs(track_adata, md, cols_to_merge)

            target_class_kwargs = {}
            if use_target_class:
                touching_col = touching_column_name(target_ct)
                if target_source == "track":
                    adata_target = ad.read_h5ad(str(self._track_adata_path(target_ct)))
                    target_class_lookup = build_target_class_lookup_from_track_adata(
                        adata_target, class_col="ClusterID",
                    )
                    time_varying = False
                else:
                    adata_target = ad.read_h5ad(str(self._state_adata_path(target_ct)))
                    state_col = (
                        FULL_STATE_COL if target_state_choice == "full_behavioral_cluster" else target_state_choice
                    )
                    target_class_lookup = build_target_class_lookup_from_state_adata(
                        adata_target, state_col=state_col,
                    )
                    time_varying = True
                target_class_kwargs = dict(
                    target_class_lookup=target_class_lookup,
                    touching_col=touching_col,
                    time_varying=time_varying,
                    target_cell_type_label=target_ct,
                )

            return save_track_contact_group_analysis(
                track_adata,
                df_timepoints,
                contact_dir,
                contact_col=contact_col,
                min_bout_length=min_bout_length,
                sample_col="sample_name",
                class_col=cluster_col,
                extra_group_cols=selected_extra_cols or None,
                group_x=group_x,
                group_y=group_y,
                group_x_levels_map=group_x_levels_map,
                group_y_levels_map=group_y_levels_map,
                verbose=True,
                **target_class_kwargs,
            )

        self._bg.run(
            fn=_run,
            desc=f"Contact-vs-no-contact analysis ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_contact_analysis],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Contact-vs-no-contact analysis done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Contact-vs-no-contact analysis failed: {e}"),
        )

    def _on_duration_comparison(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        contact_col = self.combo_contact_col.currentText()
        if not contact_col:
            QMessageBox.warning(self, "Missing selection", "Select a contact column to group tracks by.")
            return
        if not self.chk_use_target_class.isChecked():
            QMessageBox.warning(
                self, "Target classification required",
                "Enable \"Use contact cell classification\" to compare durations by class.",
            )
            return
        target_available, target_warning = self._target_class_availability()
        if not target_available:
            QMessageBox.warning(self, "Target classification unavailable", target_warning)
            return
        min_bout_length = int(self.spin_contact_min_bout.value())
        csv_path = self._track_features_csv_path(ct)
        if not csv_path or not csv_path.exists():
            QMessageBox.warning(self, "No data", "Track-features CSV not found. Run feature extraction first.")
            return
        test_mode = self.combo_duration_test_mode.currentData()
        pairing_col = None
        if test_mode == "paired":
            pairing_col = self.combo_duration_pairing_col.currentText().strip()
            if not pairing_col:
                QMessageBox.warning(self, "Missing selection", "Select a pairing column for the paired t-test.")
                return
        comparisons_per_page = int(self.spin_duration_comparisons_per_page.value())
        target_ct = self._contact_target_cell_type()
        target_source = self.combo_target_class_source.currentData()
        target_state_choice = self.combo_target_state_col.currentText()
        out = self._out_dir()
        self._log(f"▶ Running contact duration comparison for '{ct}'…")
        track_adata = self._track_adata
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None

        def _run(**kw):
            import pandas as pd
            import anndata as ad
            from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
            from behav3d.analysis.behavior.track.contact_grouping import (
                touching_column_name,
                build_target_class_lookup_from_state_adata,
                build_target_class_lookup_from_track_adata,
            )
            from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
            from behav3d.analysis.behavior.track.visualization.plots.contact_duration_report import (
                save_track_contact_duration_comparison,
            )
            from behav3d.core.metadata import merge_condition_columns_into_obs
            from behav3d.core.utils import minutes_per_frame_from_metadata
            df_timepoints = pd.read_csv(csv_path)
            contact_dir = _resolve_dtaidistance_paths(str(out) if out else "", ct)["outfolder"]

            if pairing_col and pairing_col not in track_adata.obs.columns and md is not None:
                merge_condition_columns_into_obs(track_adata, md, [pairing_col])

            touching_col = touching_column_name(target_ct)
            if target_source == "track":
                adata_target = ad.read_h5ad(str(self._track_adata_path(target_ct)))
                target_class_lookup = build_target_class_lookup_from_track_adata(
                    adata_target, class_col="ClusterID",
                )
                time_varying = False
            else:
                adata_target = ad.read_h5ad(str(self._state_adata_path(target_ct)))
                state_col = (
                    FULL_STATE_COL if target_state_choice == "full_behavioral_cluster" else target_state_choice
                )
                target_class_lookup = build_target_class_lookup_from_state_adata(
                    adata_target, state_col=state_col,
                )
                time_varying = True

            minutes_per_frame, minutes_valid = minutes_per_frame_from_metadata(md)

            return save_track_contact_duration_comparison(
                track_adata,
                df_timepoints,
                contact_dir,
                contact_col=contact_col,
                min_bout_length=min_bout_length,
                target_class_lookup=target_class_lookup,
                touching_col=touching_col,
                time_varying=time_varying,
                target_cell_type_label=target_ct,
                test_mode=test_mode,
                pairing_col=pairing_col,
                minutes_per_frame=minutes_per_frame if minutes_valid else None,
                comparisons_per_page=comparisons_per_page,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Contact duration comparison ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_duration_comparison],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Contact duration comparison done for '{ct}' ({r.get('n_comparisons')} comparisons)."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Contact duration comparison failed: {e}"),
        )

    def _on_contact_state_shift_analysis(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        contact_col = self.combo_contact_col.currentText()
        if not contact_col:
            QMessageBox.warning(self, "Missing selection", "Select a contact column to analyze.")
            return
        state_adata_path = self._state_adata_path(ct)
        if not state_adata_path or not state_adata_path.exists():
            QMessageBox.warning(
                self, "No data",
                "Behavioral states h5ad not found. Run State Classification first.",
            )
            return
        csv_path = self._track_features_csv_path(ct)
        if not csv_path or not csv_path.exists():
            QMessageBox.warning(self, "No data", "Track-features CSV not found. Run feature extraction first.")
            return
        min_bout_length = int(self.spin_contact_min_bout.value())
        window_mode = self.combo_contact_shift_window_mode.currentText().lower()
        fixed_window_length = int(self.spin_contact_shift_window_length.value())
        state_col_choice = self.combo_contact_shift_state_col.currentText()
        out = self._out_dir()
        track_adata = self._track_adata
        self._log(f"▶ Running contact state-shift analysis for '{ct}'…")

        def _run(**kw):
            import pandas as pd
            import anndata as _ad
            from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
            from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
            from behav3d.analysis.behavior.track.visualization.plots.contact_state_shift_report import (
                save_track_contact_state_shift_report,
            )
            state_col = FULL_STATE_COL if state_col_choice == "full_behavioral_cluster" else state_col_choice
            full_adata = _ad.read_h5ad(str(state_adata_path))
            df_timepoints = pd.read_csv(csv_path)
            contact_dir = _resolve_dtaidistance_paths(str(out) if out else "", ct)["outfolder"]
            return save_track_contact_state_shift_report(
                track_adata,
                df_timepoints,
                full_adata,
                contact_dir,
                contact_col=contact_col,
                min_bout_length=min_bout_length,
                state_col=state_col,
                window_mode=window_mode,
                fixed_window_length=fixed_window_length,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Contact state-shift analysis ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_contact_state_shift],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Contact state-shift analysis done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Contact state-shift analysis failed: {e}"),
        )

    def _on_track_contact_overview(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        contact_col = self.combo_contact_col.currentText()
        if not contact_col:
            QMessageBox.warning(self, "Missing selection", "Select a contact column to analyze.")
            return
        state_adata_path = self._state_adata_path(ct)
        if not state_adata_path or not state_adata_path.exists():
            QMessageBox.warning(
                self, "No data",
                "Behavioral states h5ad not found. Run State Classification first.",
            )
            return
        csv_path = self._track_features_csv_path(ct)
        if not csv_path or not csv_path.exists():
            QMessageBox.warning(self, "No data", "Track-features CSV not found. Run feature extraction first.")
            return
        min_bout_length = int(self.spin_contact_min_bout.value())
        state_col_choice = self.combo_contact_shift_state_col.currentText()
        rows_per_page = int(self.spin_track_overview_rows_per_page.value())
        out = self._out_dir()
        track_adata = self._track_adata
        self._log(f"▶ Creating track contact overview for '{ct}'…")

        def _run(**kw):
            import pandas as pd
            import anndata as _ad
            from behav3d.analysis.behavior.state.classification import FULL_STATE_COL
            from behav3d.analysis.behavior.track.utils import _resolve_dtaidistance_paths
            from behav3d.analysis.behavior.track.visualization.plots.contact_state_shift_report import (
                save_track_contact_overview_report,
            )
            state_col = FULL_STATE_COL if state_col_choice == "full_behavioral_cluster" else state_col_choice
            full_adata = _ad.read_h5ad(str(state_adata_path))
            df_timepoints = pd.read_csv(csv_path)
            contact_dir = _resolve_dtaidistance_paths(str(out) if out else "", ct)["outfolder"]
            return save_track_contact_overview_report(
                track_adata,
                df_timepoints,
                full_adata,
                contact_dir,
                contact_col=contact_col,
                min_bout_length=min_bout_length,
                state_col=state_col,
                rows_per_page=rows_per_page,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Track contact overview ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_track_contact_overview],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Track contact overview done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Track contact overview failed: {e}"),
        )

    # ── Backprojection ───────────────────────────────────────────────────

    def _on_bp_layer_display_changed(self, event=None):
        self._display_save_timer.start()

    def _persist_bp_viewer_display(self):
        _bp_save_channel_display(self.viewer, self.metadata_loader, self._out_dir)

    def _current_viewer_frame(self) -> int:
        if self.viewer is None:
            return 0
        try:
            return int(self.viewer.dims.current_step[0])
        except Exception:
            return 0

    def _teardown_track_bp_preview(self):
        if self._track_bp_dims_callback is not None and self.viewer is not None:
            try:
                self.viewer.dims.events.current_step.disconnect(self._track_bp_dims_callback)
            except Exception:
                pass
        unregister_preview_dims_listener(self.viewer, self)
        self._track_bp_dims_callback = None
        self._track_bp_preview = None

    def _connect_track_bp_dims_listener(self):
        if self._track_bp_dims_callback is not None and self.viewer is not None:
            try:
                self.viewer.dims.events.current_step.disconnect(self._track_bp_dims_callback)
            except Exception:
                pass
            unregister_preview_dims_listener(self.viewer, self)
            self._track_bp_dims_callback = None
        if self.viewer is None:
            return

        def _on_step(*_):
            self._refresh_track_bp_layer()

        try:
            self.viewer.dims.events.current_step.connect(_on_step)
            self._track_bp_dims_callback = _on_step
            register_preview_dims_listener(self.viewer, self, _on_step)
        except Exception:
            self._track_bp_dims_callback = None

    def _refresh_track_bp_layer(self):
        """Recompute the track-cluster overlay for whichever timepoint napari's
        time slider is currently on. Mirrors the Feature Backprojection tab's
        ``_refresh_feature_layer`` (see ``behav3d.napari._feature_backprojection``)
        so opening the preview never requires writing a full per-sample zarr."""
        from behav3d.analysis.behavior.track.visualization.backprojection import (
            backproject_track_cluster_at_timepoint,
        )
        from behav3d.io.images import load_image_timepoint

        preview = self._track_bp_preview
        if preview is None or self.viewer is None:
            return

        t = self._current_viewer_frame()
        try:
            labels_frame = np.asarray(load_image_timepoint(preview["tracked_path"], t))
        except Exception as exc:
            self._log(f"❌ Could not read frame {t}: {exc}")
            return
        if labels_frame.ndim == 4:
            labels_frame = labels_frame[0]

        mapped, _ids_with_value = backproject_track_cluster_at_timepoint(
            labels_frame=labels_frame,
            cluster_code_lookup=preview["code_lookup"],
            time_index=t,
        )

        name = preview["layer_name"]
        try:
            layer = self.viewer.layers[name]
            layer.data = mapped
            layer.refresh()
        except (KeyError, ValueError):
            layer = self.viewer.add_labels(mapped, name=name, opacity=preview["opacity"])
            from behav3d.analysis.behavior.state.visualization.backprojection import (
                _apply_state_code_colors_to_layer,
            )
            _apply_state_code_colors_to_layer(layer, preview.get("code_colors", {}))

    def _on_show_track_bp(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        track_path = self._track_adata_path(ct)
        if not track_path or not track_path.exists():
            QMessageBox.warning(self, "No track adata", "Run Track Clustering first.")
            return
        opacity = self.spin_track_opacity.value() / 100.0
        state_adata_path = self._state_adata_path(ct)
        if not state_adata_path or not state_adata_path.exists():
            QMessageBox.warning(
                self, "State Classification Required",
                f"Full state adata not found at:\n{state_adata_path}\n\n"
                "Run State Classification first."
            )
            return
        self._log(f"▶ Loading track backprojection for '{ct}' / sample '{sample}'…")
        try:
            import scanpy as sc
            from behav3d.analysis.behavior.state.visualization.backprojection import (
                _resolve_raw_image_path,
                _resolve_tracked_image_path,
                _build_code_map,
                _build_state_code_color_map,
                _add_mapping_dock_widget,
                _build_state_mapping_text,
                _align_labels_to_raw_shape_for_view,
            )
            from behav3d.analysis.behavior.state.utils import (
                _get_classification_state_colors,
                _get_classification_state_order,
                _normalize_label_color_map,
            )
            from behav3d.analysis.behavior.track.visualization.backprojection import (
                build_track_cluster_frame_assignments,
                prepare_track_cluster_code_lookup,
                _add_track_statebar_click_dock,
            )
            from behav3d.analysis.backprojection import filter_track_image_to_ids
            from behav3d.io.images import load_image
            out_dir = self._out_dir()
            if not out_dir:
                raise ValueError("No output directory set.")
            adata_tracks = sc.read_h5ad(str(track_path))
            self._sync_track_cluster_combo(adata_tracks)
            color_by = self.combo_track_color_by.currentText()
            adata_full = sc.read_h5ad(str(state_adata_path))
            cluster_col = color_by if color_by else "ClusterID"
            output_col = "track_behavioral_cluster"
            obs_samples = adata_tracks.obs["sample_name"].astype(str)
            sample_name = sample if sample else obs_samples.iloc[0]

            backproj_obs, _missing_windows = build_track_cluster_frame_assignments(
                adata_full=adata_full,
                adata_tracks=adata_tracks,
                cluster_col=cluster_col,
                output_col=output_col,
                sample_name=sample_name,
                verbose=False,
            )
            state_order = _get_classification_state_order(adata_tracks, cluster_col)
            code_map = _build_code_map(backproj_obs, state_col=output_col, state_order=state_order)
            if len(code_map) == 0:
                raise ValueError(f"'{cluster_col}' has no non-empty labels for sample '{sample_name}'.")
            state_colors = _normalize_label_color_map(
                code_map.keys(), colors=_get_classification_state_colors(adata_tracks, cluster_col)
            )
            code_colors = _build_state_code_color_map(code_map, state_colors=state_colors)
            label_map = {str(code): str(label) for label, code in code_map.items()}
            cluster_code_lookup = prepare_track_cluster_code_lookup(
                backproj_obs, code_map, output_col=output_col,
            )

            gui_metadata_csv_path = resolve_metadata_csv_path(self.metadata_loader)
            raw_path = _resolve_raw_image_path(
                out_dir, sample_name, verbose=False, metadata_csv_path=gui_metadata_csv_path
            )
            if raw_path is None or not Path(raw_path).exists():
                raise FileNotFoundError(f"Raw image not found for sample '{sample_name}'.")
            tracked_path = _resolve_tracked_image_path(
                out_dir, sample_name, ct, verbose=False, metadata_csv_path=gui_metadata_csv_path
            )
            if tracked_path is None or not Path(tracked_path).exists():
                raise FileNotFoundError(
                    f"Tracked image not found for sample '{sample_name}', cell_type '{ct}'."
                )
            raw_img = load_image(raw_path)
            tracked_img = load_image(tracked_path)
            tracked_view = _align_labels_to_raw_shape_for_view(tracked_img, raw_img, "TrackID", verbose=False)
            keep_ids = cluster_code_lookup["TrackID"].unique()
            tracked_view = filter_track_image_to_ids(tracked_view, keep_ids)

            # Full viewer reset — see the matching comment in
            # ``_on_show_state_bp`` above.
            self._teardown_track_bp_preview()
            stop_dim_playback(self.viewer)
            disconnect_all_preview_dims_listeners(self.viewer)
            clear_viewer_layers(self.viewer)
            saved_channels = (
                getattr(self.metadata_loader, "behav3d_parameters", {})
                .get("viewer_display", {})
                .get("channels", {})
            )
            try:
                md = self.metadata_loader.metadata
                row = md[md["sample_name"] == sample_name].iloc[0]
                dim_order = str(row.get("dimension_order", "TCZYX")).strip() or "TCZYX"
            except Exception:
                dim_order = "TCZYX"
            ch_names = _bp_add_raw_channels(self.viewer, raw_img, sample_name, saved_channels, dim_order)
            for lname in ch_names:
                try:
                    layer = self.viewer.layers[lname]
                    layer.events.contrast_limits.connect(self._on_bp_layer_display_changed)
                    layer.events.colormap.connect(self._on_bp_layer_display_changed)
                except (KeyError, IndexError):
                    pass
            self.viewer.add_labels(tracked_view, name="filtered TrackID", visible=False, opacity=opacity)

            self._track_bp_preview = {
                "tracked_path": Path(tracked_path),
                "code_lookup": cluster_code_lookup,
                "layer_name": "behavioral_state_class",
                "opacity": opacity,
                "code_colors": code_colors,
            }
            self._refresh_track_bp_layer()
            self._connect_track_bp_dims_listener()

            if self.chk_show_track_trajectories.isChecked():
                try:
                    from behav3d.analysis.behavior.track.visualization.plots.exemplar_coordinate_utils import (
                        ensure_exemplar_coordinate_columns,
                    )
                    from behav3d.analysis.behavior.track.visualization.backprojection import (
                        prepare_track_cluster_trajectory_data,
                        add_track_cluster_trajectory_layers,
                    )
                    ensure_exemplar_coordinate_columns(
                        adata_full, output_dir=out_dir, cell_type=ct, require_pixel_for_video=True,
                    )
                    trajectory_data = prepare_track_cluster_trajectory_data(
                        adata_full=adata_full,
                        backproj_obs=backproj_obs,
                        output_col=output_col,
                        sample_name=sample_name,
                    )
                    add_track_cluster_trajectory_layers(
                        self.viewer,
                        trajectory_data=trajectory_data,
                        code_colors=code_colors,
                        label_map=label_map,
                        output_col=output_col,
                        tail_length=int(tracked_img.shape[0]),
                    )
                except Exception as exc:
                    self._log(f"⚠️ Could not add trajectory layers: {exc}")

            mapping_text = _build_state_mapping_text(label_map, code_colors)
            _add_mapping_dock_widget(
                self.viewer,
                mapping_text=mapping_text,
                label_map=label_map,
                code_colors=code_colors,
                title="Track Cluster Mapping",
            )
            _add_track_statebar_click_dock(
                self.viewer,
                sample_name=sample_name,
                adata_full=adata_full,
                adata_tracks=adata_tracks,
                cluster_col=cluster_col,
                title="Track State Bar",
            )
            self._log("✅ Track backprojection loaded (current timepoint; updates as you scrub).")
        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Backprojection failed: {e}")

    def _on_export_track_bp(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        track_path = self._track_adata_path(ct)
        if not track_path or not track_path.exists():
            QMessageBox.warning(self, "No track adata", "Run Track Clustering first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        import scanpy as sc
        adata_tracks = sc.read_h5ad(str(track_path))
        self._sync_track_cluster_combo(adata_tracks)
        color_by = self.combo_track_color_by.currentText()
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Exporting track backprojection for '{ct}'…")

        def _run(**kw):
            from behav3d.analysis.behavior.track.visualization.backprojection import (
                export_track_cluster_backprojection as _track_bp_export,
            )
            state_adata_path = self._state_adata_path(ct)
            if not state_adata_path or not state_adata_path.exists():
                raise FileNotFoundError(
                    f"Full state adata not found at '{state_adata_path}'. "
                    "Run State Classification first."
                )
            adata_full = sc.read_h5ad(str(state_adata_path))
            cluster_col = color_by if color_by else "ClusterID"
            sample_name = sample if sample else None
            logger("▶ Writing track backprojection zarrs…")
            return _track_bp_export(
                adata_full=adata_full,
                adata_tracks=adata_tracks,
                output_dir=out,
                cell_type=ct,
                cluster_col=cluster_col,
                sample_name=sample_name,
                n_workers=1,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Exporting track backprojection ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_export_track_bp, self.btn_show_track_bp],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: _log_backprojection_manifest(self._log, r),
            on_failed=lambda e: self._log(f"❌ Export failed: {e}"),
        )

    # ── View helpers ─────────────────────────────────────────────────────

    def _on_view(self, kind: str, btn: Optional[QPushButton] = None):
        ct = self._cell_type()
        if not ct:
            return
        out = self._out_dir()
        traj_dir = (out / "analysis" / ct / "behavorial_trajectories") if out else None
        candidates = []
        if kind == "track_exemplars" and traj_dir:
            exemplar_dir = traj_dir / "example_tracks"
            candidates = [
                (f.relative_to(exemplar_dir).as_posix(), f)
                for f in sorted(exemplar_dir.glob("**/*.pdf"))
            ]
        elif kind == "track_diagnostics" and traj_dir:
            candidates = [
                (f.stem, f)
                for qc_dir in self._track_diagnostics_qc_dirs(traj_dir)
                for f in sorted(qc_dir.glob("*.pdf"))
            ]
        elif kind == "track_proportions" and traj_dir:
            proportions_dir = traj_dir / "behavior_proportions"
            candidates = [
                (f.stem, f)
                for f in sorted(proportions_dir.glob("*.pdf"))
            ]
        elif kind == "window_transitions" and traj_dir:
            wt_dir = traj_dir / "window_transitions"
            candidates = [(f.stem, f) for f in sorted(wt_dir.glob("*.pdf"))]
            candidates += [
                (f"per-sample/{f.stem}", f)
                for f in sorted(wt_dir.glob("sankey_pdf_pages/*.pdf"))
            ]
        elif kind == "track_condition_comparison" and traj_dir:
            comparisons_dir = traj_dir / "behavior_comparisons"
            candidates = [
                (f.stem.replace("condition_comparison_", ""), f)
                for f in sorted(comparisons_dir.glob("condition_comparison_*.pdf"))
            ]
        elif kind == "contact_analysis" and traj_dir:
            contact_dir = traj_dir / "contact_analysis"
            candidates = [
                (f"{f.parent.name}/{f.stem}", f)
                for f in sorted(contact_dir.glob("*/*.pdf"))
            ]
        elif kind == "contact_state_shift" and traj_dir:
            contact_dir = traj_dir / "contact_analysis"
            candidates = [
                (f.parent.name, f)
                for f in sorted(contact_dir.glob("*/contact_state_shift.pdf"))
            ]
        elif kind == "track_contact_overview" and traj_dir:
            contact_dir = traj_dir / "contact_analysis"
            candidates = [
                (f.parent.name, f)
                for f in sorted(contact_dir.glob("*/track_contact_overview.pdf"))
            ]
        elif kind == "contact_duration" and traj_dir:
            contact_dir = traj_dir / "contact_analysis"
            candidates = [
                (f.parent.name, f)
                for f in sorted(contact_dir.glob("*/contact_duration_comparison.pdf"))
            ]

        existing = [(lbl, p) for lbl, p in candidates if p and p.exists()]
        if not existing:
            return
        if len(existing) == 1:
            self._open_result(existing[0][1])
        else:
            menu = QMenu(self)
            for lbl, path in existing:
                act = menu.addAction(f"👁  {lbl}")
                act.triggered.connect(lambda _=False, _p=path: self._open_result(_p))
            if btn is not None:
                menu.exec_(btn.mapToGlobal(btn.rect().bottomLeft()))
            else:
                menu.exec_()

    def _open_result(self, path: Path):
        if not path.suffix.lower() == ".pdf":
            try:
                from qtpy.QtCore import QUrl
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(path.parent)))
            except Exception:
                pass
            return
        try:
            panel = self._results_panel()
            dpi = int(panel.spin_dpi.value()) if panel else 150
            reuse = getattr(panel, "chk_reuse", None)
            reuse = reuse.isChecked() if reuse else True
            open_pdf_in_napari(path, dpi=dpi, reuse=reuse)
        except Exception as e:
            traceback.print_exc()
            self._log(f"Could not open in napari: {e}")

    def _results_panel(self):
        node = self.parent()
        while node and not hasattr(node, "results_panel"):
            node = node.parent()
        return getattr(node, "results_panel", None) if node else None

    def _notify_results(self):
        panel = self._results_panel()
        if panel:
            try:
                panel.refresh()
            except Exception:
                pass


# ═══════════════════════════════════════════════════════════════════════════
# Outer SingleCellTab
# ═══════════════════════════════════════════════════════════════════════════

class SingleCellTab(QWidget):
    """Outer single-cell tab: shared cell-type + sample dropdowns + two inner sub-tabs."""

    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader

        self._init_ui()

        if metadata_loader is not None and hasattr(metadata_loader, "metadata_loaded"):
            metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

        if (
            metadata_loader is not None
            and getattr(metadata_loader, "metadata", None) is not None
        ):
            self._on_metadata_updated()

    # ── UI ───────────────────────────────────────────────────────────────

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(4)

        # ── Shared header ────────────────────────────────────────────────
        hdr_lay = QHBoxLayout()
        hdr_lay.setSpacing(6)

        hdr_lay.addWidget(QLabel("Cell type:"))
        self.cell_type_combo = QComboBox()
        self.cell_type_combo.setMinimumWidth(160)
        self.cell_type_combo.setToolTip("Immune and other cell types only (non-multicolor).")
        self.cell_type_combo.currentTextChanged.connect(self._on_cell_type_changed)
        hdr_lay.addWidget(self.cell_type_combo)
        hdr_lay.addStretch()

        self.status_lbl = QLabel("Load metadata to begin.")
        self.status_lbl.setStyleSheet("color: #888; font-size: 11px;")
        hdr_lay.addWidget(self.status_lbl)
        outer.addLayout(hdr_lay)

        # ── Data-consistency warning (h5ad vs. filtered CSV track IDs) ────
        self.data_consistency_warning_label = QLabel("")
        self.data_consistency_warning_label.setWordWrap(True)
        self.data_consistency_warning_label.setStyleSheet(
            "QLabel { background: #3d2200; color: #ffaa44; border-radius: 4px; "
            "padding: 6px 8px; font-size: 11px; }"
        )
        self.data_consistency_warning_label.hide()
        outer.addWidget(self.data_consistency_warning_label)
        self._consistency_bg = BackgroundOperation(self)

        # ── Guided overview + isolated settings pages ────────────────────
        from behav3d.napari._guided import GuidedPanel, make_back_header
        from behav3d.napari.analysis_guided_copy import (
            SINGLE_CELL_ANALYSES, GUIDED_INTRO,
        )

        self._stack = QStackedWidget()
        outer.addWidget(self._stack, stretch=1)

        # Page 0 — Guided explainers.
        self._guided_panel = GuidedPanel(
            SINGLE_CELL_ANALYSES,
            start_cb=self._on_guided_start,
            intro=GUIDED_INTRO,
        )
        self._stack.addWidget(self._guided_panel)

        # Page 1 — the existing inner sub-tabs (settings forms), with a Back
        # header above them. The tab bar is hidden so only the sub-tab
        # selected via Start is visible — this is an isolated settings view,
        # not a tab switcher between the two analyses.
        settings_page = QWidget()
        settings_lay = QVBoxLayout(settings_page)
        settings_lay.setContentsMargins(0, 0, 0, 0)
        settings_lay.setSpacing(0)
        back_bar, self._focus_title = make_back_header(
            on_back=lambda: self._stack.setCurrentIndex(0)
        )
        settings_lay.addWidget(back_bar)

        self.inner_tabs = QTabWidget()
        self.inner_tabs.tabBar().setVisible(False)
        settings_lay.addWidget(self.inner_tabs)
        self._stack.addWidget(settings_page)
        reset_scroll_on_page_change(self._stack)

        self.state_tab = StateClassificationSubTab(
            viewer=self.viewer,
            metadata_loader=self.metadata_loader,
            cell_type_getter=self._current_cell_type,
            parent=self,
        )
        # Display name only — internal SC_STATE_CLUSTER StepType / methods /
        # outputs are unchanged.
        self.inner_tabs.addTab(self.state_tab, "🔬 Behavioral State")

        self.track_tab = TrackClassificationSubTab(
            viewer=self.viewer,
            metadata_loader=self.metadata_loader,
            cell_type_getter=self._current_cell_type,
            parent=self,
        )
        # Display name only — internal SC_TRACK_CLUSTER StepType / methods /
        # outputs are unchanged.
        self.inner_tabs.addTab(self.track_tab, "🛤️ State Trajectory")

        # When switching to State Trajectory (index 1), trigger a reload
        # so the behavioral-states path auto-fills from disk if it exists.
        self.inner_tabs.currentChanged.connect(self._on_inner_tab_changed)
        reset_scroll_on_page_change(self.inner_tabs)

        # Always land on the Guided overview; a specific analysis's settings
        # are only reached via that analysis's Start button.
        self._stack.setCurrentIndex(0)

    # ── Guided overview / focused settings ──────────────────────────────────
    def _on_guided_start(self, analysis_id: str):
        """Start from a Guided card: show only that analysis's settings."""
        from behav3d.napari.analysis_guided_copy import (
            BEHAVIORAL_STATE, STATE_TRAJECTORY,
        )

        index = {"state": 0, "track": 1}.get(analysis_id)
        if index is not None:
            self.inner_tabs.setCurrentIndex(index)
        title = {
            "state": BEHAVIORAL_STATE["title"],
            "track": STATE_TRAJECTORY["title"],
        }.get(analysis_id, "")
        self._focus_title.setText(title)
        self._stack.setCurrentIndex(1)

    def _on_inner_tab_changed(self, index: int):
        """Auto-fill path fields when switching to Track Classification tab."""
        if index == 1:  # Track Classification
            self.track_tab._autofill_paths(self._current_cell_type())
            # Re-check prerequisites too: a behavioral-states h5ad may have been
            # written (or rewritten) by the State Classification tab since the
            # last time this tab was visited, and _autofill_paths alone doesn't
            # refresh the "states not found" warning/lock state.
            self.track_tab._check_prerequisites()
            self.track_tab._apply_max_trajectory_size(self._current_cell_type())
            if self.track_tab.chk_use_original.isChecked():
                self.track_tab._show_original_dtw_disclaimer()

    # ── Metadata update ──────────────────────────────────────────────────

    def _on_metadata_updated(self, *_):
        cell_types = _detect_sc_cell_types(self.metadata_loader)
        current = self.cell_type_combo.currentText()
        self.cell_type_combo.blockSignals(True)
        self.cell_type_combo.clear()
        if cell_types:
            self.cell_type_combo.addItems(cell_types)
            saved_ct = (
                getattr(self.metadata_loader, "behav3d_parameters", {})
                .get("single_cell", {})
                .get("selected_cell_type")
            )
            if saved_ct and saved_ct in cell_types:
                self.cell_type_combo.setCurrentText(saved_ct)
            elif current in cell_types:
                self.cell_type_combo.setCurrentText(current)
            self.status_lbl.setText(f"Ready — {len(cell_types)} cell type(s) available.")
        else:
            self.status_lbl.setText("No immune/other cell types detected.")
        self.cell_type_combo.blockSignals(False)

        # Populate sample combos in sub-tabs
        md = getattr(self.metadata_loader, "metadata", None) if self.metadata_loader else None
        samples = ["— All samples —"]
        if md is not None and "sample_name" in md.columns:
            samples.extend(sorted(md["sample_name"].astype(str).unique()))
            
        for combo in (self.state_tab.sample_combo, self.track_tab.sample_combo):
            cur_s = combo.currentText()
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(samples)
            if cur_s in samples:
                combo.setCurrentText(cur_s)
            combo.blockSignals(False)

        self._propagate_metadata_update()
        self._refresh_data_consistency_warning()

    def _propagate_metadata_update(self):
        self.state_tab.on_metadata_updated()
        self.track_tab.on_metadata_updated()

    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

    def _on_cell_type_changed(self, text: str):
        params = getattr(self.metadata_loader, "behav3d_parameters", None)
        if isinstance(params, dict) and text:
            params.setdefault("single_cell", {})["selected_cell_type"] = text
            _save_behav3d_params(self.metadata_loader, self._out_dir)
        self._propagate_metadata_update()
        self._refresh_data_consistency_warning()

    def _current_cell_type(self) -> str:
        return self.cell_type_combo.currentText()

    # ── Data-consistency warning ─────────────────────────────────────────

    def showEvent(self, event):
        """Fires whenever this tab becomes visible (initial show and every tab
        switch back to it) — the natural "entering Single-Cell Analysis" hook,
        so the warning reflects any h5ad/CSV drift since it was last shown."""
        super().showEvent(event)
        self._refresh_data_consistency_warning()

    def _refresh_data_consistency_warning(self):
        """Cheap background check: do the filtered track-features CSV and the
        track/behavioral-states h5ad(s) for the current cell type still describe
        the same tracks?

        Only reads the two ID columns from the CSV (not the full file) and the
        ``.obs`` table from each h5ad in ``backed="r"`` mode (not the full
        AnnData, so no X matrix is loaded) — this is far cheaper than the
        full contact/no-contact analysis pipeline and safe to run on every tab
        entry. Runs off the Qt thread so it never blocks tab switching; if a
        previous check is still in flight, this call is skipped (a subsequent
        trigger — e.g. the next tab visit — will pick it up).
        """
        ct = self._current_cell_type()
        out = self._out_dir()
        self.data_consistency_warning_label.hide()
        if not ct or not out:
            return

        csv_path = self.track_tab._track_features_csv_path(ct)
        if not csv_path or not csv_path.exists():
            return
        h5ad_sources = [
            ("track classification h5ad", self.track_tab._track_adata_path(ct)),
            ("behavioral states h5ad", self.track_tab._state_adata_path(ct)),
        ]
        h5ad_sources = [(label, p) for label, p in h5ad_sources if p and p.exists()]
        if not h5ad_sources:
            return

        if self._consistency_bg.is_running():
            return

        groupby_cols = ["sample_name", "TrackID"]

        def _key_set(df):
            from behav3d.analysis.behavior.track.contact_grouping import _normalize_id_column
            d = df[groupby_cols].copy()
            for col in groupby_cols:
                d[col] = _normalize_id_column(d[col])
            return set(map(tuple, d.drop_duplicates().to_numpy()))

        def _check():
            import pandas as pd
            import anndata as ad

            csv_keys = _key_set(pd.read_csv(csv_path, usecols=groupby_cols))
            mismatches = []
            for label, path in h5ad_sources:
                adata = ad.read_h5ad(str(path), backed="r")
                try:
                    h5ad_keys = _key_set(adata.obs[groupby_cols])
                finally:
                    if getattr(adata, "isbacked", False):
                        adata.file.close()
                missing = h5ad_keys - csv_keys
                if missing:
                    mismatches.append((label, len(missing)))
            return mismatches

        def _on_done(mismatches):
            if not mismatches:
                self.data_consistency_warning_label.hide()
                return
            details = "; ".join(f"{n} track(s) missing from {label}" for label, n in mismatches)
            self.data_consistency_warning_label.setText(
                f"⚠ Filtered track data no longer matches the behavioral-analysis h5ad output "
                f"for cell type '{ct}' ({details}). Re-run Track / State Classification to "
                f"refresh the h5ad from the current filtered data."
            )
            self.data_consistency_warning_label.show()

        def _on_failed(err):
            print(f"[BEHAV3D] Data-consistency check failed: {err}")

        self._consistency_bg.run(
            fn=_check,
            desc="Checking data consistency…",
            inject_progress=False,
            on_done=_on_done,
            on_failed=_on_failed,
        )
