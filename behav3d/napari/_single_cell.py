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

import yaml
from qtpy.QtCore import Qt, QTimer
from qtpy.QtGui import QDesktopServices
from qtpy.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
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
    QTabWidget,
    QTextEdit,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from behav3d.napari._analysis import CollapsibleSection
from behav3d.napari._background_runner import (
    BackgroundOperation,
    ProgressBarRow,
    ThreadSafeLogger,
)
from behav3d.napari._pdf_view import open_pdf_in_napari
from behav3d.napari._rename_dialog import RenameClusterDialog
from behav3d.core.qt_help import HelpButton, make_help_row


# ═══════════════════════════════════════════════════════════════════════════
# Shared helpers
# ═══════════════════════════════════════════════════════════════════════════

def _detect_sc_cell_types(metadata_loader) -> list[str]:
    """Return immune + other cell types (non-multicolor)."""
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

        return imm + oth
    except Exception:
        traceback.print_exc()
        return []


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


def _make_info_label(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setWordWrap(True)
    lbl.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
    lbl.setAlignment(Qt.AlignLeft | Qt.AlignTop)
    return lbl


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
        self._bg = BackgroundOperation(self)

        # Queue buttons (wired by _widget.py)
        self.btn_queue_state_cluster = _make_queue_btn()
        self.btn_queue_train_state = _make_queue_btn()
        self.btn_queue_apply_state = _make_queue_btn()

        self._init_ui()

    # ── UI ──────────────────────────────────────────────────────────────

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        scroll.setWidget(content)

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

        self.feat_sel_lay.addWidget(QLabel("<b>Timepoint features</b>"))
        self.timepoint_features_lay = QVBoxLayout()
        self.feat_sel_lay.addLayout(self.timepoint_features_lay)

        self.feat_sel_lay.addWidget(QLabel("<b>Window features</b>"))
        win_feat_form = QFormLayout()
        self.spin_window_size = QSpinBox()
        self.spin_window_size.setRange(1, 500)
        self.spin_window_size.setValue(5)
        win_feat_form.addRow("Window size:", make_help_row(
            self.spin_window_size, "Window size",
            "Number of timepoints used for computing rolling window features "
            "(e.g. net displacement, straightness over a sliding window)."
        ))
        self.feat_sel_lay.addLayout(win_feat_form)

        self.chk_net_disp = QCheckBox("net_displacement")
        self.chk_straight = QCheckBox("straightness")
        self.chk_msd = QCheckBox("mean_square_displacement")
        self.feat_sel_lay.addLayout(_make_chk_help_row(self.chk_net_disp,
            "net_displacement", "Include net displacement window feature."))
        self.feat_sel_lay.addLayout(_make_chk_help_row(self.chk_straight,
            "straightness", "Include straightness window feature."))
        self.feat_sel_lay.addLayout(_make_chk_help_row(self.chk_msd,
            "mean_square_displacement", "Include mean square displacement window feature."))

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

        self.feat_proc_lay.addWidget(QLabel("<b>Log scaling</b>"))
        self.log_scale_lay = QVBoxLayout()
        self.feat_proc_lay.addLayout(self.log_scale_lay)

        proc_form = QFormLayout()
        self.spin_hmm_feature_smoothing_window = QSpinBox()
        self.spin_hmm_feature_smoothing_window.setRange(1, 100)
        self.spin_hmm_feature_smoothing_window.setValue(1)
        proc_form.addRow("Smooth window:", make_help_row(
            self.spin_hmm_feature_smoothing_window, "Feature smoothing window",
            "Size of the rolling smoothing window applied to features before HMM fitting. "
            "1 = no smoothing."
        ))

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
        self.spin_quant_hi.setValue(0.99)
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

        # Advanced Settings
        adv_sec = CollapsibleSection("⚙ Advanced Configuration", expanded=False)
        adv_form = QFormLayout()

        self.combo_hmm_n_states_mode = QComboBox()
        self.combo_hmm_n_states_mode.addItems(["fixed", "auto"])
        adv_form.addRow("State selection mode:", make_help_row(
            self.combo_hmm_n_states_mode, "State selection mode",
            "'fixed' uses n_states directly. 'auto' searches between k_min and k_max "
            "and selects the best number by BIC/AIC."
        ))

        self.spin_hmm_k_min = QSpinBox()
        self.spin_hmm_k_min.setRange(2, 50)
        self.spin_hmm_k_min.setValue(2)
        self.row_k_min = QLabel("k_min (Auto mode):")
        adv_form.addRow(self.row_k_min, make_help_row(
            self.spin_hmm_k_min, "k_min",
            "Minimum number of states to test in auto mode."
        ))

        self.spin_hmm_k_max = QSpinBox()
        self.spin_hmm_k_max.setRange(2, 50)
        self.spin_hmm_k_max.setValue(8)
        self.row_k_max = QLabel("k_max (Auto mode):")
        adv_form.addRow(self.row_k_max, make_help_row(
            self.spin_hmm_k_max, "k_max",
            "Maximum number of states to test in auto mode."
        ))

        self.spin_hmm_start_offset = QSpinBox()
        self.spin_hmm_start_offset.setRange(0, 100000)
        self.spin_hmm_start_offset.setValue(0)
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
        adv_form.addRow("", make_help_row(
            self.chk_hmm_sticky, "Sticky HMM",
            "Adds a self-transition bias (kappa) to make the model prefer staying "
            "in the same state, producing longer, more stable behavioral bouts."
        ))

        self.spin_hmm_stickiness_kappa = QDoubleSpinBox()
        self.spin_hmm_stickiness_kappa.setRange(0.0, 100.0)
        self.spin_hmm_stickiness_kappa.setValue(8.0)
        self.row_kappa = QLabel("kappa (Sticky):")
        adv_form.addRow(self.row_kappa, make_help_row(
            self.spin_hmm_stickiness_kappa, "Stickiness kappa",
            "Self-transition bias strength. Higher = longer bouts before switching states."
        ))

        self.spin_hmm_transmat_alpha = QDoubleSpinBox()
        self.spin_hmm_transmat_alpha.setRange(0.0, 100.0)
        self.spin_hmm_transmat_alpha.setValue(1.0)
        self.row_alpha = QLabel("alpha (Sticky):")
        adv_form.addRow(self.row_alpha, make_help_row(
            self.spin_hmm_transmat_alpha, "Transition matrix alpha",
            "Dirichlet prior concentration for the transition matrix in sticky HMM mode."
        ))

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

        rename_intrinsic_row = QHBoxLayout()
        self.btn_rename_intrinsic = QPushButton("✏  Rename Primary Dynamic State Clusters")
        _style_rename(self.btn_rename_intrinsic)
        self.btn_rename_intrinsic.setEnabled(False)
        rename_intrinsic_row.addWidget(self.btn_rename_intrinsic, stretch=1)
        self.btn_view_rename_intrinsic = _make_view_btn()
        rename_intrinsic_row.addWidget(self.btn_view_rename_intrinsic)
        g2.addLayout(rename_intrinsic_row)

        rename_full_row = QHBoxLayout()
        self.btn_rename_full = QPushButton("✏  Rename Full Behavioral Clusters (Binary Groups)")
        _style_rename(self.btn_rename_full)
        self.btn_rename_full.setEnabled(False)
        rename_full_row.addWidget(self.btn_rename_full, stretch=1)
        self.btn_view_rename_full = _make_view_btn()
        rename_full_row.addWidget(self.btn_view_rename_full)
        g2.addLayout(rename_full_row)
        lay.addWidget(self.grp2)

        # ── Step 3: Reports ──────────────────────────────────────────────
        self.grp3 = QGroupBox("Step 3 — Reports")
        g3 = QVBoxLayout(self.grp3)
        g3.setSpacing(4)

        comp_row = QHBoxLayout()
        self.btn_state_composition = QPushButton("▶ State Composition Report")
        _style_secondary(self.btn_state_composition)
        comp_row.addWidget(self.btn_state_composition, stretch=1)
        self.btn_view_composition = _make_view_btn()
        comp_row.addWidget(self.btn_view_composition)
        g3.addLayout(comp_row)

        g3.addWidget(QLabel("Group composition plots by (Ctrl/Cmd click for multiple):"))
        self.list_composition_group_cols = QListWidget()
        self.list_composition_group_cols.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.list_composition_group_cols.setMaximumHeight(80)
        g3.addWidget(self.list_composition_group_cols)

        trans_row = QHBoxLayout()
        self.btn_state_transition = QPushButton("▶ State Transition Report")
        _style_secondary(self.btn_state_transition)
        trans_row.addWidget(self.btn_state_transition, stretch=1)
        self.btn_view_transition = _make_view_btn()
        trans_row.addWidget(self.btn_view_transition)
        g3.addLayout(trans_row)
        lay.addWidget(self.grp3)

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
            ["full_behavioral_cluster", "intrinsic_behavioral_cluster"]
        )
        self.combo_state_color_by.setMinimumWidth(220)
        bp_form.addRow("Color by:", make_help_row(
            self.combo_state_color_by, "Color by",
            "Which clustering label to use for coloring the backprojection overlay. "
            "'full_behavioral_cluster' includes binary grouping; "
            "'intrinsic_behavioral_cluster' shows only the HMM states."
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
        g_view.addLayout(bp_form)

        view_row = QHBoxLayout()
        self.btn_show_state_bp = QPushButton("▶ Show State Backprojection in Napari")
        _style_primary(self.btn_show_state_bp)
        view_row.addWidget(self.btn_show_state_bp, stretch=1)
        g_view.addLayout(view_row)
        g_bp.addWidget(grp_view)

        # Export options (collapsed)
        export_sec = CollapsibleSection("⚙ Export Options", expanded=False)
        export_form = QFormLayout()
        self.spin_state_export_dpi = QSpinBox()
        self.spin_state_export_dpi.setRange(50, 600)
        self.spin_state_export_dpi.setValue(150)
        self.spin_state_export_dpi.setSuffix(" dpi")
        self.spin_state_export_dpi.setMaximumWidth(110)
        export_form.addRow("DPI (PDF):", make_help_row(
            self.spin_state_export_dpi, "Export DPI",
            "Resolution for exported PDF images."
        ))
        export_sec.addLayout(export_form)

        export_opts = QHBoxLayout()
        self.chk_state_export_pdf = QCheckBox("PDF")
        self.chk_state_export_pdf.setChecked(True)
        export_opts.addWidget(self.chk_state_export_pdf)
        self.chk_state_export_mp4 = QCheckBox("MP4")
        self.chk_state_export_mp4.setChecked(False)
        export_opts.addWidget(self.chk_state_export_mp4)
        export_opts.addStretch()
        export_sec.addLayout(export_opts)

        export_row = QHBoxLayout()
        self.btn_export_state_bp = QPushButton("▶ Export State Backprojection")
        _style_secondary(self.btn_export_state_bp)
        export_row.addWidget(self.btn_export_state_bp, stretch=1)
        export_sec.addLayout(export_row)
        g_bp.addWidget(export_sec)

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
        self.btn_view_state.clicked.connect(lambda: self._on_view("state_qc"))
        self.btn_rename_intrinsic.clicked.connect(self._on_rename_intrinsic)
        self.btn_rename_full.clicked.connect(self._on_rename_full)
        self.btn_state_composition.clicked.connect(self._on_state_composition)
        self.btn_view_composition.clicked.connect(lambda: self._on_view("state_composition"))
        self.btn_state_transition.clicked.connect(self._on_state_transition)
        self.btn_view_transition.clicked.connect(lambda: self._on_view("state_transition"))
        self.btn_browse_hmm.clicked.connect(self._browse_hmm_artifact)
        self.btn_apply_hmm.clicked.connect(self._on_apply_hmm)
        self.btn_show_state_bp.clicked.connect(self._on_show_state_bp)
        self.btn_export_state_bp.clicked.connect(self._on_export_state_bp)

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

        csv_path = out / "analysis" / ct / "track_features" / f"{ct}_track_features.csv"
        if not csv_path.exists():
            self.warning_label.setText(
                f"⚠ Feature data not found for cell type '{ct}'.\n"
                "Run Feature Extraction first, then return here.\n"
                # f"Expected: {csv_path}" ## Removed: Too long for display
            )
            self.warning_label.show()
            for grp in [self.grp_train, self.grp2, self.grp3, self.grp_bp]:
                grp.setEnabled(False)
            return False
        else:
            self.warning_label.hide()
            for grp in [self.grp_train, self.grp2, self.grp3, self.grp_bp]:
                grp.setEnabled(True)
            return True

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
        self.spin_hmm_k_max.setVisible(auto)
        self.row_k_max.setVisible(auto)

    def _toggle_sticky_hmm(self, checked):
        self.spin_hmm_stickiness_kappa.setVisible(checked)
        self.row_kappa.setVisible(checked)
        self.spin_hmm_transmat_alpha.setVisible(checked)
        self.row_alpha.setVisible(checked)

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
            from behav3d.analysis.behavior.state.utils import _resolve_state_paths
            paths = _resolve_state_paths(out, ct)
            apply_hmm_deployment_artifact_to_full_dataset(
                deployment_artifact_path=Path(pkl_path),
                input_adata_path=paths.prepared_input_adata_path,
                output_adata_path=paths.full_output_adata_path,
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

    def _reload(self):
        try:
            ct = self._cell_type()
            if not ct:
                self._model_adata = None
                self._refresh_buttons()
                return

            self._populate_dynamic_features(ct)
            path = self._model_adata_path(ct)
            if path and path.exists():
                self._load_model_adata(path)
            else:
                self._model_adata = None

            self._check_prerequisites()
            self._refresh_buttons()
            self._update_view_buttons()
            self._update_bp_buttons()
        except Exception:
            traceback.print_exc()

    def _populate_dynamic_features(self, ct):
        from behav3d.widgets.utils import behav3d_calculated_features
        from behav3d.core.utils import expand_column_patterns
        import pandas as pd
        out = self._out_dir()
        if not out:
            return

        csv_path = out / "analysis" / ct / "track_features" / f"{ct}_track_features.csv"
        self._timepoint_checkboxes = {}
        self._logscale_checkboxes = {}
        self._bingrp_checkboxes = {}

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
            return

        df = pd.read_csv(csv_path, nrows=0)
        cols = list(df.columns)
        excluded = [
            "TrackID", "position_t", "position_x", "position_y", "position_z",
            "frame", "file", "index", "id", "sample_name", "Condition", "Timepoint"
        ]
        usable_cols = [c for c in cols if c not in excluded]
        bin_cols, feat_cols = [], []
        for c in usable_cols:
            dtype = df.dtypes.get(c)
            if (dtype == 'object' or dtype == 'bool'
                    or pd.api.types.is_string_dtype(dtype)
                    or pd.api.types.is_bool_dtype(dtype)):
                bin_cols.append(c)
            else:
                feat_cols.append(c)

        from copy import deepcopy
        base_groups = deepcopy(behav3d_calculated_features)
        matched = set()
        for gname, patterns in base_groups.items():
            vals = []
            for pat in patterns:
                vals.extend(expand_column_patterns(pat, feat_cols))
            clean_vals = sorted({x for x in vals if x in feat_cols})
            if clean_vals:
                self.timepoint_features_lay.addWidget(QLabel(f"<b>{gname}</b>"))
                grid = QGridLayout()
                for i, f in enumerate(clean_vals):
                    cb = QCheckBox(f)
                    self._timepoint_checkboxes[f] = cb
                    grid.addWidget(cb, i // 3, i % 3)
                self.timepoint_features_lay.addLayout(grid)
                matched.update(clean_vals)

        other = sorted([c for c in feat_cols if c not in matched])
        if other:
            self.timepoint_features_lay.addWidget(QLabel("<b>other</b>"))
            grid = QGridLayout()
            for i, f in enumerate(other):
                cb = QCheckBox(f)
                self._timepoint_checkboxes[f] = cb
                grid.addWidget(cb, i // 3, i % 3)
            self.timepoint_features_lay.addLayout(grid)

        grid_log = QGridLayout()
        for i, f in enumerate(feat_cols):
            cb = QCheckBox(f)
            self._logscale_checkboxes[f] = cb
            grid_log.addWidget(cb, i // 3, i % 3)
        self.log_scale_lay.addLayout(grid_log)

        if not bin_cols:
            self.bin_grp_lay.addWidget(QLabel("<i>No binary columns detected yet.</i>"))
        else:
            for b in bin_cols:
                cb = QCheckBox(b)
                self._bingrp_checkboxes[b] = cb
                self.bin_grp_lay.addWidget(cb)

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
        return (
            out / "analysis" / cell_type / "behavioral_states" / "results"
            / f"{report}_{cell_type}.pdf"
        )

    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

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

    def _update_view_buttons(self):
        ct = self._cell_type()
        if not ct:
            for btn in (self.btn_view_state, self.btn_view_composition, self.btn_view_transition):
                btn.setEnabled(False)
            return
        model_path = self._model_adata_path(ct)
        self.btn_view_state.setEnabled(bool(model_path and model_path.exists()))
        comp = self._report_path(ct, "state_composition_report")
        self.btn_view_composition.setEnabled(bool(comp and comp.exists()))
        trans = self._report_path(ct, "state_transition_report")
        self.btn_view_transition.setEnabled(bool(trans and trans.exists()))

    def _update_bp_buttons(self):
        ct = self._cell_type()
        state_path = self._full_adata_path(ct) if ct else None
        has_state = bool(state_path and state_path.exists())
        self.btn_show_state_bp.setEnabled(has_state)
        self.btn_export_state_bp.setEnabled(has_state)

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
            res = run_hmm_state_clustering(**params, verbose=True)
            self._hmm_model = res["hmm_model"]
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

    def _on_state_done(self, result):
        ct = self._cell_type()
        self._log(f"✅ State classification complete for '{ct}'.")
        path = self._model_adata_path(ct)
        if path and path.exists():
            self._load_model_adata(path)
        self._refresh_buttons()
        self._update_view_buttons()
        self._update_bp_buttons()
        self._notify_results()

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
            res = run_hmm_state_clustering(**params, verbose=True)
            self._hmm_model = res["hmm_model"]
            return res["model_adata"]

        def _done(result):
            self._on_state_done(result)
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
        dlg.exec_()
        self._refresh_buttons()
        self._update_view_buttons()

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
        dlg.exec_()
        self._refresh_buttons()
        self._update_view_buttons()

    def _on_state_composition(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        if getattr(self, "_model_adata", None) is None:
            QMessageBox.warning(self, "No data", "Run state classification first.")
            return
        out = self._out_dir()
        self._log(f"▶ Generating state composition report for '{ct}'…")
        selected_cols = [item.text() for item in self.list_composition_group_cols.selectedItems()]

        def _run(**kw):
            from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
                save_state_composition_report,
            )
            composition_dir = out / "state_composition"
            composition_dir.mkdir(parents=True, exist_ok=True)
            return save_state_composition_report(
                adata=self._model_adata,
                output_pdf_path=composition_dir / "state_composition_report.pdf",
                output_csv_path=composition_dir / "state_composition_report.csv",
                time_col="position_t",
                state_col="ClusterID",
                sample_col="sample_name",
                include_pooled_summary=True,
                group_cols=selected_cols,
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
        if getattr(self, "_model_adata", None) is None:
            QMessageBox.warning(self, "No data", "Run state classification first.")
            return
        out = self._out_dir()
        self._log(f"▶ Generating state transition report for '{ct}'…")

        def _run(**kw):
            from behav3d.analysis.behavior.state.visualization.plots.state_transitions import (
                save_state_transition_report,
            )
            transition_dir = out / "state_transitions"
            transition_dir.mkdir(parents=True, exist_ok=True)
            return save_state_transition_report(
                adata=self._model_adata,
                output_dir=transition_dir,
                state_col="ClusterID",
                time_col="position_t",
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

    # ── Backprojection ───────────────────────────────────────────────────

    def _on_show_state_bp(self):
        ct = self._cell_type()
        sample = self._sample()
        if not ct:
            QMessageBox.warning(self, "No cell type", "Select a cell type first.")
            return
        state_path = self._full_adata_path(ct)
        if not state_path or not state_path.exists():
            QMessageBox.warning(
                self, "No state adata",
                f"State adata not found:\n{state_path}\n\nRun State Classification first."
            )
            return
        color_by = self.combo_state_color_by.currentText()
        opacity = self.spin_state_opacity.value() / 100.0
        self._log(f"▶ Loading state backprojection for '{ct}' / sample '{sample}'…")
        try:
            from behav3d.napari.backprojection import show_behavioral_state_backprojection
            show_behavioral_state_backprojection(
                viewer=self.viewer,
                adata_path=str(state_path),
                cell_type=ct,
                sample=sample if sample else None,
                color_by=color_by,
                opacity=opacity,
                metadata=getattr(self.metadata_loader, "metadata", None),
                output_dir=str(self._out_dir()) if self._out_dir() else "",
            )
            self._log("✅ State backprojection loaded.")
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
        dpi = int(self.spin_state_export_dpi.value())
        make_pdf = self.chk_state_export_pdf.isChecked()
        make_mp4 = self.chk_state_export_mp4.isChecked()
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Exporting state backprojection for '{ct}'…")

        def _run(**kw):
            from behav3d.napari.backprojection import export_behavioral_state_backprojection
            return export_behavioral_state_backprojection(
                adata_path=str(state_path),
                cell_type=ct,
                sample=sample if sample else None,
                color_by=color_by,
                output_dir=str(out) if out else "",
                metadata=getattr(self.metadata_loader, "metadata", None),
                dpi=dpi,
                make_pdf=make_pdf,
                make_mp4=make_mp4,
                log_callback=logger,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Exporting state backprojection ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_export_state_bp, self.btn_show_state_bp],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: self._log("✅ State backprojection export done."),
            on_failed=lambda e: self._log(f"❌ Export failed: {e}"),
        )

    # ── View helpers ─────────────────────────────────────────────────────

    def _on_view(self, kind: str):
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
            menu.exec_()

    def _get_view_candidates(self, kind: str, ct: str):
        out = self._out_dir()
        if not out:
            return []
        if kind == "state_qc":
            p = self._model_adata_path(ct)
            return [(f"Model adata ({ct})", p)] if p else []
        if kind == "state_composition":
            p = self._report_path(ct, "state_composition_report")
            return [(f"Composition report ({ct})", p)] if p else []
        if kind == "state_transition":
            p = self._report_path(ct, "state_transition_report")
            return [(f"Transition report ({ct})", p)] if p else []
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
        self._bg = BackgroundOperation(self)

        # Queue buttons (wired by _widget.py)
        self.btn_queue_track_cluster = _make_queue_btn()
        self.btn_queue_train_track = _make_queue_btn()
        self.btn_queue_apply_track = _make_queue_btn()

        self._init_ui()

    # ── UI ───────────────────────────────────────────────────────────────

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)

        content = QWidget()
        lay = QVBoxLayout(content)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        scroll.setWidget(content)

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
        basic_form.addRow("Trajectory size:", make_help_row(
            self.spin_traj_size, "Trajectory size",
            "Number of timepoints to include per track (trims from the end by default). "
            "Shorter trajectories are discarded. Match the value used in your notebook."
        ))
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
        dtw_form = QFormLayout()
        dtw_form.setSpacing(3)

        self.combo_linkage = QComboBox()
        self.combo_linkage.addItems(["average", "complete", "single"])
        self.combo_linkage.setMaximumWidth(130)
        dtw_form.addRow("Linkage:", make_help_row(
            self.combo_linkage, "Linkage",
            "Agglomerative clustering linkage method. 'average' is the most commonly used; "
            "'complete' produces more compact clusters."
        ))

        self.combo_trim = QComboBox()
        self.combo_trim.addItems(["last", "first"])
        self.combo_trim.setMaximumWidth(130)
        dtw_form.addRow("Trim mode:", make_help_row(
            self.combo_trim, "Trim mode",
            "How to trim trajectories to behavioral_trajectory_size: "
            "'last' removes trailing timepoints; 'first' removes leading ones."
        ))
        self.adv1.addLayout(dtw_form)

        self.chk_parallel = QCheckBox("Parallel computation")
        self.chk_parallel.setChecked(True)
        self.adv1.addLayout(_make_chk_help_row(
            self.chk_parallel, "Parallel computation",
            "Use parallel computing for DTW distance matrix calculation. "
            "Faster on multi-core machines."
        ))

        self.chk_save_dist = QCheckBox("Save distance matrix CSV")
        self.chk_save_dist.setChecked(False)
        self.adv1.addLayout(_make_chk_help_row(
            self.chk_save_dist, "Save distance matrix",
            "Save the full DTW pairwise distance matrix to a CSV file. "
            "Can be large for many tracks."
        ))

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
        self.btn_view_track = _make_view_btn()
        run_row.addWidget(self.btn_view_track)
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

        rename_row = QHBoxLayout()
        self.btn_rename_track = QPushButton("✏  Rename Track Clusters")
        _style_rename(self.btn_rename_track)
        self.btn_rename_track.setEnabled(False)
        rename_row.addWidget(self.btn_rename_track, stretch=1)
        self.btn_view_rename_track = _make_view_btn()
        rename_row.addWidget(self.btn_view_rename_track)
        g2.addLayout(rename_row)
        lay.addWidget(self.grp2)

        # ── Step 3: Classify Tracks ──────────────────────────────────────
        self.grp3 = QGroupBox("Step 3 — Classify Tracks")
        g3 = QVBoxLayout(self.grp3)
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
        self.combo_track_max_features.addItems(["sqrt", "log2", "auto"])
        self.combo_track_max_features.setMaximumWidth(100)
        adv3_form.addRow("Max features:", make_help_row(
            self.combo_track_max_features, "Max features",
            "Number of features to consider at each split: 'sqrt' (recommended), "
            "'log2', or 'auto' (same as sqrt)."
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
        self.btn_view_train_track = _make_view_btn()
        train_rf_run_row.addWidget(self.btn_view_train_track)
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
        self.btn_view_apply_track = _make_view_btn()
        apply_clf_run_row.addWidget(self.btn_view_apply_track)
        g_apply_clf.addLayout(apply_clf_run_row)
        g3.addWidget(grp_apply_clf)

        lay.addWidget(self.grp3)

        # ── Step 4: Create Plots ─────────────────────────────────────────
        self.grp4 = QGroupBox("Step 4 — Create Plots")
        g4 = QVBoxLayout(self.grp4)
        g4.setSpacing(6)

        # Diagnostic sub-section
        grp_diag = QGroupBox("Diagnostic")
        g_diag = QVBoxLayout(grp_diag)
        g_diag.setSpacing(4)

        diag_row = QHBoxLayout()
        self.btn_diagnostics = QPushButton("▶ Create Diagnostics")
        _style_secondary(self.btn_diagnostics)
        diag_row.addWidget(self.btn_diagnostics, stretch=1)
        self.btn_view_diagnostics = _make_view_btn()
        diag_row.addWidget(self.btn_view_diagnostics)
        g_diag.addLayout(diag_row)
        g4.addWidget(grp_diag)

        # Exemplar Tracks sub-section
        grp_exemplar = QGroupBox("Exemplar Tracks")
        g_exemplar = QVBoxLayout(grp_exemplar)
        g_exemplar.setSpacing(4)

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
        g4.addWidget(grp_exemplar)
        lay.addWidget(self.grp4)

        # ── Step 5: Backprojection ───────────────────────────────────────
        self.grp_bp = QGroupBox("Step 5 — Backprojection")
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
        self.combo_track_color_by.addItems(
            ["behavioral_trajectory_cluster", "dtaidistance_cluster", "track_cluster"]
        )
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
        g_track_view.addLayout(tbp_form)

        track_view_row = QHBoxLayout()
        self.btn_show_track_bp = QPushButton("▶ Show Track Backprojection in Napari")
        _style_primary(self.btn_show_track_bp)
        track_view_row.addWidget(self.btn_show_track_bp, stretch=1)
        g_track_view.addLayout(track_view_row)
        g_bp5.addWidget(grp_track_view)

        # Export options (collapsed)
        track_export_sec = CollapsibleSection("⚙ Export Options", expanded=False)
        track_export_form = QFormLayout()
        self.spin_track_export_dpi = QSpinBox()
        self.spin_track_export_dpi.setRange(50, 600)
        self.spin_track_export_dpi.setValue(150)
        self.spin_track_export_dpi.setSuffix(" dpi")
        self.spin_track_export_dpi.setMaximumWidth(110)
        track_export_form.addRow("DPI (PDF):", make_help_row(
            self.spin_track_export_dpi, "Export DPI",
            "Resolution for exported PDF images."
        ))
        track_export_sec.addLayout(track_export_form)

        track_export_opts = QHBoxLayout()
        self.chk_track_pdf = QCheckBox("PDF")
        self.chk_track_pdf.setChecked(True)
        track_export_opts.addWidget(self.chk_track_pdf)
        self.chk_track_mp4 = QCheckBox("MP4")
        self.chk_track_mp4.setChecked(False)
        track_export_opts.addWidget(self.chk_track_mp4)
        track_export_opts.addStretch()
        track_export_sec.addLayout(track_export_opts)

        track_export_run_row = QHBoxLayout()
        self.btn_export_track_bp = QPushButton("▶ Export Track Backprojection")
        _style_secondary(self.btn_export_track_bp)
        track_export_run_row.addWidget(self.btn_export_track_bp, stretch=1)
        track_export_sec.addLayout(track_export_run_row)
        g_bp5.addWidget(track_export_sec)
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
        self.chk_apply_pretrained.toggled.connect(self._toggle_pretrained_mode)
        self.chk_use_original.toggled.connect(self._on_original_toggled_from_adv)
        self.chk_use_original_top.toggled.connect(self._on_original_toggled_from_top)
        self.btn_run_track.clicked.connect(self._on_run_cluster)
        self.btn_view_track.clicked.connect(lambda: self._on_view("track_adata"))
        self.btn_rename_track.clicked.connect(self._on_rename_track)
        self.btn_train_track.clicked.connect(self._on_train_track)
        self.btn_view_train_track.clicked.connect(lambda: self._on_view("track_classifier"))
        self.btn_apply_track.clicked.connect(self._on_apply_track)
        self.btn_view_apply_track.clicked.connect(lambda: self._on_view("track_applied"))
        self.btn_exemplars.clicked.connect(self._on_exemplars)
        self.btn_view_exemplars.clicked.connect(lambda: self._on_view("track_exemplars"))
        self.btn_diagnostics.clicked.connect(self._on_diagnostics)
        self.btn_view_diagnostics.clicked.connect(lambda: self._on_view("track_diagnostics"))
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
            self.grp1.setEnabled(True)
            for grp in [self.grp2, self.grp3, self.grp4, self.grp_bp]:
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
            for grp in [self.grp1, self.grp2, self.grp3, self.grp4, self.grp_bp]:
                grp.setEnabled(True)

            # Re-enable the checkboxes and revert to standard mode automatically.
            for chk in (self.chk_use_original, self.chk_use_original_top):
                chk.blockSignals(True)
                chk.setChecked(False)
                chk.setEnabled(True)
                chk.blockSignals(False)
            self._apply_original_mode(False)

            return True

    # ── Toggle helpers ───────────────────────────────────────────────────

    def _toggle_pretrained_mode(self, checked: bool):
        """Hide Steps 1–3 when apply-pretrained mode is active."""
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

    def _on_original_toggled_from_top(self, checked: bool):
        """Top checkbox toggled (user unchecks from top) → sync adv checkbox + apply mode."""
        self.chk_use_original.blockSignals(True)
        self.chk_use_original.setChecked(checked)
        self.chk_use_original.blockSignals(False)
        self._apply_original_mode(checked)

    def _apply_original_mode(self, checked: bool):
        """Update visibility of UI sections for original vs dtaidistance mode."""
        self.chk_use_original_top.setVisible(checked)
        self.adv1.setVisible(not checked)
        self._umap_frame.setVisible(checked)
        if checked:
            self.btn_run_track.setText("▶ Run Original BEHAV3D DTW")
        else:
            self.btn_run_track.setText("▶ Run Track Clustering")

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

    # ── Metadata / reload ────────────────────────────────────────────────

    def on_metadata_updated(self):
        self._reload()

    def _reload(self):
        ct = self._cell_type()
        if not ct:
            self._track_adata = None
            self._refresh_buttons()
            return
        path = self._track_adata_path(ct)
        if path and path.exists():
            try:
                import anndata as ad
                self._track_adata = ad.read_h5ad(str(path))
            except Exception:
                traceback.print_exc()
                self._track_adata = None
        else:
            self._track_adata = None

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

        self._check_prerequisites()
        self._refresh_buttons()
        self._update_view_buttons()
        self._update_bp_buttons()

    # ── Path helpers ─────────────────────────────────────────────────────

    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

    def _track_adata_path(self, ct: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        d = out / "analysis" / ct / "behavorial_trajectories"
        if d.exists():
            for f in sorted(d.glob("*.h5ad")):
                return f
        return d / f"BEHAV3D_{ct}_track_trajectories.h5ad"

    def _track_classifier_path(self, ct: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return out / "analysis" / ct / "behavorial_trajectories" / f"classifier_{ct}.pkl"

    def _track_applied_path(self, ct: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return (
            out / "analysis" / ct / "behavorial_trajectories"
            / f"BEHAV3D_{ct}_track_clusters.csv"
        )

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
        if has_adata:
            self.rename_track_status.setText(
                f"✅ Track adata loaded: {self._track_adata.n_obs} rows."
            )
        else:
            self.rename_track_status.setText("ℹ Run clustering first to enable renaming.")

    def _update_view_buttons(self):
        ct = self._cell_type()
        if not ct:
            for btn in (
                self.btn_view_track, self.btn_view_train_track,
                self.btn_view_apply_track, self.btn_view_exemplars,
                self.btn_view_diagnostics,
            ):
                btn.setEnabled(False)
            return
        adata_path = self._track_adata_path(ct)
        self.btn_view_track.setEnabled(bool(adata_path and adata_path.exists()))
        clf_path = self._track_classifier_path(ct)
        self.btn_view_train_track.setEnabled(bool(clf_path and clf_path.exists()))
        app_path = self._track_applied_path(ct)
        self.btn_view_apply_track.setEnabled(bool(app_path and app_path.exists()))
        out = self._out_dir()
        traj_dir = (out / "analysis" / ct / "behavorial_trajectories") if out else None
        self.btn_view_exemplars.setEnabled(
            bool(traj_dir and any(traj_dir.glob("exemplar_tracks*.pdf")))
        )
        self.btn_view_diagnostics.setEnabled(
            bool(traj_dir and any(traj_dir.glob("diagnostics*.pdf")))
        )

    def _update_bp_buttons(self):
        ct = self._cell_type()
        track_path = self._track_adata_path(ct) if ct else None
        has_track = bool(track_path and track_path.exists())
        self.btn_show_track_bp.setEnabled(has_track)
        self.btn_export_track_bp.setEnabled(has_track)

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
        params = {
            "output_dir": str(out) if out else "",
            "cell_type": ct,
            "behavioral_trajectory_size": int(self.spin_traj_size.value()),
            "n_clusters": int(self.spin_n_clusters.value()),
            "random_state": int(self.spin_seed.value()),
            "linkage": self.combo_linkage.currentText(),
            "trajectory_trim_mode": self.combo_trim.currentText(),
            "parallel": self.chk_parallel.isChecked(),
            "save_distance_matrix": self.chk_save_dist.isChecked(),
        }

        on_done_ext = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_ext = extra_callbacks.get("on_failed") if extra_callbacks else None

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import (
                run_categorical_dtaidistance_trajectory_clustering,
            )
            return run_categorical_dtaidistance_trajectory_clustering(**params, verbose=True)

        def _done(r):
            self._log(f"✅ Track clustering done for '{ct}'.")
            self._reload()
            self._notify_results()
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
            )

        self._bg.run(
            fn=_run,
            desc=f"Original BEHAV3D DTW ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_run_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Original BEHAV3D clustering done for '{ct}'."),
                self._reload(),
                self._notify_results(),
            ),
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
        dlg.exec_()
        self._refresh_buttons()
        self._update_view_buttons()

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
        logger = ThreadSafeLogger(self._log)
        n_per = int(self.spin_n_per_cluster.value())
        make_overview = self.chk_overview_statebars.isChecked()
        make_bp_pdf = self.chk_backproj_pdf.isChecked()
        make_bp_mp4 = self.chk_backproj_mp4.isChecked()
        self._log(f"▶ Creating exemplar PDFs for '{ct}'…")
        track_adata = self._track_adata

        def _run(**kw):
            from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
                save_exemplar_statebar_backprojection_pdf,
            )
            return save_exemplar_statebar_backprojection_pdf(
                adata=track_adata,
                output_dir=str(out) if out else "",
                cell_type=ct,
                metadata=getattr(self.metadata_loader, "metadata", None),
                n_per_cluster=n_per,
                make_overview_statebars=make_overview,
                make_backprojection_pdf=make_bp_pdf,
                make_backprojection_mp4=make_bp_mp4,
                log_callback=logger,
                verbose=True,
            )

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
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Creating diagnostics for '{ct}'…")
        track_adata = self._track_adata

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import save_dtaidistance_diagnostics
            return save_dtaidistance_diagnostics(
                adata=track_adata,
                output_dir=str(out) if out else "",
                cell_type=ct,
                log_callback=logger,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Track diagnostics ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_diagnostics],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Diagnostics done for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Diagnostics failed: {e}"),
        )

    # ── Backprojection ───────────────────────────────────────────────────

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
        color_by = self.combo_track_color_by.currentText()
        opacity = self.spin_track_opacity.value() / 100.0
        self._log(f"▶ Loading track backprojection for '{ct}' / sample '{sample}'…")
        try:
            from behav3d.napari.backprojection import show_track_cluster_backprojection
            show_track_cluster_backprojection(
                viewer=self.viewer,
                adata_path=str(track_path),
                cell_type=ct,
                sample=sample if sample else None,
                color_by=color_by,
                opacity=opacity,
                metadata=getattr(self.metadata_loader, "metadata", None),
                output_dir=str(self._out_dir()) if self._out_dir() else "",
            )
            self._log("✅ Track backprojection loaded.")
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
        color_by = self.combo_track_color_by.currentText()
        dpi = int(self.spin_track_export_dpi.value())
        make_pdf = self.chk_track_pdf.isChecked()
        make_mp4 = self.chk_track_mp4.isChecked()
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Exporting track backprojection for '{ct}'…")

        def _run(**kw):
            from behav3d.napari.backprojection import export_track_cluster_backprojection
            return export_track_cluster_backprojection(
                adata_path=str(track_path),
                cell_type=ct,
                sample=sample if sample else None,
                color_by=color_by,
                output_dir=str(out) if out else "",
                metadata=getattr(self.metadata_loader, "metadata", None),
                dpi=dpi,
                make_pdf=make_pdf,
                make_mp4=make_mp4,
                log_callback=logger,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Exporting track backprojection ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_export_track_bp, self.btn_show_track_bp],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: self._log("✅ Track backprojection export done."),
            on_failed=lambda e: self._log(f"❌ Export failed: {e}"),
        )

    # ── View helpers ─────────────────────────────────────────────────────

    def _on_view(self, kind: str):
        ct = self._cell_type()
        if not ct:
            return
        out = self._out_dir()
        traj_dir = (out / "analysis" / ct / "behavorial_trajectories") if out else None
        candidates = []
        if kind == "track_adata":
            p = self._track_adata_path(ct)
            if p:
                candidates = [(f"Track adata ({ct})", p)]
        elif kind == "track_classifier":
            p = self._track_classifier_path(ct)
            if p:
                candidates = [(f"Track classifier ({ct})", p)]
        elif kind == "track_applied":
            p = self._track_applied_path(ct)
            if p:
                candidates = [(f"Track clusters CSV ({ct})", p)]
        elif kind == "track_exemplars" and traj_dir:
            candidates = [(f.stem, f) for f in sorted(traj_dir.glob("exemplar_tracks*.pdf"))]
        elif kind == "track_diagnostics" and traj_dir:
            candidates = [(f.stem, f) for f in sorted(traj_dir.glob("diagnostics*.pdf"))]

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

        # ── Inner sub-tabs ───────────────────────────────────────────────
        self.inner_tabs = QTabWidget()
        outer.addWidget(self.inner_tabs, stretch=1)

        self.state_tab = StateClassificationSubTab(
            viewer=self.viewer,
            metadata_loader=self.metadata_loader,
            cell_type_getter=self._current_cell_type,
            parent=self,
        )
        self.inner_tabs.addTab(self.state_tab, "🔬 State Classification")

        self.track_tab = TrackClassificationSubTab(
            viewer=self.viewer,
            metadata_loader=self.metadata_loader,
            cell_type_getter=self._current_cell_type,
            parent=self,
        )
        self.inner_tabs.addTab(self.track_tab, "🛤️ Track Classification")

        # When switching to Track Classification (index 1), trigger a reload
        # so the behavioral-states path auto-fills from disk if it exists.
        self.inner_tabs.currentChanged.connect(self._on_inner_tab_changed)

    def _on_inner_tab_changed(self, index: int):
        """Reload track tab when switched to, so states path auto-fills from disk."""
        if index == 1:  # Track Classification
            self.track_tab.on_metadata_updated()

    # ── Metadata update ──────────────────────────────────────────────────

    def _on_metadata_updated(self, *_):
        cell_types = _detect_sc_cell_types(self.metadata_loader)
        current = self.cell_type_combo.currentText()
        self.cell_type_combo.blockSignals(True)
        self.cell_type_combo.clear()
        if cell_types:
            self.cell_type_combo.addItems(cell_types)
            if current in cell_types:
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

    def _propagate_metadata_update(self):
        self.state_tab.on_metadata_updated()
        self.track_tab.on_metadata_updated()

    def _on_cell_type_changed(self, _text: str):
        self._propagate_metadata_update()

    def _current_cell_type(self) -> str:
        return self.cell_type_combo.currentText()
