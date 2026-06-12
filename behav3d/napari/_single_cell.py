"""
BEHAV3D napari plugin – Single Cell Analysis Tab.

Provides two inner sub-tabs:

- **🔬 State Classification** — HMM-based behavioral state clustering,
  intrinsic + full cluster renaming, random-forest train / apply,
  and composition / transition reports.
- **🛤️ Track Classification** — DTW-based trajectory clustering (dtaidistance
  one-hot or original feature-based BEHAV3D mode), track cluster renaming,
  random-forest train / apply, and exemplar / diagnostics plots.

Cell-type scope: immune + other types only (non-multicolor), mirroring the
rest of the pipeline from Feature Extraction onwards.

Background execution: heavy steps use :class:`BackgroundOperation` from
``_background_runner.py``.  Computationally intensive steps expose a
``+🛒`` queue button (wired by ``_widget.py``).  Rename steps open a modal
:class:`RenameClusterDialog` popup.  Result-view (👁) buttons are enabled
once the expected output file exists on disk.
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
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
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
    QFrame,
)

from behav3d.napari._analysis import CollapsibleSection
from behav3d.napari._background_runner import (
    BackgroundOperation,
    ProgressBarRow,
    ThreadSafeLogger,
)
from behav3d.napari._pdf_view import open_pdf_in_napari
from behav3d.napari._rename_dialog import RenameClusterDialog
from behav3d.core.qt_help import make_help_row


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


# ═══════════════════════════════════════════════════════════════════════════
# State Classification sub-tab
# ═══════════════════════════════════════════════════════════════════════════

class StateClassificationSubTab(QWidget):
    """Steps 1-4 for HMM-based behavioral state classification."""

    def __init__(self, viewer=None, metadata_loader=None, cell_type_getter=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._get_cell_type = cell_type_getter  # callable → str
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

        # ── Step 1: Clustering ───────────────────────────────────────
        grp1 = QGroupBox("Step 1 — State Clustering")
        g1 = QVBoxLayout(grp1)
        g1.setSpacing(4)

        row1 = QHBoxLayout()
        self.btn_run_state = QPushButton("▶ Run State Classification")
        _style_primary(self.btn_run_state)
        self.btn_run_state.clicked.connect(self._on_run_state)
        row1.addWidget(self.btn_run_state, stretch=1)
        row1.addWidget(self.btn_queue_state_cluster)
        self.btn_view_state = _make_view_btn()
        self.btn_view_state.clicked.connect(lambda: self._on_view("state_qc"))
        row1.addWidget(self.btn_view_state)
        g1.addLayout(row1)

        # Advanced config
        adv1 = CollapsibleSection("⚙  Advanced Configuration", expanded=False)

        # Feature Selection sub-section
        feat_frame = QFrame()
        feat_frame.setFrameShape(QFrame.StyledPanel)
        feat_frame.setStyleSheet("QFrame { background: #222233; border: 1px solid #333; border-radius: 2px; }")
        feat_lay = QVBoxLayout(feat_frame)
        feat_lay.setContentsMargins(6, 4, 6, 4)
        feat_lay.setSpacing(2)
        feat_lay.addWidget(QLabel("Feature groups  (leave all unchecked to use defaults from params file)"))

        self._feat_group_checks: dict[str, QCheckBox] = {}
        feat_groups_row = QHBoxLayout()
        for gname in ("speed", "displacement", "morphology", "contact"):
            chk = QCheckBox(gname)
            chk.setChecked(False)
            self._feat_group_checks[gname] = chk
            feat_groups_row.addWidget(chk)
        feat_groups_row.addStretch()
        feat_lay.addLayout(feat_groups_row)

        self.chk_use_binary_groups = QCheckBox("Include binary group labels as features")
        self.chk_use_binary_groups.setChecked(True)
        feat_lay.addWidget(self.chk_use_binary_groups)
        adv1.addWidget(feat_frame)

        # Windowed processing sub-section
        win_form = QFormLayout()
        win_form.setSpacing(3)
        self.spin_window_size = QSpinBox()
        self.spin_window_size.setRange(1, 500)
        self.spin_window_size.setValue(5)
        self.spin_window_size.setMaximumWidth(80)
        win_form.addRow("Trailing window size:", make_help_row(
            self.spin_window_size, "Trailing window size", "Number of prior timepoints to consider for computing sliding window statistics."
        ))

        self.combo_window_partial = QComboBox()
        self.combo_window_partial.addItems(["partial", "full", "ignore"])
        self.combo_window_partial.setMaximumWidth(120)
        win_form.addRow("Partial window policy:", make_help_row(
            self.combo_window_partial, "Partial window policy", "How to handle frames where a full trailing window is not available (e.g., at track start)."
        ))

        quant_row = QHBoxLayout()
        self.spin_quant_lo = QDoubleSpinBox()
        self.spin_quant_lo.setRange(0.0, 0.5)
        self.spin_quant_lo.setSingleStep(0.01)
        self.spin_quant_lo.setDecimals(3)
        self.spin_quant_lo.setValue(0.0)
        self.spin_quant_lo.setMaximumWidth(90)
        quant_row.addWidget(QLabel("Quantile cap lower:"))
        quant_row.addWidget(self.spin_quant_lo)
        self.spin_quant_hi = QDoubleSpinBox()
        self.spin_quant_hi.setRange(0.5, 1.0)
        self.spin_quant_hi.setSingleStep(0.01)
        self.spin_quant_hi.setDecimals(3)
        self.spin_quant_hi.setValue(0.99)
        self.spin_quant_hi.setMaximumWidth(90)
        quant_row.addWidget(QLabel("upper:"))
        quant_row.addWidget(self.spin_quant_hi)
        quant_row.addStretch()
        adv1.addLayout(win_form)
        adv1.addLayout(quant_row)

        # Clustering params
        clust_form = QFormLayout()
        clust_form.setSpacing(3)
        self.combo_clust_method = QComboBox()
        self.combo_clust_method.addItems(["leiden", "kmeans"])
        self.combo_clust_method.setMaximumWidth(120)
        clust_form.addRow("Clustering method:", self.combo_clust_method)

        self.spin_resolution = QDoubleSpinBox()
        self.spin_resolution.setRange(0.01, 10.0)
        self.spin_resolution.setSingleStep(0.05)
        self.spin_resolution.setDecimals(2)
        self.spin_resolution.setValue(0.2)
        self.spin_resolution.setMaximumWidth(90)
        clust_form.addRow("Resolution (leiden):", make_help_row(
            self.spin_resolution, "Resolution", "Leiden resolution parameter. Higher values yield more clusters."
        ))

        self.spin_neighbors = QSpinBox()
        self.spin_neighbors.setRange(2, 500)
        self.spin_neighbors.setValue(60)
        self.spin_neighbors.setMaximumWidth(80)
        clust_form.addRow("Neighbors:", make_help_row(
            self.spin_neighbors, "Neighbors", "Number of neighbors to use for kNN graph construction."
        ))

        self.spin_umap_min_dist = QDoubleSpinBox()
        self.spin_umap_min_dist.setRange(0.001, 1.0)
        self.spin_umap_min_dist.setSingleStep(0.05)
        self.spin_umap_min_dist.setDecimals(3)
        self.spin_umap_min_dist.setValue(0.1)
        self.spin_umap_min_dist.setMaximumWidth(90)
        clust_form.addRow("UMAP min_dist:", self.spin_umap_min_dist)

        self.spin_pca_var = QDoubleSpinBox()
        self.spin_pca_var.setRange(0.5, 1.0)
        self.spin_pca_var.setSingleStep(0.01)
        self.spin_pca_var.setDecimals(2)
        self.spin_pca_var.setValue(0.95)
        self.spin_pca_var.setMaximumWidth(90)
        clust_form.addRow("PCA variance:", self.spin_pca_var)

        self.spin_seed = QSpinBox()
        self.spin_seed.setRange(0, 99999)
        self.spin_seed.setValue(123)
        self.spin_seed.setMaximumWidth(90)
        clust_form.addRow("Random seed:", self.spin_seed)

        self.chk_reuse_prepared = QCheckBox("Reuse already prepared dataset (skip PCA/UMAP)")
        self.chk_reuse_prepared.setChecked(False)
        adv1.addLayout(clust_form)
        adv1.addWidget(self.chk_reuse_prepared)

        g1.addWidget(adv1)
        lay.addWidget(grp1)

        # ── Step 2: Rename Clusters ──────────────────────────────────
        grp2 = QGroupBox("Step 2 — Rename Clusters")
        g2 = QVBoxLayout(grp2)
        g2.setSpacing(4)

        self.rename_status_lbl = QLabel(
            "ℹ Run state classification first to enable renaming."
        )
        self.rename_status_lbl.setStyleSheet("color: #999; font-size: 11px;")
        self.rename_status_lbl.setWordWrap(True)
        g2.addWidget(self.rename_status_lbl)

        rename_intrinsic_row = QHBoxLayout()
        self.btn_rename_intrinsic = QPushButton("✏  Rename Primary Dynamic State Clusters")
        _style_rename(self.btn_rename_intrinsic)
        self.btn_rename_intrinsic.setEnabled(False)
        self.btn_rename_intrinsic.clicked.connect(self._on_rename_intrinsic)
        rename_intrinsic_row.addWidget(self.btn_rename_intrinsic, stretch=1)
        self.btn_view_rename_intrinsic = _make_view_btn()
        rename_intrinsic_row.addWidget(self.btn_view_rename_intrinsic)
        g2.addLayout(rename_intrinsic_row)

        rename_full_row = QHBoxLayout()
        self.btn_rename_full = QPushButton("✏  Rename Full Behavioral Clusters (Binary Groups)")
        _style_rename(self.btn_rename_full)
        self.btn_rename_full.setEnabled(False)
        self.btn_rename_full.clicked.connect(self._on_rename_full)
        rename_full_row.addWidget(self.btn_rename_full, stretch=1)
        self.btn_view_rename_full = _make_view_btn()
        rename_full_row.addWidget(self.btn_view_rename_full)
        g2.addLayout(rename_full_row)

        lay.addWidget(grp2)

        # ── Step 3: Train & Apply ────────────────────────────────────
        grp3 = QGroupBox("Step 3 — Train & Apply Classifier")
        g3 = QVBoxLayout(grp3)
        g3.setSpacing(4)

        train_row = QHBoxLayout()
        self.btn_train = QPushButton("▶ Train State Classifier")
        _style_primary(self.btn_train)
        self.btn_train.setEnabled(False)
        self.btn_train.clicked.connect(self._on_train)
        train_row.addWidget(self.btn_train, stretch=1)
        train_row.addWidget(self.btn_queue_train_state)
        self.btn_view_train = _make_view_btn()
        self.btn_view_train.clicked.connect(lambda: self._on_view("state_classifier"))
        train_row.addWidget(self.btn_view_train)
        g3.addLayout(train_row)

        apply_row = QHBoxLayout()
        self.btn_apply = QPushButton("▶ Apply Classifier to Full Dataset")
        _style_secondary(self.btn_apply)
        self.btn_apply.setEnabled(False)
        self.btn_apply.clicked.connect(self._on_apply)
        apply_row.addWidget(self.btn_apply, stretch=1)
        apply_row.addWidget(self.btn_queue_apply_state)
        self.btn_view_apply = _make_view_btn()
        self.btn_view_apply.clicked.connect(lambda: self._on_view("state_full"))
        apply_row.addWidget(self.btn_view_apply)
        g3.addLayout(apply_row)

        adv3 = CollapsibleSection("⚙  Advanced Configuration (Classifier)", expanded=False)
        clf_form = QFormLayout()
        clf_form.setSpacing(3)
        self.spin_clf_n_est = QSpinBox()
        self.spin_clf_n_est.setRange(10, 2000)
        self.spin_clf_n_est.setValue(100)
        self.spin_clf_n_est.setMaximumWidth(90)
        clf_form.addRow("n_estimators:", self.spin_clf_n_est)

        self.spin_clf_test_size = QDoubleSpinBox()
        self.spin_clf_test_size.setRange(0.05, 0.5)
        self.spin_clf_test_size.setSingleStep(0.05)
        self.spin_clf_test_size.setDecimals(2)
        self.spin_clf_test_size.setValue(0.2)
        self.spin_clf_test_size.setMaximumWidth(90)
        clf_form.addRow("Test holdout:", self.spin_clf_test_size)

        self.spin_clf_n_jobs = QSpinBox()
        self.spin_clf_n_jobs.setRange(-1, 64)
        self.spin_clf_n_jobs.setValue(-1)
        self.spin_clf_n_jobs.setMaximumWidth(80)
        clf_form.addRow("n_jobs:", self.spin_clf_n_jobs)

        adv3.addLayout(clf_form)
        g3.addWidget(adv3)
        lay.addWidget(grp3)

        # ── Step 4: Reports ──────────────────────────────────────────
        grp4 = QGroupBox("Step 4 — Reports")
        g4 = QVBoxLayout(grp4)
        g4.setSpacing(4)

        comp_row = QHBoxLayout()
        self.btn_state_composition = QPushButton("▶ State Composition Report")
        _style_secondary(self.btn_state_composition)
        self.btn_state_composition.clicked.connect(self._on_state_composition)
        comp_row.addWidget(self.btn_state_composition, stretch=1)
        self.btn_view_composition = _make_view_btn()
        self.btn_view_composition.clicked.connect(lambda: self._on_view("state_composition"))
        comp_row.addWidget(self.btn_view_composition)
        g4.addLayout(comp_row)

        trans_row = QHBoxLayout()
        self.btn_state_transition = QPushButton("▶ State Transition Report")
        _style_secondary(self.btn_state_transition)
        self.btn_state_transition.clicked.connect(self._on_state_transition)
        trans_row.addWidget(self.btn_state_transition, stretch=1)
        self.btn_view_transition = _make_view_btn()
        self.btn_view_transition.clicked.connect(lambda: self._on_view("state_transition"))
        trans_row.addWidget(self.btn_view_transition)
        g4.addLayout(trans_row)

        lay.addWidget(grp4)

        # ── Progress + Log ───────────────────────────────────────────
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

    # ── Metadata update ─────────────────────────────────────────────────
    def on_metadata_updated(self):
        self._reload()

    def _reload(self):
        """Re-check disk for existing model adata and refresh button states."""
        try:
            ct = self._get_cell_type() if self._get_cell_type else ""
            if not ct:
                self._model_adata = None
                self._refresh_buttons()
                return
            path = self._model_adata_path(ct)
            if path and path.exists():
                self._load_model_adata(path)
            else:
                self._model_adata = None
            self._refresh_buttons()
            self._update_view_buttons()
        except Exception:
            traceback.print_exc()

    def _model_adata_path(self, cell_type: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return (
            out
            / "analysis"
            / cell_type
            / "behavioral_states"
            / "processing"
            / f"BEHAV3D_{cell_type}_behavioral_states_modeldata.h5ad"
        )

    def _full_adata_path(self, cell_type: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return (
            out
            / "analysis"
            / cell_type
            / "behavioral_states"
            / f"BEHAV3D_{cell_type}_behavioral_states.h5ad"
        )

    def _classifier_path(self, cell_type: str, kind="full") -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        return (
            out
            / "analysis"
            / cell_type
            / "behavioral_states"
            / "processing"
            / "full_behavioral_classification"
            / f"state_classification_random_forest_{kind}_{cell_type}.pkl"
        )

    def _report_path(self, cell_type: str, report: str) -> Optional[Path]:
        """Return expected PDF path for composition / transition report."""
        out = self._out_dir()
        if not out:
            return None
        return (
            out
            / "analysis"
            / cell_type
            / "behavioral_states"
            / "results"
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

    # ── Button state management ─────────────────────────────────────────
    def _refresh_buttons(self):
        has_model = self._model_adata is not None
        has_intrinsic = has_model and "intrinsic_behavioral_cluster" in self._model_adata.obs.columns
        has_full = has_model and "full_behavioral_cluster" in self._model_adata.obs.columns

        self.btn_rename_intrinsic.setEnabled(has_intrinsic)
        self.btn_rename_full.setEnabled(has_full)
        self.btn_train.setEnabled(has_model)
        self.btn_apply.setEnabled(has_model)

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
            self.rename_status_lbl.setText(
                "ℹ Run state classification first to enable renaming."
            )

    def _update_view_buttons(self):
        ct = self._get_cell_type() if self._get_cell_type else ""
        if not ct:
            for btn in (
                self.btn_view_state,
                self.btn_view_train,
                self.btn_view_apply,
                self.btn_view_composition,
                self.btn_view_transition,
            ):
                btn.setEnabled(False)
            return

        # State QC: check for model adata
        model_path = self._model_adata_path(ct)
        self.btn_view_state.setEnabled(bool(model_path and model_path.exists()))

        # Classifier
        clf_path = self._classifier_path(ct, "full")
        self.btn_view_train.setEnabled(bool(clf_path and clf_path.exists()))

        # Full adata
        full_path = self._full_adata_path(ct)
        self.btn_view_apply.setEnabled(bool(full_path and full_path.exists()))

        # Reports
        comp = self._report_path(ct, "state_composition_report")
        self.btn_view_composition.setEnabled(bool(comp and comp.exists()))

        trans = self._report_path(ct, "state_transition_report")
        self.btn_view_transition.setEnabled(bool(trans and trans.exists()))

    # ── Click handlers ──────────────────────────────────────────────────

    def _cell_type(self) -> str:
        return self._get_cell_type() if self._get_cell_type else ""

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
        logger = ThreadSafeLogger(self._log)

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import (
                run_hmm_state_clustering,
            )
            res = run_hmm_state_clustering(
                **params,
                verbose=True,
            )
            self._hmm_model = res["hmm_model"]
            return res["model_adata"]

        self._bg.run(
            fn=_run,
            desc=f"State classification ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_run_state, self.btn_train, self.btn_apply],
            viewer=self.viewer,
            inject_progress=False,
            on_done=self._on_state_done,
            on_failed=lambda e: self._log(f"❌ State classification failed: {e}"),
        )

    def _collect_state_params(self, ct: str) -> dict:
        out = self._out_dir()
        params_dict = getattr(self.metadata_loader, "behav3d_parameters", {})
        if not isinstance(params_dict, dict):
            params_dict = {}
        section = params_dict.get("behavioral_state_classification", {})
        defaults = section.get("defaults", {})
        cell_cfg = section.get(ct, {})
        eff_cfg = {**defaults, **cell_cfg}
        
        features = eff_cfg.get("features", [])
        if not features:
            features = ["speed", "displacement", "direction_autocorr"]
            
        lower_cap = float(self.spin_quant_lo.value())
        upper_cap = float(self.spin_quant_hi.value())
        
        return {
            "output_dir": str(out) if out else "",
            "cell_type": ct,
            "features": features,
            "binary_features_to_group": eff_cfg.get("binary_features_to_group", []),
            "additional_window_features": eff_cfg.get("additional_window_features", []),
            "log_scale_features": eff_cfg.get("log_scale_features", []),
            "n_states": "auto",
            "k_min": int(eff_cfg.get("hmm_k_min", 1)),
            "k_max": int(eff_cfg.get("hmm_k_max", 8)),
            "covariance_type": str(eff_cfg.get("hmm_covariance_type", "full")),
            "n_iter": int(eff_cfg.get("hmm_n_iter", 200)),
            "tol": float(eff_cfg.get("hmm_tol", 1e-3)),
            "sticky": bool(eff_cfg.get("hmm_sticky", False)),
            "stickiness_kappa": float(eff_cfg.get("hmm_stickiness_kappa", 8.0)),
            "transmat_alpha": float(eff_cfg.get("hmm_transmat_alpha", 1.0)),
            "min_covar": float(eff_cfg.get("hmm_min_covar", 1e-3)),
            "feature_smoothing_window": int(eff_cfg.get("hmm_feature_smoothing_window", 1)),
            "smoothing_min_periods": int(eff_cfg.get("hmm_smoothing_min_periods", 1)),
            "start_offset": int(eff_cfg.get("hmm_start_offset", 0)),
            "start_offset_fill_mode": str(eff_cfg.get("hmm_start_offset_fill_mode", "backfill")),
            "window_features_window": int(self.spin_window_size.value()),
            "lower_quantile_cap": lower_cap if lower_cap > 0 else None,
            "upper_quantile_cap": upper_cap if upper_cap < 1.0 else None,
            "random_state": int(self.spin_seed.value()),
            "return_details": True,
        }

    def _on_state_done(self, result):
        ct = self._cell_type()
        self._log(f"✅ State classification complete for '{ct}'.")
        # Reload the model adata from disk
        path = self._model_adata_path(ct)
        if path and path.exists():
            self._load_model_adata(path)
        self._refresh_buttons()
        self._update_view_buttons()
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
        logger = ThreadSafeLogger(self._log)

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import (
                run_hmm_state_clustering,
            )
            res = run_hmm_state_clustering(
                **params,
                verbose=True,
            )
            self._hmm_model = res["hmm_model"]
            return res["model_adata"]

        on_done_cb = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_cb = extra_callbacks.get("on_failed") if extra_callbacks else None

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
            lambda mapping: self._log(
                f"✅ Intrinsic clusters renamed: {mapping}"
            )
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
                "Run intrinsic cluster renaming first (full clusters are "
                "created when intrinsic states are combined with binary groups)."
            )
            return
        dlg = RenameClusterDialog(
            mode="full",
            model_adata=self._model_adata,
            adata_path=self._model_adata_path(ct),
            parent=self,
        )
        dlg.clusters_renamed.connect(
            lambda mapping: self._log(
                f"✅ Full behavioral clusters renamed: {mapping}"
            )
        )
        dlg.exec_()
        self._refresh_buttons()
        self._update_view_buttons()

    def _on_train(self):
        ct = self._cell_type()
        if self._model_adata is None:
            QMessageBox.warning(self, "No model", "Run state classification first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return

        self._log(f"▶ Training state classifier for '{ct}'…")
        out = self._out_dir()
        n_est = int(self.spin_clf_n_est.value())
        test_size = float(self.spin_clf_test_size.value())
        n_jobs = int(self.spin_clf_n_jobs.value())
        logger = ThreadSafeLogger(self._log)

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import (
                save_hmm_deployment_artifact,
            )
            clf_path = self._classifier_path(ct, "full")
            if not clf_path:
                raise ValueError("No classifier path.")
            
            hmm_model = getattr(self, "_hmm_model", None)
            if hmm_model is None:
                raise ValueError("HMM model not found. Re-run state clustering first.")
                
            return save_hmm_deployment_artifact(
                output_path=str(clf_path),
                model_adata=self._model_adata,
                hmm_model=hmm_model,
                output_dir=str(out) if out else "",
                cell_type=ct,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Training state classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_train, self.btn_apply],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ State classifier trained for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Train failed: {e}"),
        )

    def run_train_state(self, interactive=True, extra_callbacks=None):
        """Called from queue runner."""
        ct = self._cell_type()
        if not ct or self._model_adata is None:
            err = "No model adata available — run state classification first."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        on_done_cb = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_cb = extra_callbacks.get("on_failed") if extra_callbacks else None
        out = self._out_dir()
        n_est = int(self.spin_clf_n_est.value())
        test_size = float(self.spin_clf_test_size.value())
        n_jobs = int(self.spin_clf_n_jobs.value())
        logger = ThreadSafeLogger(self._log)

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import (
                save_hmm_deployment_artifact,
            )
            clf_path = self._classifier_path(ct, "full")
            if not clf_path:
                raise ValueError("No classifier path.")
            
            hmm_model = getattr(self, "_hmm_model", None)
            if hmm_model is None:
                raise ValueError("HMM model not found. Re-run state clustering first.")
                
            return save_hmm_deployment_artifact(
                output_path=str(clf_path),
                model_adata=self._model_adata,
                hmm_model=hmm_model,
                output_dir=str(out) if out else "",
                cell_type=ct,
                verbose=True,
            )

        def _done(r):
            self._log(f"✅ State classifier trained for '{ct}'.")
            self._update_view_buttons()
            self._notify_results()
            if on_done_cb:
                on_done_cb(r)

        def _fail(e):
            self._log(f"❌ Train failed: {e}")
            if on_fail_cb:
                on_fail_cb(e)

        self._bg.run(
            fn=_run,
            desc=f"Training state classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_train],
            viewer=self.viewer,
            inject_progress=False,
            on_done=_done,
            on_failed=_fail,
        )

    def _on_apply(self):
        ct = self._cell_type()
        if self._model_adata is None:
            QMessageBox.warning(self, "No model", "Run state classification first.")
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        # Check classifier exists
        clf_path = self._classifier_path(ct, "full")
        if not clf_path or not clf_path.exists():
            QMessageBox.warning(
                self, "No classifier",
                "Train the state classifier first before applying."
            )
            return

        self._log(f"▶ Applying state classifier to full dataset for '{ct}'…")
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import (
                apply_hmm_deployment_artifact_to_full_dataset,
                load_hmm_deployment_artifact,
            )
            artifact = load_hmm_deployment_artifact(str(clf_path))
            return apply_hmm_deployment_artifact_to_full_dataset(
                output_dir=str(out) if out else "",
                cell_type=ct,
                hmm_deployment_artifact=artifact,
                start_offset=0,
                start_offset_fill_mode="backfill",
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Applying state classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_apply],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Classifier applied to full dataset for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Apply failed: {e}"),
        )

    def run_apply_state(self, interactive=True, extra_callbacks=None):
        """Called from queue runner."""
        ct = self._cell_type()
        if not ct or self._model_adata is None:
            err = "No model adata available — run state classification first."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        clf_path = self._classifier_path(ct, "full")
        if not clf_path or not clf_path.exists():
            err = "State classifier not found — train first."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return

        on_done_cb = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_cb = extra_callbacks.get("on_failed") if extra_callbacks else None
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)

        def _run(**kw):
            from behav3d.analysis.behavior.state.classification import (
                apply_hmm_deployment_artifact_to_full_dataset,
                load_hmm_deployment_artifact,
            )
            artifact = load_hmm_deployment_artifact(str(clf_path))
            return apply_hmm_deployment_artifact_to_full_dataset(
                output_dir=str(out) if out else "",
                cell_type=ct,
                hmm_deployment_artifact=artifact,
                start_offset=0,
                start_offset_fill_mode="backfill",
                verbose=True,
            )

        def _done(r):
            self._log(f"✅ Classifier applied for '{ct}'.")
            self._update_view_buttons()
            self._notify_results()
            if on_done_cb:
                on_done_cb(r)

        def _fail(e):
            self._log(f"❌ Apply failed: {e}")
            if on_fail_cb:
                on_fail_cb(e)

        self._bg.run(
            fn=_run,
            desc=f"Applying state classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_apply],
            viewer=self.viewer,
            inject_progress=False,
            on_done=_done,
            on_failed=_fail,
        )

    def _on_state_composition(self):
        ct = self._cell_type()
        if not ct:
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Generating state composition report for '{ct}'…")

        def _run(**kw):
            from behav3d.analysis.behavior.state.visualization.plots.state_composition import (
                save_state_composition_report,
            )
            return save_state_composition_report(
                output_dir=str(out) if out else "",
                cell_type=ct,
                metadata=getattr(self.metadata_loader, "metadata", None),
                log_callback=logger,
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
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Generating state transition report for '{ct}'…")

        def _run(**kw):
            from behav3d.analysis.behavior.state.visualization.plots.state_transitions import (
                save_state_transition_report,
            )
            return save_state_transition_report(
                output_dir=str(out) if out else "",
                cell_type=ct,
                metadata=getattr(self.metadata_loader, "metadata", None),
                log_callback=logger,
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

    # ── View buttons ────────────────────────────────────────────────────

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
        if kind == "state_classifier":
            p = self._classifier_path(ct, "full")
            return [(f"Classifier ({ct})", p)] if p else []
        if kind == "state_full":
            p = self._full_adata_path(ct)
            return [(f"Full adata ({ct})", p)] if p else []
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
    """Steps 1-4 for DTW-based track trajectory classification."""

    def __init__(self, viewer=None, metadata_loader=None, cell_type_getter=None, parent=None):
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

        # ── Step 1: Clustering ───────────────────────────────────────
        grp1 = QGroupBox("Step 1 — Track Clustering")
        g1 = QVBoxLayout(grp1)
        g1.setSpacing(4)

        basic_form = QFormLayout()
        basic_form.setSpacing(3)
        self.spin_traj_size = QSpinBox()
        self.spin_traj_size.setRange(1, 9999)
        self.spin_traj_size.setValue(100)
        self.spin_traj_size.setSpecialValueText("Variable-length")
        self.spin_traj_size.setMaximumWidth(120)
        basic_form.addRow("Trajectory size:", self.spin_traj_size)

        self.spin_n_clusters = QSpinBox()
        self.spin_n_clusters.setRange(2, 200)
        self.spin_n_clusters.setValue(6)
        self.spin_n_clusters.setMaximumWidth(90)
        basic_form.addRow("N clusters:", self.spin_n_clusters)

        self.spin_seed = QSpinBox()
        self.spin_seed.setRange(0, 99999)
        self.spin_seed.setValue(123)
        self.spin_seed.setMaximumWidth(90)
        basic_form.addRow("Random seed:", self.spin_seed)

        g1.addLayout(basic_form)

        run_row = QHBoxLayout()
        self.btn_run_track = QPushButton("▶ Run Track Clustering")
        _style_primary(self.btn_run_track)
        self.btn_run_track.clicked.connect(self._on_run_track)
        run_row.addWidget(self.btn_run_track, stretch=1)
        run_row.addWidget(self.btn_queue_track_cluster)
        self.btn_view_track = _make_view_btn()
        self.btn_view_track.clicked.connect(lambda: self._on_view("track_adata"))
        run_row.addWidget(self.btn_view_track)
        g1.addLayout(run_row)

        # Advanced config
        adv1 = CollapsibleSection("⚙  Advanced Configuration (DTW)", expanded=False)
        dtw_form = QFormLayout()
        dtw_form.setSpacing(3)

        self.edit_window = QLineEdit()
        self.edit_window.setPlaceholderText("blank = unconstrained")
        self.edit_window.setMaximumWidth(130)
        dtw_form.addRow("Window:", self.edit_window)

        self.edit_max_dist = QLineEdit()
        self.edit_max_dist.setPlaceholderText("blank = off")
        self.edit_max_dist.setMaximumWidth(130)
        dtw_form.addRow("Max dist:", self.edit_max_dist)

        self.edit_penalty = QLineEdit()
        self.edit_penalty.setPlaceholderText("blank = off")
        self.edit_penalty.setMaximumWidth(130)
        dtw_form.addRow("Penalty:", self.edit_penalty)

        self.edit_psi = QLineEdit()
        self.edit_psi.setPlaceholderText("blank = off")
        self.edit_psi.setMaximumWidth(130)
        dtw_form.addRow("Psi:", self.edit_psi)

        self.combo_linkage = QComboBox()
        self.combo_linkage.addItems(["average", "complete", "single"])
        self.combo_linkage.setMaximumWidth(130)
        dtw_form.addRow("Linkage:", self.combo_linkage)

        self.combo_trim = QComboBox()
        self.combo_trim.addItems(["last", "first"])
        self.combo_trim.setMaximumWidth(130)
        dtw_form.addRow("Trim mode:", self.combo_trim)

        self.combo_missing = QComboBox()
        self.combo_missing.addItem("Keep as category", "keep")
        self.combo_missing.addItem("Drop missing timepoints", "drop")
        self.combo_missing.setMaximumWidth(200)
        dtw_form.addRow("Missing policy:", self.combo_missing)

        adv1.addLayout(dtw_form)

        self.chk_parallel = QCheckBox("Parallel DTW computation")
        self.chk_parallel.setChecked(True)
        adv1.addWidget(self.chk_parallel)

        self.chk_save_dist = QCheckBox("Save distance matrix CSV")
        self.chk_save_dist.setChecked(False)
        adv1.addWidget(self.chk_save_dist)

        # Original BEHAV3D mode
        orig_sep = QFrame()
        orig_sep.setFrameShape(QFrame.HLine)
        orig_sep.setStyleSheet("color: #444;")
        adv1.addWidget(orig_sep)

        self.chk_use_original = QCheckBox("Use original feature-based BEHAV3D DTW clustering")
        self.chk_use_original.setChecked(False)
        self.chk_use_original.toggled.connect(self._on_original_toggled)
        adv1.addWidget(self.chk_use_original)

        self._original_frame = QFrame()
        self._original_frame.setVisible(False)
        orig_form = QFormLayout(self._original_frame)
        orig_form.setSpacing(3)

        self.spin_orig_traj_size = QSpinBox()
        self.spin_orig_traj_size.setRange(1, 9999)
        self.spin_orig_traj_size.setValue(100)
        self.spin_orig_traj_size.setMaximumWidth(90)
        orig_form.addRow("Trajectory size:", self.spin_orig_traj_size)

        self.spin_orig_n_clusters = QSpinBox()
        self.spin_orig_n_clusters.setRange(2, 200)
        self.spin_orig_n_clusters.setValue(5)
        self.spin_orig_n_clusters.setMaximumWidth(90)
        orig_form.addRow("N clusters:", self.spin_orig_n_clusters)

        self.spin_orig_umap_neighbors = QSpinBox()
        self.spin_orig_umap_neighbors.setRange(2, 200)
        self.spin_orig_umap_neighbors.setValue(15)
        self.spin_orig_umap_neighbors.setMaximumWidth(90)
        orig_form.addRow("UMAP n_neighbors:", self.spin_orig_umap_neighbors)

        self.spin_orig_umap_min_dist = QDoubleSpinBox()
        self.spin_orig_umap_min_dist.setRange(0.001, 1.0)
        self.spin_orig_umap_min_dist.setSingleStep(0.05)
        self.spin_orig_umap_min_dist.setDecimals(3)
        self.spin_orig_umap_min_dist.setValue(0.1)
        self.spin_orig_umap_min_dist.setMaximumWidth(90)
        orig_form.addRow("UMAP min_dist:", self.spin_orig_umap_min_dist)

        orig_run_row = QHBoxLayout()
        self.btn_run_original = QPushButton("▶ Run Original BEHAV3D DTW")
        _style_secondary(self.btn_run_original)
        self.btn_run_original.clicked.connect(self._on_run_original)
        orig_run_row.addWidget(self.btn_run_original)
        orig_run_row.addStretch()
        self._original_frame.layout().addRow("", QWidget())  # spacer row trick

        adv1.addWidget(self._original_frame)
        adv1.addLayout(orig_run_row)

        g1.addWidget(adv1)
        lay.addWidget(grp1)

        # ── Step 2: Rename ───────────────────────────────────────────
        grp2 = QGroupBox("Step 2 — Rename Track Clusters")
        g2 = QVBoxLayout(grp2)
        g2.setSpacing(4)

        self.rename_track_status = QLabel("ℹ Run clustering first to enable renaming.")
        self.rename_track_status.setStyleSheet("color: #999; font-size: 11px;")
        self.rename_track_status.setWordWrap(True)
        g2.addWidget(self.rename_track_status)

        rename_row = QHBoxLayout()
        self.btn_rename_track = QPushButton("✏  Rename Track Clusters")
        _style_rename(self.btn_rename_track)
        self.btn_rename_track.setEnabled(False)
        self.btn_rename_track.clicked.connect(self._on_rename_track)
        rename_row.addWidget(self.btn_rename_track, stretch=1)
        self.btn_view_rename_track = _make_view_btn()
        rename_row.addWidget(self.btn_view_rename_track)
        g2.addLayout(rename_row)

        lay.addWidget(grp2)

        # ── Step 3: Train & Apply ────────────────────────────────────
        grp3 = QGroupBox("Step 3 — Train & Apply Track Classifier")
        g3 = QVBoxLayout(grp3)
        g3.setSpacing(4)

        train_row = QHBoxLayout()
        self.btn_train_track = QPushButton("▶ Train Track Classifier")
        _style_primary(self.btn_train_track)
        self.btn_train_track.setEnabled(False)
        self.btn_train_track.clicked.connect(self._on_train_track)
        train_row.addWidget(self.btn_train_track, stretch=1)
        train_row.addWidget(self.btn_queue_train_track)
        self.btn_view_train_track = _make_view_btn()
        self.btn_view_train_track.clicked.connect(lambda: self._on_view("track_classifier"))
        train_row.addWidget(self.btn_view_train_track)
        g3.addLayout(train_row)

        apply_row = QHBoxLayout()
        self.btn_apply_track = QPushButton("▶ Apply Track Classifier")
        _style_secondary(self.btn_apply_track)
        self.btn_apply_track.setEnabled(False)
        self.btn_apply_track.clicked.connect(self._on_apply_track)
        apply_row.addWidget(self.btn_apply_track, stretch=1)
        apply_row.addWidget(self.btn_queue_apply_track)
        self.btn_view_apply_track = _make_view_btn()
        self.btn_view_apply_track.clicked.connect(lambda: self._on_view("track_applied"))
        apply_row.addWidget(self.btn_view_apply_track)
        g3.addLayout(apply_row)

        adv3 = CollapsibleSection("⚙  Advanced Configuration (Track Classifier)", expanded=False)
        clf_form = QFormLayout()
        clf_form.setSpacing(3)
        self.spin_track_n_est = QSpinBox()
        self.spin_track_n_est.setRange(10, 2000)
        self.spin_track_n_est.setValue(100)
        self.spin_track_n_est.setMaximumWidth(90)
        clf_form.addRow("n_estimators:", self.spin_track_n_est)

        self.spin_track_test_size = QDoubleSpinBox()
        self.spin_track_test_size.setRange(0.05, 0.5)
        self.spin_track_test_size.setSingleStep(0.05)
        self.spin_track_test_size.setDecimals(2)
        self.spin_track_test_size.setValue(0.2)
        self.spin_track_test_size.setMaximumWidth(90)
        clf_form.addRow("Test holdout:", self.spin_track_test_size)

        adv3.addLayout(clf_form)
        g3.addWidget(adv3)
        lay.addWidget(grp3)

        # ── Step 4: Exemplar Plots ───────────────────────────────────
        grp4 = QGroupBox("Step 4 — Exemplar Plots & Diagnostics")
        g4 = QVBoxLayout(grp4)
        g4.setSpacing(4)

        ex_form = QFormLayout()
        ex_form.setSpacing(3)
        self.spin_n_per_cluster = QSpinBox()
        self.spin_n_per_cluster.setRange(1, 100)
        self.spin_n_per_cluster.setValue(10)
        self.spin_n_per_cluster.setMaximumWidth(80)
        ex_form.addRow("Exemplars / cluster:", self.spin_n_per_cluster)
        g4.addLayout(ex_form)

        opts_row = QHBoxLayout()
        self.chk_overview_statebars = QCheckBox("Overview statebars")
        self.chk_overview_statebars.setChecked(True)
        opts_row.addWidget(self.chk_overview_statebars)
        self.chk_backproj_pdf = QCheckBox("Backprojection PDFs")
        self.chk_backproj_pdf.setChecked(False)
        opts_row.addWidget(self.chk_backproj_pdf)
        self.chk_backproj_mp4 = QCheckBox("Backprojection MP4")
        self.chk_backproj_mp4.setChecked(False)
        opts_row.addWidget(self.chk_backproj_mp4)
        opts_row.addStretch()
        g4.addLayout(opts_row)

        exemplar_row = QHBoxLayout()
        self.btn_exemplars = QPushButton("▶ Create Exemplar PDFs")
        _style_secondary(self.btn_exemplars)
        self.btn_exemplars.clicked.connect(self._on_exemplars)
        exemplar_row.addWidget(self.btn_exemplars, stretch=1)
        self.btn_view_exemplars = _make_view_btn()
        self.btn_view_exemplars.clicked.connect(lambda: self._on_view("track_exemplars"))
        exemplar_row.addWidget(self.btn_view_exemplars)
        g4.addLayout(exemplar_row)

        diag_row = QHBoxLayout()
        self.btn_diagnostics = QPushButton("▶ Create Diagnostics")
        _style_secondary(self.btn_diagnostics)
        self.btn_diagnostics.clicked.connect(self._on_diagnostics)
        diag_row.addWidget(self.btn_diagnostics, stretch=1)
        self.btn_view_diagnostics = _make_view_btn()
        self.btn_view_diagnostics.clicked.connect(lambda: self._on_view("track_diagnostics"))
        diag_row.addWidget(self.btn_view_diagnostics)
        g4.addLayout(diag_row)

        lay.addWidget(grp4)

        # ── Progress + Log ───────────────────────────────────────────
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

    # ── Helpers ─────────────────────────────────────────────────────────

    def _on_original_toggled(self, checked: bool):
        self._original_frame.setVisible(checked)

    def on_metadata_updated(self):
        self._reload()

    def _reload(self):
        ct = self._get_cell_type() if self._get_cell_type else ""
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
        self._refresh_buttons()
        self._update_view_buttons()

    def _out_dir(self) -> Optional[Path]:
        od = getattr(self.metadata_loader, "output_dir", None) if self.metadata_loader else None
        return Path(str(od)).expanduser() if od else None

    def _track_adata_path(self, ct: str) -> Optional[Path]:
        out = self._out_dir()
        if not out:
            return None
        # dtaidistance output
        d = out / "analysis" / ct / "behavorial_trajectories"
        # pick first h5ad found
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
        return out / "analysis" / ct / "behavorial_trajectories" / f"BEHAV3D_{ct}_track_clusters.csv"

    def _cell_type(self) -> str:
        return self._get_cell_type() if self._get_cell_type else ""

    def _log(self, msg: str):
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        formatted = f"[{ts}] {msg}"
        self.log_edit.append(formatted)
        print(formatted)

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
                self.btn_view_track,
                self.btn_view_train_track,
                self.btn_view_apply_track,
                self.btn_view_exemplars,
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

        # Exemplar / diagnostics: look for PDFs
        out = self._out_dir()
        traj_dir = (out / "analysis" / ct / "behavorial_trajectories") if out else None
        has_exemplar = bool(
            traj_dir and any(traj_dir.glob("exemplar_tracks*.pdf"))
        )
        has_diag = bool(
            traj_dir and any(traj_dir.glob("diagnostics*.pdf"))
        )
        self.btn_view_exemplars.setEnabled(has_exemplar)
        self.btn_view_diagnostics.setEnabled(has_diag)

    # ── Click handlers ──────────────────────────────────────────────────

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
        logger = ThreadSafeLogger(self._log)
        window_raw = self.edit_window.text().strip()
        max_dist_raw = self.edit_max_dist.text().strip()
        penalty_raw = self.edit_penalty.text().strip()
        psi_raw = self.edit_psi.text().strip()

        def _parse_opt_int(s):
            try:
                return int(s) if s else None
            except ValueError:
                return None

        def _parse_opt_float(s):
            try:
                return float(s) if s else None
            except ValueError:
                return None

        params = {
            "output_dir": str(out) if out else "",
            "cell_type": ct,
            "metadata": getattr(self.metadata_loader, "metadata", None),
            "behavioral_trajectory_size": int(self.spin_traj_size.value()),
            "n_clusters": int(self.spin_n_clusters.value()),
            "random_state": int(self.spin_seed.value()),
            "window": _parse_opt_int(window_raw),
            "max_dist": _parse_opt_float(max_dist_raw),
            "penalty": _parse_opt_float(penalty_raw),
            "psi": _parse_opt_int(psi_raw),
            "linkage": self.combo_linkage.currentText(),
            "trajectory_trim_mode": self.combo_trim.currentText(),
            "missing_policy": self.combo_missing.currentData() or "keep",
            "parallel": self.chk_parallel.isChecked(),
            "save_distance_matrix": self.chk_save_dist.isChecked(),
        }

        on_done_ext = extra_callbacks.get("on_done") if extra_callbacks else None
        on_fail_ext = extra_callbacks.get("on_failed") if extra_callbacks else None

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import (
                run_categorical_dtaidistance_trajectory_clustering,
            )
            return run_categorical_dtaidistance_trajectory_clustering(
                **params,
                log_callback=logger,
                verbose=True,
            )

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
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Running original BEHAV3D feature-based DTW for '{ct}'…")

        n_clust = int(self.spin_orig_n_clusters.value())
        traj_size = int(self.spin_orig_traj_size.value())
        n_neigh = int(self.spin_orig_umap_neighbors.value())
        min_dist = float(self.spin_orig_umap_min_dist.value())

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
            buttons=[self.btn_run_original, self.btn_run_track],
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
        test_size = float(self.spin_track_test_size.value())
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Training track classifier for '{ct}'…")

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import (
                train_dtaidistance_trajectory_classifier,
            )
            return train_dtaidistance_trajectory_classifier(
                adata=self._track_adata,
                output_dir=str(out) if out else "",
                cell_type=ct,
                n_estimators=n_est,
                test_size=test_size,
                log_callback=logger,
                verbose=True,
            )

        self._bg.run(
            fn=_run,
            desc=f"Training track classifier ({ct})…",
            progress_row=self.progress_row,
            buttons=[self.btn_train_track, self.btn_apply_track],
            viewer=self.viewer,
            inject_progress=False,
            on_done=lambda r: (
                self._log(f"✅ Track classifier trained for '{ct}'."),
                self._update_view_buttons(),
                self._notify_results(),
            ),
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
        test_size = float(self.spin_track_test_size.value())
        logger = ThreadSafeLogger(self._log)

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import (
                train_dtaidistance_trajectory_classifier,
            )
            return train_dtaidistance_trajectory_classifier(
                adata=self._track_adata,
                output_dir=str(out) if out else "",
                cell_type=ct,
                n_estimators=n_est,
                test_size=test_size,
                log_callback=logger,
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
        if self._track_adata is None:
            QMessageBox.warning(self, "No data", "Run track clustering first.")
            return
        clf_path = self._track_classifier_path(ct)
        if not clf_path or not clf_path.exists():
            QMessageBox.warning(
                self, "No classifier",
                "Train the track classifier first before applying."
            )
            return
        if self._bg.is_running():
            QMessageBox.warning(self, "Busy", "Another operation is running.")
            return
        out = self._out_dir()
        logger = ThreadSafeLogger(self._log)
        self._log(f"▶ Applying track classifier for '{ct}'…")

        def _run(**kw):
            from behav3d.analysis.behavior.track.bouts import (
                apply_track_classifier_to_subtracks,
            )
            return apply_track_classifier_to_subtracks(
                adata=self._track_adata,
                output_dir=str(out) if out else "",
                cell_type=ct,
                classifier_path=str(clf_path),
                log_callback=logger,
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
                self._update_view_buttons(),
                self._notify_results(),
            ),
            on_failed=lambda e: self._log(f"❌ Apply track failed: {e}"),
        )

    def run_apply_track(self, interactive=True, extra_callbacks=None):
        """Called from queue runner."""
        ct = self._cell_type()
        if not ct or self._track_adata is None:
            err = "No track adata available."
            if extra_callbacks and extra_callbacks.get("on_failed"):
                extra_callbacks["on_failed"](err)
            return
        clf_path = self._track_classifier_path(ct)
        if not clf_path or not clf_path.exists():
            err = "Track classifier not found — train first."
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
        logger = ThreadSafeLogger(self._log)

        def _run(**kw):
            from behav3d.analysis.behavior.track.bouts import (
                apply_track_classifier_to_subtracks,
            )
            return apply_track_classifier_to_subtracks(
                adata=self._track_adata,
                output_dir=str(out) if out else "",
                cell_type=ct,
                classifier_path=str(clf_path),
                log_callback=logger,
                verbose=True,
            )

        def _done(r):
            self._log(f"✅ Track classifier applied for '{ct}'.")
            self._update_view_buttons()
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

        def _run(**kw):
            from behav3d.analysis.behavior.track.visualization.plots.exemplar_track_per_cluster import (
                save_exemplar_statebar_backprojection_pdf,
            )
            return save_exemplar_statebar_backprojection_pdf(
                adata=self._track_adata,
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

        def _run(**kw):
            from behav3d.analysis.behavior.track.state_dtw import (
                save_dtaidistance_diagnostics,
            )
            return save_dtaidistance_diagnostics(
                adata=self._track_adata,
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

    # ── View helpers ────────────────────────────────────────────────────

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
            candidates = [
                (f.stem, f) for f in sorted(traj_dir.glob("exemplar_tracks*.pdf"))
            ]
        elif kind == "track_diagnostics" and traj_dir:
            candidates = [
                (f.stem, f) for f in sorted(traj_dir.glob("diagnostics*.pdf"))
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
            menu.exec_()

    def _open_result(self, path: Path):
        if not path.suffix.lower() == ".pdf":
            # For non-PDF (h5ad, csv) open folder in OS file manager
            try:
                QDesktopServices.openUrl(
                    __import__("qtpy.QtCore", fromlist=["QUrl"]).QUrl.fromLocalFile(
                        str(path.parent)
                    )
                )
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
    """Outer single-cell tab: shared cell-type dropdown + two inner sub-tabs."""

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

    # ── UI ──────────────────────────────────────────────────────────────

    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(4)

        # ── Shared cell-type header ──────────────────────────────────
        hdr_lay = QHBoxLayout()
        hdr_lay.setSpacing(6)
        hdr_lay.addWidget(QLabel("Cell type:"))
        self.cell_type_combo = QComboBox()
        self.cell_type_combo.setMinimumWidth(160)
        self.cell_type_combo.setToolTip(
            "Immune and other cell types only (non-multicolor)."
        )
        self.cell_type_combo.currentTextChanged.connect(self._on_cell_type_changed)
        hdr_lay.addWidget(self.cell_type_combo)

        self.btn_refresh = QPushButton("🔄 Refresh")
        self.btn_refresh.setFixedHeight(26)
        self.btn_refresh.setStyleSheet(
            "QPushButton { background: #333; color: #ddd; border-radius: 3px; padding: 2px 8px; }"
            "QPushButton:hover { background: #555; }"
        )
        self.btn_refresh.clicked.connect(self._on_refresh)
        hdr_lay.addWidget(self.btn_refresh)
        hdr_lay.addStretch()

        self.status_lbl = QLabel("Load metadata to begin.")
        self.status_lbl.setStyleSheet("color: #888; font-size: 11px;")
        hdr_lay.addWidget(self.status_lbl)

        outer.addLayout(hdr_lay)

        # ── Inner sub-tabs ───────────────────────────────────────────
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

    # ── Metadata update ─────────────────────────────────────────────────

    def _on_metadata_updated(self, *_):
        cell_types = _detect_sc_cell_types(self.metadata_loader)
        current = self.cell_type_combo.currentText()
        self.cell_type_combo.blockSignals(True)
        self.cell_type_combo.clear()
        if cell_types:
            self.cell_type_combo.addItems(cell_types)
            if current in cell_types:
                self.cell_type_combo.setCurrentText(current)
            self.status_lbl.setText(
                f"Ready — {len(cell_types)} cell type(s) available."
            )
        else:
            self.status_lbl.setText("No immune/other cell types detected.")
        self.cell_type_combo.blockSignals(False)
        self._propagate_metadata_update()

    def _propagate_metadata_update(self):
        self.state_tab.on_metadata_updated()
        self.track_tab.on_metadata_updated()

    def _on_cell_type_changed(self, _text: str):
        self._propagate_metadata_update()

    def _on_refresh(self):
        self._on_metadata_updated()

    def _current_cell_type(self) -> str:
        return self.cell_type_combo.currentText()
