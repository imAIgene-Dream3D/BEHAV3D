"""
BEHAV3D napari plugin – Tracking Tab.

Provides per-cell-type sub-tabs with method selection (LAP / TrackPy / Propagation),
method-specific parameters, and batch-tracking options.
"""
import sys
import traceback
from pathlib import Path

import yaml
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout,
    QTabWidget, QGroupBox, QLabel, QPushButton,
    QComboBox, QSpinBox, QDoubleSpinBox, QCheckBox,
    QScrollArea, QTextEdit, QStackedWidget,
    QLineEdit, QFileDialog,
)
from qtpy.QtCore import Qt

from behav3d.napari._widgets import make_help_row


def _cfg_get(cfg: dict, dotted_key: str, default=None):
    cur = cfg
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


# ═══════════════════════════════════════════════════════════════════════════
# CellTypeTrackingPanel — one per cell type
# ═══════════════════════════════════════════════════════════════════════════
class CellTypeTrackingPanel(QWidget):
    """Parameter panel for tracking a single cell type."""

    def __init__(self, cell_type: str, category: str, metadata_loader,
                 all_cell_types: list, category_types: list,
                 log_callback=None, viewer=None, parent=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.category = category          # "organoid" / "immune" / "other"
        self.metadata_loader = metadata_loader
        self.all_cell_types = all_cell_types
        self.category_types = category_types   # all types in same category
        self.log = log_callback or print
        self.viewer = viewer

        # Determine defaults
        if category == "organoid":
            def_method = "propagation"
        else:
            def_method = "lap"

        # Load saved config
        tcfg = _cfg_get(self.metadata_loader.behav3d_parameters,
                        f"tracking.{cell_type}", {})
        if not isinstance(tcfg, dict):
            tcfg = {}

        self._build_ui(tcfg, def_method)

    # ------------------------------------------------------------------
    def _build_ui(self, tcfg, def_method):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(6)

        # ── Method selector ──────────────────────────────────────────
        method_group = QGroupBox("Tracking Method")
        method_layout = QHBoxLayout()
        method_layout.setContentsMargins(6, 4, 6, 4)
        self.combo_method = QComboBox()
        self.combo_method.addItems([
            "LAP (laptrack)", "TrackPy", "Propagation",
            "btrack (Bayesian)", "Import tracking",
        ])
        saved_method = tcfg.get("method", def_method)
        idx_map = {"lap": 0, "trackpy": 1, "propagation": 2, "btrack": 3, "import": 4}
        self.combo_method.setCurrentIndex(idx_map.get(saved_method, 0))
        method_layout.addWidget(QLabel("Method:"))
        method_layout.addWidget(self.combo_method)
        method_group.setLayout(method_layout)

        layout.addWidget(method_group)

        # ── Stacked method params ────────────────────────────────────
        self.param_stack = QStackedWidget()
        self.combo_method.currentIndexChanged.connect(self.param_stack.setCurrentIndex)

        # Page 0 — LAP params
        lap_cfg = tcfg.get("lap", {})
        lap_page = QWidget()
        lap_form = QFormLayout(lap_page)
        lap_form.setContentsMargins(6, 4, 6, 4)
        lap_form.setSpacing(3)
        lap_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.lap_track_cost = QSpinBox()
        self.lap_track_cost.setRange(1, 999)
        self.lap_track_cost.setValue(int(lap_cfg.get("track_cost_px", 45)))
        self.lap_track_cost.setMaximumWidth(80)
        lap_form.addRow("Track cost (px):", make_help_row(
            self.lap_track_cost,
            "Track Cost (pixels)",
            "Maximum pixel distance a cell can travel between two "
            "consecutive frames to be linked as the same track.\n\n"
            "Increase if cells move fast; decrease to avoid false links."
        ))

        self.lap_gap_cost = QSpinBox()
        self.lap_gap_cost.setRange(1, 999)
        self.lap_gap_cost.setValue(int(lap_cfg.get("gap_close_cost_px", 60)))
        self.lap_gap_cost.setMaximumWidth(80)
        lap_form.addRow("Gap close cost (px):", make_help_row(
            self.lap_gap_cost,
            "Gap Closing Cost (pixels)",
            "Maximum distance (in pixels) allowed when reconnecting a "
            "track that was temporarily lost for one or more frames.\n\n"
            "Should be >= Track cost. Increase if cells disappear "
            "briefly due to segmentation gaps."
        ))

        self.lap_gap_frames = QSpinBox()
        self.lap_gap_frames.setRange(0, 100)
        self.lap_gap_frames.setValue(int(lap_cfg.get("gap_close_max_frames", 3)))
        self.lap_gap_frames.setMaximumWidth(60)
        lap_form.addRow("Gap close max frames:", make_help_row(
            self.lap_gap_frames,
            "Gap Closing Max Frames",
            "Maximum number of consecutive frames a cell can be missing "
            "before the gap is too large to close.\n\n"
            "Higher values recover longer disappearances but may "
            "introduce false links."
        ))

        self.lap_merge_cost = QSpinBox()
        self.lap_merge_cost.setRange(0, 999)
        self.lap_merge_cost.setValue(int(lap_cfg.get("merging_cost_px", 0)))
        self.lap_merge_cost.setMaximumWidth(80)
        lap_form.addRow("Merging cost (px):", make_help_row(
            self.lap_merge_cost,
            "Merging Cost (pixels)",
            "Maximum distance for detecting merge events, where two "
            "tracks combine into one object.\n\n"
            "Set to 0 to disable merging detection.\n"
            "Useful when cells fuse or cluster together."
        ))

        self.lap_split_cost = QSpinBox()
        self.lap_split_cost.setRange(0, 999)
        self.lap_split_cost.setValue(int(lap_cfg.get("splitting_cost_px", 0)))
        self.lap_split_cost.setMaximumWidth(80)
        lap_form.addRow("Splitting cost (px):", make_help_row(
            self.lap_split_cost,
            "Splitting Cost (pixels)",
            "Maximum distance for detecting split events, where one "
            "object divides into two tracks.\n\n"
            "Set to 0 to disable splitting detection.\n"
            "Useful for cell division or organoid fragmentation."
        ))

        self.param_stack.addWidget(lap_page)

        # Page 1 — TrackPy params
        tp_cfg = tcfg.get("trackpy", {})
        tp_page = QWidget()
        tp_form = QFormLayout(tp_page)
        tp_form.setContentsMargins(6, 4, 6, 4)
        tp_form.setSpacing(3)
        tp_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.tp_search_range = QSpinBox()
        self.tp_search_range.setRange(1, 999)
        self.tp_search_range.setValue(int(tp_cfg.get("search_range_px", 31)))
        self.tp_search_range.setMaximumWidth(80)
        tp_form.addRow("Search range (px):", make_help_row(
            self.tp_search_range,
            "Search Range (pixels)",
            "Maximum pixel distance to look for a cell in the next frame.\n\n"
            "Should be large enough to cover the fastest-moving cells."
        ))

        self.tp_memory = QSpinBox()
        self.tp_memory.setRange(0, 100)
        self.tp_memory.setValue(int(tp_cfg.get("memory_frames", 2)))
        self.tp_memory.setMaximumWidth(60)
        tp_form.addRow("Memory (frames):", make_help_row(
            self.tp_memory,
            "Memory (frames)",
            "Number of frames a cell can disappear and still be "
            "reconnected to its previous track.\n\n"
            "Similar to 'gap closing' in LAP."
        ))

        self.tp_adaptive_stop = QDoubleSpinBox()
        self.tp_adaptive_stop.setRange(0.0, 100.0)
        self.tp_adaptive_stop.setSingleStep(0.5)
        self.tp_adaptive_stop.setValue(float(tp_cfg.get("adaptive_stop", 10.0)))
        self.tp_adaptive_stop.setMaximumWidth(80)
        tp_form.addRow("Adaptive stop:", make_help_row(
            self.tp_adaptive_stop,
            "Adaptive Stop",
            "Factor that limits how much the search range can shrink "
            "adaptively.\n\n"
            "Higher = more conservative shrinking.\n"
            "Leave at default unless tracking quality is poor."
        ))

        self.tp_adaptive_step = QDoubleSpinBox()
        self.tp_adaptive_step.setRange(0.01, 1.0)
        self.tp_adaptive_step.setSingleStep(0.01)
        self.tp_adaptive_step.setDecimals(3)
        self.tp_adaptive_step.setValue(float(tp_cfg.get("adaptive_step", 0.95)))
        self.tp_adaptive_step.setMaximumWidth(80)
        tp_form.addRow("Adaptive step:", make_help_row(
            self.tp_adaptive_step,
            "Adaptive Step",
            "Multiplier applied to the search range at each iteration "
            "of the adaptive search.\n\n"
            "Values close to 1.0 = slow reduction.\n"
            "Values close to 0.5 = aggressive reduction.\n\n"
            "Default 0.95 works well in most cases."
        ))

        self.param_stack.addWidget(tp_page)

        # Page 2 — Propagation params (Simplified)
        prop_page = QWidget()
        prop_lay = QVBoxLayout(prop_page)
        prop_notice = QLabel("No tunable parameters for\nPropagation tracking method.")
        prop_notice.setWordWrap(True)
        prop_notice.setAlignment(Qt.AlignCenter)
        prop_notice.setStyleSheet("color: #666; font-style: italic; padding: 10px;")
        prop_lay.addWidget(prop_notice)
        prop_lay.addStretch()

        self.param_stack.addWidget(prop_page)

        # Page 3 — btrack (Bayesian tracking)
        bt_cfg = tcfg.get("btrack", {})
        btrack_page = QWidget()
        btrack_lay = QVBoxLayout(btrack_page)
        btrack_lay.setContentsMargins(6, 4, 6, 4)
        btrack_lay.setSpacing(4)

        # Info banner
        _bt_info = QLabel(
            "<b>btrack</b> \u2014 Bayesian multi-object tracker.<br>"
            "<b>Step 1</b> (Kalman filter): always runs \u2014 tune and validate first.<br>"
            "<b>Step 2</b> (optimizer): opt-in \u2014 enable only once Step 1 looks correct.<br>"
            "See <tt>models/README_btrack.md</tt> for full parameter reference."
        )
        _bt_info.setWordWrap(True)
        _bt_info.setStyleSheet("color: #888; font-size: 10px; padding: 2px 0 6px 0;")
        btrack_lay.addWidget(_bt_info)

        # ── Sub-group A: Core Tracking (Step 1) ─────────────────
        step1_group = QGroupBox("Step 1 \u2014 Kalman Filter Tracking")
        step1_form = QFormLayout(step1_group)
        step1_form.setContentsMargins(6, 4, 6, 4)
        step1_form.setSpacing(3)
        step1_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.bt_config_preset = QComboBox()
        self.bt_config_preset.addItems([
            "Cell (bundled)", "Particle (bundled)", "Custom JSON…",
        ])
        preset_val = bt_cfg.get("config_preset", "cell")
        preset_idx = {"cell": 0, "particle": 1}.get(str(preset_val).lower(), 2)
        self.bt_config_preset.setCurrentIndex(preset_idx)
        step1_form.addRow("Config preset:", make_help_row(
            self.bt_config_preset,
            "Configuration Preset",
            "Select a bundled motion/hypothesis model, or provide\n"
            "your own JSON configuration file.\n\n"
            "• Cell — tuned for biological cell tracking\n"
            "• Particle — tuned for small particle tracking\n"
            "• Custom JSON — browse to your own config file"
        ))

        config_path_row = QHBoxLayout()
        self.bt_config_path = QLineEdit()
        self.bt_config_path.setPlaceholderText("Path to custom btrack JSON…")
        self.bt_config_path.setText(bt_cfg.get("config_path", ""))
        self.bt_config_path.setEnabled(preset_idx == 2)
        self.bt_browse_btn = QPushButton("Browse")
        self.bt_browse_btn.setFixedWidth(60)
        self.bt_browse_btn.setEnabled(preset_idx == 2)
        self.bt_browse_btn.clicked.connect(self._bt_browse_config)
        config_path_row.addWidget(self.bt_config_path)
        config_path_row.addWidget(self.bt_browse_btn)
        config_path_widget = QWidget()
        config_path_widget.setLayout(config_path_row)
        step1_form.addRow("", config_path_widget)

        def _on_preset_changed(idx):
            is_custom = idx == 2
            self.bt_config_path.setEnabled(is_custom)
            self.bt_browse_btn.setEnabled(is_custom)
        self.bt_config_preset.currentIndexChanged.connect(_on_preset_changed)

        self.bt_max_search_radius = QSpinBox()
        self.bt_max_search_radius.setRange(1, 9999)
        self.bt_max_search_radius.setValue(int(bt_cfg.get("max_search_radius", 100)))
        self.bt_max_search_radius.setMaximumWidth(80)
        step1_form.addRow("Max search radius (px):", make_help_row(
            self.bt_max_search_radius,
            "Max Search Radius (pixels)",
            "Maximum isotropic distance (pixels) to search for\n"
            "linking objects between frames.\n\n"
            "Increase for fast-moving cells; decrease to\n"
            "prevent long-range false links."
        ))

        self.bt_update_method = QComboBox()
        self.bt_update_method.addItems(["EXACT (recommended)", "APPROXIMATE (large datasets)"])
        self.bt_update_method.setCurrentIndex(
            1 if bt_cfg.get("update_method", "EXACT").upper() == "APPROXIMATE" else 0
        )
        step1_form.addRow("Update method:", make_help_row(
            self.bt_update_method,
            "Update Method",
            "EXACT — full Bayesian belief matrix (accurate, slower).\n"
            "APPROXIMATE — local spatial search (faster, for >1000\n"
            "cells per frame).\n\n"
            "Use EXACT unless tracking is too slow."
        ))

        self.bt_step_size = QSpinBox()
        self.bt_step_size.setRange(1, 10000)
        self.bt_step_size.setValue(int(bt_cfg.get("step_size", 100)))
        self.bt_step_size.setMaximumWidth(80)
        step1_form.addRow("Step size (frames):", make_help_row(
            self.bt_step_size,
            "Step Size",
            "Number of frames processed per tracking iteration.\n\n"
            "Lower values use less RAM but may be slightly slower."
        ))

        btrack_lay.addWidget(step1_group)

        # ── Sub-group B: Global Optimizer (Step 2 — opt-in) ─────
        step2_group = QGroupBox("Step 2 — Global Hypothesis Optimizer")
        step2_lay = QVBoxLayout(step2_group)
        step2_lay.setContentsMargins(6, 4, 6, 4)
        step2_lay.setSpacing(3)

        self.bt_use_optimize = QCheckBox("Enable global track optimization")
        self.bt_use_optimize.setChecked(bt_cfg.get("use_optimize", False))
        self.bt_use_optimize.setToolTip(
            "Step 2 runs a global integer-programming optimizer\n"
            "to resolve track merges, splits, and false positives.\n\n"
            "Recommended: first tune Step 1 without optimization,\n"
            "then enable this for refinement."
        )
        step2_lay.addWidget(self.bt_use_optimize)

        # Hypotheses checkboxes
        hyp_group = QGroupBox("Hypotheses")
        hyp_lay = QVBoxLayout(hyp_group)
        hyp_lay.setContentsMargins(4, 2, 4, 2)
        hyp_lay.setSpacing(1)
        saved_hyps = bt_cfg.get("hypotheses",
                                ["P_FP", "P_init", "P_term", "P_link"])
        self.bt_hyp_checks = {}
        for hyp_name, hyp_desc, default_on in [
            ("P_FP",     "False positive",    True),
            ("P_init",   "Track initialization", True),
            ("P_term",   "Track termination",  True),
            ("P_link",   "Track linking",      True),
            ("P_branch", "Track branching",    False),
            ("P_dead",   "Cell death",         False),
            ("P_merge",  "Track merging",      False),
        ]:
            cb = QCheckBox(f"{hyp_name} — {hyp_desc}")
            is_on = hyp_name in saved_hyps if saved_hyps else default_on
            cb.setChecked(is_on)
            if hyp_name == "P_FP":
                cb.setEnabled(False)  # P_FP is always required
            hyp_lay.addWidget(cb)
            self.bt_hyp_checks[hyp_name] = cb
        step2_lay.addWidget(hyp_group)

        step2_form = QFormLayout()
        step2_form.setContentsMargins(0, 0, 0, 0)
        step2_form.setSpacing(3)
        step2_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.bt_dist_thresh = QSpinBox()
        self.bt_dist_thresh.setRange(1, 9999)
        self.bt_dist_thresh.setValue(int(bt_cfg.get("dist_thresh", 60)))
        self.bt_dist_thresh.setMaximumWidth(80)
        step2_form.addRow("Distance threshold:", make_help_row(
            self.bt_dist_thresh,
            "Distance Threshold",
            "Maximum distance (pixels) for generating\n"
            "link/branch hypotheses in the optimizer."
        ))

        self.bt_time_thresh = QSpinBox()
        self.bt_time_thresh.setRange(1, 999)
        self.bt_time_thresh.setValue(int(bt_cfg.get("time_thresh", 3)))
        self.bt_time_thresh.setMaximumWidth(60)
        step2_form.addRow("Time threshold:", make_help_row(
            self.bt_time_thresh,
            "Time Threshold",
            "Maximum frame gap for generating link\n"
            "hypotheses in the optimizer."
        ))
        step2_lay.addLayout(step2_form)

        btrack_lay.addWidget(step2_group)

        # Toggle Step 2 widgets enabled state
        self._bt_step2_widgets = [
            hyp_group, self.bt_dist_thresh, self.bt_time_thresh,
        ]
        def _on_optimize_toggled(checked):
            for w in self._bt_step2_widgets:
                w.setEnabled(checked)
        self.bt_use_optimize.toggled.connect(_on_optimize_toggled)
        _on_optimize_toggled(self.bt_use_optimize.isChecked())

        btrack_lay.addStretch()
        self.param_stack.addWidget(btrack_page)

        # Page 4 — Import tracking (Coming soon)
        import_page = QWidget()
        import_lay = QVBoxLayout(import_page)
        import_notice = QLabel(
            "Import tracking will be\navailable in a future release."
        )
        import_notice.setWordWrap(True)
        import_notice.setAlignment(Qt.AlignCenter)
        import_notice.setStyleSheet("color: #999; font-style: italic; padding: 10px;")
        import_lay.addWidget(import_notice)
        import_lay.addStretch()
        self.param_stack.addWidget(import_page)

        # Set active page
        self.param_stack.setCurrentIndex(self.combo_method.currentIndex())
        layout.addWidget(self.param_stack)

        # ── Apply settings buttons ──────────────────────────────────
        cat_label = self.category.capitalize() + "s" if self.category != "other" else "Other types"
        sync_row = QHBoxLayout()
        self.btn_apply_cat = QPushButton(f"Apply to all {cat_label}")
        self.btn_apply_cat.clicked.connect(lambda: self._apply_to_others(category_only=True))
        self.btn_apply_all = QPushButton("Apply to all")
        self.btn_apply_all.clicked.connect(lambda: self._apply_to_others(category_only=False))
        
        sync_row.addWidget(self.btn_apply_cat)
        sync_row.addWidget(self.btn_apply_all)
        layout.addLayout(sync_row)

        # ── Run button ───────────────────────────────────────────────
        self.btn_run = QPushButton(f"Run {self.cell_type.capitalize()} Tracking")
        self.btn_run.setStyleSheet(
            "background-color: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px;"
        )
        self.btn_run.clicked.connect(self._on_run_clicked)
        layout.addWidget(self.btn_run)

        # Disable run button for coming-soon methods (only Import at index 4)
        def _on_method_idx_changed(idx):
            is_coming_soon = idx >= 4
            self.btn_run.setEnabled(not is_coming_soon)
            if is_coming_soon:
                self.btn_run.setToolTip("This tracking method is not yet available.")
            else:
                self.btn_run.setToolTip("")
        self.combo_method.currentIndexChanged.connect(_on_method_idx_changed)
        _on_method_idx_changed(self.combo_method.currentIndex())

        layout.addStretch()

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------
    def _get_method_key(self):
        return ["lap", "trackpy", "propagation", "btrack", "import"][
            self.combo_method.currentIndex()
        ]

    def _bt_browse_config(self):
        """Open file dialog for custom btrack JSON config."""
        path, _ = QFileDialog.getOpenFileName(
            self, "Select btrack config JSON", "",
            "JSON files (*.json);;All files (*)",
        )
        if path:
            self.bt_config_path.setText(path)

    def _bt_get_config_preset(self):
        """Return the resolved config preset value."""
        idx = self.bt_config_preset.currentIndex()
        if idx == 0:
            return "cell"
        elif idx == 1:
            return "particle"
        else:
            return self.bt_config_path.text().strip()

    def _bt_get_hypotheses(self):
        """Return list of checked hypothesis names."""
        return [name for name, cb in self.bt_hyp_checks.items() if cb.isChecked()]

    def _collect_params(self) -> dict:
        """Collect current widget values into a dict."""
        return {
            "method": self._get_method_key(),
            "lap": {
                "track_cost_px": int(self.lap_track_cost.value()),
                "gap_close_cost_px": int(self.lap_gap_cost.value()),
                "gap_close_max_frames": int(self.lap_gap_frames.value()),
                "merging_cost_px": int(self.lap_merge_cost.value()),
                "splitting_cost_px": int(self.lap_split_cost.value()),
            },
            "trackpy": {
                "search_range_px": int(self.tp_search_range.value()),
                "memory_frames": int(self.tp_memory.value()),
                "adaptive_stop": float(self.tp_adaptive_stop.value()),
                "adaptive_step": float(self.tp_adaptive_step.value()),
            },
            "propagation": {
                # Notice: no tunable params currently exposed
            },
            "btrack": {
                "config_preset": self._bt_get_config_preset(),
                "config_path": self.bt_config_path.text().strip(),
                "max_search_radius": int(self.bt_max_search_radius.value()),
                "update_method": "APPROXIMATE" if self.bt_update_method.currentIndex() == 1 else "EXACT",
                "step_size": int(self.bt_step_size.value()),
                "use_optimize": self.bt_use_optimize.isChecked(),
                "hypotheses": self._bt_get_hypotheses(),
                "dist_thresh": int(self.bt_dist_thresh.value()),
                "time_thresh": int(self.bt_time_thresh.value()),
            },
        }

    def _apply_to_others(self, category_only=False):
        """Copy current settings to other panels."""
        parent_tab = self.parent()
        while parent_tab and not hasattr(parent_tab, 'panels'):
            parent_tab = parent_tab.parent()
        
        if not parent_tab:
            return

        settings = self._collect_params()
        targets = self.category_types if category_only else self.all_cell_types
        
        count = 0
        for ct in targets:
            if ct == self.cell_type:
                continue
            if ct in parent_tab.panels:
                panel = parent_tab.panels[ct]
                # Apply settings
                idx_map = {"lap": 0, "trackpy": 1, "propagation": 2, "btrack": 3}
                panel.combo_method.setCurrentIndex(idx_map.get(settings["method"], 0))
                
                # LAP
                panel.lap_track_cost.setValue(settings["lap"]["track_cost_px"])
                panel.lap_gap_cost.setValue(settings["lap"]["gap_close_cost_px"])
                panel.lap_gap_frames.setValue(settings["lap"]["gap_close_max_frames"])
                panel.lap_merge_cost.setValue(settings["lap"]["merging_cost_px"])
                panel.lap_split_cost.setValue(settings["lap"]["splitting_cost_px"])
                
                # TrackPy
                panel.tp_search_range.setValue(settings["trackpy"]["search_range_px"])
                panel.tp_memory.setValue(settings["trackpy"]["memory_frames"])
                panel.tp_adaptive_stop.setValue(settings["trackpy"]["adaptive_stop"])
                panel.tp_adaptive_step.setValue(settings["trackpy"]["adaptive_step"])
                
                # btrack
                bt = settings.get("btrack", {})
                preset = bt.get("config_preset", "cell")
                preset_idx = {"cell": 0, "particle": 1}.get(
                    str(preset).lower(), 2
                )
                panel.bt_config_preset.setCurrentIndex(preset_idx)
                panel.bt_config_path.setText(bt.get("config_path", ""))
                panel.bt_max_search_radius.setValue(bt.get("max_search_radius", 100))
                panel.bt_update_method.setCurrentIndex(
                    1 if bt.get("update_method", "EXACT").upper() == "APPROXIMATE" else 0
                )
                panel.bt_step_size.setValue(bt.get("step_size", 100))
                panel.bt_use_optimize.setChecked(bt.get("use_optimize", False))
                for hyp_name, cb in panel.bt_hyp_checks.items():
                    if hyp_name == "P_FP":
                        continue
                    cb.setChecked(hyp_name in bt.get("hypotheses", []))
                panel.bt_dist_thresh.setValue(bt.get("dist_thresh", 60))
                panel.bt_time_thresh.setValue(bt.get("time_thresh", 3))
                
                count += 1
        
        scope = "category" if category_only else "all"
        self.log(f"Applied settings to {count} other cell types ({scope}).")

    def _persist(self):
        """Save this cell type's tracking params to YAML."""
        params = self.metadata_loader.behav3d_parameters
        tracking = params.setdefault("tracking", {})
        tracking[self.cell_type] = self._collect_params()

        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception as e:
                self.log(f"Warning: Could not save parameters: {e}")

    # ------------------------------------------------------------------
    # Running
    # ------------------------------------------------------------------
    def _determine_targets(self) -> list:
        """Which cell types to track based on batch checkboxes."""
        if self.check_batch_all.isChecked():
            return list(self.all_cell_types)
        if self.check_batch_category.isChecked():
            return list(self.category_types)
        return [self.cell_type]

    def _run_tracking_for(self, cell_type: str, overwrite: bool = False):
        """Run tracking for a single cell type using current panel settings."""
        method = self._get_method_key()
        metadata = self.metadata_loader.metadata
        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        print(f"\n{'='*50}", file=sys.stderr)
        print(f"  Tracking: {cell_type} ({method.upper()})", file=sys.stderr)
        print(f"{'='*50}", file=sys.stderr)

        if method == "lap":
            from behav3d.preprocessing.tracking.laptracking import run_laptracking
            tc = int(self.lap_track_cost.value()) ** 2
            gc = int(self.lap_gap_cost.value()) ** 2
            mc = int(self.lap_merge_cost.value())
            sc = int(self.lap_split_cost.value())
            new_md = run_laptracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                track_cost_cutoff=tc,
                gap_closing_cost_cutoff=gc,
                gap_closing_max_frame_count=int(self.lap_gap_frames.value()),
                merging_cost_cutoff=(mc ** 2 if mc > 0 else False),
                splitting_cost_cutoff=(sc ** 2 if sc > 0 else False),
                overwrite=overwrite,
            )

        elif method == "trackpy":
            from behav3d.preprocessing.tracking.trackpy_tracking import run_trackpy_tracking_generic
            new_md = run_trackpy_tracking_generic(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                search_range=int(self.tp_search_range.value()),
                memory=int(self.tp_memory.value()),
                adaptive_stop=float(self.tp_adaptive_stop.value()),
                adaptive_step=float(self.tp_adaptive_step.value()),
                log_callback=self.log,
            )

        elif method == "btrack":
            from behav3d.preprocessing.tracking.btrack_tracking import run_btracking
            config_preset = self._bt_get_config_preset()
            update_method = "APPROXIMATE" if self.bt_update_method.currentIndex() == 1 else "EXACT"
            use_opt = self.bt_use_optimize.isChecked()
            new_md = run_btracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite,
                config_preset=config_preset,
                max_search_radius=int(self.bt_max_search_radius.value()),
                update_method=update_method,
                step_size=int(self.bt_step_size.value()),
                use_optimize=use_opt,
                hypotheses=self._bt_get_hypotheses() if use_opt else None,
                dist_thresh=int(self.bt_dist_thresh.value()) if use_opt else None,
                time_thresh=int(self.bt_time_thresh.value()) if use_opt else None,
                log_callback=self.log,
            )

        else:  # propagation
            from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking
            new_md = run_propagation_tracking(
                metadata=metadata, output_dir=out_dir, cell_type=cell_type,
                overwrite=overwrite
            )

        print(f"\u2705 {cell_type} tracking finished.", file=sys.stderr)

        # Update shared metadata
        self.metadata_loader.metadata = new_md
        # Save updated metadata CSV
        csv_path = self.metadata_loader.behav3d_parameters.get("paths", {}).get("metadata_csv")
        if csv_path:
            new_md.to_csv(csv_path, sep=",", index=False)

    def _check_existing_tracking(self, cell_types: list) -> list:
        """Return descriptions of existing tracking data that would be overwritten."""
        warnings = []
        md = self.metadata_loader.metadata
        if md is None:
            return warnings
        out_dir = Path(self.metadata_loader.output_dir)
        for ct in cell_types:
            for _, sample in md.iterrows():
                sn = sample.get("sample_name", "unknown")
                zarr_path = out_dir / "images" / sn / f"{sn}_{ct}_tracked.zarr"
                csv_dir = out_dir / "trackdata" / sn / ct
                if zarr_path.exists() or (csv_dir.exists() and list(csv_dir.glob("*.csv"))):
                    warnings.append(f"{ct} tracking data for {sn}")
        return warnings

    def _on_run_clicked(self):
        self._persist()
        method = self._get_method_key()
        self.log(f"Running {method.upper()} tracking for: {self.cell_type}")

        # Check for existing data across all samples
        existing = self._check_existing_tracking([self.cell_type])

        overwrite = True
        if existing:
            from qtpy.QtWidgets import QMessageBox
            details = "\n".join(f"  \u2022 {w}" for w in existing)
            box = QMessageBox(self)
            box.setWindowTitle("Overwrite Existing Tracking?")
            box.setText(
                f"The following tracking data already exists:\n\n{details}\n\n"
                "What do you want to do?"
            )
            btn_overwrite = box.addButton("Overwrite", QMessageBox.DestructiveRole)
            btn_skip = box.addButton("Skip", QMessageBox.AcceptRole)
            btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(btn_cancel)
            box.exec_()
            clicked = box.clickedButton()
            if clicked != btn_overwrite:
                self.log(f"Tracking for {self.cell_type} cancelled.")
                return
            overwrite = True

        try:
            self._run_tracking_for(self.cell_type, overwrite=overwrite)
            self.log(f"\u2705 {self.cell_type} tracking finished.")
            
            # Show visualizer prompt
            from qtpy.QtWidgets import QMessageBox
            res = QMessageBox.question(
                self, "Tracking Finished",
                f"Tracking for {self.cell_type} finished! \n\nDo you want to switch to the Visualization Tab and see the tracks?",
                QMessageBox.Yes | QMessageBox.No
            )
            if res == QMessageBox.Yes:
                self._switch_to_viz_and_show_tracks()

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during tracking: {e}")

    def _switch_to_viz_and_show_tracks(self):
        """Utility to switch to visualization tab and enable track layers."""
        # Stop any running napari animation to avoid stale slider references
        # (prevents RuntimeError: wrapped C/C++ object has been deleted)
        try:
            qt_dims = self.viewer.window._qt_viewer.dims
            if hasattr(qt_dims, '_animation_thread') and qt_dims._animation_thread is not None:
                qt_dims._animation_thread.quit()
                qt_dims._animation_thread.wait()
        except Exception:
            pass

        parent = self.parent()
        while parent and not hasattr(parent, 'tabs'):
            parent = parent.parent()
        
        if parent and hasattr(parent, 'tabs'):
            parent.tabs.setCurrentIndex(1) # Visualization Tab
            if hasattr(parent, 'visualization_tab'):
                parent.visualization_tab.sample_combo.setCurrentIndex(0)
                parent.visualization_tab._on_load_dataset()
                
                # Make 'Tracks' layers visible
                for layer in self.viewer.layers:
                    if "Tracks" in layer.name:
                        layer.visible = True


# ═══════════════════════════════════════════════════════════════════════════
# TrackingTab — main tab with per-cell-type sub-tabs
# ═══════════════════════════════════════════════════════════════════════════
class TrackingTab(QWidget):
    def __init__(self, viewer, metadata_loader, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.panels: dict[str, CellTypeTrackingPanel] = {}

        self._init_ui()

        # Listen to metadata changes
        if hasattr(self.metadata_loader, "metadata_loaded"):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

    # ------------------------------------------------------------------
    def _init_ui(self):
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        outer.addWidget(scroll)
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)
        scroll.setWidget(content)

        # Sub-tab widget (West position = left tabs)
        self.cell_tabs = QTabWidget()
        self.cell_tabs.setTabPosition(QTabWidget.West)
        layout.addWidget(self.cell_tabs)

        # Global Run Button + Queue button
        self.btn_run_batch = QPushButton("Run Batch Tracking (All Cell Types)")
        self.btn_run_batch.setStyleSheet(
            "background-color: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 10px; font-size: 14px;"
        )
        self.btn_run_batch.clicked.connect(self._on_run_batch_clicked)

        self.btn_queue_track = QPushButton("+🛒")
        self.btn_queue_track.setFixedSize(36, 32)
        self.btn_queue_track.setToolTip("Add Batch Tracking to Processing Queue")
        self.btn_queue_track.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )

        batch_btn_row = QHBoxLayout()
        batch_btn_row.setSpacing(4)
        batch_btn_row.addWidget(self.btn_run_batch, stretch=1)
        batch_btn_row.addWidget(self.btn_queue_track)
        layout.addLayout(batch_btn_row)

        # Hidden until metadata is loaded
        self.btn_run_batch.setVisible(False)
        self.btn_queue_track.setVisible(False)

        # Log window
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(120)
        self.log_box.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log_box)

        # Placeholder when no metadata
        self._placeholder = QLabel("Load metadata in the Data Preparation tab to see tracking options.")
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._placeholder.setWordWrap(True)
        self._placeholder.setStyleSheet("color: #888; font-style: italic;")
        self.cell_tabs.addTab(self._placeholder, "—")

    def _log(self, msg):
        import datetime
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_box.append(f"[{ts}] {msg}")
        self.log_box.verticalScrollBar().setValue(
            self.log_box.verticalScrollBar().maximum()
        )

    # ------------------------------------------------------------------
    def _on_metadata_updated(self):
        self._log("Metadata updated — refreshing tracking tabs…")
        self._rebuild_tabs()

    def _detect_cell_types(self):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
        )
        md = self.metadata_loader.metadata
        if md is None:
            return [], [], []
        return (
            detect_organoid_types_from_metadata(md),
            detect_immune_cell_types_from_metadata(md),
            detect_other_cell_types_from_metadata(md),
        )

    def _rebuild_tabs(self):
        # Clear old tabs
        self.cell_tabs.clear()
        self.panels.clear()

        org, imm, oth = self._detect_cell_types()
        all_types = org + imm + oth

        if not all_types:
            self.cell_tabs.addTab(self._placeholder, "—")
            self.btn_run_batch.setVisible(False)
            self.btn_queue_track.setVisible(False)
            return

        self.btn_run_batch.setVisible(True)
        self.btn_queue_track.setVisible(True)

        for ct in org:
            panel = CellTypeTrackingPanel(
                cell_type=ct, category="organoid",
                metadata_loader=self.metadata_loader,
                all_cell_types=all_types, category_types=org,
                log_callback=self._log,
                viewer=self.viewer
            )
            self.panels[ct] = panel
            self.cell_tabs.addTab(panel, f"🟣 {ct.capitalize()}")

        for ct in imm:
            panel = CellTypeTrackingPanel(
                cell_type=ct, category="immune",
                metadata_loader=self.metadata_loader,
                all_cell_types=all_types, category_types=imm,
                log_callback=self._log,
                viewer=self.viewer
            )
            self.panels[ct] = panel
            self.cell_tabs.addTab(panel, f"🔵 {ct.capitalize()}")

        for ct in oth:
            panel = CellTypeTrackingPanel(
                cell_type=ct, category="other",
                metadata_loader=self.metadata_loader,
                all_cell_types=all_types, category_types=oth,
                log_callback=self._log,
                viewer=self.viewer
            )
            self.panels[ct] = panel
            self.cell_tabs.addTab(panel, f"🟡 {ct.capitalize()}")

    def _on_run_batch_clicked(self):
        self.run_batch_tracking(interactive=True)

    def run_batch_tracking(self, interactive=True, skip_existing=False):
        """Sequential run for all configured cell type panels.
        When interactive=False, skips overwrite and visualization dialogs."""
        if not self.panels:
            self._log("No cell type panels to track.")
            return

        total = len(self.panels)
        self._log(f"Starting batch tracking for {total} cell types...")

        print(f"\n{'='*60}", file=sys.stderr)
        print(f"  Running batch tracking for {total} cell types", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

        # Check for existing tracking data across all cell types
        all_cts = list(self.panels.keys())
        existing = []
        existing_cts = set()
        out_dir = Path(self.metadata_loader.output_dir)
        md = self.metadata_loader.metadata
        for ct in all_cts:
            for _, sample in md.iterrows():
                sn = sample.get("sample_name", "unknown")
                zarr_path = out_dir / "images" / sn / f"{sn}_{ct}_tracked.zarr"
                csv_dir = out_dir / "trackdata" / sn / ct
                if zarr_path.exists() or (csv_dir.exists() and list(csv_dir.glob("*.csv"))):
                    existing.append(f"{ct} tracking data for {sn}")
                    existing_cts.add(ct)

        skip_existing_flag = skip_existing
        overwrite = not skip_existing
        if existing and interactive:
            from qtpy.QtWidgets import QMessageBox
            details = "\n".join(f"  \u2022 {w}" for w in existing)
            box = QMessageBox(self)
            box.setWindowTitle("Overwrite Existing Tracking?")
            box.setText(
                f"The following tracking data already exists:\n\n{details}\n\n"
                "What do you want to do?"
            )
            btn_overwrite = box.addButton("Overwrite All", QMessageBox.DestructiveRole)
            btn_skip = box.addButton("Skip Existing", QMessageBox.AcceptRole)
            btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
            box.setDefaultButton(btn_cancel)
            box.exec_()
            clicked = box.clickedButton()
            if clicked == btn_cancel:
                self._log("Batch tracking cancelled by user.")
                return
            skip_existing_flag = (clicked == btn_skip)
            overwrite = not skip_existing_flag

        # Persist all panel params before starting so config is saved even if an error occurs
        for ct, panel in self.panels.items():
            panel._persist()

        try:
            for i, (ct, panel) in enumerate(self.panels.items(), 1):
                if skip_existing_flag and ct in existing_cts:
                    self._log(f"--- [{i}/{total}] Skipping {ct} (existing data) ---")
                    continue
                print(f"\n\u25b6 [{i}/{total}] Tracking {ct}...", file=sys.stderr)
                self._log(f"--- [{i}/{total}] Tracking {ct} ---")
                panel._run_tracking_for(ct, overwrite=overwrite)
                self._log(f"Done: {ct}")
            
            self._log("\u2705 Batch tracking finished.")
            print(f"\n{'='*60}", file=sys.stderr)
            print(f"  \u2705 Batch tracking complete", file=sys.stderr)
            print(f"{'='*60}\n", file=sys.stderr)

            if interactive:
                from qtpy.QtWidgets import QMessageBox
                res = QMessageBox.question(
                    self, "Batch Tracking Finished",
                    "Successfully tracked all cell types! \n\nDo you want to switch to the Visualization Tab and see the tracks?",
                    QMessageBox.Yes | QMessageBox.No
                )
                if res == QMessageBox.Yes:
                    first_panel = next(iter(self.panels.values()))
                    first_panel._switch_to_viz_and_show_tracks()

        except Exception as e:
            traceback.print_exc()
            self._log(f"\u274c Batch tracking error: {e}")
