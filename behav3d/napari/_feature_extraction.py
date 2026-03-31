"""
BEHAV3D napari plugin – Feature Extraction Tab.

Provides per-cell-type sub-tabs with feature checkboxes (movement, intensity,
morphology, contact, death), thresholds, workers, and batch-run capability.
Mirrors the architecture of _tracking.py.
"""
import sys
import os
import traceback
from pathlib import Path

import yaml
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QFormLayout, QLabel,
    QPushButton, QTabWidget, QTextEdit, QCheckBox,
    QDoubleSpinBox, QSpinBox, QGroupBox, QMessageBox, QScrollArea,
    QSlider,
)
from qtpy.QtCore import Qt

from behav3d.napari._segmentation import make_help_row


# ═══════════════════════════════════════════════════════════════════════════
# Per-cell-type panel
# ═══════════════════════════════════════════════════════════════════════════
class CellTypeFeaturePanel(QWidget):
    """Feature extraction controls for one cell type."""

    def __init__(
        self,
        cell_type: str,
        category: str,
        metadata_loader,
        all_cell_types: list,
        category_types: list,
        log_callback=None,
        parent=None,
    ):
        super().__init__(parent)
        self.cell_type = cell_type
        self.category = category
        self.metadata_loader = metadata_loader
        self.all_cell_types = all_cell_types
        self.category_types = category_types
        self.log = log_callback or (lambda m: None)
        self.viewer = None
        self._preview_connected = False  # track whether spinner signal is connected

        # Read saved config
        params = self.metadata_loader.behav3d_parameters
        fcfg = params.get("features", {}).get(self.cell_type, {}) or {}

        # Detect dead channel
        md = metadata_loader.metadata
        self._has_dead = False
        if md is not None:
            self._has_dead = (
                "dead_channel" in md.columns and md["dead_channel"].notna().any()
            )

        self._init_ui(fcfg)

    # ---- UI ---------------------------------------------------------------
    def _init_ui(self, fcfg):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        # ── Feature checkboxes ────────────────────────────────────────
        feat_group = QGroupBox("Features to extract")
        feat_lay = QVBoxLayout(feat_group)
        feat_lay.setSpacing(2)

        all_feats = ["movement", "intensity", "morphology", "contact", "death"]
        default_feats = fcfg.get("features_choice", all_feats)
        if not isinstance(default_feats, (list, tuple)):
            default_feats = all_feats

        self.feature_checks: dict[str, QCheckBox] = {}
        for f in all_feats:
            if f == "death" and not self._has_dead:
                continue
            cb = QCheckBox(f.capitalize())
            cb.setChecked(f in default_feats)
            self.feature_checks[f] = cb

        # Contact threshold (only visible when Contact checked)
        self.contact_threshold = QDoubleSpinBox()
        self.contact_threshold.setRange(0.0, 500.0)
        self.contact_threshold.setSingleStep(0.5)
        self.contact_threshold.setDecimals(2)
        self.contact_threshold.setValue(float(fcfg.get("contact_threshold", 0.0)))
        self.contact_threshold.setMaximumWidth(90)
        self.contact_threshold.setSuffix(" µm")

        contact_row = QHBoxLayout()
        contact_row.setSpacing(6)
        contact_row.addWidget(self.feature_checks.get("contact", QCheckBox()))
        
        contact_help = make_help_row(
            self.contact_threshold,
            "Contact Threshold (µm)",
            "Distance threshold (in micrometers) for detecting intensity transfer "
            "between neighboring cells (contact features).\n\n"
            "Segments closer than this will be checked for contact signal."
        )
        contact_row.addLayout(contact_help)
        contact_row.addStretch()

        def _toggle_ct(state=None):
            vis = self.feature_checks.get("contact", QCheckBox()).isChecked()
            self.contact_threshold.setEnabled(vis)

        _toggle_ct()
        if "contact" in self.feature_checks:
            self.feature_checks["contact"].stateChanged.connect(_toggle_ct)

        for f in all_feats:
            if f == "contact":
                feat_lay.addLayout(contact_row)
            elif f in self.feature_checks:
                feat_lay.addWidget(self.feature_checks[f])

        layout.addWidget(feat_group)

        # ── Dead threshold ────────────────────────────────────────
        dead_group = QGroupBox("Death Threshold")
        dead_lay = QVBoxLayout(dead_group)
        dead_lay.setSpacing(4)

        dead_desc = QLabel(
            "Set the percentage of dead-mask pixels overlapping a segment\n"
            "above which the cell is classified as 'dead'. The 'dead' column\n"
            "is created during feature extraction when this threshold is set."
        )
        dead_desc.setWordWrap(True)
        dead_desc.setStyleSheet("color: #888; font-size: 10px;")
        dead_lay.addWidget(dead_desc)

        dead_form = QFormLayout()
        dead_form.setContentsMargins(0, 0, 0, 0)
        dead_form.setFieldGrowthPolicy(QFormLayout.FieldsStayAtSizeHint)

        self.spin_dead_threshold = QDoubleSpinBox()
        self.spin_dead_threshold.setRange(0.0, 100.0)
        self.spin_dead_threshold.setSingleStep(1.0)
        self.spin_dead_threshold.setDecimals(1)
        self.spin_dead_threshold.setValue(float(fcfg.get("dead_mask_percentage_threshold", 0.0)))
        self.spin_dead_threshold.setSuffix(" %")
        self.spin_dead_threshold.setMaximumWidth(100)
        dead_form.addRow("Dead mask % threshold:", make_help_row(
            self.spin_dead_threshold,
            "Dead Mask Percentage Threshold",
            "Percentage of dead-mask pixels overlapping a segment's\n"
            "volume required to classify the cell as dead.\n\n"
            "Set to 0 to skip dead classification.\n"
            "Typical range: 5–30%."
        ))
        dead_lay.addLayout(dead_form)

        # Preview button
        self.btn_preview_dead = QPushButton("\U0001F441 Preview Dead Threshold in Viewer")
        self.btn_preview_dead.setStyleSheet(
            "QPushButton { background: #37474F; color: white; padding: 5px 10px; "
            "border-radius: 3px; font-size: 11px; } "
            "QPushButton:hover { background: #546E7A; }"
        )
        self.btn_preview_dead.setToolTip(
            "Load the first sample's segments + dead mask data and show\n"
            "a dead/alive overlay in the napari viewer. Adjust the threshold\n"
            "slider to see how the classification changes."
        )
        self.btn_preview_dead.clicked.connect(self._on_preview_dead_clicked)
        dead_lay.addWidget(self.btn_preview_dead)

        dead_group.setVisible(self._has_dead)
        layout.addWidget(dead_group)
        self._dead_group = dead_group

        # Toggle dead group visibility with Death feature checkbox
        def _toggle_dead_group(state=None):
            cb = self.feature_checks.get("death")
            self._dead_group.setVisible(self._has_dead and (cb is not None and cb.isChecked()))
        if "death" in self.feature_checks:
            self.feature_checks["death"].stateChanged.connect(_toggle_dead_group)
        _toggle_dead_group()



        # ── Workers ────────────────────────────────────────────────────
        n_cores = os.cpu_count() or 4
        max_allowed = max(1, n_cores - 1)
        self.spin_workers = QSpinBox()
        self.spin_workers.setRange(1, max_allowed)
        default_workers = min(int(fcfg.get("n_workers", max(4, n_cores // 2))), max_allowed)
        self.spin_workers.setValue(default_workers)
        self.spin_workers.setMaximumWidth(60)
        self.spin_workers.valueChanged.connect(self._on_workers_changed)

        workers_form = QFormLayout()
        workers_form.setContentsMargins(0, 0, 0, 0)
        workers_form.addRow("Workers:", make_help_row(
            self.spin_workers,
            "Number of Workers",
            f"Number of CPU cores to use for parallel processing.\n\n"
            f"Your machine has {n_cores} cores.\n"
            f"Recommendation: Use at most {max(1, n_cores - 1)} cores to keep the system responsive."
        ))
        layout.addLayout(workers_form)

        # ── Apply settings buttons ────────────────────────────────────
        cat_label = self.category.capitalize() + "s" if self.category != "other" else "Other types"
        sync_row = QHBoxLayout()
        btn_apply_cat = QPushButton(f"Apply to all {cat_label}")
        btn_apply_cat.clicked.connect(lambda: self._apply_to_others(category_only=True))
        btn_apply_all = QPushButton("Apply to all")
        btn_apply_all.clicked.connect(lambda: self._apply_to_others(category_only=False))
        sync_row.addWidget(btn_apply_cat)
        sync_row.addWidget(btn_apply_all)
        layout.addLayout(sync_row)

        # ── Run button ────────────────────────────────────────────────
        self.btn_run = QPushButton(f"Run {self.cell_type.capitalize()} Feature Extraction")
        self.btn_run.setStyleSheet(
            "background-color: #28a745; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 6px;"
        )
        self.btn_run.clicked.connect(self._on_run_clicked)
        layout.addWidget(self.btn_run)

        layout.addStretch()

    # ---- Helpers ----------------------------------------------------------
    def _on_workers_changed(self, value):
        n_cores = os.cpu_count() or 4
        max_allowed = max(1, n_cores - 1)
        if value > max_allowed:
            self.spin_workers.setValue(max_allowed)
            self.log(
                f"⚠️ Workers capped to {max_allowed} (system has {n_cores} cores). "
                f"Using all cores can freeze the system."
            )

    def _selected_features(self) -> list:
        return [f for f, cb in self.feature_checks.items() if cb.isChecked()]

    def _collect_params(self) -> dict:
        thr = float(self.spin_dead_threshold.value())
        return {
            "features_choice": self._selected_features(),
            "contact_threshold": float(self.contact_threshold.value()),
            "dead_mask_percentage_threshold": thr if thr > 0 else None,
            "n_workers": int(self.spin_workers.value()),
        }

    def _apply_to_others(self, category_only=False):
        parent_tab = self.parent()
        while parent_tab and not hasattr(parent_tab, "panels"):
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
                p = parent_tab.panels[ct]
                for f, cb in p.feature_checks.items():
                    cb.setChecked(f in settings["features_choice"])
                p.contact_threshold.setValue(settings["contact_threshold"])
                p.spin_workers.setValue(settings["n_workers"])
                if hasattr(p, 'spin_dead_threshold'):
                    thr_val = settings.get("dead_mask_percentage_threshold")
                    p.spin_dead_threshold.setValue(thr_val if thr_val is not None else 0.0)
                count += 1
        scope = "category" if category_only else "all"
        self.log(f"Applied feature settings to {count} other cell types ({scope}).")

    def _persist(self):
        params = self.metadata_loader.behav3d_parameters
        features = params.setdefault("features", {})
        features[self.cell_type] = self._collect_params()

        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(params, f, sort_keys=False)
            except Exception as e:
                self.log(f"Warning: Could not save parameters: {e}")

    def _run_feature_extraction_for(self, cell_type: str, overwrite: bool = False):
        """Run feature extraction for a single cell type."""
        from behav3d.features.timepoint_features import run_feature_extraction

        print(f"\n{'='*50}", file=sys.stderr)
        print(f"  Feature Extraction: {cell_type}", file=sys.stderr)
        print(f"{'='*50}", file=sys.stderr)

        out_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        # Dead threshold
        dead_thr = float(self.spin_dead_threshold.value())
        dead_mask_pct = dead_thr if dead_thr > 0 else None

        run_feature_extraction(
            metadata=self.metadata_loader.metadata,
            output_dir=out_dir,
            cell_type=cell_type,
            features_choice=self._selected_features(),
            contact_threshold=float(self.contact_threshold.value()),
            dead_mask_percentage_threshold=dead_mask_pct,
            n_workers=int(self.spin_workers.value()),
            overwrite=overwrite,
        )

        print(f"✅ {cell_type} feature extraction finished.", file=sys.stderr)

    def _check_existing_features(self, cell_types: list) -> list:
        """Return descriptions of existing feature data."""
        warnings = []
        out_dir = Path(self.metadata_loader.output_dir)
        for ct in cell_types:
            feat_dir = out_dir / "analysis" / ct / "track_features"
            combined = feat_dir / f"BEHAV3D_{ct}_combined_track_features.csv"
            if combined.exists():
                warnings.append(f"{ct} feature data ({combined.name})")
        return warnings

    # ---- Click handler ----------------------------------------------------
    def _on_run_clicked(self, interactive=True):
        self._persist()
        self.log(f"Running feature extraction for: {self.cell_type}")

        overwrite = False
        existing = self._check_existing_features([self.cell_type])
        if existing:
            if interactive:
                details = "\n".join(f"  • {w}" for w in existing)
                box = QMessageBox(self)
                box.setWindowTitle("Overwrite Existing Features?")
                box.setText(
                    f"The following feature data already exists:\n\n{details}\n\n"
                    "What do you want to do?"
                )
                btn_overwrite = box.addButton("Overwrite", QMessageBox.DestructiveRole)
                btn_skip = box.addButton("Skip", QMessageBox.AcceptRole)
                btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
                box.setDefaultButton(btn_cancel)
                box.exec_()
                clicked = box.clickedButton()
                if clicked != btn_overwrite:
                    self.log(f"Feature extraction for {self.cell_type} cancelled.")
                    return
                overwrite = True
            else:
                # Default for queue (assuming dependency check was already done or user wants latest)
                overwrite = True

        try:
            self._run_feature_extraction_for(self.cell_type, overwrite=overwrite)
            self.log(f"✅ {self.cell_type} feature extraction finished.")
        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during feature extraction: {e}")

    # ---- Dead Threshold Preview -------------------------------------------
    def _on_preview_dead_clicked(self):
        """Load segments + dead mask for 1st sample & show dead/alive overlay."""
        try:
            import numpy as np
            from behav3d.io.images import load_image
            import pandas as pd

            viewer = self.viewer
            if viewer is None:
                self.log("⚠️ No viewer available for dead threshold preview.")
                return

            md = self.metadata_loader.metadata
            if md is None or md.empty:
                self.log("⚠️ No metadata loaded.")
                return

            # Find first sample with both segments & dead mask
            sample_row = None
            seg_path = None
            dead_mask_path = None
            for _, row in md.iterrows():
                # Find segments path for this cell type
                for prefix in ['or', 'im', 'ot']:
                    col = f"{prefix}_{self.cell_type}_segments_image_path"
                    if col in row.index and pd.notna(row.get(col)):
                        seg_path = row[col]
                        break
                # Find dead mask path
                dead_col = f"dead_mask_path"
                if dead_col in row.index and pd.notna(row.get(dead_col)):
                    dead_mask_path = row[dead_col]
                if seg_path and dead_mask_path:
                    sample_row = row
                    break

            if sample_row is None or seg_path is None or dead_mask_path is None:
                self.log("⚠️ Could not find a sample with both segments and dead mask paths.")
                return

            sample_name = sample_row.get('sample_name', 'Unknown')
            self.log(f"Loading preview for sample: {sample_name}")

            # Load segments and dead mask (use first timepoint)
            segments = load_image(seg_path)
            dead_mask = load_image(dead_mask_path)
            t = 0  # Preview at first timepoint
            seg_t = np.asarray(segments[t])
            dead_t = np.asarray(dead_mask[t])

            # Also try to load raw image for context (dead channel)
            raw_path = sample_row.get('raw_image_path')
            dead_ch_idx = sample_row.get('dead_channel')

            # Remove old preview layers
            for name in ['Dead Mask Preview', 'Dead/Alive Preview', 'Segments Preview', 'Dead Channel']:
                try:
                    viewer.layers.remove(name)
                except (KeyError, ValueError):
                    pass

            # Add dead channel if available
            if raw_path and pd.notna(dead_ch_idx):
                try:
                    raw = load_image(raw_path)
                    raw_t = np.asarray(raw[t])
                    ch_idx = int(dead_ch_idx)
                    if raw_t.ndim >= 2 and ch_idx < raw_t.shape[0]:
                        dead_ch_data = raw_t[ch_idx]
                        viewer.add_image(
                            dead_ch_data, name='Dead Channel',
                            colormap='magenta', blending='additive',
                            opacity=0.5
                        )
                except Exception:
                    pass  # Skip if raw loading fails

            # Add segments preview
            viewer.add_labels(seg_t, name='Segments Preview', opacity=0.3)

            # Add dead mask
            viewer.add_image(
                dead_t.astype(float), name='Dead Mask Preview',
                colormap='magenta', blending='additive', opacity=0.4
            )

            # Compute and add dead/alive overlay
            threshold = float(self.spin_dead_threshold.value())
            self._update_dead_alive_overlay(viewer, seg_t, dead_t, threshold)

            # Connect spinner to live-update if not already connected
            if not self._preview_connected:
                self._preview_seg_t = seg_t
                self._preview_dead_t = dead_t
                self.spin_dead_threshold.valueChanged.connect(self._on_threshold_spin_changed)
                self._preview_connected = True
            else:
                # Update cached data
                self._preview_seg_t = seg_t
                self._preview_dead_t = dead_t

            self.log(f"✅ Dead threshold preview loaded. Adjust threshold to see changes.")

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error loading dead threshold preview: {e}")

    def _on_threshold_spin_changed(self, value):
        """Live-update the dead/alive overlay when threshold changes."""
        viewer = self.viewer
        if viewer is None:
            return
        if not hasattr(self, '_preview_seg_t') or not hasattr(self, '_preview_dead_t'):
            return
        self._update_dead_alive_overlay(
            viewer, self._preview_seg_t, self._preview_dead_t, float(value)
        )

    @staticmethod
    def _update_dead_alive_overlay(viewer, seg_t, dead_t, threshold_pct):
        """Compute dead/alive per-segment overlay and update napari layer."""
        import numpy as np
        from skimage.measure import regionprops

        # Create a colored label image: green=alive (1), red=dead (2)
        overlay = np.zeros_like(seg_t, dtype=np.uint8)

        # Binarize dead mask
        dead_binary = dead_t > 0

        for region in regionprops(seg_t):
            label_val = region.label
            mask = seg_t == label_val
            total_pixels = mask.sum()
            if total_pixels == 0:
                continue
            dead_pixels = (mask & dead_binary).sum()
            pct = (dead_pixels / total_pixels) * 100.0

            if pct >= threshold_pct and threshold_pct > 0:
                overlay[mask] = 2  # Dead → red
            else:
                overlay[mask] = 1  # Alive → green

        # Define custom colormap: 0=transparent, 1=green, 2=red
        from napari.utils.colormaps import DirectLabelColormap
        try:
            cmap = DirectLabelColormap(
                color_dict={0: [0, 0, 0, 0], 1: [0, 0.8, 0, 0.5], 2: [0.9, 0, 0, 0.5]}
            )
        except Exception:
            cmap = None

        layer_name = 'Dead/Alive Preview'
        try:
            layer = viewer.layers[layer_name]
            layer.data = overlay
        except (KeyError, ValueError):
            if cmap is not None:
                viewer.add_labels(
                    overlay, name=layer_name, opacity=0.6,
                    colormap=cmap,
                )
            else:
                viewer.add_labels(overlay, name=layer_name, opacity=0.6)


# ═══════════════════════════════════════════════════════════════════════════
# FeatureExtractionTab — main tab with per-cell-type sub-tabs
# ═══════════════════════════════════════════════════════════════════════════
class FeatureExtractionTab(QWidget):
    def __init__(self, viewer=None, metadata_loader=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.panels: dict[str, CellTypeFeaturePanel] = {}

        self._init_ui()

        if hasattr(self.metadata_loader, "metadata_loaded"):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

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

        # Sub-tab widget
        self.cell_tabs = QTabWidget()
        self.cell_tabs.setTabPosition(QTabWidget.West)
        layout.addWidget(self.cell_tabs)

        # Global Run + Queue buttons
        self.btn_run_batch = QPushButton("Run Batch Feature Extraction (All Cell Types)")
        self.btn_run_batch.setStyleSheet(
            "background-color: #007bff; color: white; font-weight: bold; "
            "border-radius: 4px; padding: 10px; font-size: 14px;"
        )
        self.btn_run_batch.clicked.connect(self._on_run_batch_clicked)

        self.btn_queue_feature = QPushButton("+🛒")
        self.btn_queue_feature.setFixedSize(36, 32)
        self.btn_queue_feature.setToolTip("Add Feature Extraction to Processing Queue")
        self.btn_queue_feature.setStyleSheet(
            "QPushButton { background: #1a1a2e; color: #ffc107; border: 1px solid #ffc107; "
            "border-radius: 4px; font-size: 11px; font-weight: bold; }"
            "QPushButton:hover { background: #ffc107; color: #1a1a2e; }"
        )

        batch_row = QHBoxLayout()
        batch_row.setSpacing(4)
        batch_row.addWidget(self.btn_run_batch, stretch=1)
        batch_row.addWidget(self.btn_queue_feature)
        layout.addLayout(batch_row)

        # Hidden until metadata is loaded
        self.btn_run_batch.setVisible(False)
        self.btn_queue_feature.setVisible(False)

        # Log
        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(120)
        self.log_box.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log_box)

        # Placeholder
        self._placeholder = QLabel("Load metadata in the Data Preparation tab to see feature extraction options.")
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._placeholder.setStyleSheet("color: #888; font-style: italic;")
        self.cell_tabs.addTab(self._placeholder, "—")

    def _log(self, msg):
        import datetime
        ts = datetime.datetime.now().strftime("%H:%M:%S")
        self.log_box.append(f"[{ts}] {msg}")
        self.log_box.verticalScrollBar().setValue(
            self.log_box.verticalScrollBar().maximum()
        )

    def _on_metadata_updated(self):
        self._log("Metadata updated — refreshing feature extraction tabs…")
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
        self.cell_tabs.clear()
        self.panels.clear()

        org, imm, oth = self._detect_cell_types()
        all_types = org + imm + oth

        if not all_types:
            self.cell_tabs.addTab(self._placeholder, "—")
            self.btn_run_batch.setVisible(False)
            self.btn_queue_feature.setVisible(False)
            return

        self.btn_run_batch.setVisible(True)
        self.btn_queue_feature.setVisible(True)

        color_map = {"organoid": "🟣", "immune": "🔵", "other": "🟡"}

        for ct in org:
            self._add_panel(ct, "organoid", all_types, org, color_map)
        for ct in imm:
            self._add_panel(ct, "immune", all_types, imm, color_map)
        for ct in oth:
            self._add_panel(ct, "other", all_types, oth, color_map)

    def _add_panel(self, ct, category, all_types, cat_types, color_map):
        panel = CellTypeFeaturePanel(
            cell_type=ct, category=category,
            metadata_loader=self.metadata_loader,
            all_cell_types=all_types, category_types=cat_types,
            log_callback=self._log,
        )
        panel.viewer = self.viewer  # Pass viewer for dead threshold preview
        self.panels[ct] = panel
        icon = color_map.get(category, "")
        self.cell_tabs.addTab(panel, f"{icon} {ct}")

    # ---- Batch run --------------------------------------------------------
    def _on_run_batch_clicked(self):
        self.run_batch_feature_extraction(interactive=True)

    def run_batch_feature_extraction(self, interactive=True, skip_existing=False):
        """Run feature extraction for all cell types sequentially."""
        if not self.panels:
            self._log("No cell type panels available.")
            return

        total = len(self.panels)
        self._log(f"Starting batch feature extraction for {total} cell types…")

        print(f"\n{'='*60}", file=sys.stderr)
        print(f"  Running batch feature extraction for {total} cell types", file=sys.stderr)
        print(f"{'='*60}", file=sys.stderr)

        # Overwrite check
        all_cts = list(self.panels.keys())
        existing = []
        existing_cts = set()
        out_dir = Path(self.metadata_loader.output_dir)
        for ct in all_cts:
            feat_dir = out_dir / "analysis" / ct / "track_features"
            combined = feat_dir / f"BEHAV3D_{ct}_combined_track_features.csv"
            if combined.exists():
                existing.append(f"{ct} feature data ({combined.name})")
                existing_cts.add(ct)

        skip_existing_flag = skip_existing
        overwrite = not skip_existing
        if existing:
            if interactive:
                details = "\n".join(f"  • {w}" for w in existing)
                box = QMessageBox(self)
                box.setWindowTitle("Overwrite Existing Features?")
                box.setText(
                    f"The following feature data already exists:\n\n{details}\n\n"
                    "What do you want to do?"
                )
                btn_overwrite = box.addButton("Overwrite All", QMessageBox.DestructiveRole)
                btn_skip = box.addButton("Skip Existing", QMessageBox.AcceptRole)
                btn_cancel = box.addButton("Cancel", QMessageBox.RejectRole)
                box.setDefaultButton(btn_cancel)
                box.exec_()
                clicked = box.clickedButton()
                if clicked == btn_cancel:
                    self._log("Batch feature extraction cancelled.")
                    return
                skip_existing_flag = (clicked == btn_skip)
                overwrite = not skip_existing_flag
            else:
                overwrite = True

        try:
            for i, (ct, panel) in enumerate(self.panels.items(), 1):
                if skip_existing_flag and ct in existing_cts:
                    self._log(f"--- [{i}/{total}] Skipping {ct} (existing data) ---")
                    continue
                print(f"\n▶ [{i}/{total}] Feature extraction: {ct}…", file=sys.stderr)
                self._log(f"--- [{i}/{total}] Feature extraction: {ct} ---")
                panel._persist()
                panel._run_feature_extraction_for(ct, overwrite=overwrite)
                self._log(f"Done: {ct}")

            self._log("✅ Batch feature extraction finished.")
            print(f"\n{'='*60}", file=sys.stderr)
            print(f"  ✅ Batch feature extraction complete", file=sys.stderr)
            print(f"{'='*60}\n", file=sys.stderr)

        except Exception as e:
            traceback.print_exc()
            self._log(f"❌ Batch feature extraction error: {e}")
