import ipywidgets as widgets
from pathlib import Path
import pandas as pd
import yaml
import traceback
import fnmatch
from copy import deepcopy
from .utils import (
    _cfg_get, 
    _DEFAULT_CONFIG,
    spinning_loader
)
from behav3d.analysis.backprojection import (
    backproject_mean_features_behav3d,
    backproject_time_features_behav3d
)

# Need this for the default groups
try:
    from behav3d.defaults import behav3d_calculated_features
except ImportError:
    behav3d_calculated_features = {
        "movement": ["speed", "msd", "track_duration", "dist_to_origin"],
        "intensity": ["mean_intensity_*"],
        "morphology": ["volume", "sphericity", "surface_area"],
        "contact": ["*_contact", "*_distance"],
        "death": ["mean_dead_dye", "percentage_dead_mask"],
        "active_killing": ["is_active_killing", "killing_efficiency"]
    }

class BackprojectionPanel:
    """
    UI for backprojecting features to images (napari).
    """
    def __init__(self, metadata_loader, cell_type=None):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        from behav3d.io.images import load_image, load_zarr, save_as_zarr
        from behav3d.core.metadata import (
            detect_immune_cell_types_from_metadata,
            detect_organoid_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel
        )
        md = self.metadata_loader.metadata
        if md is None: raise RuntimeError("Metadata not loaded.")
        
        self.all_cell_types = detect_immune_cell_types_from_metadata(md) + \
                              detect_organoid_types_from_metadata(md) + \
                              detect_other_cell_types_from_metadata(md)
        if not self.all_cell_types: raise RuntimeError("No cell types detected.")
        if cell_type is None: cell_type = self.all_cell_types[0]
        
        # Load and prepare feature groups
        self._base_groups = deepcopy(behav3d_calculated_features)
        if not has_dead_channel(md):
            if "death" in self._base_groups: del self._base_groups["death"]
            if "intensity" in self._base_groups:
                self._base_groups["intensity"] = [f for f in self._base_groups["intensity"] if f != "mean_dead_dye"]
        
        for ct in self.all_cell_types:
            if f"{ct}_contact" not in self._base_groups["contact"]:
                self._base_groups["contact"].extend([f"{ct}_contact", f"{ct}_contact_pixels"])

        params = dict(self.metadata_loader.behav3d_parameters or {})
        self._cfg = params.setdefault("backprojection", self._default_backproj_config(cell_type))
        self._params = params

        # UI components
        sample_list = sorted(map(str, md["sample_name"].unique().tolist()))
        self.sample_dd = widgets.Dropdown(description="Sample", options=sample_list, value=self._cfg.get("last_sample") or (sample_list[0] if sample_list else None), layout=widgets.Layout(width="360px"))
        
        self._celltype_map = {ct.replace('_', ' ').title(): ct for ct in self.all_cell_types}
        ct_inv = {v: k for k, v in self._celltype_map.items()}
        self.celltype_dd = widgets.Dropdown(options=list(self._celltype_map.keys()), value=ct_inv.get(str(self._cfg.get("cell_type", cell_type)), list(self._celltype_map.keys())[0]), description="Cell type", layout=widgets.Layout(width="360px"))
        
        self._mode_map = {"Mean features": "mean", "Time features": "time"}
        self.mode_tb = widgets.ToggleButtons(options=list(self._mode_map.keys()), value="Mean features" if self._cfg.get("mode") == "mean" else "Time features", description="Mode")
        
        self.save_cb = widgets.Checkbox(description="Save .zarr to disk", value=bool(self._cfg.get("save", False)))
        self.refresh_btn = widgets.Button(description="Refresh Features", icon="refresh", layout=widgets.Layout(width="150px"))
        self.refresh_btn.on_click(self._on_refresh_clicked)
        
        # Feature Selection Dropdown
        self.feature_dd = widgets.Dropdown(description="Feature:", options=["Behavioral States Only"], value="Behavioral States Only", layout=widgets.Layout(width="360px"), style={'description_width': '100px'})
        
        self.status_html = widgets.HTML("<i>Loading available features...</i>")
        
        self.selection_container = widgets.VBox([])
        self._mean_group_rows = {}
        self._time_group_rows = {}
        self._time_selection_box = self._build_time_selection_box()
        self._rebuild_mean_selection_from_file()
        
        self.btn_run = widgets.Button(description="Run backprojection", button_style="success")
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()
        
        self.mode_tb.observe(lambda _: self._swap_ui(), names="value")
        self.celltype_dd.observe(lambda _: self._on_celltype_changed(), names="value")
        
        self.ui = widgets.VBox([
            widgets.HTML('<div style="font-size:22px;font-weight:700;">Backprojection</div>'),
            widgets.HTML('<div style="color:#666;font-size:12px;margin-bottom:10px;">Visualize behavioral states and individual features in Napari.</div>'),
            self.sample_dd, self.celltype_dd, self.mode_tb, self.feature_dd, self.save_cb,
            widgets.HBox([self.status_html, self.refresh_btn]),
            # self.selection_container is now less critical as we use the dropdown
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out
        ])
        self._on_celltype_changed()

    def _default_backproj_config(self, cell_type):
        return {"cell_type": cell_type, "mode": "mean", "save": False, "columns_input_time": [], "columns_input_mean": []}

    def _build_time_selection_box(self):
        # Implementation simplified for migration
        return widgets.HTML("<i>Time selection UI placeholder</i>")

    def _rebuild_mean_selection_from_file(self):
        # Implementation simplified for migration
        self.selection_container.children = [widgets.HTML("<i>Mean selection UI placeholder</i>")]

    def _swap_ui(self):
        m = self._mode_map[self.mode_tb.value]
        self.selection_container.children = [self._time_selection_box] if m == "time" else [widgets.HTML("<i>Mean selection placeholder</i>")]

    def _on_celltype_changed(self):
        self._update_available_features()

    def _update_available_features(self):
        """Populate the dropdown with available numeric features."""
        ct = self._celltype_map[self.celltype_dd.value]
        fdir = Path(self.output_dir, "analysis", ct, "track_features")
        
        # Check for filtered features first, then advanced
        p_filtered = fdir / f"BEHAV3D_{ct}_combined_track_features_filtered.csv"
        p_advanced = Path(self.output_dir, "analysis", ct, "active_killing", f"BEHAV3D_{ct}_advanced_track_features.csv")
        
        available_cols = []
        if p_advanced.exists():
            try: available_cols = pd.read_csv(p_advanced, nrows=0).columns.tolist()
            except Exception: pass
        elif p_filtered.exists():
            try: available_cols = pd.read_csv(p_filtered, nrows=0).columns.tolist()
            except Exception: pass
        
        if not available_cols:
            self.feature_dd.options = ["Behavioral States Only"]
            self.status_html.value = '<span style="color:#c00;">⚠️ No features found. Run extraction first.</span>'
            return

        # Filter numeric and exclude boilerplate
        exclude = ["TrackID", "sample_name", "position_t", "position_x", "position_y", "position_z", "ClusterID", "UMAP1", "UMAP2"]
        numeric_features = [c for c in available_cols if c not in exclude]
        
        # Sort and add default
        self.feature_dd.options = ["Behavioral States (Clusters)"] + sorted(numeric_features)
        self.status_html.value = f'<span style="color:green;">✓ Found {len(numeric_features)} features</span>'

    def _on_refresh_clicked(self, *_):
        self._update_available_features()

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True; self.spinner_html.layout.display = None; self.out.clear_output()
        with self.out:
            try:
                m = self._mode_map[self.mode_tb.value]
                ct = self._celltype_map[self.celltype_dd.value]
                
                # Collect columns: Behavioral States (ClusterID) + Selected Feature
                cols = ["ClusterID"]
                selected = self.feature_dd.value
                if selected != "Behavioral States (Clusters)":
                    cols.append(selected)
                
                fn = backproject_mean_features_behav3d if m == "mean" else backproject_time_features_behav3d
                fn(
                    metadata=self.metadata_loader.metadata, 
                    sample_name=self.sample_dd.value, 
                    output_dir=self.output_dir, 
                    cell_type=ct, 
                    columns=cols, 
                    save=bool(self.save_cb.value)
                )
                print("✅ Backprojection finished.")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False
