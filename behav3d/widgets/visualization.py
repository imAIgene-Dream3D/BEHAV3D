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
        from behav3d.metadata import (
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
        self.status_html = widgets.HTML("")
        self.refresh_btn = widgets.Button(description="Refresh", icon="refresh")
        self.refresh_btn.on_click(self._on_refresh_clicked)
        
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
            self.sample_dd, self.celltype_dd, self.mode_tb, self.save_cb,
            widgets.HBox([self.status_html, self.refresh_btn]),
            self.selection_container,
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out
        ])
        self._swap_ui()

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

    def _on_celltype_changed(self): self._rebuild_mean_selection_from_file()
    def _on_refresh_clicked(self, *_): self._rebuild_mean_selection_from_file()

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True; self.spinner_html.layout.display = None
        with self.out:
            try:
                m = self._mode_map[self.mode_tb.value]
                fn = backproject_mean_features_behav3d if m == "mean" else backproject_time_features_behav3d
                fn(metadata=self.metadata_loader.metadata, sample_name=self.sample_dd.value, output_dir=self.output_dir, cell_type=self._celltype_map[self.celltype_dd.value], columns=[], save=bool(self.save_cb.value))
                print("✅ Backprojection finished.")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False
