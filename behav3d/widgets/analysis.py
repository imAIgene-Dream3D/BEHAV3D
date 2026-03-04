import ipywidgets as widgets
from pathlib import Path
import os
import yaml
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from io import BytesIO
from PIL import Image as PILImage, ImageDraw, ImageFont
from copy import deepcopy
from .utils import (
    _cfg_get, 
    spinning_loader,
    detect_cell_type_category,
    _DEFAULT_CONFIG
)
from behav3d.core.utils import expand_column_patterns
from behav3d.io.formats.zarr import load_zarr
from behav3d.io.images import load_image
from behav3d.features.timepoint_features import run_feature_extraction
from behav3d.analysis.tcell_analysis import (
    filter_cell_tracks,
    run_tcell_analysis
)
from behav3d.analysis.organoid_analysis import (
    filter_organoid_tracks,
    run_organoid_analysis,
    plot_multi_organoid_death_dynamics,
    run_organoid_morphology_dead_analysis
)

from behav3d.io.images import load_zarr
from behav3d.analysis import summarize_track_features
from behav3d.features.advanced_timepoint_features import run_active_killing_analysis
from behav3d.analysis.interaction_analysis import run_interaction_analysis

# Fallback for default features
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

class FeatureExtractionPanel:
    """
    Minimal UI for feature extraction (per-cell-type).
    """
    def __init__(self, metadata_loader, cell_type="tcell"):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()

        fcfg = _cfg_get(self.metadata_loader.behav3d_parameters, f"features.{self.cell_type}", {}) or {}

        has_dead_channel = False
        if hasattr(metadata_loader, 'metadata') and metadata_loader.metadata is not None:
            has_dead_channel = 'dead_channel' in metadata_loader.metadata.columns and \
                              metadata_loader.metadata['dead_channel'].notna().any()
        
        self.dead_mask_threshold = widgets.FloatText(
            description="Dead mask % thr",
            value=float(fcfg.get("dead_mask_percentage_threshold", 0.05)),
            style={'description_width':'160px'}
        )
        if not has_dead_channel: self.dead_mask_threshold.layout.display = "none"

        self._all_features = ["movement", "intensity", "morphology", "contact", "death"]
        default_feats = fcfg.get("features_choice", self._all_features)
        if not isinstance(default_feats, (list, tuple)): default_feats = self._all_features
        default_feats = [f for f in default_feats if f in self._all_features] or self._all_features

        self.feature_checks = {
            f: widgets.Checkbox(description=f.capitalize(), value=(f in default_feats), indent=False)
            for f in self._all_features
        }

        self.contact_threshold = widgets.FloatText(
            description="Contact Threshold (µm)",
            value=float(fcfg.get("contact_threshold", 0.0)),
            layout=widgets.Layout(width="240px"),
            style={'description_width': '180px'}
        )

        contact_row = widgets.HBox([self.feature_checks["contact"], self.contact_threshold], layout=widgets.Layout(align_items="center", gap="12px"))

        def _toggle_contact_threshold(change=None):
            self.contact_threshold.layout.display = None if self.feature_checks["contact"].value else "none"

        _toggle_contact_threshold()
        self.feature_checks["contact"].observe(_toggle_contact_threshold, names="value")

        feat_rows = []
        for f in self._all_features:
            if f == "contact": feat_rows.append(contact_row)
            elif f == "death" and not has_dead_channel: continue
            else: feat_rows.append(self.feature_checks[f])
                
        self.features_box = widgets.VBox([widgets.HTML("<b>Features</b>")] + feat_rows)

        self.n_workers = widgets.IntText(
            description="Workers",
            value=int(fcfg.get("n_workers", max(8, (os.cpu_count() or 8)))),
            max=max(8, (os.cpu_count() or 8)),
            style={'description_width':'160px'}
        )

        self.overwrite = widgets.Checkbox(description="Overwrite existing", value=bool(fcfg.get("overwrite", False)))

        self.btn_run = widgets.Button(description=f"Run {cell_type} feature extraction", button_style="success", layout=widgets.Layout(width="fit-content", flex="0 0 auto"))
        self.btn_run.on_click(self._on_run_clicked)
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"

        self.run_row = widgets.HBox([self.btn_run, self.spinner_html], layout=widgets.Layout(align_items="center", gap="10px"))
        self.out = widgets.Output()

        self.ui = widgets.VBox([
            widgets.HTML(f"<b> {self.cell_type} Feature Extraction</b>"),
            widgets.HBox([self.dead_mask_threshold, self.n_workers]),
            self.features_box,
            self.overwrite,
            self.run_row,
            widgets.HTML("<hr>"),
            self.out
        ])

    def _selected_features(self):
        return [f for f in self._all_features if self.feature_checks[f].value]

    def _persist_params(self):
        params = self.metadata_loader.behav3d_parameters
        prof = params.setdefault("features", {}).setdefault(self.cell_type, {})
        prof.update({
            "dead_mask_percentage_threshold": float(self.dead_mask_threshold.value),
            "features_choice": self._selected_features(),
            "n_workers": int(self.n_workers.value),
            "overwrite": bool(self.overwrite.value),
            "contact_threshold": float(self.contact_threshold.value)
        })
        if getattr(self.metadata_loader, "metadata_csv_path", None):
            params.setdefault("paths", {})["metadata_csv"] = str(Path(self.metadata_loader.metadata_csv_path).expanduser())
        if getattr(self.metadata_loader, "output_dir", None):
            params.setdefault("paths", {})["output_dir"] = str(Path(self.metadata_loader.output_dir).expanduser())

        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(params, f, sort_keys=False)

    def _lock(self, state: bool):
        for w in [self.dead_mask_threshold, self.n_workers, self.overwrite, self.btn_run, *self.feature_checks.values()]:
            if hasattr(w, "disabled"): w.disabled = state

    def _on_run_clicked(self, _):
        self._lock(True); self.spinner_html.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                out_dir = Path(self.metadata_loader.output_dir).expanduser()
                out_dir.mkdir(parents=True, exist_ok=True)
                run_feature_extraction(
                    dead_mask_percentage_threshold=float(self.dead_mask_threshold.value),
                    contact_threshold=float(self.contact_threshold.value),
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(out_dir),
                    features_choice=self._selected_features(),
                    cell_type=self.cell_type,
                    n_workers=int(self.n_workers.value),
                    overwrite=bool(self.overwrite.value)
                )
                print("✅ Feature extraction finished.")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self._lock(False)

class TrackFilterPanel:
    """
    Generic track filtering panel that works for ANY cell type.
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        self.category = detect_cell_type_category(self.cell_type, metadata_loader.metadata)
        
        active_killing_dir = Path(self.output_dir, "analysis", self.cell_type, "active_killing")
        self._advanced_features_path = Path(active_killing_dir, f"BEHAV3D_{self.cell_type}_advanced_track_features.csv")
        self._use_advanced_features = self._advanced_features_path.exists()
        
        params = self.metadata_loader.behav3d_parameters
        cfg = params.setdefault("track_filtering", {}).setdefault(self.cell_type, {})
        
        self.exp_duration = widgets.IntText(description=f"Max timepoints ({cfg.get('time_type','frames')})", value=int(cfg.get("exp_duration", 350)), style={'description_width': '180px'})
        self.en_exp_duration = widgets.Checkbox(description="Trim down full time series", value=bool(cfg.get("exp_duration_enabled", True)), indent=False)
        self.row_exp = widgets.HBox([self.exp_duration], layout=widgets.Layout(display=(None if self.en_exp_duration.value else "none")))
        self.en_exp_duration.observe(lambda c: setattr(self.row_exp.layout, 'display', None if c['new'] else 'none'), names="value")

        self.min_track_length = widgets.IntText(description=f"Minimal length ({cfg.get('time_type','frames')})", value=int(cfg.get("min_track_length", 30)), style={'description_width': '180px'})
        self.en_min_length = widgets.Checkbox(description="Select only tracks with minimal length", value=bool(cfg.get("min_length_enabled", True)), indent=False)
        self.row_min = widgets.HBox([self.min_track_length], layout=widgets.Layout(display=(None if self.en_min_length.value else "none")))
        self.en_min_length.observe(lambda c: setattr(self.row_min.layout, 'display', None if c['new'] else 'none'), names="value")

        self.max_track_length = widgets.IntText(description=f"Maximal length ({cfg.get('time_type','frames')})", value=int(cfg.get("max_track_length", 30)), style={'description_width': '180px'})
        self.en_max_length = widgets.Checkbox(description="Trim down tracks to supplied length", value=bool(cfg.get("max_length_enabled", True)), indent=False)
        self.row_max = widgets.HBox([self.max_track_length], layout=widgets.Layout(display=(None if self.en_max_length.value else "none")))
        self.en_max_length.observe(lambda c: setattr(self.row_max.layout, 'display', None if c['new'] else 'none'), names="value")

        from behav3d.core.metadata import has_dead_channel
        self.has_dead = has_dead_channel(self.metadata_loader.metadata)
        
        if self.category == "organoid":
            self.filter_t0_dead = widgets.Checkbox(description="Filter by minimal size at t=1", value=bool(cfg.get("filter_min_size_t1", True)), indent=False)
            self.min_size_t1 = widgets.IntText(description="Minimal size (px @ t=1)", value=int(cfg.get("min_size_t1", 1000)), style={'description_width': '180px'})
            self.row_size_t1 = widgets.HBox([self.min_size_t1], layout=widgets.Layout(display=(None if self.filter_t0_dead.value else "none")))
            self.filter_t0_dead.observe(lambda c: setattr(self.row_size_t1.layout, 'display', None if c['new'] else 'none'), names="value")
        else:
            self.filter_t0_dead = widgets.Checkbox(description="Filter tracks that are dead at t=0", value=bool(cfg.get("filter_t0_dead", True)), indent=False)
            if not self.has_dead: self.filter_t0_dead.layout.display = "none"
            self.min_size_t1 = self.row_size_t1 = None
        
        self.time_type = widgets.ToggleButtons(options=["frames", "hours"], value=cfg.get("time_type", "frames"))
        def _on_time_unit_change(change):
            if change['name'] == 'value':
                u = change['new']
                self.exp_duration.description = f"Max timepoints ({u})"
                self.min_track_length.description = f"Minimal length ({u})"
                self.max_track_length.description = f"Maximal length ({u})"
        self.time_type.observe(_on_time_unit_change, names='value')

        self.btn_run = widgets.Button(description=f"Filter {self.cell_type} tracks & summarize", button_style="success", layout=widgets.Layout(width="fit-content"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out_run = widgets.Output()
        
        ui_elements = [widgets.HTML(f'<div style="font-size:22px;font-weight:700;">{self.cell_type.capitalize()} Track Filtering</div>'), self.en_exp_duration, self.row_exp, self.en_min_length, self.row_min, self.en_max_length, self.row_max, self.filter_t0_dead]
        if self.row_size_t1: ui_elements.append(self.row_size_t1)
        if self.category != "organoid": 
             ui_elements.extend([widgets.HTML('<b>Unit for time-based filters:</b>'), self.time_type])
        ui_elements.extend([widgets.HBox([self.btn_run, self.spinner_html]), self.out_run])
        self.ui = widgets.VBox(ui_elements)
        self._filter_fn = filter_organoid_tracks if self.category == "organoid" else filter_cell_tracks

    def _lock(self, locked):
        w_list = [self.en_exp_duration, self.exp_duration, self.en_min_length, self.min_track_length, self.en_max_length, self.max_track_length, self.filter_t0_dead, self.time_type, self.btn_run]
        if self.min_size_t1: w_list.append(self.min_size_t1)
        for w in w_list: w.disabled = locked

    def _on_run_clicked(self, *_):
        self._lock(True); self.out_run.clear_output(); self.spinner_html.layout.display = None
        with self.out_run:
            try:
                cfg = self.metadata_loader.behav3d_parameters["track_filtering"][self.cell_type]
                cfg.update({"exp_duration_enabled": self.en_exp_duration.value, "exp_duration": int(self.exp_duration.value), "min_length_enabled": self.en_min_length.value, "min_track_length": int(self.min_track_length.value), "max_length_enabled": self.en_max_length.value, "max_track_length": int(self.max_track_length.value), "time_type": str(self.time_type.value)})
                if self.category == "organoid":
                    cfg.update({"filter_min_size_t1": self.filter_t0_dead.value, "min_size_t1": int(self.min_size_t1.value)})
                else: cfg["filter_t0_dead"] = bool(self.filter_t0_dead.value)
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f: yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)

                df_input_path = str(self._advanced_features_path) if self._use_advanced_features else None
                filter_kwargs = {"metadata": self.metadata_loader.metadata, "output_dir": self.output_dir, "exp_duration": (int(self.exp_duration.value) if self.en_exp_duration.value else None), "min_track_length": (int(self.min_track_length.value) if self.en_min_length.value else None), "max_track_length": (int(self.max_track_length.value) if self.en_max_length.value else None), "time_type": str(self.time_type.value), "cell_type": self.cell_type, "df_input_path": df_input_path}
                if self.category == "organoid": filter_kwargs["min_size"] = (int(self.min_size_t1.value) if self.filter_t0_dead.value else None)
                else: filter_kwargs.update({"filter_t0_dead": (bool(self.filter_t0_dead.value) if self.has_dead else False), "plot_results": True})
                
                self._filter_fn(**filter_kwargs)
                print("✅ Filtering done. Summarizing tracks...")
                summarize_track_features(output_dir=self.output_dir, cell_type=self.cell_type)
                print(f"✅ {self.cell_type} filtering complete!")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self._lock(False)

class ActiveKillingPanel:
    """
    Advanced feature extraction panel for Active Killing Analysis.
    """
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        self.killing_results = None
        
        from behav3d.core.metadata import (
            detect_immune_cell_types_from_metadata,
            detect_organoid_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel
        )
        md = self.metadata_loader.metadata
        if md is None: raise RuntimeError("Metadata not loaded.")
        
        self.immune_types = detect_immune_cell_types_from_metadata(md)
        self.organoid_types = detect_organoid_types_from_metadata(md)
        self.other_types = detect_other_cell_types_from_metadata(md)
        self.potential_immune = self.immune_types + self.other_types
        self.target_types = self.organoid_types
        
        params = dict(self.metadata_loader.behav3d_parameters or {})
        self._cfg = params.setdefault("active_killing", deepcopy(_DEFAULT_CONFIG.get("active_killing", {})))
        
        # UI
        immune_options = self.potential_immune if self.potential_immune else ["(none detected)"]
        self.immune_dd = widgets.Dropdown(options=immune_options, value=immune_options[0] if immune_options else None, description="Immune cell:", style={'description_width': '120px'}, layout=widgets.Layout(width="280px"))
        
        target_info = ", ".join(self.target_types) if self.target_types else "(none detected)"
        self.target_info_html = widgets.HTML(f'<div style="padding:5px;background:#f0f0f0;border-radius:4px;"><b>Target cell types:</b> {target_info}</div>')
        
        self.observation_window = widgets.IntText(description="Observation window:", value=int(self._cfg.get("observation_window", 5)), style={'description_width': '150px'}, layout=widgets.Layout(width="220px"))
        self.death_signal_dd = widgets.Dropdown(options=["mean_dead_dye", "percentage_dead_mask", "nr_dead_mask_pixels"], value=self._cfg.get("death_signal_column", "mean_dead_dye"), description="Death signal:", style={'description_width': '150px'}, layout=widgets.Layout(width="300px"))
        self.killing_threshold = widgets.FloatText(description="Killing threshold:", value=float(self._cfg.get("killing_threshold_multiplier", 1.5)), style={'description_width': '150px'}, layout=widgets.Layout(width="220px"))
        self.min_contact_duration = widgets.IntText(description="Min contact duration:", value=int(self._cfg.get("min_contact_duration", 1)), style={'description_width': '150px'}, layout=widgets.Layout(width="220px"))
        
        # Absolute threshold option
        self.use_absolute_threshold = widgets.Checkbox(description="Use absolute threshold", value=bool(self._cfg.get("use_absolute_threshold", False)), indent=False, layout=widgets.Layout(width="200px"))
        self.absolute_threshold = widgets.FloatText(description="Absolute threshold:", value=float(self._cfg.get("absolute_killing_threshold", 0.0) or 0.0), style={'description_width': '150px'}, layout=widgets.Layout(width="220px"), disabled=not self._cfg.get("use_absolute_threshold", False))
        self.use_absolute_threshold.observe(self._on_absolute_threshold_toggle, names="value")
        
        self.save_results = widgets.Checkbox(description="Save results to CSV", value=bool(self._cfg.get("save_results", True)), indent=False)
        self.gallery_item_count = widgets.IntText(description="Immune cells in gallery:", value=5, style={'description_width': '180px'}, layout=widgets.Layout(width="280px"))
        
        self.btn_run = widgets.Button(description="Run Active Killing Analysis", button_style="danger", icon="bolt", layout=widgets.Layout(width="260px"))
        self.btn_run.on_click(self._on_run_clicked)
        
        self.btn_show_results = widgets.Button(description="Display Existing Results", button_style="success", icon="eye", layout=widgets.Layout(width="260px"))
        self.btn_show_results.on_click(lambda _: self._update_gallery_from_file(self._existing_results_path))
        self.btn_show_results.layout.display = "none"
        self._existing_results_path = None
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
        self.out = widgets.Output()
        self.insights_out = widgets.Output()
        self.gallery_out = widgets.Output()
        
        self.insights_accordion = widgets.Accordion(children=[widgets.VBox([self.insights_out, self.gallery_out])], selected_index=None)
        self.insights_accordion.set_title(0, "Active Killing Insights & Gallery")
        self.insights_accordion.layout.display = "none"
        self.validation_html = widgets.HTML("")
        self._validate_inputs()
        self.immune_dd.observe(lambda _: self._validate_inputs(), names="value")
        
        self.ui = widgets.VBox([
            widgets.HTML('<div style="font-size:22px;font-weight:700;">Active Killing Analysis</div>'),
            widgets.HTML('<div style="color:#555;font-size:13px;margin-bottom:10px;">Detects functional killing events. <b>Targets:</b> Auto-detected organoids.</div>'),
            self.immune_dd, self.target_info_html, self.validation_html, widgets.HTML("<hr>"),
            widgets.HBox([self.observation_window, widgets.HTML('<span style="color:#666;font-size:12px;">timepoints after contact</span>'), widgets.HTML("&nbsp;&nbsp;&nbsp;"), self.killing_threshold, widgets.HTML('<span style="color:#666;font-size:12px;">× background rate</span>')], layout=widgets.Layout(align_items="center")),
            widgets.HBox([self.min_contact_duration, widgets.HTML('<span style="color:#666;font-size:12px;">timepoints</span>'), widgets.HTML("&nbsp;&nbsp;&nbsp;"), self.death_signal_dd], layout=widgets.Layout(align_items="center")),
            widgets.HBox([self.use_absolute_threshold, self.absolute_threshold, widgets.HTML('<span style="color:#666;font-size:12px;">death signal increase (bypasses multiplier)</span>')], layout=widgets.Layout(align_items="center")),
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html, self.save_results, self.gallery_item_count, widgets.HTML('<span style="color:#666;font-size:12px;">x sample</span>')], layout=widgets.Layout(align_items="center", gap="15px")),
            self.btn_show_results,
            self.insights_accordion,
            self.out
        ])
        
        # Check if results already exist to enable visualize button
        self._check_existing_results()
    
    def _validate_inputs(self):
        immune = self.immune_dd.value
        messages = []
        valid = True
        if immune == "(none detected)": messages.append("⚠️ No immune cell types detected."); valid = False
        if not self.target_types: messages.append("⚠️ No organoid types detected."); valid = False
        if valid:
            p = Path(self.output_dir, "analysis", immune, "track_features", f"BEHAV3D_{immune}_combined_track_features.csv")
            if not p.with_name(f"BEHAV3D_{immune}_combined_track_features_filtered.csv").exists() and not p.exists():
                messages.append(f"⚠️ {immune} tracks not found. Run feature extraction first."); valid = False
        self.validation_html.value = '<span style="color:green;">✓ Ready</span>' if valid else '<br>'.join([f'<span style="color:#c00;">{m}</span>' for m in messages])
        self.btn_run.disabled = not valid

    def _on_absolute_threshold_toggle(self, change):
        """Enable/disable absolute threshold input based on checkbox."""
        self.absolute_threshold.disabled = not change["new"]
        # Disable multiplier display hint when using absolute
        if change["new"]:
            self.killing_threshold.disabled = True
        else:
            self.killing_threshold.disabled = False

    def _check_existing_results(self):
        """Check if analysis results already exist to enable the display button."""
        immune = self.immune_dd.value
        if immune == "(none detected)": return
        
        results_dir = Path(self.output_dir, "analysis", immune, "active_killing")
        advanced_path = results_dir / f"BEHAV3D_{immune}_advanced_track_features.csv"
        if advanced_path.exists():
            self._existing_results_path = advanced_path
            self.btn_show_results.layout.display = None
        else:
            self.btn_show_results.layout.display = "none"

    def _update_gallery_from_file(self, advanced_path):
        """Load results from CSV and display gallery."""
        self.gallery_out.clear_output()
        with self.gallery_out:
            display(widgets.HTML('<b style="color:#007bff;">⏳ Loading gallery insights...</b>'))
            try:
                df_killing = pd.read_csv(advanced_path)
                # Normalize columns immediately
                if "TrackID" in df_killing.columns and "immune_track_id" not in df_killing.columns:
                    df_killing["immune_track_id"] = df_killing["TrackID"]
                
                self.killing_results = df_killing # Store for fallback
                
                # Coordinate normalization
                coord_map = {"centroid-0": "position_z", "centroid-1": "position_y", "centroid-2": "position_x"}
                for old_col, new_col in coord_map.items():
                    if old_col in df_killing.columns and new_col not in df_killing.columns:
                        df_killing[new_col] = df_killing[old_col]
                
                # Ensure visibility
                self.insights_accordion.layout.display = None
                self.insights_accordion.selected_index = 0
                
                self._display_insights(df_killing)
                self._display_gallery(df_killing)
            except Exception:
                traceback.print_exc()

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()
        self.insights_out.clear_output()
        self.gallery_out.clear_output()
        self.insights_accordion.selected_index = None # Collapse accordion
        self.insights_accordion.layout.display = "none"
        
        immune = self.immune_dd.value
        if immune == "(none detected)":
            with self.out: print("Please select an immune cell type.")
            self.btn_run.disabled = False
            self.spinner_html.layout.display = "none"
            return

        with self.out:
            try:
                print(f"Starting Active Killing Analysis for {immune}...")
                
                # Save current parameters
                self._cfg["observation_window"] = self.observation_window.value
                self._cfg["death_signal_column"] = self.death_signal_dd.value
                self._cfg["killing_threshold_multiplier"] = self.killing_threshold.value
                self._cfg["min_contact_duration"] = self.min_contact_duration.value
                self._cfg["use_absolute_threshold"] = self.use_absolute_threshold.value
                self._cfg["absolute_killing_threshold"] = self.absolute_threshold.value
                self._cfg["save_results"] = self.save_results.value
                
                # Ensure the 'active_killing' key exists in behav3d_parameters
                if "active_killing" not in self.metadata_loader.behav3d_parameters:
                    self.metadata_loader.behav3d_parameters["active_killing"] = {}
                self.metadata_loader.behav3d_parameters["active_killing"].update(self._cfg)
                
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
                
                # Run analysis
                results_dir = Path(self.output_dir, "analysis", immune, "active_killing")
                results_dir.mkdir(parents=True, exist_ok=True)
                
                # Determine absolute threshold value
                abs_thresh = float(self.absolute_threshold.value) if self.use_absolute_threshold.value else None
                
                df_killing, _, _ = run_active_killing_analysis(
                    metadata=self.metadata_loader.metadata,
                    output_dir=self.output_dir,
                    immune_cell_type=immune,
                    target_cell_types=self.target_types,
                    observation_window=self.observation_window.value,
                    death_signal_column=self.death_signal_dd.value,
                    killing_threshold_multiplier=float(self.killing_threshold.value),
                    min_contact_duration=int(self.min_contact_duration.value),
                    absolute_killing_threshold=abs_thresh,
                    save_results=bool(self.save_results.value)
                )
                
                print("✅ Active Killing Analysis complete! Loading gallery...")
                # Use the automated loader to ensure enriched columns (coordinates) are present
                results_dir = Path(self.output_dir, "analysis", immune, "active_killing")
                advanced_path = results_dir / f"BEHAV3D_{immune}_advanced_track_features.csv"
                
                if advanced_path.exists():
                    self._update_gallery_from_file(advanced_path)
                else:
                    # Fallback if file wasn't created for some reason
                    self.killing_results = df_killing
                    self.insights_accordion.layout.display = None
                    self.insights_accordion.selected_index = 0
                    self._display_insights(df_killing)
                    self._display_gallery(df_killing)

            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False

    def _display_insights(self, df_killing):
        """Display graphical insights in the accordion."""
        import matplotlib.pyplot as plt
        import seaborn as sns
        
        self.insights_out.clear_output()
        if df_killing.empty: return
        
        # Normalize columns: enriched CSV uses TrackID, analysis uses immune_track_id
        if "TrackID" in df_killing.columns and "immune_track_id" not in df_killing.columns:
            df_killing["immune_track_id"] = df_killing["TrackID"]
        
        with self.insights_out:
            # 1. Per-sample kinetics (4 plots each)
            samples = df_killing["sample_name"].unique()
            for sample_name in samples:
                df_sample_all = df_killing[df_killing["sample_name"] == sample_name]
                df_sample_active = df_sample_all[df_sample_all["is_active_killing"]].copy()
                
                display(widgets.HTML(f"<h3 style='margin-bottom:5px;'>Sample: {sample_name}</h3>"))
                
                fig, axes = plt.subplots(2, 2, figsize=(16, 10))
                axes = axes.flatten()
                
                # Plot 1: Efficiency Distribution (Per Sample)
                if not df_sample_active.empty:
                    sns.histplot(df_sample_active["killing_efficiency"], kde=True, ax=axes[0], color='red')
                axes[0].set_title("1. Killing Efficiency Distribution", fontsize=14, fontweight='bold')
                axes[0].set_xlabel("Efficiency Score (signal increase / expected background)")
                axes[0].set_ylabel("Active Killing Events")
                
                # Plot 2: Smoothed Kinetics
                temp_counts = df_sample_active.groupby("position_t").size()
                if not temp_counts.empty:
                    max_t = int(df_sample_all["position_t"].max())
                    full_t = pd.Series(0, index=range(max_t + 1))
                    full_t.update(temp_counts)
                    window = max(5, max_t // 20)
                    smoothed = full_t.rolling(window=window, center=True).mean()
                    
                    axes[1].plot(full_t.index, full_t.values, color='darkred', alpha=0.2, label='Raw Counts')
                    axes[1].plot(smoothed.index, smoothed.values, color='red', linewidth=2, label='Kinetics Trend')
                    axes[1].fill_between(smoothed.index, 0, smoothed.values, color='red', alpha=0.1)
                    axes[1].set_title("2. Killing Intensity (Smoothed)", fontsize=14, fontweight='bold')
                    axes[1].set_xlabel("Timepoint")
                    axes[1].set_ylabel("Events / Timepoint")
                    axes[1].legend()

                # Plot 3: Cumulative Progress
                if not temp_counts.empty:
                    cumulative = full_t.cumsum()
                    axes[2].plot(cumulative.index, cumulative.values, color='darkblue', linewidth=3)
                    axes[2].fill_between(cumulative.index, 0, cumulative.values, color='blue', alpha=0.1)
                    axes[2].set_title("3. Cumulative Killing Progress", fontsize=14, fontweight='bold')
                    axes[2].set_xlabel("Timepoint")
                    axes[2].set_ylabel("Cumulative Sum of Active Killing Events")
                    axes[2].grid(True, linestyle='--', alpha=0.6)

                # Plot 4: Raster
                if not df_sample_active.empty:
                    df_sample_active['hitter_label'] = df_sample_active['immune_track_id'].astype(str)
                    hitters = df_sample_active.groupby('hitter_label').size().sort_values(ascending=False)
                    hitter_map = {name: i for i, name in enumerate(hitters.index)}
                    df_sample_active['y_idx'] = df_sample_active['hitter_label'].map(hitter_map)
                    
                    scatter = axes[3].scatter(df_sample_active["position_t"], df_sample_active["y_idx"], 
                                              c=df_sample_active["killing_efficiency"], cmap='YlOrRd', 
                                              s=50, edgecolors='black', alpha=0.8)
                    axes[3].set_title("4. Killing Event Raster (By T-cell)", fontsize=14, fontweight='bold')
                    axes[3].set_xlabel("Timepoint")
                    axes[3].invert_yaxis()
                    
                    if len(hitters) <= 25:
                        axes[3].set_yticks(range(len(hitters)))
                        axes[3].set_yticklabels([str(n) for n in hitters.index], fontsize=8)
                        axes[3].set_ylabel("T-cell TrackID (Ranked)")
                    else:
                        axes[3].set_yticks(range(15))
                        axes[3].set_yticklabels([str(n) for n in hitters.index[:15]], fontsize=8)
                        axes[3].set_ylabel("T-cell TrackID (Top 15 shown, ranked)")
                    plt.colorbar(scatter, ax=axes[3], label='Efficiency')

                plt.tight_layout()
                
                # Save per-sample plot
                sample_plot_dir = Path(self.output_dir, "analysis", self.immune_dd.value, "active_killing", "plots", sample_name)
                sample_plot_dir.mkdir(parents=True, exist_ok=True)
                plot_path = sample_plot_dir / f"killing_kinetics_summary_{sample_name}.png"
                plt.savefig(plot_path, dpi=150, bbox_inches='tight')
                plt.show()
                display(widgets.HTML(f"<div style='color:green;font-size:11px;margin-top:-10px;margin-bottom:20px;'>📊 <b>Kinetics saved to:</b> {plot_path.relative_to(Path(self.output_dir))}</div>"))

            # 2. Combined Killing Efficiency Distribution (ONLY that plot)
            df_active = df_killing[df_killing["is_active_killing"]].copy()
            if not df_active.empty:
                display(widgets.HTML("<hr><h3 style='margin-bottom:5px;'>Combined Samples</h3>"))
                fig, ax = plt.subplots(figsize=(10, 6))
                sns.histplot(df_active["killing_efficiency"], kde=True, ax=ax, color='purple')
                ax.set_title("Combined Killing Efficiency Distribution", fontsize=16, fontweight='bold')
                ax.set_xlabel("Efficiency Score (signal increase / expected background)")
                ax.set_ylabel("Active Killing Events")
                
                # Save combined plot
                combined_plot_dir = Path(self.output_dir, "analysis", self.immune_dd.value, "active_killing", "plots")
                combined_plot_dir.mkdir(parents=True, exist_ok=True)
                combined_path = combined_plot_dir / "combined_killing_efficiency_distribution.png"
                plt.savefig(combined_path, dpi=150, bbox_inches='tight')
                plt.show()
                display(widgets.HTML(f"<div style='color:green;font-size:11px;margin-top:-10px;'>📊 <b>Combined distribution saved to:</b> {combined_path.relative_to(Path(self.output_dir))}</div>"))

            # Summary Table
            n_items = int(self.gallery_item_count.value)
            top_hitters = df_killing[df_killing["is_active_killing"]].groupby(["sample_name", "immune_track_id"]).size().sort_values(ascending=False).head(n_items)
            if not top_hitters.empty:
                display(widgets.HTML("<br><b>Top Killing T cells across all samples (by events):</b>"))
                display(top_hitters.to_frame("Count"))

    def _display_gallery(self, df_killing):
        """Display the Active Killing Events Gallery organized by sample."""
        n_items_per_sample = int(self.gallery_item_count.value)
        
        display(widgets.HTML(f'<hr><b style="font-size:18px;">Active Killing Events Gallery (Top {n_items_per_sample} per sample)</b>'))
        display(widgets.HTML('<p style="font-size:12px;color:#666;">Showing Maximum Intensity Projections (MIP) using raw image data.<br>'
                             '<b>Legend:</b> <span style="color:#666;">Grayscale (Raw Structure)</span>, '
                             '<span style="color:red;">Translucent Red (Dead Mask)</span>, '
                             '<span style="color:purple;">Translucent Purple (ACTIVE T-cell)</span></p>'))
        
        df_active = df_killing[df_killing["is_active_killing"]].copy()
        if not df_active.empty:
            samples = df_active["sample_name"].unique()
            
            for sample_name in samples:
                display(widgets.HTML(f"<h3 style='margin-top:20px;background:#f0faf0;padding:5px;border-left:5px solid #28a745;'>Sample: {sample_name}</h3>"))
                
                df_sample = df_active[df_active["sample_name"] == sample_name].copy()
                all_killers = df_sample.groupby("immune_track_id").size().sort_values(ascending=False)
                max_available = len(all_killers)
                
                current_n = min(n_items_per_sample, max_available)
                if n_items_per_sample > max_available:
                     display(widgets.HTML(f'<small style="color:orange;">Only {max_available} active killers found in this sample.</small>'))
                
                top_hitters_idx = all_killers.head(current_n).index
                
                gallery_items = []
                for tid in top_hitters_idx:
                    df_hitter = df_sample[df_sample["immune_track_id"] == tid]
                    best_event = df_hitter.sort_values("killing_efficiency", ascending=False).iloc[0]
                    
                    display(widgets.HTML(f"<i>📸 Processing T-cell {tid}...</i>"))
                    item = self._generate_killing_gallery_item(best_event, df_killing)
                    if item: gallery_items.append(item)
                
                if gallery_items:
                    # No easy way to clear previous "Processing" messages without clearing everything, 
                    # but we can wrap them in a VBox to show them all at once at the end of sample processing.
                    # Or just display them. Let's just display them.
                    display(widgets.VBox(gallery_items))
                else:
                    display(widgets.HTML("<i>Could not generate gallery items for this sample.</i>"))
        else:
            display(widgets.HTML("<i>No active killing events found to display in gallery.</i>"))

    def _generate_killing_gallery_item(self, event_row, df_killing=None):
        """Generate a widget representing one killing event with images."""
        try:
            # Redundant imports removed as they are at the top
            sample_name = event_row["sample_name"]
            t_id = int(event_row["immune_track_id"])
            o_id = int(event_row["targeted_track_id"]) if pd.notna(event_row["targeted_track_id"]) else -1
            t_start = int(event_row["position_t"])
            observation_window = int(self.observation_window.value)
            t_end = t_start + (2 * observation_window) # Extend window for context
            
            # Metadata lookup
            md_match = self.metadata_loader.metadata[self.metadata_loader.metadata["sample_name"] == sample_name]
            if md_match.empty:
                print(f"  ❌ Sample '{sample_name}' not found in metadata.")
                return None
            md_row = md_match.iloc[0]
            
            # Paths - Robust column lookup (handles im_ and or_ prefixes)
            def get_col_path(suffix):
                cols = [c for c in self.metadata_loader.metadata.columns if c.endswith(suffix)]
                if not cols: return None
                val = md_row[cols[0]]
                if pd.isna(val) or str(val).strip() == "": return None
                return Path(val)

            immune_tracks_path = get_col_path(f"{self.immune_dd.value}_tracks_image_path")
            if immune_tracks_path is None:
                print(f"  ❌ Metadata column not found for: {self.immune_dd.value}_tracks_image_path")
                return None
            if not immune_tracks_path.exists():
                print(f"  ❌ Immune tracks file not found at: {immune_tracks_path}")
                return None
            
            # Auto-detect organoid type
            org_type = "organoid"
            for ct in self.metadata_loader.metadata.columns:
                if "_tracks_image_path" in ct and f"im_{self.immune_dd.value}" not in ct and f"_{self.immune_dd.value}" not in ct:
                    org_type = ct.replace("_tracks_image_path", "").replace("or_", "")
                    break
            
            raw_img_path = get_col_path("raw_image_path")
            if raw_img_path is None:
                 raw_img_path = get_col_path("image_path") # Typical fallback

            org_tracks_path = get_col_path(f"{org_type}_tracks_image_path")
            if org_tracks_path is None:
                print(f"  ❌ Metadata column not found for: {org_type}_tracks_image_path")
                return None
            if not org_tracks_path.exists():
                print(f"  ❌ Target tracks file not found at: {org_tracks_path}")
                return None
                
            dead_mask_path = Path(self.output_dir, "images", sample_name, f"{sample_name}_mask_dead.zarr")
            
            # Helper to get MIP crop
            def get_mip_crop(t, track_id, organoid_id, center_coords=None):
                try:
                    res_xy = float(md_row.get("pixel_distance_xy", 1.0))
                    res_z = float(md_row.get("pixel_distance_z", 1.0))
                    
                    # Determine coordinates: Use provided center or look up from data
                    if center_coords is not None:
                        cz, cy, cx = center_coords
                    else:
                        # Try to get coordinates from provided dataframe first
                        row_t = pd.DataFrame()
                        if df_killing is not None:
                            row_t = df_killing[(df_killing["sample_name"] == sample_name) & 
                                               (df_killing["immune_track_id"] == track_id) & 
                                               (df_killing["position_t"] == t)]
                        if row_t.empty:
                            # Fallback to self.killing_results
                            if hasattr(self, "killing_results") and self.killing_results is not None:
                                 row_t = self.killing_results[(self.killing_results["sample_name"] == sample_name) & 
                                                        (self.killing_results["immune_track_id"] == track_id) & 
                                                        (self.killing_results["position_t"] == t)]
                        
                        if row_t.empty: return None
                        cz = int(row_t["position_z"].iloc[0] / res_z)
                        cy = int(row_t["position_y"].iloc[0] / res_xy)
                        cx = int(row_t["position_x"].iloc[0] / res_xy)
                    
                    win = 60 # Slightly larger window for raw data
                    
                    def load_tp(path, timepoint):
                        if not path or not path.exists(): return None
                        try:
                            if str(path).endswith(".zarr") or str(path).endswith(".zarr.zip"):
                                img = load_zarr(path)
                                if timepoint >= img.shape[0]: return None
                                return img[timepoint]
                            img = load_image(path)
                            if timepoint >= img.shape[0]: return None
                            return img[timepoint]
                        except: return None

                    def get_crop(img):
                        if img is None: return None
                        try:
                            if img.ndim == 4: # C,Z,Y,X
                                c, z, y, x = img.shape
                                y1, y2 = max(0, cy-win), min(y, cy+win)
                                x1, x2 = max(0, cx-win), min(x, cx+win)
                                return img[:, :, y1:y2, x1:x2]
                            else: # Z,Y,X
                                z, y, x = img.shape
                                y1, y2 = max(0, cy-win), min(y, cy+win)
                                x1, x2 = max(0, cx-win), min(x, cx+win)
                                return img[:, y1:y2, x1:x2]
                        except: return None

                    # Load raw data and masks
                    raw_tp = load_tp(raw_img_path, t)
                    m_t_full = load_tp(immune_tracks_path, t)
                    m_d_full = load_tp(dead_mask_path, t)

                    c_raw = get_crop(raw_tp)
                    c_mt = get_crop(m_t_full)
                    c_md = get_crop(m_d_full)

                    if c_raw is None: return None
                    
                    # Ensure numpy arrays for computation
                    if c_mt is not None:
                        if hasattr(c_mt, 'compute'): c_mt = c_mt.compute()
                        c_mt = np.asarray(c_mt)
                    if c_raw is not None:
                        if hasattr(c_raw, 'compute'): c_raw = c_raw.compute()
                        c_raw = np.asarray(c_raw)
                    if c_md is not None:
                        if hasattr(c_md, 'compute'): c_md = c_md.compute()
                        c_md = np.asarray(c_md)

                    # 1. Grayscale Base (all channels)
                    if c_raw is not None and c_raw.ndim == 4:
                        raw_gray = np.max(c_raw, axis=0) # Z,Y,X
                        mip_gray = np.max(raw_gray, axis=0).astype(np.float32)
                        p1, p99 = np.percentile(mip_gray, [1, 99])
                        mip_gray = np.clip((mip_gray - p1) / (p99 - p1 + 1e-10), 0, 1)
                    else:
                        mip_gray = np.zeros(c_raw.shape[1:] if c_raw is not None else (120, 120), dtype=np.float32)

                    # 2. Tracks & Color Blending
                    mip_dead = np.max(c_md > 0, axis=0) if c_md is not None else np.zeros_like(mip_gray)
                    mip_active = np.max(c_mt == track_id, axis=0) if c_mt is not None else np.zeros_like(mip_gray)

                    rgb = np.zeros((mip_gray.shape[0], mip_gray.shape[1], 3), dtype=np.float32)
                    rgb[:, :, 0] = mip_gray; rgb[:, :, 1] = mip_gray; rgb[:, :, 2] = mip_gray
                    
                    alpha_red = 0.3
                    alpha_purple = 0.5
                    
                    rgb[mip_dead > 0, 0] = rgb[mip_dead > 0, 0] * (1-alpha_red) + 1.0 * alpha_red
                    rgb[mip_dead > 0, 1] = rgb[mip_dead > 0, 1] * (1-alpha_red)
                    rgb[mip_dead > 0, 2] = rgb[mip_dead > 0, 2] * (1-alpha_red)
                    
                    rgb[mip_active > 0, 0] = rgb[mip_active > 0, 0] * (1-alpha_purple) + 1.0 * alpha_purple
                    rgb[mip_active > 0, 1] = rgb[mip_active > 0, 1] * (1-alpha_purple)
                    rgb[mip_active > 0, 2] = rgb[mip_active > 0, 2] * (1-alpha_purple) + 1.0 * alpha_purple
                    
                    return np.clip(rgb * 255, 0, 255).astype(np.uint8)
                except: return None

            # Max time calculation
            if "timelapse_duration" in md_row:
                max_t = int(md_row["timelapse_duration"]) - 1
            elif df_killing is not None:
                max_t = int(df_killing[df_killing["sample_name"] == sample_name]["position_t"].max())
            else:
                max_t = t_end  # Fallback
            
            # Determine fixed center coordinates at event start
            res_xy = float(md_row.get("pixel_distance_xy", 1.0))
            res_z = float(md_row.get("pixel_distance_z", 1.0))
            fixed_coords = (
                int(event_row["position_z"] / res_z),
                int(event_row["position_y"] / res_xy),
                int(event_row["position_x"] / res_xy)
            )

            frames = []
            for t in range(t_start, min(t_end + 1, max_t + 1)):
                img = get_mip_crop(t, t_id, o_id, center_coords=fixed_coords)
                if img is not None:
                    pil_img = PILImage.fromarray(img)
                    draw = ImageDraw.Draw(pil_img)
                    text = f"T={t}"
                    # Simple text overlay at top-right
                    # Using default font and drawing a small shadow for visibility
                    w, h = pil_img.size
                    pos = (w - 50, 5) 
                    draw.text((pos[0]+1, pos[1]+1), text, fill="black") # shadow
                    draw.text(pos, text, fill="white")
                    frames.append(pil_img)
            
            if not frames:
                return None
            
            # Save GIF to output directory
            gallery_dir = Path(self.output_dir, "analysis", self.immune_dd.value, "active_killing", "gallery", sample_name)
            gallery_dir.mkdir(parents=True, exist_ok=True)
            gif_name = f"killing_event_{sample_name}_T{t_id}_start{t_start}.gif"
            gif_path = gallery_dir / gif_name
            
            # Save with Pillow
            frames[0].save(
                gif_path,
                save_all=True,
                append_images=frames[1:],
                duration=200, # 200ms per frame
                loop=0
            )
            
            with open(gif_path, "rb") as f:
                gif_data = f.read()
                
            img_widget = widgets.Image(value=gif_data, format='gif')
            
            killing_window_end = min(t_start + observation_window, max_t)
            viz_window_end = min(t_end, max_t)
            
            info_html = widgets.HTML(f"""
                <div style='border: 1px solid #ddd; padding: 10px; margin-bottom: 5px; background: #f9f9f9; width: 600px;'>
                    <b>Event Details:</b> Sample: {sample_name} | T-cell ID: {t_id} | Target ID: {o_id if o_id != -1 else 'N/A'}<br>
                    <b>Efficiency:</b> {event_row['killing_efficiency']:.2f}<br>
                    <div style='font-size: 11px; color: #444; margin-top: 5px;'>
                        <b>🧪 Killing Window:</b> T={t_start} to T={killing_window_end} (Analysis period)<br>
                        <b>🎬 Visualization Window:</b> T={t_start} to T={viz_window_end} (Forensic context)
                    </div>
                    <small style='color: #666;'>Saved to: {gif_path.relative_to(Path(self.output_dir))}</small>
                </div>
            """)
            return widgets.VBox([info_html, img_widget])
            
        except Exception as e:
            with self.gallery_out:
                print(f"  ❌ Error generating gallery item for {sample_name}:")
                traceback.print_exc()
            return None


class DeathThresholdPreview:
    """
    Interactive visual preview for tuning the organoid death threshold.
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = Path(self.metadata_loader.output_dir).expanduser()
        
        # Paths
        feature_outdir = self.output_dir / "analysis" / self.cell_type / "track_features"
        #self.df_tracks_path = feature_outdir / f"BEHAV3D_{self.cell_type}_combined_track_features.csv"
        self.df_tracks_path = feature_outdir / f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv"

        #if not self.df_tracks_path.exists():
            #self.df_tracks_path = feature_outdir / f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv"

        if not self.df_tracks_path.exists():
            self.df_tracks = None
        else:
            self.df_tracks = pd.read_csv(self.df_tracks_path)
            self.df_tracks['TrackID'] = self.df_tracks['TrackID'].astype(str)

        self.samples = self.metadata_loader.metadata['sample_name'].tolist()
        
        # UI Elements
        self.sample_dd = widgets.Dropdown(options=self.samples, description="Sample:", layout={'width': '250px'})
        self.time_slider = widgets.IntSlider(description="Time:", continuous_update=False, layout={'width': '400px'})
        #self.threshold_slider = widgets.FloatSlider(value=0.02, min=0, max=0.5, step=0.005, description="Threshold:", continuous_update=False, layout={'width': '400px'})
        self.threshold_input = widgets.BoundedFloatText(
            value=0.02, min=0, max=1, step=0.0005,
            description="Threshold:", style={'description_width': 'initial'},
            layout={'width': '120px'}
            )
        self.org_ch_input = widgets.IntText(value=0, description=f"{self.cell_type} Channel:", style={'description_width': 'initial'}, layout={'width': '200px'})
        self.dead_ch_input = widgets.IntText(value=1, description="Dead Channel:", style={'description_width': 'initial'}, layout={'width': '200px'})

        self.out = widgets.Output(layout={'overflow': 'visible', 'height': 'auto'})
        
        # Observers
        self.sample_dd.observe(self._on_sample_change, names='value')
        self.time_slider.observe(self._update, names='value')
        self.threshold_input.observe(self._update, names='value')
        self.org_ch_input.observe(self._update, names='value')
        self.dead_ch_input.observe(self._update, names='value')

        self.ui = widgets.VBox([
            widgets.HTML(f"<b>Visual Threshold Preview</b> <i>(Max Intensity Projection)</i>"),
            widgets.HBox([self.sample_dd, self.threshold_input]),
            widgets.HBox([
                widgets.HTML("<b>Channels: </b>"),
                self.org_ch_input,
                self.dead_ch_input,
                widgets.HTML("<i>(Indices starting at 0)</i>")
            ]),
            self.time_slider,
            self.out,
        ])

        # Initial Load
        if self.samples:
            self._on_sample_change({'new': self.samples[0]})

    def _on_sample_change(self, change):
        sample = change['new']
        row = self.metadata_loader.metadata[self.metadata_loader.metadata['sample_name'] == sample].iloc[0]
        
        sample_dir = self.output_dir / "images" / sample
        self.raw_path = sample_dir / f"{sample}.zarr"
        
        # Segmentation/tracked image path from metadata
        self.label_path = None
        for prefix in ['or', 'im', 'ot']:
            col = f"{prefix}_{self.cell_type}_tracks_image_path"
            if col in row and pd.notna(row[col]) and str(row[col]).strip():
                p = Path(row[col])
                if p.exists():
                    self.label_path = p
                    break
        
        # Fallback:tracked.zarr
        if self.label_path is None or not self.label_path.exists():
            tracked_fallback = sample_dir / f"{sample}_{self.cell_type}_tracked.zarr"
            #segments_fallback = sample_dir / f"{sample}_{self.cell_type}_segments.zarr"
            if tracked_fallback.exists():
                self.label_path = tracked_fallback
            #elif segments_fallback.exists():
            #    self.label_path = segments_fallback
            else:
                self.label_path = tracked_fallback  # will fail with informative error
            
        self.dead_mask_path = sample_dir / f"{sample}_mask_dead.zarr"
        if 'dead_mask_path' in row and pd.notna(row['dead_mask_path']):
            dp = Path(row['dead_mask_path'])
            if dp.exists():
                self.dead_mask_path = dp
        
        try:
            raw = load_zarr(self.raw_path)
            # T, C, Z, Y, X or T, Z, Y, X
            self.time_slider.max = raw.shape[0] - 1
            self.time_slider.value = raw.shape[0] // 2
        except: pass
        
        self._update(None)

    @staticmethod
    def _contrast_limits(img):
        """Compute display vmin/vmax using 1st–99th percentile of non-zero pixels."""
        pos = img[img > 0]
        if pos.size == 0:
            return 0, max(1, img.max())
        return float(np.percentile(pos, 1)), float(np.percentile(pos, 99))

    def _update(self, _):
        sample = self.sample_dd.value
        t = self.time_slider.value
        thr = self.threshold_input.value
        
        self.out.clear_output(wait=False)
        with self.out:
            try:
                org_ch = self.org_ch_input.value
                dead_ch = self.dead_ch_input.value

                # Load Z stacks and apply Maximum Intensity Projection per channel
                raw_zarr = load_zarr(self.raw_path)
                has_channels = (raw_zarr.ndim == 5)
                if has_channels:
                    n_channels = raw_zarr.shape[1]
                    if org_ch >= n_channels:
                        print(f"{self.cell_type} Channel {org_ch} out of range (0–{n_channels-1}). Select a valid channel.")
                        return
                    if dead_ch >= n_channels:
                        print(f"Dead Channel {dead_ch} out of range (0–{n_channels-1}). Select a valid channel.")
                        return
                    raw_img = np.asarray(raw_zarr[t, org_ch]).max(axis=0)
                    dead_raw = np.asarray(raw_zarr[t, dead_ch]).max(axis=0)
                else:
                    raw_img = np.asarray(raw_zarr[t]).max(axis=0)
                    dead_raw = raw_img

                # Segmentation: first-nonzero projection (avoids label mixing from MIP)
                label_zarr = load_zarr(self.label_path)
                label_vol = np.asarray(label_zarr[t])
                if label_vol.ndim == 3:
                    first_nz_idx = np.argmax(label_vol > 0, axis=0)
                    label_img = np.take_along_axis(label_vol, first_nz_idx[np.newaxis], axis=0)[0]
                else:
                    label_img = label_vol

                # Dead mask: MIP over Z (binary OR: correct for 0/1 masks)
                dead_mask_zarr = load_zarr(self.dead_mask_path)
                dead_vol = np.asarray(dead_mask_zarr[t])
                dead_mask = dead_vol.max(axis=0) if dead_vol.ndim == 3 else dead_vol

                # Classification image (red=dead, green=alive, gray=unknown)
                class_img = np.zeros((*label_img.shape, 3), dtype=float)
                if self.df_tracks is not None:
                    df_t = self.df_tracks[(self.df_tracks['sample_name'] == sample) & (self.df_tracks['position_t'] == t)]
                    for label_id in np.unique(label_img):
                        if label_id == 0: continue
                        row = df_t[df_t['TrackID'] == str(label_id)]
                        if not row.empty:
                            perc = row['percentage_dead_mask'].values[0]
                            color = [1, 0, 0] if perc >= thr else [0, 1, 0]
                            class_img[label_img == label_id] = color
                        else:
                            class_img[label_img == label_id] = [0.5, 0.5, 0.5]

                # Contrast limits
                vmin_r, vmax_r = self._contrast_limits(raw_img)
                vmin_d, vmax_d = self._contrast_limits(dead_raw)

                # Masks for overlay
                seg_display = np.ma.masked_where(label_img == 0, label_img)
                dead_mask_overlay = np.ma.masked_where(dead_mask == 0, dead_mask)
                red_cmap = ListedColormap(['#FF2020'])

                #fig, axes = plt.subplots(1, 5, figsize=(25, 5))
                fig, axes = plt.subplots(1, 5, figsize=(30, 6), dpi=150)
                # 1) Raw Organoid
                axes[0].imshow(raw_img, cmap='gray', vmin=vmin_r, vmax=vmax_r)
                axes[0].set_title(f"Raw {self.cell_type} (Ch {org_ch})")
                # 2) Segmentation
                axes[1].imshow(np.zeros_like(label_img), cmap='gray', vmin=0, vmax=1)
                axes[1].imshow(seg_display, cmap='nipy_spectral', interpolation='nearest')
                axes[1].set_title(f"Segments ({Path(self.label_path).stem})")
                # 3) Raw Dead channel
                axes[2].imshow(dead_raw, cmap='gray', vmin=vmin_d, vmax=vmax_d)
                axes[2].set_title(f"Raw Dead (Ch {dead_ch})")
                # 4) Dead mask (red) on raw organoid
                axes[3].imshow(raw_img, cmap='gray', vmin=vmin_r, vmax=vmax_r)
                axes[3].imshow(dead_mask_overlay, cmap=red_cmap, alpha=0.55)
                axes[3].set_title(f"Dead Mask on {self.cell_type}")
                # 5) Classification
                axes[4].imshow(class_img)
                axes[4].set_title(f"Classification (Thr: {thr:.3f})\nRed=Dead, Green=Alive")
                for ax in axes: ax.axis('off')
                plt.tight_layout()
                plt.show()
                
                ch_note = f"Raw: {raw_zarr.shape} ({n_channels} ch)" if has_channels else f"Raw: {raw_zarr.shape} (no channel dim)"
                info_html = f"""
                <div style='padding: 5px; border: 1px solid #ccc; background-color: #f9f9f9; font-size: 12px;'>
                    <b>Legend:</b>
                    <span style='color: green;'>●</span> Alive |
                    <span style='color: red;'>●</span> Dead |
                    <span style='color: gray;'>●</span> Track missing/filtered |
                </div>
                """
                display(widgets.HTML(info_html))
                
                if self.df_tracks is None:
                    print("Please run Feature Extraction and track filtering first.")
            except Exception as e:
                import traceback; traceback.print_exc()
                print(f"Preview unavailable: {e}")

class DeathDynamicsPanel:
    """
    Death dynamics analysis panel (organoid-specific).
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        feature_outdir = Path(self.output_dir, "analysis", self.cell_type, "track_features")
        p = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features.csv")
        self.has_death_features = False
        if p.exists():
            try: self.has_death_features = bool({'mean_dead_dye', 'percentage_dead_mask', 'nr_dead_mask_pixels'}.intersection(set(pd.read_csv(p, nrows=0).columns)))
            except Exception: pass
        
        params = dict(self.metadata_loader.behav3d_parameters or {})
        cfg = params.setdefault("death_dynamics", {}).setdefault(self.cell_type, deepcopy(_DEFAULT_CONFIG.get("death_dynamics", {}).get("organoid", {})))
        self._params = params
        self._panel_cfg = cfg
        self.metadata_loader.behav3d_parameters = self._params
        
        self.dead_perc_threshold = widgets.FloatText(description="Dead % threshold", value=float(cfg.get("dead_perc_threshold", 0.02)), style={'description_width': '160px'}, layout=widgets.Layout(width="220px"))
        self.btn_run = widgets.Button(description=f"Run {cell_type} death dynamics", button_style="warning", layout=widgets.Layout(width="300px"))
        self.btn_run.on_click(self._on_run_clicked)
        
        self.btn_preview = widgets.ToggleButton(description="Preview Threshold", button_style="info", layout=widgets.Layout(width="200px"))
        self.btn_preview.observe(self._on_preview_toggled, names='value')
        self.preview_container = widgets.VBox([])

        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()
        
        if not self.has_death_features:
            self.ui = widgets.VBox([widgets.HTML(f'<b>{self.cell_type} Death Dynamics</b>'), widgets.HTML('<div style="color:#b00;">No death features found. Run feature extraction first.</div>')])
        else:
            self.ui = widgets.VBox([
                widgets.HTML(f'<b>{self.cell_type} Death Dynamics</b>'), 
                self.btn_preview,
                self.preview_container,
                widgets.HTML("<hr>"),
                widgets.HBox([self.dead_perc_threshold, widgets.HTML('<div style="font-size:12px;color:#666;">(Adjust value manually or use preview slider)</div>')], layout=widgets.Layout(align_items="center")),
                widgets.HBox([self.btn_run, self.spinner_html]), 
                self.out
            ])

    def _on_preview_toggled(self, change):
        if change['new']:
            try:
                preview = DeathThresholdPreview(self.metadata_loader, self.cell_type)
                widgets.jslink((self.dead_perc_threshold, 'value'), (preview.threshold_input, 'value'))
                self.preview_container.children = [preview.ui]
            except Exception as e:
                with self.out: print(f"Error launching preview: {e}")
        else:
            self.preview_container.children = []

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True; self.spinner_html.layout.display = None; self.out.clear_output()
        with self.out:
            try:
                self._panel_cfg["dead_perc_threshold"] = float(self.dead_perc_threshold.value)
                self.metadata_loader.behav3d_parameters = self._params
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f: yaml.safe_dump(self._params, f, sort_keys=False)
                run_organoid_analysis(dead_perc_threshold=float(self.dead_perc_threshold.value), output_dir=self.output_dir, df_tracks_path=None, org_type=self.cell_type, metadata=self.metadata_loader.metadata)
                print(f"✅ {self.cell_type} death dynamics complete!")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False


class MultiOrganoidDeathDynamicsPanel:
    """
    Combined death dynamics comparison panel for multiple organoid types.
    Shows all organoid types together on the same plot for easy comparison.
    Only appears when 2+ organoid types are detected.
    """
    def __init__(self, metadata_loader, organoid_types):
        self.metadata_loader = metadata_loader
        self.organoid_types = organoid_types  # List of organoid type names
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        # Check which organoid types have death dynamics data available
        self.available_data = {}
        for org_type in self.organoid_types:
            csv_path = Path(self.output_dir, "analysis", org_type, "results", f"combined_general_{org_type}_dynamics_analysis.csv")
            if csv_path.exists():
                self.available_data[org_type] = csv_path
        
        self.btn_refresh = widgets.Button(description="🔄 Refresh", button_style="info", layout=widgets.Layout(width="100px"))
        self.btn_refresh.on_click(self._on_refresh_clicked)
        self.btn_run = widgets.Button(description="Generate Comparison Plot", button_style="warning", layout=widgets.Layout(width="250px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.status_html = widgets.HTML("")
        self.out = widgets.Output()
        
        self._refresh_status()
        
        self.ui = widgets.VBox([
            widgets.HTML('<b>Multi-Organoid Death Dynamics Comparison</b>'),
            widgets.HTML('<div style="font-size:12px;color:#666;">Compare death dynamics across all organoid types on the same plot.</div>'),
            widgets.HBox([self.status_html, self.btn_refresh], layout=widgets.Layout(align_items="center", gap="10px")),
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out
        ])
    
    def _refresh_status(self):
        """Check which organoid types have death dynamics data available."""
        self.available_data = {}
        for org_type in self.organoid_types:
            csv_path = Path(self.output_dir, "analysis", org_type, "results", f"combined_general_{org_type}_dynamics_analysis.csv")
            if csv_path.exists():
                self.available_data[org_type] = csv_path
        
        if len(self.available_data) < 2:
            missing = [ot for ot in self.organoid_types if ot not in self.available_data]
            self.status_html.value = f'<div style="color:#b00;">Waiting for death dynamics data from: {", ".join(missing)}</div>'
            self.btn_run.disabled = True
        else:
            self.status_html.value = f'<div style="color:#080;">Ready: {", ".join(self.available_data.keys())}</div>'
            self.btn_run.disabled = False
    
    def _on_refresh_clicked(self, *_):
        self._refresh_status()
    
    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()
        
        with self.out:
            try:
                # Call the analysis function from organoid_analysis module
                plot_multi_organoid_death_dynamics(
                    output_dir=self.output_dir,
                    organoid_types=self.organoid_types
                )
                print("✅ Multi-organoid comparison complete!")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False


class MultiOrganoidMorphologyDeathPanel:
    """
    Morphology vs Death analysis for all organoid types combined.
    """
    def __init__(self, metadata_loader, organoid_types):
        self.metadata_loader = metadata_loader
        self.organoid_types = organoid_types
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        params = dict(self.metadata_loader.behav3d_parameters or {})
        cfg = params.setdefault("morpho_dead", {})

        self.window_start = widgets.IntText(
            description="Window Start (t >= N)",
            value=int(cfg.get("window_start", 0)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.window_end = widgets.IntText(
            description="Window End (t <= N)",
            value=int(cfg.get("window_end", 5)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.umap_n_neighbors = widgets.IntText(
            description="UMAP neighbors",
            value=int(cfg.get("umap_n_neighbors", 15)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.umap_n_neighbors.tooltip = "Number of local neighbors used for UMAP manifold construction. Controls local vs global structure."
        
        self.umap_min_dist = widgets.FloatText(
            description="UMAP min_dist",
            value=float(cfg.get("umap_min_dist", 0.1)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.umap_min_dist.tooltip = "Minimum distance between points in UMAP embedding. Controls how tightly points are packed."
        
        self.umap_spread = widgets.FloatText(
            description="UMAP spread",
            value=float(cfg.get("umap_spread", 1.0)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.umap_spread.tooltip = "The effective scale of embedded points."
        
        self.distance_metric = widgets.Text(
            value=cfg.get("distance_metric", "euclidean"),
            description="UMAP metric",
            placeholder="euclidean",
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.distance_metric.tooltip = "Distance metric for UMAP. See umap-learn docs for full options."
                
        self.leiden_resolution = widgets.FloatText(
            description="Leiden resolution",
            value=float(cfg.get("leiden_resolution", 1.0)),
            style={'description_width': '180px'},
            layout=widgets.Layout(width="260px")
        )
        self.leiden_resolution.tooltip = "Controls cluster granularity. Higher values lead to more clusters."

        self.status_html = widgets.HTML("")
        self.btn_refresh = widgets.Button(description="Refresh", button_style="info", layout=widgets.Layout(width="100px"))
        self.btn_refresh.on_click(self._on_refresh_clicked)
        self.btn_run = widgets.Button(description="Run Initial Morphology vs Death", button_style="warning", layout=widgets.Layout(width="280px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()

        self._refresh_status()

        self.ui = widgets.VBox([
            widgets.HTML('<b>Initial Morphology Related to Death</b>'),
            widgets.HTML('<div style="font-size:12px;color:#666;">Uses baseline morphology averaged over the selected window, excludes organoids dead in that window, and labels dead if any death occurs after.</div>'),
            widgets.HTML('<div style="font-size:12px;color:#666;">Features: 4 key morphology descriptors (volume, sphericity, solidity, surface-to-volume ratio).</div>'),
            widgets.HBox([self.status_html, self.btn_refresh], layout=widgets.Layout(align_items="center", gap="10px")),
            widgets.HTML("<hr>"),
            widgets.HBox([self.window_start, self.window_end], layout=widgets.Layout(gap="12px")),
            widgets.HBox([self.umap_n_neighbors, self.umap_min_dist], layout=widgets.Layout(gap="12px")),
            widgets.HBox([self.distance_metric, self.leiden_resolution], layout=widgets.Layout(gap="12px")),
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html]),
            self.out
        ])

        # Check for organoid presence to show/hide the entire panel
        md = self.metadata_loader.metadata
        has_organoids = any(col.startswith('or_') for col in md.columns)
        if not has_organoids:
            self.ui.layout.display = "none"

    def _refresh_status(self):
        self.available_data = {}
        missing = []
        for org_type in self.organoid_types:
            csv_path = Path(self.output_dir, "analysis", org_type, "track_features", f"BEHAV3D_{org_type}_combined_track_features_filtered.csv")
            if csv_path.exists():
                self.available_data[org_type] = csv_path
            else:
                missing.append(org_type)

        if not self.available_data:
            self.status_html.value = '<div style="color:#b00;">No organoid track features found. Run track filtering first.</div>'
            self.btn_run.disabled = True
        elif missing:
            self.status_html.value = (
                f'<div style="color:#b80;">Missing data for: {", ".join(missing)}. '
                f'Will run with available: {", ".join(self.available_data.keys())}</div>'
            )
            self.btn_run.disabled = False
        else:
            self.status_html.value = f'<div style="color:#080;">Ready: {", ".join(self.available_data.keys())}</div>'
            self.btn_run.disabled = False

    def _on_refresh_clicked(self, *_):
        self._refresh_status()

    def _persist_params(self):
        params = self.metadata_loader.behav3d_parameters
        cfg = params.setdefault("morpho_dead", {})
        cfg.update({
            "window_start": int(self.window_start.value),
            "window_end": int(self.window_end.value),
            "umap_n_neighbors": int(self.umap_n_neighbors.value),
            "umap_min_dist": float(self.umap_min_dist.value),
            #"umap_spread": float(self.umap_spread.value),
            "distance_metric": str(self.distance_metric.value).strip() or "euclidean",
            "leiden_resolution": float(self.leiden_resolution.value),
        })
        if hasattr(self.metadata_loader, "behav3d_parameters_path"):
            yaml.safe_dump(
                params,
                self.metadata_loader.behav3d_parameters_path.open("w"),
                sort_keys=False
            )

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True
        self.spinner_html.layout.display = None
        self.out.clear_output()
        with self.out:
            try:
                self._persist_params()
                run_organoid_morphology_dead_analysis(
                    output_dir=self.output_dir,
                    organoid_types=self.organoid_types,
                    initial_window=(int(self.window_start.value), int(self.window_end.value)),
                    umap_n_neighbors=int(self.umap_n_neighbors.value),
                    umap_min_dist=float(self.umap_min_dist.value),
                    #umap_spread=float(self.umap_spread.value),
                    distance_metric=str(self.distance_metric.value).strip() or "euclidean",
                    leiden_resolution=float(self.leiden_resolution.value),
                    metadata=self.metadata_loader.metadata,
                )
                print("Initial morphology vs death analysis complete!")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self.btn_run.disabled = False



class InteractionAnalysisPanel:
    """
    Interaction analysis panel for organoids.
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        from behav3d.core.metadata import detect_immune_cell_types_from_metadata, detect_other_cell_types_from_metadata
        md = self.metadata_loader.metadata
        self.available_types = detect_immune_cell_types_from_metadata(md) + detect_other_cell_types_from_metadata(md)
        self.df_tracks_path = Path(self.output_dir, "analysis", self.cell_type, "track_features", f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")
        
        self.status_html = widgets.HTML("")
        self.cell_type_box = widgets.VBox([])
        self.cell_type_checkboxes = {}
        self.btn_refresh = widgets.Button(description="🔄 Refresh", button_style="info", layout=widgets.Layout(width="100px"))
        self.btn_refresh.on_click(self._on_refresh_clicked)
        self.btn_run = widgets.Button(description="Run Interaction Analysis", button_style="warning", layout=widgets.Layout(width="250px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()
        
        self.ui = widgets.VBox([widgets.HTML(f'<b>{self.cell_type} Interaction Analysis</b>'), widgets.HBox([self.status_html, self.btn_refresh], layout=widgets.Layout(align_items="center", gap="10px")), widgets.HTML('<b>Select cell types:</b>'), self.cell_type_box, widgets.HTML("<hr>"), widgets.HBox([self.btn_run, self.spinner_html]), self.out])
        self._refresh_data_status()

    def _refresh_data_status(self):
        self.has_data = self.df_tracks_path.exists()
        contact_types = []
        if self.has_data:
            cols = pd.read_csv(self.df_tracks_path, nrows=0).columns
            contact_types = [ct for ct in self.available_types if f"{ct}_contact" in cols]
        
        if not self.has_data or not contact_types:
            self.status_html.value = '<div style="color:#b00;">⚠️ Waiting for data/contacts.</div>'
            self.btn_run.disabled = True
        else:
            self.status_html.value = f'<div style="color:#080;">✅ Ready: {", ".join(contact_types)}</div>'
            self.btn_run.disabled = False
            self.cell_type_checkboxes = {ct: widgets.Checkbox(value=True, description=ct, indent=False) for ct in contact_types}
            self.cell_type_box.children = list(self.cell_type_checkboxes.values())

    def _on_refresh_clicked(self, *_): self._refresh_data_status()

    def _on_run_clicked(self, *_):
        sel = [ct for ct, cb in self.cell_type_checkboxes.items() if cb.value]
        if not sel: return
        self.btn_run.disabled = True; self.spinner_html.layout.display = None; self.out.clear_output()
        with self.out:
            try:
                thr = _cfg_get(self.metadata_loader.behav3d_parameters, f"death_dynamics.{self.cell_type}.dead_perc_threshold", 0.02)
                run_interaction_analysis(output_dir=self.output_dir, cell_type=self.cell_type, interacting_cell_types=sel, dead_threshold=thr, df_tracks_path=str(self.df_tracks_path), show_plots=True)
                print("✅ Interaction Analysis complete!")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False

class MotileCellAnalysisPanel:
    """
    Generic behavioral analysis panel (DTW/UMAP/clustering) for ANY cell type.
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        self.category = detect_cell_type_category(self.cell_type, metadata_loader.metadata)
        
        params = dict(self.metadata_loader.behav3d_parameters or {})
        self._panel_cfg = params.setdefault("analysis", {}).setdefault(self.cell_type, deepcopy(_DEFAULT_CONFIG["analysis"].get(self.category, {})))
        self._params = params
        
        groups = deepcopy(behav3d_calculated_features)
        feature_outdir = Path(self.output_dir, "analysis", self.cell_type, "track_features")
        df_p = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")
        if not df_p.exists(): df_p = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features.csv")
        
        cols = []
        if df_p.exists():
            try: cols = pd.read_csv(df_p, nrows=0).columns.tolist()
            except Exception: pass
        
        # Filter groups based on actual columns
        if cols:
            expanded = {}
            for g, patterns in groups.items():
                matches = []
                for p in patterns:
                    if '*' in p or '?' in p: matches.extend(expand_column_patterns(p, cols))
                    elif p in cols: matches.append(p)
                if matches: expanded[g] = sorted(set(matches))
            groups = expanded

        self.seed_widget = widgets.IntText(description="Seed", value=int(self._panel_cfg.get("seed", 42)), style={"description_width": "80px"})
        
        self._group_rows = {}
        sel_boxes = []
        preset_set = set(self._panel_cfg.get("dtw_features_input", []))
        
        for g, feats in groups.items():
            child_cbs = [widgets.Checkbox(value=(f in preset_set), description=f, indent=True) for f in feats]
            gcb = widgets.Checkbox(value=any(cb.value for cb in child_cbs), indent=False)
            header = widgets.HBox([gcb, widgets.HTML(f"<b>{g}</b>")])
            grid = widgets.GridBox(child_cbs, layout=widgets.Layout(grid_template_columns="repeat(3, max-content)", grid_gap="4px 12px", margin="0 0 0 24px"))
            self._group_rows[g] = {"group_cb": gcb, "child_cbs": child_cbs, "container": widgets.VBox([header, grid])}
            sel_boxes.append(self._group_rows[g]["container"])
            
            def make_h(gn):
                def _h(change):
                    for cb in self._group_rows[gn]["child_cbs"]: cb.value = change["new"]
                return _h
            gcb.observe(make_h(g), names="value")

        self.umap_dist = widgets.FloatText(description="UMAP min_dist", value=float(self._panel_cfg.get("umap_min_dist", 0.1)))
        self.umap_neigh = widgets.IntText(description="UMAP n_neighbors", value=int(self._panel_cfg.get("umap_n_neighbors", 15)))
        self.clusters = widgets.IntText(description="# clusters", value=int(self._panel_cfg.get("nr_of_clusters", 5)))
        
        self.btn_run = widgets.Button(description=f"Run {cell_type} analysis", button_style="success", layout=widgets.Layout(width="300px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()
        
        self.ui = widgets.VBox([
            widgets.HTML(f'<div style="font-size:22px;font-weight:700;">{self.cell_type} behavioral analysis</div>'),
            self.seed_widget, widgets.HTML('<b>Select features for DTW:</b>'),
            widgets.VBox(sel_boxes),
            widgets.HBox([self.umap_dist, self.umap_neigh, self.clusters]),
            widgets.HBox([self.btn_run, self.spinner_html]), self.out
        ])

    def _on_run_clicked(self, *_):
        sel = [cb.description for r in self._group_rows.values() for cb in r["child_cbs"] if cb.value]
        if not sel: return
        self.btn_run.disabled = True; self.spinner_html.layout.display = None; self.out.clear_output()
        with self.out:
            try:
                self._panel_cfg.update({"seed": int(self.seed_widget.value), "umap_min_dist": float(self.umap_dist.value), "umap_n_neighbors": int(self.umap_neigh.value), "nr_of_clusters": int(self.clusters.value), "dtw_features_input": sel})
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f: yaml.safe_dump(self._params, f, sort_keys=False)
                run_tcell_analysis(cell_type=self.cell_type, output_dir=self.output_dir, df_tracks_path=str(Path(self.output_dir, "analysis", self.cell_type, "track_features", f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")), columns_to_use=sel, columns_to_normalize=sel, umap_minimal_distance=float(self.umap_dist.value), umap_n_neighbors=int(self.umap_neigh.value), nr_of_clusters=int(self.clusters.value), plot_results=True, seed=int(self.seed_widget.value))
                print(f"✅ {self.cell_type} analysis complete!")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False
