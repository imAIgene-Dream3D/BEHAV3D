import ipywidgets as widgets
from pathlib import Path
import os
import yaml
import traceback
import pandas as pd
from copy import deepcopy
from .utils import (
    _cfg_get, 
    spinning_loader,
    detect_cell_type_category,
    _DEFAULT_CONFIG
)
from behav3d.core.utils import expand_column_patterns
from behav3d.features.timepoint_features import run_feature_extraction
from behav3d.analysis.tcell_analysis import (
    filter_cell_tracks,
    run_tcell_analysis
)
from behav3d.analysis.organoid_analysis import (
    filter_organoid_tracks,
    run_organoid_analysis
)
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
        self.save_results = widgets.Checkbox(description="Save results to CSV", value=bool(self._cfg.get("save_results", True)), indent=False)
        
        self.btn_run = widgets.Button(description="Run Active Killing Analysis", button_style="danger", icon="bolt", layout=widgets.Layout(width="260px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()
        self.validation_html = widgets.HTML("")
        self._validate_inputs()
        self.immune_dd.observe(lambda _: self._validate_inputs(), names="value")
        
        self.ui = widgets.VBox([
            widgets.HTML('<div style="font-size:22px;font-weight:700;">Active Killing Analysis</div>'),
            widgets.HTML('<div style="color:#555;font-size:13px;margin-bottom:10px;">Detects functional killing events. <b>Targets:</b> Auto-detected organoids.</div>'),
            self.immune_dd, self.target_info_html, self.validation_html, widgets.HTML("<hr>"),
            widgets.HBox([self.observation_window, widgets.HTML('<span style="color:#666;font-size:12px;">timepoints after contact</span>'), widgets.HTML("&nbsp;&nbsp;&nbsp;"), self.killing_threshold, widgets.HTML('<span style="color:#666;font-size:12px;">× background rate</span>')], layout=widgets.Layout(align_items="center")),
            widgets.HBox([self.min_contact_duration, widgets.HTML('<span style="color:#666;font-size:12px;">timepoints</span>'), widgets.HTML("&nbsp;&nbsp;&nbsp;"), self.death_signal_dd], layout=widgets.Layout(align_items="center")),
            widgets.HTML("<hr>"),
            widgets.HBox([self.btn_run, self.spinner_html, self.save_results], layout=widgets.Layout(align_items="center", gap="15px")),
            self.out
        ])
    
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

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True; self.spinner_html.layout.display = None; self.out.clear_output()
        with self.out:
            try:
                self._cfg.update({"observation_window": int(self.observation_window.value), "death_signal_column": str(self.death_signal_dd.value), "killing_threshold_multiplier": float(self.killing_threshold.value), "min_contact_duration": int(self.min_contact_duration.value), "save_results": bool(self.save_results.value), "last_immune_cell": str(self.immune_dd.value)})
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f: yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
                run_active_killing_analysis(metadata=self.metadata_loader.metadata, output_dir=self.output_dir, immune_cell_type=str(self.immune_dd.value), target_cell_types=None, observation_window=int(self.observation_window.value), death_signal_column=str(self.death_signal_dd.value), min_contact_duration=int(self.min_contact_duration.value), killing_threshold_multiplier=float(self.killing_threshold.value), save_results=bool(self.save_results.value))
                print(f"✅ Killing analysis finished.")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False

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
        
        # Initialize death_dynamics structure in the actual parameters dictionary
        cfg = self.metadata_loader.behav3d_parameters.setdefault("death_dynamics", {}).setdefault(
            self.cell_type, 
            deepcopy(_DEFAULT_CONFIG.get("death_dynamics", {}).get("organoid", {}))
        )
        
        self.dead_perc_threshold = widgets.FloatText(description="Dead % threshold", value=float(cfg.get("dead_perc_threshold", 0.02)), style={'description_width': '160px'}, layout=widgets.Layout(width="220px"))
        self.btn_run = widgets.Button(description=f"Run {cell_type} death dynamics", button_style="warning", layout=widgets.Layout(width="300px"))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        self.out = widgets.Output()
        
        if not self.has_death_features:
            self.ui = widgets.VBox([widgets.HTML(f'<b>{self.cell_type} Death Dynamics</b>'), widgets.HTML('<div style="color:#b00;">⚠️ No death features found.</div>')])
        else:
            self.ui = widgets.VBox([widgets.HTML(f'<b>{self.cell_type} Death Dynamics</b>'), widgets.HBox([self.dead_perc_threshold, widgets.HTML('<div style="font-size:12px;color:#666;">(threshold for classification)</div>')], layout=widgets.Layout(align_items="center")), widgets.HTML("<hr>"), widgets.HBox([self.btn_run, self.spinner_html]), self.out])

    def _on_run_clicked(self, *_):
        self.btn_run.disabled = True; self.spinner_html.layout.display = None; self.out.clear_output()
        with self.out:
            try:
                self.metadata_loader.behav3d_parameters["death_dynamics"][self.cell_type]["dead_perc_threshold"] = float(self.dead_perc_threshold.value)
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f: yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
                run_organoid_analysis(dead_perc_threshold=float(self.dead_perc_threshold.value), output_dir=self.output_dir, df_tracks_path=None, org_type=self.cell_type, metadata=self.metadata_loader.metadata)
                print(f"✅ {self.cell_type} death dynamics complete!")
            except Exception: traceback.print_exc()
            finally: self.spinner_html.layout.display = "none"; self.btn_run.disabled = False

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
