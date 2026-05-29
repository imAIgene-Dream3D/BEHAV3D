import ipywidgets as widgets
from pathlib import Path
import pandas as pd
import traceback
import yaml
from copy import deepcopy
from .utils import (
    _cfg_get, 
    spinning_loader,
    detect_cell_type_category,
    ExternalImageImporter,
)
from behav3d.preprocessing.tracking.laptracking import run_laptracking
from behav3d.preprocessing.tracking.trackpy_tracking import run_trackpy_tracking_generic
from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking
from behav3d.preprocessing.tracking.btrack_tracking import run_btracking
from behav3d.preprocessing.tracking import visualize_tracks, combine_multicolor_tracked_outputs_for_base
from behav3d.core.metadata import (
    detect_merged_cell_types_from_metadata,
    is_multicolor_celltype,
    multicolor_base_name,
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
)


_TRACKING_PANEL_REGISTRY = {}


def _multicolor_suffix_index(cell_type):
    import re
    match = re.search(r'_(\d+)_multicolor$', str(cell_type))
    return int(match.group(1)) if match else 0


def _sorted_multicolor_family(cell_types):
    return sorted({str(ct).strip() for ct in cell_types}, key=_multicolor_suffix_index)


_TRACKING_PROFILE_DEFAULTS = {
    "immune": {
        "method": "trackpy",
        "overwrite": False,
        "lap": {
            "track_cost_px": 60,
            "gap_close_cost_px": 45,
            "gap_close_max_frames": 5,
            "merging_cost_px": 0,
            "splitting_cost_px": 0,
        },
        "trackpy": {
            "search_range_px": 31,
            "memory_frames": 5,
            "adaptive_stop": 5.0,
            "adaptive_step": 0.95,
        },
        "btrack": {
            "config_preset": "cell",
            "config_path": "",
            "use_visual_features": False,
            "max_search_radius": 100,
            "update_method": "EXACT",
            "step_size": 100,
            "n_workers": 16,
            "use_optimize": False,
            "hypotheses": ["P_FP", "P_init", "P_term", "P_link"],
            "dist_thresh": 60,
            "time_thresh": 3,
        },
    },
    "other": {
        "method": "trackpy",
        "overwrite": False,
        "lap": {
            "track_cost_px": 60,
            "gap_close_cost_px": 45,
            "gap_close_max_frames": 5,
            "merging_cost_px": 0,
            "splitting_cost_px": 0,
        },
        "trackpy": {
            "search_range_px": 31,
            "memory_frames": 5,
            "adaptive_stop": 5.0,
            "adaptive_step": 0.95,
        },
        "btrack": {
            "config_preset": "cell",
            "config_path": "",
            "use_visual_features": False,
            "max_search_radius": 100,
            "update_method": "EXACT",
            "step_size": 100,
            "n_workers": 8,
            "use_optimize": False,
            "hypotheses": ["P_FP", "P_init", "P_term", "P_link"],
            "dist_thresh": 60,
            "time_thresh": 3,
        },
    },
    "organoid": {
        "method": "propagation",
        "overwrite": False,
        "lap": {
            "track_cost_px": 60,
            "gap_close_cost_px": 80,
            "gap_close_max_frames": 3,
            "merging_cost_px": 0,
            "splitting_cost_px": 0,
        },
        "trackpy": {
            "search_range_px": 35,
            "memory_frames": 2,
            "adaptive_stop": 10.0,
            "adaptive_step": 0.95,
        },
        "btrack": {
            "config_preset": "cell",
            "config_path": "",
            "use_visual_features": False,
            "max_search_radius": 100,
            "update_method": "EXACT",
            "step_size": 100,
            "n_workers": 8,
            "use_optimize": False,
            "hypotheses": ["P_FP", "P_init", "P_term", "P_link"],
            "dist_thresh": 60,
            "time_thresh": 3,
        },
    },
    "all_organoids": {
        "method": "propagation_all_organoids",
        "overwrite": False,
    },
}


def _tracking_profile_defaults(category):
    return deepcopy(_TRACKING_PROFILE_DEFAULTS[category])


class TrackingPanel:
    """
    Minimal UI for tracking (LAP, TrackPy, or Propagation) for any cell type.
    - Automatically detects cell type category (organoid/immune/other)
    - Loads appropriate default configuration based on category
    - Uses metadata_loader.output_dir and metadata_loader.metadata_csv_path directly
    - Updates metadata_loader.metadata and always saves CSV
    - cell_type is supplied in __init__
    """
    def __init__(self, metadata_loader, cell_type, organoid_group_context=None, organoid_panel_index=None, multicolor_group_context=None):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()

        self._detect_all_cell_types()
        self._detect_category()
        self._ensure_cfg_skeleton()
        self.organoid_panel_index = organoid_panel_index
        self.organoid_group_context = self._prepare_organoid_group_context(organoid_group_context)
        self._syncing_organoid_toggle = False
        self.multicolor_group_context = self._prepare_multicolor_group_context(multicolor_group_context)
        self._syncing_multicolor_toggle = False
        self._family_group_widgets = None
        self._family_group_run_button = None
        
        params = dict(self.metadata_loader.behav3d_parameters or {})
        tcfg = self._get_config_for_cell_type(params)
        individual_method_options = [
            ("LAP (laptrack)", "lap"),
            ("TrackPy", "trackpy"),
            ("Propagation", "propagation"),
        ]
        if self.category == "organoid" and not self._has_linked_organoid_mode():
            individual_method_options.append(("Propagation (all organoids)", "propagation_all_organoids"))
        individual_method_options.append(("btrack (Bayesian)", "btrack"))
        self._individual_method_options = individual_method_options
        self._combined_method_options = [("Propagation (all organoids)", "propagation_all_organoids")]
        default_method = str(tcfg.get("method", "propagation" if self.category == "organoid" else "trackpy"))
        valid_methods = {value for _, value in individual_method_options}
        if default_method not in valid_methods:
            default_method = "propagation" if self.category == "organoid" else "trackpy"

        self.method = widgets.Dropdown(
            description="Tracking method",
            options=individual_method_options,
            value=default_method,
            style={'description_width': '160px'}
        )

        self.overwrite = widgets.Checkbox(
            description="Overwrite existing",
            value=bool(tcfg.get("overwrite", False))
        )

        lap = tcfg.get("lap", {})
        self.track_cost_dist = widgets.IntText(
            description="Track cost (px)",
            value=int(lap.get("track_cost_px", 45)),
            style={'description_width':'160px'}
        )
        self.gap_cost_dist = widgets.IntText(
            description="Gap close cost (px)",
            value=int(lap.get("gap_close_cost_px", 60)),
            style={'description_width':'160px'}
        )
        self.gap_max_frames = widgets.IntText(
            description="Gap close max frames",
            value=int(lap.get("gap_close_max_frames", 3)),
            style={'description_width':'160px'}
        )
        self.merging_cost_dist = widgets.IntText(
            description="Merging cost (px)",
            value=int(lap.get("merging_cost_px", 0)),
            style={'description_width':'160px'}
        )
        self.splitting_cost_dist = widgets.IntText(
            description="Splitting cost (px)",
            value=int(lap.get("splitting_cost_px", 0)),
            style={'description_width':'160px'}
        )
        self.lap_params = widgets.VBox([
            widgets.HTML("<b>LAP parameters</b>"),
            self.track_cost_dist, self.gap_cost_dist, self.gap_max_frames,
            self.merging_cost_dist, self.splitting_cost_dist
        ])

        tpy = tcfg.get("trackpy", {})
        self.tp_search_range = widgets.IntText(
            description="Search range (px)",
            value=int(tpy.get("search_range_px", 31)),
            style={'description_width':'160px'}
        )
        self.tp_memory = widgets.IntText(
            description="Memory (frames)",
            value=int(tpy.get("memory_frames", 2)),
            style={'description_width':'160px'}
        )
        self.tp_adaptive_stop = widgets.FloatText(
            description="Adaptive stop",
            value=float(tpy.get("adaptive_stop", 10.0)),
            style={'description_width':'160px'}
        )
        self.tp_adaptive_step = widgets.FloatText(
            description="Adaptive step",
            value=float(tpy.get("adaptive_step", 0.95)),
            style={'description_width':'160px'}
        )
        self.tp_params = widgets.VBox([
            widgets.HTML("<b>TrackPy parameters</b>"),
            self.tp_search_range, self.tp_memory, self.tp_adaptive_stop, self.tp_adaptive_step
        ])

        self.prop_params = widgets.VBox([
            widgets.HTML("<b>Propagation tracking</b>"),
            widgets.HTML("<i>No tunable parameters.</i>")
        ])

        # btrack params
        bt = tcfg.get("btrack", {})
        bt_preset_value, bt_config_path_value = self._normalize_btrack_config(bt)
        self.bt_config_preset = widgets.Dropdown(
            description="Config preset",
            options=[("Cell (bundled)", "cell"),
                     ("Particle (bundled)", "particle"),
                     ("Custom JSON", "custom")],
            value=bt_preset_value,
            style={'description_width': '160px'}
        )
        self.bt_config_path = widgets.Text(
            description="JSON path",
            value=bt_config_path_value,
            placeholder="Path to custom btrack JSON…",
            style={'description_width': '160px'},
            layout=widgets.Layout(width='400px')
        )
        self.bt_config_path.disabled = (self.bt_config_preset.value != "custom")
        def _bt_preset_changed(change):
            self.bt_config_path.disabled = (change['new'] != "custom")
        self.bt_config_preset.observe(_bt_preset_changed, names='value')
        self.bt_use_visual_features = widgets.Checkbox(
            description="Use visual features",
            value=bool(bt.get("use_visual_features", False)),
            indent=False
        )
        self.bt_visual_help = widgets.HTML(
            "<span style='color:#666'>Requires <code>raw_image_path</code> and uses "
            "raw-image-derived features during linking.</span>"
        )

        self.bt_max_search_radius = widgets.IntText(
            description="Max search radius (px)",
            value=int(bt.get("max_search_radius", 100)),
            style={'description_width': '160px'}
        )
        self.bt_update_method = widgets.Dropdown(
            description="Update method",
            options=[("EXACT (recommended)", "EXACT"),
                     ("APPROXIMATE (large datasets)", "APPROXIMATE")],
            value=str(bt.get("update_method", "EXACT")),
            style={'description_width': '160px'}
        )
        self.bt_step_size = widgets.IntText(
            description="Step size (frames)",
            value=int(bt.get("step_size", 100)),
            style={'description_width': '160px'}
        )
        self.bt_n_workers = widgets.IntText(
            description="Workers",
            value=int(bt.get("n_workers", 8)),
            style={'description_width': '160px'}
        )
        self.bt_use_optimize = widgets.Checkbox(
            description="Enable global track optimization",
            value=bool(bt.get("use_optimize", False)),
            indent=False
        )
        saved_hyps = bt.get("hypotheses", ["P_FP", "P_init", "P_term", "P_link"])
        self.bt_hyp_checks = {}
        hyp_checkboxes = []
        for hyp_name, hyp_desc, default_on in [
            ("P_FP",     "False positive",    True),
            ("P_init",   "Track initialization", True),
            ("P_term",   "Track termination",  True),
            ("P_link",   "Track linking",      True),
            ("P_branch", "Track branching",    False),
            ("P_dead",   "Cell death",         False),
            ("P_merge",  "Track merging",      False),
        ]:
            is_on = hyp_name in saved_hyps if saved_hyps else default_on
            cb = widgets.Checkbox(
                description=f"{hyp_name} — {hyp_desc}",
                value=is_on,
                indent=False,
                disabled=(hyp_name == "P_FP"),
            )
            self.bt_hyp_checks[hyp_name] = cb
            hyp_checkboxes.append(cb)

        self.bt_dist_thresh = widgets.IntText(
            description="Distance threshold",
            value=int(bt.get("dist_thresh", 60)),
            style={'description_width': '160px'}
        )
        self.bt_time_thresh = widgets.IntText(
            description="Time threshold",
            value=int(bt.get("time_thresh", 3)),
            style={'description_width': '160px'}
        )

        self.bt_optimizer_box = widgets.VBox([
            widgets.HTML("<b>Hypotheses</b>"),
            *hyp_checkboxes,
            self.bt_dist_thresh,
            self.bt_time_thresh,
        ])
        # Toggle optimizer widgets visibility
        self.bt_optimizer_box.layout.display = "flex" if self.bt_use_optimize.value else "none"
        def _bt_optimize_toggled(change):
            self.bt_optimizer_box.layout.display = "flex" if change['new'] else "none"
        self.bt_use_optimize.observe(_bt_optimize_toggled, names='value')

        self.bt_params = widgets.VBox([
            widgets.HTML(
                "<b>btrack — Bayesian multi-object tracker</b><br>"
                "<span style='color:#666'>"
                "<b>Step 1</b> (Kalman filter): always runs — tune and validate first "
                "by running <i>without</i> optimization.<br>"
                "<b>Step 2</b> (optimizer): opt-in — enable only once Step 1 looks correct. "
                "See <code>models/README_btrack.md</code> for full parameter reference."
                "</span>"
            ),
            self.bt_config_preset,
            self.bt_config_path,
            self.bt_use_visual_features,
            self.bt_visual_help,
            self.bt_max_search_radius,
            self.bt_update_method,
            self.bt_step_size,
            self.bt_n_workers,
            widgets.HTML("<hr>"),
            self.bt_use_optimize,
            self.bt_optimizer_box,
        ])

        self.param_container = widgets.VBox()
        self.method.observe(self._toggle_param_groups, names='value')
        self._toggle_param_groups()

        self.title_html = widgets.HTML(f"<b>{cell_type} Tracking</b>")
        self.track_organoids_together = widgets.Checkbox(
            description="Track organoids together",
            value=self._organoids_together_enabled(),
            indent=False,
        )
        self.track_organoids_together.observe(self._on_track_organoids_together_changed, names='value')
        self.track_organoids_together_box = widgets.HBox([self.track_organoids_together])
        self.track_organoids_together_box.layout.display = "flex" if self._has_linked_organoid_mode() else "none"
        self.organoid_group_info = widgets.HTML("")
        self.organoid_group_info.layout.display = "none"

        self.btn_run = widgets.Button(
            description=f"Run {cell_type} tracking",
            button_style="success",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        self.spinner_run = widgets.HTML(value=spinning_loader)
        self.spinner_run.layout.display = "none"
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_run],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        self.btn_run.on_click(self._on_run_clicked)
        self.out = widgets.Output()

        self._build_import_section()

        _TRACKING_PANEL_REGISTRY[self.cell_type] = self

        self.ui = widgets.VBox([
            self.title_html,
            self.track_organoids_together_box,
            self.organoid_group_info,
            self._import_accordion,
            self.method,
            widgets.HBox([self.overwrite]),
            self.param_container,
            self.run_row,
            widgets.HTML("<hr>"),
            self.out
        ])

        self._load_individual_profile_into_widgets(tcfg)
        if self._has_linked_organoid_mode():
            self._register_organoid_group_panel()
            self._apply_organoid_group_mode(load_from_config=True)
        else:
            self._render_individual_panel(load_from_config=True)

        self._register_multicolor_group_panel()
        self._apply_multicolor_group_mode(load_from_config=True)

    def display(self):
        widgets.display(self.ui)

    def _build_import_section(self):
        """Build external tracking image import widget for this cell type."""
        ct = self.cell_type
        cat = self.category
        prefix = {"organoid": "or", "immune": "im", "other": "ot"}[cat]
        img_col = f"{prefix}_{ct}_tracks_image_path"
        csv_col = f"{prefix}_{ct}_tracks_csv_path"
        metadata = self.metadata_loader.metadata
        sample_names = list(metadata["sample_name"].astype(str))

        out_dir = Path(self.metadata_loader.output_dir) / "images"
        trackdata_dir = Path(self.metadata_loader.output_dir) / "trackdata"
        # Note: per-sample subdirs are created inside the callback based on sample_name

        # Build a callback that creates per-sample trackdata dirs
        def _tracking_callback_factory():
            from .utils import make_tracking_metadata_callback as _make_cb
            def _on_converted(sample_name, src, zarr_path):
                sample_trackdata = trackdata_dir / sample_name / ct
                cb = _make_cb(self.metadata_loader, img_col, csv_col, str(sample_trackdata))
                cb(sample_name, src, zarr_path)
            return _on_converted

        self._track_importer = ExternalImageImporter(
            label=f"{ct} tracking",
            output_dir=str(out_dir),
            on_converted=_tracking_callback_factory(),
            sample_names=sample_names,
            output_filename="{sample_name}_" + ct + "_tracked.zarr",
            is_label_image=True,
        )
        self._import_accordion = widgets.Accordion(
            children=[self._track_importer],
        )
        self._import_accordion.set_title(0, "Import external tracking")
        self._import_accordion.selected_index = None

    def _prepare_multicolor_group_context(self, multicolor_group_context):
        if not is_multicolor_celltype(self.cell_type):
            return None

        base_cell_type = multicolor_base_name(self.cell_type)
        context = dict(multicolor_group_context or {})
        context.setdefault("base_cell_type", base_cell_type)
        context.setdefault("enabled", True)
        context.setdefault("panels", [])

        members = list(context.get("members") or [])
        if not members:
            candidates = list(self.organoid_types) + list(self.immune_types) + list(self.other_types)
            members = [ct for ct in candidates if is_multicolor_celltype(ct) and multicolor_base_name(ct) == base_cell_type]
        context["members"] = _sorted_multicolor_family(members)
        return context

    def _register_multicolor_group_panel(self):
        if self.multicolor_group_context is None:
            return
        panels = self.multicolor_group_context.setdefault("panels", [])
        if self not in panels:
            panels.append(self)

    def _has_linked_multicolor_mode(self):
        if self.multicolor_group_context is None:
            return False
        return len(self.multicolor_group_context.get("members", [])) >= 2

    def _is_first_multicolor_panel(self):
        if not self._has_linked_multicolor_mode():
            return False
        members = self.multicolor_group_context.get("members", [])
        return bool(members) and self.cell_type == members[0]

    def _multicolor_group_panels(self):
        if not self._has_linked_multicolor_mode():
            return []
        return [
            _TRACKING_PANEL_REGISTRY.get(ct)
            for ct in self.multicolor_group_context.get("members", [])
        ]

    def _apply_multicolor_group_mode(self, load_from_config=True):
        if not self._has_linked_multicolor_mode():
            return
        
        if not self._is_first_multicolor_panel():
            self.ui.layout.display = "none"
            return
        
        # Hide individual run button and replace with consolidated multicolor run button
        self.run_row.layout.display = "none"
        
        if self._family_group_widgets is None:
            self._family_group_widgets = widgets.Checkbox(
                description=f"Track {self.multicolor_group_context['base_cell_type']} multicolor together",
                value=True,
                disabled=True,
                indent=False,
            )
            self._family_group_run_button = widgets.Button(
                description=f"Run {self.multicolor_group_context['base_cell_type']} multicolor tracking",
                button_style="warning",
                layout=widgets.Layout(width="fit-content", flex="0 0 auto")
            )
            self._family_group_run_button.on_click(self._on_multicolor_group_run_clicked)

            group_info = widgets.HTML(
                "<div style='color:#555;font-size:13px;'>"
                f"Multicolor family: <b>{', '.join(self.multicolor_group_context.get('members', []))}</b>"
                "</div>"
            )
            group_row = widgets.HBox([
                self._family_group_widgets,
                self._family_group_run_button,
            ], layout=widgets.Layout(align_items="center", gap="8px"))
            group_box = widgets.VBox([group_info, group_row])

            children = list(self.ui.children)
            children.insert(1, group_box)
            self.ui.children = tuple(children)
            self.ui.layout.display = "flex"

    def _on_multicolor_group_run_clicked(self, _):
        if not self._has_linked_multicolor_mode() or not self._is_first_multicolor_panel():
            return

        members = self.multicolor_group_context.get("members", [])
        base_cell_type = self.multicolor_group_context.get("base_cell_type", multicolor_base_name(self.cell_type))
        panels = [panel for panel in self._multicolor_group_panels() if panel is not None]

        self._lock(True)
        self.spinner_run.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                self._force_commit_pending_changes()
                for panel in panels:
                    panel._run_tracking_internal(merge_multicolor=False)
                self.metadata_loader.metadata = combine_multicolor_tracked_outputs_for_base(
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(Path(self.metadata_loader.output_dir).expanduser()),
                    base_cell_type=base_cell_type,
                    n_channels=len(members),
                    overwrite=bool(self.overwrite.value),
                    n_workers=1,
                )
                csv_path = Path(self.metadata_loader.metadata_csv_path).expanduser()
                self.metadata_loader.metadata.to_csv(csv_path, sep=",", index=False)
                print("✅ Multicolor tracking finished.")
            except FileNotFoundError:
                print("⚠️ Multicolor merge skipped: expected output files were not found.")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_run.layout.display = "none"
                self._lock(False)

    def _detect_all_cell_types(self):
        from behav3d.io.images import load_image, load_zarr, save_as_zarr
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
        )
        metadata = self.metadata_loader.metadata
        if metadata is None: raise ValueError("Metadata not loaded.")
        self.organoid_types = detect_organoid_types_from_metadata(metadata)
        self.immune_types = detect_immune_cell_types_from_metadata(metadata)
        self.other_types = detect_other_cell_types_from_metadata(metadata)
        self.merged_types = detect_merged_cell_types_from_metadata(metadata)

    def _detect_category(self):
        self.category = detect_cell_type_category(self.cell_type, self.metadata_loader.metadata)
    
    def _get_config_for_cell_type(self, params):
        tracking_cfg = params.get("tracking", {})
        if self.cell_type in tracking_cfg: return tracking_cfg[self.cell_type]
        if self.category == 'organoid': return tracking_cfg.get('organoid', {})
        if self.category == 'immune': return tracking_cfg.get('immune', {})
        return tracking_cfg.get('other', {})

    def _normalize_btrack_config(self, bt):
        valid_presets = {"cell", "particle", "custom"}

        preset_raw = bt.get("config_preset", "")
        path_raw = bt.get("config_path", "")

        preset = str(preset_raw).strip() if preset_raw is not None else ""
        config_path = str(path_raw).strip() if path_raw is not None else ""

        if preset in valid_presets:
            return preset, config_path
        if preset:
            return "custom", preset
        if config_path:
            return "custom", config_path
        return "cell", ""

    def _prepare_organoid_group_context(self, organoid_group_context):
        if self.category != "organoid" or organoid_group_context is None:
            return None
        if not isinstance(organoid_group_context, dict):
            raise TypeError("organoid_group_context must be a dict or None.")

        tracking_cfg = dict((self.metadata_loader.behav3d_parameters or {}).get("tracking", {}))
        organoid_group_context["organoid_types"] = list(
            organoid_group_context.get("organoid_types") or self.organoid_types
        )
        organoid_group_context.setdefault(
            "enabled",
            bool(tracking_cfg.get("track_organoids_together", False)),
        )
        organoid_group_context.setdefault("panels", [])
        return organoid_group_context

    def _register_organoid_group_panel(self):
        if self.organoid_group_context is None:
            return
        panels = self.organoid_group_context.setdefault("panels", [])
        if self not in panels:
            panels.append(self)

    def _has_linked_organoid_mode(self):
        if self.category != "organoid" or self.organoid_group_context is None:
            return False
        return len(self.organoid_group_context.get("organoid_types", [])) >= 2

    def _organoids_together_enabled(self):
        return self._has_linked_organoid_mode() and bool(self.organoid_group_context.get("enabled", False))

    def _is_first_organoid_panel(self):
        if not self._has_linked_organoid_mode():
            return False
        if self.organoid_panel_index is not None:
            return int(self.organoid_panel_index) == 0
        organoid_types = self.organoid_group_context.get("organoid_types", [])
        return bool(organoid_types) and self.cell_type == organoid_types[0]

    def _write_behav3d_parameters(self):
        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)

    def _get_all_organoids_config(self, params):
        tracking_cfg = params.get("tracking", {})
        return tracking_cfg.get("all_organoids", {})

    def _load_individual_profile_into_widgets(self, profile):
        profile = dict(profile or {})
        self.overwrite.value = bool(profile.get("overwrite", False))

        lap = profile.get("lap", {})
        self.track_cost_dist.value = int(lap.get("track_cost_px", 45))
        self.gap_cost_dist.value = int(lap.get("gap_close_cost_px", 60))
        self.gap_max_frames.value = int(lap.get("gap_close_max_frames", 3))
        self.merging_cost_dist.value = int(lap.get("merging_cost_px", 0))
        self.splitting_cost_dist.value = int(lap.get("splitting_cost_px", 0))

        tpy = profile.get("trackpy", {})
        self.tp_search_range.value = int(tpy.get("search_range_px", 31))
        self.tp_memory.value = int(tpy.get("memory_frames", 2))
        self.tp_adaptive_stop.value = float(tpy.get("adaptive_stop", 10.0))
        self.tp_adaptive_step.value = float(tpy.get("adaptive_step", 0.95))

        bt = profile.get("btrack", {})
        bt_preset_value, bt_config_path_value = self._normalize_btrack_config(bt)
        self.bt_config_preset.value = bt_preset_value
        self.bt_config_path.value = bt_config_path_value
        self.bt_use_visual_features.value = bool(bt.get("use_visual_features", False))
        self.bt_max_search_radius.value = int(bt.get("max_search_radius", 100))
        self.bt_update_method.value = str(bt.get("update_method", "EXACT"))
        self.bt_step_size.value = int(bt.get("step_size", 100))
        self.bt_n_workers.value = int(bt.get("n_workers", 8))
        self.bt_use_optimize.value = bool(bt.get("use_optimize", False))

        saved_hyps = bt.get("hypotheses", ["P_FP", "P_init", "P_term", "P_link"])
        default_hyps = {
            "P_FP": True,
            "P_init": True,
            "P_term": True,
            "P_link": True,
            "P_branch": False,
            "P_dead": False,
            "P_merge": False,
        }
        for hyp_name, cb in self.bt_hyp_checks.items():
            if hyp_name == "P_FP":
                cb.value = True
            else:
                cb.value = hyp_name in saved_hyps if saved_hyps else default_hyps.get(hyp_name, False)

        self.bt_dist_thresh.value = int(bt.get("dist_thresh", 60))
        self.bt_time_thresh.value = int(bt.get("time_thresh", 3))
        self._toggle_param_groups()

    def _set_method_options_and_value(self, options, value):
        self.method.options = options
        valid_values = {option_value for _, option_value in options}
        if value not in valid_values:
            value = options[0][1]
        self.method.value = value

    def _render_individual_panel(self, load_from_config=True):
        self.ui.layout.display = "flex"
        self.title_html.value = f"<b>{self.cell_type} Tracking</b>"
        self.organoid_group_info.layout.display = "none"
        self.btn_run.description = f"Run {self.cell_type} tracking"

        if load_from_config:
            params = dict(self.metadata_loader.behav3d_parameters or {})
            profile = self._get_config_for_cell_type(params)
            default_method = "propagation" if self.category == "organoid" else "trackpy"
            method_value = str(profile.get("method", default_method))
            self._set_method_options_and_value(self._individual_method_options, method_value)
            self._load_individual_profile_into_widgets(profile)
        else:
            self._set_method_options_and_value(self._individual_method_options, self.method.value)
            self._toggle_param_groups()

    def _render_all_organoids_panel(self, load_from_config=True):
        self.ui.layout.display = "flex"
        organoid_types = self.organoid_group_context.get("organoid_types", [])
        organoid_text = ", ".join(organoid_types) if organoid_types else "No organoid types detected"
        self.title_html.value = "<b>All Organoids Tracking</b>"
        self.organoid_group_info.value = (
            "<div style='color:#555;font-size:13px;'>"
            f"Tracking together: <b>{organoid_text}</b>"
            "</div>"
        )
        self.organoid_group_info.layout.display = "flex"
        self.btn_run.description = "Run all organoid tracking"

        params = dict(self.metadata_loader.behav3d_parameters or {})
        profile = self._get_all_organoids_config(params) if load_from_config else {}
        method_value = str(profile.get("method", "propagation_all_organoids"))
        self._set_method_options_and_value(self._combined_method_options, method_value)
        self.overwrite.value = bool(profile.get("overwrite", False))
        self._toggle_param_groups()

    def _apply_organoid_group_mode(self, load_from_config=True):
        if not self._has_linked_organoid_mode():
            self.ui.layout.display = "flex"
            return

        self._syncing_organoid_toggle = True
        try:
            self.track_organoids_together.value = self._organoids_together_enabled()
        finally:
            self._syncing_organoid_toggle = False

        if self._organoids_together_enabled():
            if self._is_first_organoid_panel():
                self._render_all_organoids_panel(load_from_config=load_from_config)
            else:
                self.ui.layout.display = "none"
        else:
            self._render_individual_panel(load_from_config=load_from_config)

    def _on_track_organoids_together_changed(self, change):
        if self._syncing_organoid_toggle or not self._has_linked_organoid_mode():
            return

        panels = list(self.organoid_group_context.get("panels", []))
        for panel in panels:
            panel._persist_current_mode(write_file=False)

        self.organoid_group_context["enabled"] = bool(change["new"])
        self._ensure_cfg_skeleton()
        params = dict(self.metadata_loader.behav3d_parameters or {})
        params.setdefault("tracking", {})["track_organoids_together"] = bool(change["new"])
        self.metadata_loader.behav3d_parameters = params
        self._write_behav3d_parameters()

        for panel in panels:
            panel._apply_organoid_group_mode(load_from_config=True)

    def _toggle_param_groups(self, change=None):
        if self.method.value == "lap": self.param_container.children = (self.lap_params,)
        elif self.method.value == "trackpy": self.param_container.children = (self.tp_params,)
        elif self.method.value == "btrack": self.param_container.children = (self.bt_params,)
        else: self.param_container.children = (self.prop_params,)

    def _lock(self, state: bool):
        for w in [self.method, self.overwrite, self.track_cost_dist, self.gap_cost_dist, self.gap_max_frames,
                  self.merging_cost_dist, self.splitting_cost_dist, self.tp_search_range, self.tp_memory,
                  self.tp_adaptive_stop, self.tp_adaptive_step,
                  self.bt_config_preset, self.bt_config_path, self.bt_use_visual_features,
                  self.bt_max_search_radius,
                  self.bt_update_method, self.bt_step_size, self.bt_n_workers, self.bt_use_optimize,
                self.bt_dist_thresh, self.bt_time_thresh, self.track_organoids_together,
                self._family_group_run_button,
                  self.btn_run]:
            if hasattr(w, "disabled"): w.disabled = state
        for cb in self.bt_hyp_checks.values():
            if cb is not self.bt_hyp_checks.get("P_FP"):
                cb.disabled = state           

    def _force_commit_pending_changes(self):
        try:
            from IPython.display import Javascript, display
            widgets.display(Javascript("document.activeElement && document.activeElement.blur();"))
            import time as _t; _t.sleep(0.08)
        except Exception: pass

    def _ensure_cfg_skeleton(self):
        p = dict(self.metadata_loader.behav3d_parameters or {})
        tracking_cfg = p.setdefault("tracking", {})
        tracking_cfg.setdefault("immune", _tracking_profile_defaults("immune"))
        tracking_cfg.setdefault("other", _tracking_profile_defaults("other"))
        tracking_cfg.setdefault("organoid", _tracking_profile_defaults("organoid"))
        tracking_cfg.setdefault("all_organoids", _tracking_profile_defaults("all_organoids"))
        tracking_cfg.setdefault("track_organoids_together", False)
        self.metadata_loader.behav3d_parameters = p

    def _persist_current_mode(self, write_file=True):
        self._ensure_cfg_skeleton()
        params = dict(self.metadata_loader.behav3d_parameters or {})
        tracking_cfg = params.setdefault("tracking", {})
        tracking_cfg["track_organoids_together"] = bool(
            self.organoid_group_context.get("enabled", False)
        ) if self._has_linked_organoid_mode() else bool(tracking_cfg.get("track_organoids_together", False))

        if self._organoids_together_enabled() and not self._is_first_organoid_panel():
            self.metadata_loader.behav3d_parameters = params
            if write_file:
                self._write_behav3d_parameters()
            return

        if self._organoids_together_enabled() and self._is_first_organoid_panel():
            prof = tracking_cfg.setdefault("all_organoids", _tracking_profile_defaults("all_organoids"))
            prof["method"] = str(self.method.value)
            prof["overwrite"] = bool(self.overwrite.value)
        else:
            if self.cell_type not in tracking_cfg:
                base_key = self.category if self.category in _TRACKING_PROFILE_DEFAULTS else "other"
                tracking_cfg[self.cell_type] = _tracking_profile_defaults(base_key)
            prof = tracking_cfg[self.cell_type]
            prof["method"] = str(self.method.value)
            prof["overwrite"] = bool(self.overwrite.value)
            prof.setdefault("lap", {}).update({
                "track_cost_px": int(self.track_cost_dist.value),
                "gap_close_cost_px": int(self.gap_cost_dist.value),
                "gap_close_max_frames": int(self.gap_max_frames.value),
                "merging_cost_px": int(self.merging_cost_dist.value),
                "splitting_cost_px": int(self.splitting_cost_dist.value),
            })
            prof.setdefault("trackpy", {}).update({
                "search_range_px": int(self.tp_search_range.value),
                "memory_frames": int(self.tp_memory.value),
                "adaptive_stop": float(self.tp_adaptive_stop.value),
                "adaptive_step": float(self.tp_adaptive_step.value),
            })
            bt = prof.setdefault("btrack", {})
            preset_val = str(self.bt_config_preset.value)
            bt["config_preset"] = "custom" if preset_val == "custom" else preset_val
            bt["config_path"] = self.bt_config_path.value.strip() if preset_val == "custom" else ""
            bt["use_visual_features"] = bool(self.bt_use_visual_features.value)
            bt["max_search_radius"] = int(self.bt_max_search_radius.value)
            bt["update_method"] = str(self.bt_update_method.value)
            bt["step_size"] = int(self.bt_step_size.value)
            bt["n_workers"] = max(1, int(self.bt_n_workers.value))
            bt["use_optimize"] = bool(self.bt_use_optimize.value)
            bt["hypotheses"] = [name for name, cb in self.bt_hyp_checks.items() if cb.value]
            bt["dist_thresh"] = int(self.bt_dist_thresh.value)
            bt["time_thresh"] = int(self.bt_time_thresh.value)

        self.metadata_loader.behav3d_parameters = params
        if write_file:
            self._write_behav3d_parameters()

    def _persist_params(self):
        self._persist_current_mode(write_file=True)

    def _run_tracking_internal(self, merge_multicolor=True):
        self._lock(True)
        self.spinner_run.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                self._force_commit_pending_changes()
                self._persist_params()
                out_dir = Path(self.metadata_loader.output_dir).expanduser()
                out_dir.mkdir(parents=True, exist_ok=True)
                csv_path = Path(self.metadata_loader.metadata_csv_path).expanduser()

                method = self.method.value
                run_all_organoids = (
                    (self._organoids_together_enabled() and self._is_first_organoid_panel())
                    or (method == "propagation_all_organoids")
                )
                if method == "lap":
                    tc, gc, mc, sc = int(self.track_cost_dist.value)**2, int(self.gap_cost_dist.value)**2, int(self.merging_cost_dist.value), int(self.splitting_cost_dist.value)
                    new_md = run_laptracking(metadata=self.metadata_loader.metadata, output_dir=str(out_dir), cell_type=self.cell_type, track_cost_cutoff=tc, gap_closing_cost_cutoff=gc, gap_closing_max_frame_count=int(self.gap_max_frames.value), merging_cost_cutoff=(mc**2 if mc>0 else False), splitting_cost_cutoff=(sc**2 if sc>0 else False), overwrite=bool(self.overwrite.value))
                elif method == "trackpy":
                    new_md = run_trackpy_tracking_generic(metadata=self.metadata_loader.metadata, output_dir=str(out_dir), cell_type=self.cell_type, overwrite=bool(self.overwrite.value), search_range=int(self.tp_search_range.value), memory=int(self.tp_memory.value), adaptive_stop=float(self.tp_adaptive_stop.value), adaptive_step=float(self.tp_adaptive_step.value))
                elif method == "btrack":
                    preset_val = str(self.bt_config_preset.value)
                    config_preset = self.bt_config_path.value.strip() if preset_val == "custom" else preset_val
                    use_opt = bool(self.bt_use_optimize.value)
                    hyps = [name for name, cb in self.bt_hyp_checks.items() if cb.value] if use_opt else None
                    new_md = run_btracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        cell_type=self.cell_type,
                        overwrite=bool(self.overwrite.value),
                        config_preset=config_preset,
                        use_visual_features=bool(self.bt_use_visual_features.value),
                        max_search_radius=int(self.bt_max_search_radius.value),
                        update_method=str(self.bt_update_method.value),
                        step_size=int(self.bt_step_size.value),
                        n_workers=max(1, int(self.bt_n_workers.value)),
                        use_optimize=use_opt,
                        hypotheses=hyps,
                        dist_thresh=int(self.bt_dist_thresh.value) if use_opt else None,
                        time_thresh=int(self.bt_time_thresh.value) if use_opt else None,
                    )
                else:
                    new_md = run_propagation_tracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        cell_type=self.cell_type,
                        overwrite=bool(self.overwrite.value),
                        all_organoids=run_all_organoids,
                    )

                self.metadata_loader.metadata = new_md

                if merge_multicolor and is_multicolor_celltype(self.cell_type):
                    base_cell_type = multicolor_base_name(self.cell_type)
                    family_cell_types = [
                        ct for ct in (self.organoid_types + self.immune_types + self.other_types)
                        if is_multicolor_celltype(ct) and multicolor_base_name(ct) == base_cell_type
                    ]
                    if family_cell_types:
                        parsed = []
                        for ct in family_cell_types:
                            match = None
                            import re
                            match = re.search(r'_(\d+)_multicolor$', ct)
                            if match:
                                parsed.append((int(match.group(1)), ct))
                        family_cell_types = [ct for _, ct in sorted(parsed, key=lambda item: item[0])]
                        try:
                            self.metadata_loader.metadata = combine_multicolor_tracked_outputs_for_base(
                                metadata=self.metadata_loader.metadata,
                                output_dir=str(out_dir),
                                base_cell_type=base_cell_type,
                                n_channels=len(family_cell_types),
                                overwrite=bool(self.overwrite.value),
                                n_workers=1,
                            )
                        except FileNotFoundError:
                            pass

                self.metadata_loader.metadata.to_csv(csv_path, sep=",", index=False)
                print("✅ Tracking finished.")
            except Exception: traceback.print_exc()
            finally:
                self.spinner_run.layout.display = "none"
                self._lock(False)

    def _on_run_clicked(self, _):
        self._run_tracking_internal(merge_multicolor=True)

class TrackingVisualizationPanel:
    def __init__(self, metadata_loader, default_time_range=None, channel_colors=None):
        self.metadata_loader = metadata_loader
        self.channel_colors = tuple(channel_colors or _cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.channel_colors", ["cyan", "yellow", "red", "green", "magenta", "blue"]))
        self._viewer = None
        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(description="Refresh", tooltip="Build/Update selector from metadata_loader.metadata")
        self._refresh_btn.on_click(self._on_refresh_clicked)
        self.sample_dropdown = widgets.Dropdown(options=[], value=None, description="Sample:", layout=widgets.Layout(width="350px"), disabled=True)
        self.use_range = widgets.Checkbox(description="Use custom time range", value=bool(_cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.use_range", False)), indent=False)

        if isinstance(default_time_range, (tuple, list)) and len(default_time_range) == 2:
            _start_default, _end_default = map(int, default_time_range)
        else:
            _start_default = int(_cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.start_t", 0))
            _end_default   = int(_cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.end_t", 0))

        self.start_t = widgets.IntText(description="Start T:", value=_start_default, layout=widgets.Layout(width="180px"))
        self.end_t   = widgets.IntText(description="End T:", value=_end_default, layout=widgets.Layout(width="180px"))
        self.range_box = widgets.HBox([self.start_t, self.end_t])
        self.range_box.layout.display = "flex" if self.use_range.value else "none"
        self.use_range.observe(lambda ch: setattr(self.range_box.layout, 'display', 'flex' if ch['new'] else 'none'), names="value")

        self.open_button = widgets.Button(description="Visualize segmentation and tracking results", button_style="primary", tooltip="Launch Napari for the selected sample", icon="eye", layout=widgets.Layout(width="300px"), disabled=True)
        self.close_button = widgets.Button(description="Close viewer", button_style="danger", icon="stop", tooltip="Close the active Napari viewer", layout=widgets.Layout(width="200px", display="none"))
        self.msg = widgets.Output()
        self.open_button.on_click(self._on_open_clicked)
        self.close_button.on_click(self._on_close_clicked)
        self._panel = widgets.VBox([widgets.HBox([self._status, self._refresh_btn]), self.sample_dropdown, self.use_range, self.range_box, widgets.HBox([self.open_button, self.close_button]), self.msg])
        self._maybe_build_from_loader()

    def display(self):
        widgets.display(self._panel)

    def get_selected_row(self) -> pd.Series:
        name = self.sample_dropdown.value
        df = self.metadata_loader.metadata
        return df[df["sample_name"].astype(str) == str(name)].iloc[0]

    def open_viewer(self):
        row = self.get_selected_row()
        time_range = (int(self.start_t.value), int(self.end_t.value)) if self.use_range.value else None
        with self.msg:
            self.msg.clear_output()
            try:
                if self._viewer is not None:
                    try:
                        self._viewer.close()
                    except Exception:
                        pass
                    finally:
                        self._viewer = None

                viewer = visualize_tracks(
                    metadata_row=row,
                    timepoint_range=time_range,
                    channel_colors=self.channel_colors,
                )
                self._viewer = viewer
            except Exception as e:
                self._viewer = None
                print(f"Error: {e}")
            self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"

    def _maybe_build_from_loader(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if isinstance(df, pd.DataFrame) and not df.empty and ("sample_name" in df.columns):
            sample_names = df["sample_name"].astype(str).unique().tolist()
            desired = _cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.sample_name", None)
            self.sample_dropdown.options = sample_names
            self.sample_dropdown.value = desired if (desired in sample_names) else sample_names[0]
            self.sample_dropdown.disabled = False
            self.open_button.disabled = False
            self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
        else:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            self.sample_dropdown.disabled = True
            self.open_button.disabled = True

    def _on_open_clicked(self, _): self.open_viewer()
    def _on_close_clicked(self, _):
        with self.msg:
            try:
                if self._viewer: self._viewer.close(); self._viewer = None
            finally: self.close_button.layout.display = "none"
    def _on_refresh_clicked(self, _):
        with self.msg: self.msg.clear_output(); self._maybe_build_from_loader()
