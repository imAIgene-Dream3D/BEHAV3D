import ipywidgets as widgets
from pathlib import Path
import os
import yaml
import traceback
from functools import partial
from .utils import (
    _cfg_get, 
    _mk_timepoint_range, 
    spinning_loader,
    PathPicker
)
from behav3d.preprocessing.segmentation.napari_pixelclassifier import (
    train_pixel_classifier,
    run_pixel_classifier_segmentation
)
import napari

class PixelClassifierPanel:
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        
        self._detect_cell_types()
        
        print(f"Detected organoid types: {self.organoid_types}")
        print(f"Detected immune cell types: {self.immune_types}")
        print(f"Detected other cell types: {self.other_types}")
        print(f"Dead channel present: {self.has_death}")

        self.examples_per_sample = widgets.IntText(
            description="Examples per sample",
            value=int(pc.get("examples_per_sample", 3))
        )
        self.sample_specific_classifier = widgets.Checkbox(
            description="Sample-specific classifier",
            value=bool(pc.get("sample_specific_classifier", False))
        )
        
        self.n_workers = widgets.IntText(
            description="Workers",
            value=int(pc.get("workers", (os.cpu_count() or 8))),
            max=max(8, (os.cpu_count() or 8))
        )

        self.manual_clf_paths = widgets.Checkbox(
            description="Manually supply classifiers",
            value=False,  # Always start unchecked, user must opt-in
            indent=False
        )

        self.overwrite_existing = widgets.Checkbox(
            description="Overwrite existing",
            value=bool(pc.get("overwrite_existing", False))
        )

        self.clf_dir = PathPicker(
            mode='dir',
            description='Classifier dir:',
            default="",
            description_width='160px',
            width='100%',
        )
        
        self.clf_paths = {}
        self._create_clf_path_pickers()
        
        self.clf_death_path = PathPicker(
            mode='file',
            description='Death clf:',
            default="",
            filter_pattern='*.joblib',
            description_width='160px',
            width='100%',
        )
        
        if self.manual_clf_paths.value:
            self.clf_dir.value = str(pc.get("clf_dir", "") or "")
            self.clf_death_path.value = str(pc.get("clf_death_path", "") or "")
            for cell_type in self.all_cell_types:
                if cell_type in self.clf_paths:
                    saved_path = pc.get(f"clf_{cell_type}_path", "")
                    if saved_path:
                        self.clf_paths[cell_type].value = str(saved_path)

        self.edt_thresholds = {}
        self._create_edt_threshold_inputs()
        self._create_postprocessing_inputs()
        self.use_all_timepoints = widgets.Checkbox(
            description="Process ALL timepoints",
            value=bool(pc.get("use_all_timepoints", True))
        )
        self.tp_start = widgets.IntText(description="Start t", value=int(pc.get("tp_start", 0)))
        self.tp_end   = widgets.IntText(description="End t",   value=int(pc.get("tp_end", 0)))

        self.use_all_timepoints.observe(self._toggle_timepoint_inputs, names='value')
        self._toggle_timepoint_inputs()

        self.btn_train = widgets.Button(
            description="Train pixel classifier",
            button_style="primary",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        self.close_button = widgets.Button(
            description="Close viewer",
            button_style="danger",
            icon="stop",
            tooltip="Close the active Napari viewer",
            layout=widgets.Layout(width="200px", display="none")
        )

        self.spinner_train = widgets.HTML(value=spinning_loader)
        self.spinner_train.layout.display = "none"

        self.train_row = widgets.HBox(
            [self.btn_train, self.close_button, self.spinner_train],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        self.btn_run = widgets.Button(
            description="Run segmentation",
            button_style="success",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        self.btn_resegment = widgets.Button(
            description="Only resegment",
            button_style="warning",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        
        self.spinner_apply = widgets.HTML(value=spinning_loader)
        self.spinner_apply.layout.display = "none"
        self.apply_row = widgets.HBox(
            [self.btn_run, self.btn_resegment, self.spinner_apply],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        self.btn_train.on_click(self._on_train_clicked)
        self.close_button.on_click(self._on_close_clicked)
        self.btn_run.on_click(partial(self._on_apply_clicked, only_segment=False))
        self.btn_resegment.on_click(partial(self._on_apply_clicked, only_segment=True))

        self.out = widgets.Output()

        self.clf_paths_box = widgets.VBox()
        self._build_clf_paths_box()
        
        self.manual_clf_paths.observe(self._toggle_clf_path_section, names='value')
        self._toggle_clf_path_section()

        def _dir_changed(change):
            if not self.manual_clf_paths.value:
                return
            newd = (change.get('new') or '').strip()
            if newd:
                self._apply_dir_to_clf_paths(newd)
        self.clf_dir.text.observe(_dir_changed, names='value')

        # Build segmentation parameter widgets organized by cell type
        segmentation_widgets = [widgets.HTML("<b>Segmentation parameters per cell type</b>")]
        
        if self.organoid_types:
            segmentation_widgets.append(widgets.HTML("<i>Organoids:</i>"))
            for org_type in self.organoid_types:
                if org_type in self.edt_thresholds:
                    segmentation_widgets.append(widgets.HBox([
                        self.edt_thresholds[org_type],
                        self.segment_size_mins[org_type],
                    ]))
                    segmentation_widgets.append(widgets.HBox([
                        self.opening_nr_pixels[org_type],
                        self.fill_holes[org_type],
                    ]))
        if self.immune_types:
            segmentation_widgets.append(widgets.HTML("<i>Immune cells:</i>"))
            for immune_type in self.immune_types:
                if immune_type in self.edt_thresholds:
                    segmentation_widgets.append(widgets.HBox([
                        self.edt_thresholds[immune_type],
                        self.segment_size_mins[immune_type],
                    ]))
                    segmentation_widgets.append(widgets.HBox([
                        self.opening_nr_pixels[immune_type],
                        self.fill_holes[immune_type],
                    ]))
        if self.other_types:
            segmentation_widgets.append(widgets.HTML("<i>Other cells:</i>"))
            for other_type in self.other_types:
                if other_type in self.edt_thresholds:
                    segmentation_widgets.append(widgets.HBox([
                        self.edt_thresholds[other_type],
                        self.segment_size_mins[other_type],
                    ]))
                    segmentation_widgets.append(widgets.HBox([
                        self.opening_nr_pixels[other_type],
                        self.fill_holes[other_type],
                    ]))

        self.ui = widgets.VBox([
            widgets.VBox([
                widgets.HTML("<b>Train pixel classifier</b>"),
                widgets.HBox([self.examples_per_sample, self.n_workers]),
                self.sample_specific_classifier,
                self.train_row,
            ]),
            widgets.HTML("<hr>"),
            widgets.VBox([
                widgets.HTML("<b>Apply segmentation</b>"),
                self.manual_clf_paths,
                self.clf_paths_box,
                *segmentation_widgets,
                self.overwrite_existing,
                widgets.HBox([self.use_all_timepoints, self.tp_start, self.tp_end]),
                self.apply_row,
            ]),
            widgets.HTML("<hr>"),
            self.out
        ])

        self._viewer = None

    def _detect_cell_types(self):
        from behav3d.io.images import load_image, load_zarr, save_as_zarr
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel
        )
        metadata = self.metadata_loader.metadata
        self.organoid_types = detect_organoid_types_from_metadata(metadata)
        self.immune_types = detect_immune_cell_types_from_metadata(metadata)
        self.other_types = detect_other_cell_types_from_metadata(metadata)
        self.has_death = has_dead_channel(metadata)
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types
    
    def _create_clf_path_pickers(self):
        for cell_type in self.all_cell_types:
            self.clf_paths[cell_type] = PathPicker(
                mode='file',
                description=f'{cell_type.capitalize()} clf:',
                default="",
                filter_pattern='*.joblib',
                description_width='160px',
                width='100%',
            )
    
    def _create_edt_threshold_inputs(self):
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        for cell_type in self.all_cell_types:
            if cell_type in self.organoid_types: default_threshold = 12.0
            elif cell_type in self.immune_types: default_threshold = 2.5
            else: default_threshold = 1.0
            saved_threshold = pc.get(f"{cell_type}_edt_threshold", default_threshold)
            self.edt_thresholds[cell_type] = widgets.FloatText(
                description=f"{cell_type.capitalize()} EDT:",
                value=float(saved_threshold),
                style={'description_width': '160px'},
            )
    
    def _create_postprocessing_inputs(self):
        """Create inputs for segment_size_min, opening_nr_pixels, fill_holes per cell type"""
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        
        self.segment_size_mins = {}
        self.opening_nr_pixels = {}
        self.fill_holes = {}
        
        for cell_type in self.all_cell_types:
            # Defaults based on cell type category
            if cell_type in self.organoid_types:
                default_size = 1000
                default_opening = 3
                default_fill = True
            elif cell_type in self.immune_types:
                default_size = 10
                default_opening = 0
                default_fill = True
            else:
                default_size = 10
                default_opening = 0
                default_fill = True
            
            # Segment size min
            saved_size = pc.get(f"{cell_type}_segment_size_min", default_size)
            self.segment_size_mins[cell_type] = widgets.IntText(
                description=f"{cell_type.capitalize()} min size:",
                value=int(saved_size),
                style={'description_width': '160px'},
            )
            
            # Opening nr pixels
            saved_opening = pc.get(f"{cell_type}_opening_nr_pixels", default_opening)
            self.opening_nr_pixels[cell_type] = widgets.IntText(
                description=f"{cell_type.capitalize()} opening px:",
                value=int(saved_opening),
                style={'description_width': '160px'},
            )
            
            # Fill holes
            saved_fill = pc.get(f"{cell_type}_fill_holes", default_fill)
            self.fill_holes[cell_type] = widgets.Checkbox(
                description=f"{cell_type.capitalize()} fill holes",
                value=bool(saved_fill),
                indent=False,
            )
    
    def _persist_params(self):
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc["examples_per_sample"] = int(self.examples_per_sample.value)
        pc["sample_specific_classifier"] = bool(self.sample_specific_classifier.value)
        pc["workers"] = int(self.n_workers.value)
        pc["use_all_timepoints"] = bool(self.use_all_timepoints.value)
        pc["tp_start"] = int(self.tp_start.value)
        pc["tp_end"]   = int(self.tp_end.value)
        pc["manual_clf_paths"] = bool(self.manual_clf_paths.value)
        pc["overwrite_existing"] = bool(self.overwrite_existing.value)
        
        for cell_type, threshold_widget in self.edt_thresholds.items():
            pc[f"{cell_type}_edt_threshold"] = float(threshold_widget.value)
        
        for cell_type in self.all_cell_types:
            if cell_type in self.segment_size_mins:
                pc[f"{cell_type}_segment_size_min"] = int(self.segment_size_mins[cell_type].value)
            if cell_type in self.opening_nr_pixels:
                pc[f"{cell_type}_opening_nr_pixels"] = int(self.opening_nr_pixels[cell_type].value)
            if cell_type in self.fill_holes:
                pc[f"{cell_type}_fill_holes"] = bool(self.fill_holes[cell_type].value)
        
        if self.manual_clf_paths.value:
            pc["clf_dir"] = str(self.clf_dir.value or "")
            pc["clf_death_path"] = str(self.clf_death_path.value or "")
            for cell_type, path_picker in self.clf_paths.items():
                pc[f"clf_{cell_type}_path"] = str(path_picker.value or "")

        if hasattr(self.metadata_loader, "behav3d_parameters_path"):
            yaml.safe_dump(
                self.metadata_loader.behav3d_parameters,
                self.metadata_loader.behav3d_parameters_path.open("w"),
                sort_keys=False
            )

    def _toggle_timepoint_inputs(self, change=None):
        show = not self.use_all_timepoints.value
        disp = None if show else 'none'
        self.tp_start.layout.display = disp
        self.tp_end.layout.display   = disp
        self.tp_start.disabled = not show
        self.tp_end.disabled   = not show

    def _default_clf_paths_for_dir(self, d: str):
        if not d: return {}
        p = Path(d).expanduser()
        paths = {}
        for cell_type in self.all_cell_types:
            paths[cell_type] = str(p / f'PixelClassifier_{cell_type.capitalize()}.joblib')
        paths['death'] = str(p / 'PixelClassifier_Death.joblib')
        return paths

    def _apply_dir_to_clf_paths(self, d: str):
        paths = self._default_clf_paths_for_dir(d)
        for cell_type, path in paths.items():
            if cell_type == 'death': self.clf_death_path.value = path
            elif cell_type in self.clf_paths: self.clf_paths[cell_type].value = path

    def _toggle_clf_path_section(self, change=None):
        manual = bool(self.manual_clf_paths.value)
        self.clf_paths_box.layout.display = (None if manual else 'none')
        if not manual:
            self.clf_dir.value = ""
            self.clf_death_path.value = ""
            for picker in self.clf_paths.values(): picker.value = ""

    def _build_clf_paths_box(self):
        children = [widgets.HTML("<b>Classifier paths</b>"), self.clf_dir]
        if self.organoid_types:
            children.append(widgets.HTML("<i>Organoids:</i>"))
            for org_type in self.organoid_types:
                if org_type in self.clf_paths: children.append(self.clf_paths[org_type])
        if self.immune_types:
            children.append(widgets.HTML("<i>Immune cells:</i>"))
            for immune_type in self.immune_types:
                if immune_type in self.clf_paths: children.append(self.clf_paths[immune_type])
        if self.other_types:
            children.append(widgets.HTML("<i>Other cells:</i>"))
            for other_type in self.other_types:
                if other_type in self.clf_paths: children.append(self.clf_paths[other_type])
        if self.has_death:
            children.append(widgets.HTML("<i>Death mask:</i>"))
            children.append(self.clf_death_path)
        self.clf_paths_box.children = children

    def display(self):
        widgets.display(self.ui)

    def _lock(self, state: bool):
        to_lock = [
            self.btn_train, self.btn_run, self.btn_resegment,
            self.examples_per_sample, self.sample_specific_classifier, self.n_workers,
            self.use_all_timepoints, self.tp_start, self.tp_end,
            self.manual_clf_paths, self.overwrite_existing,
        ]
        to_lock.extend(self.edt_thresholds.values())
        to_lock.extend(self.segment_size_mins.values())
        to_lock.extend(self.opening_nr_pixels.values())
        to_lock.extend(self.fill_holes.values())
        for w in to_lock: w.disabled = state
        try:
            self.clf_dir.text.disabled = state
            self.clf_dir.btn.disabled = state
            if self.has_death:
                self.clf_death_path.text.disabled = state
                self.clf_death_path.btn.disabled = state
            for picker in self.clf_paths.values():
                picker.text.disabled = state
                picker.btn.disabled = state
        except Exception: pass

    def _get_current_params(self):
        """Collect current parameter values from widgets into a dict"""
        params = {}
        for cell_type in self.all_cell_types:
            if cell_type in self.edt_thresholds:
                params[f"{cell_type}_edt_threshold"] = float(self.edt_thresholds[cell_type].value)
            if cell_type in self.segment_size_mins:
                params[f"{cell_type}_segment_size_min"] = int(self.segment_size_mins[cell_type].value)
            if cell_type in self.opening_nr_pixels:
                params[f"{cell_type}_opening_nr_pixels"] = int(self.opening_nr_pixels[cell_type].value)
            if cell_type in self.fill_holes:
                params[f"{cell_type}_fill_holes"] = bool(self.fill_holes[cell_type].value)
        return params
    
    def _update_widgets_from_params(self, params):
        """Update widget values from a params dict and persist to config"""
        for cell_type in self.all_cell_types:
            if f"{cell_type}_edt_threshold" in params and cell_type in self.edt_thresholds:
                self.edt_thresholds[cell_type].value = float(params[f"{cell_type}_edt_threshold"])
            if f"{cell_type}_segment_size_min" in params and cell_type in self.segment_size_mins:
                self.segment_size_mins[cell_type].value = int(params[f"{cell_type}_segment_size_min"])
            if f"{cell_type}_opening_nr_pixels" in params and cell_type in self.opening_nr_pixels:
                self.opening_nr_pixels[cell_type].value = int(params[f"{cell_type}_opening_nr_pixels"])
            if f"{cell_type}_fill_holes" in params and cell_type in self.fill_holes:
                self.fill_holes[cell_type].value = bool(params[f"{cell_type}_fill_holes"])
        # Persist updated values to config file
        self._persist_params()
        print("Parameters updated from training and saved to config")
    
    def _on_train_clicked(self, _):
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                odir = Path(self.metadata_loader.output_dir).expanduser()
                odir.mkdir(parents=True, exist_ok=True)
                self._lock(True)
                self.spinner_train.layout.display = None
                
                # Collect initial parameters from widget
                initial_params = self._get_current_params()
                
                print("Starting pixel classifier training...")
                print("Napari viewer will open. Close it when done labeling.")
                
                ret = train_pixel_classifier(
                    output_dir=str(odir),
                    metadata=self.metadata_loader.metadata,
                    examples_per_sample=int(self.examples_per_sample.value),
                    sample_specific_classifier=bool(self.sample_specific_classifier.value),
                    n_workers=int(self.n_workers.value),
                    organoid_types=self.organoid_types,
                    immune_types=self.immune_types,
                    other_types=self.other_types,
                    initial_params=initial_params,
                    on_params_changed=self._update_widgets_from_params,
                )
                self._viewer = ret if (ret is not None and hasattr(ret, "close")) else None
                self.spinner_train.layout.display = "none"
                self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"
                if self._viewer is None:
                    self.out.clear_output()
                    print("✅ Training finished.")
            except Exception: traceback.print_exc()
            finally: self._lock(False)

    def _on_close_clicked(self, _):
        with self.out:
            try:
                if self._viewer is not None:
                    self._viewer.close()
                    self._viewer = None
            finally:
                self.close_button.layout.display = "none"
                self.out.clear_output()
                print("✅ Training finished.")

    def _on_apply_clicked(self, _, only_segment=False):
        self._lock(True)
        with self.out:
            self.out.clear_output()
            self._persist_params()
            try:
                odir = Path(self.metadata_loader.output_dir).expanduser()
                tpr = _mk_timepoint_range(bool(self.use_all_timepoints.value), int(self.tp_start.value), int(self.tp_end.value))
                
                clf_organoid_paths = {ct: str(self.clf_paths[ct].value) for ct in self.organoid_types if ct in self.clf_paths and self.clf_paths[ct].value} if self.manual_clf_paths.value else None
                clf_immune_paths = {ct: str(self.clf_paths[ct].value) for ct in self.immune_types if ct in self.clf_paths and self.clf_paths[ct].value} if self.manual_clf_paths.value else None
                clf_other_paths = {ct: str(self.clf_paths[ct].value) for ct in self.other_types if ct in self.clf_paths and self.clf_paths[ct].value} if self.manual_clf_paths.value else None
                clf_death_path = str(self.clf_death_path.value) if (self.manual_clf_paths.value and self.has_death and self.clf_death_path.value) else None

                self.spinner_apply.layout.display = None
                new_md = run_pixel_classifier_segmentation(
                    output_dir=str(odir),
                    metadata=self.metadata_loader.metadata,
                    organoid_edt_thresholds={ct: float(self.edt_thresholds[ct].value) for ct in self.organoid_types if ct in self.edt_thresholds},
                    immune_edt_thresholds={ct: float(self.edt_thresholds[ct].value) for ct in self.immune_types if ct in self.edt_thresholds},
                    other_edt_thresholds={ct: float(self.edt_thresholds[ct].value) for ct in self.other_types if ct in self.edt_thresholds},
                    organoid_segment_size_mins={ct: int(self.segment_size_mins[ct].value) for ct in self.organoid_types if ct in self.segment_size_mins},
                    immune_segment_size_mins={ct: int(self.segment_size_mins[ct].value) for ct in self.immune_types if ct in self.segment_size_mins},
                    other_segment_size_mins={ct: int(self.segment_size_mins[ct].value) for ct in self.other_types if ct in self.segment_size_mins},
                    organoid_opening_nr_pixels={ct: int(self.opening_nr_pixels[ct].value) for ct in self.organoid_types if ct in self.opening_nr_pixels},
                    immune_opening_nr_pixels={ct: int(self.opening_nr_pixels[ct].value) for ct in self.immune_types if ct in self.opening_nr_pixels},
                    other_opening_nr_pixels={ct: int(self.opening_nr_pixels[ct].value) for ct in self.other_types if ct in self.opening_nr_pixels},
                    organoid_fill_holes={ct: bool(self.fill_holes[ct].value) for ct in self.organoid_types if ct in self.fill_holes},
                    immune_fill_holes={ct: bool(self.fill_holes[ct].value) for ct in self.immune_types if ct in self.fill_holes},
                    other_fill_holes={ct: bool(self.fill_holes[ct].value) for ct in self.other_types if ct in self.fill_holes},
                    timepoint_range=tpr,
                    clf_organoid_paths=clf_organoid_paths,
                    clf_immune_paths=clf_immune_paths,
                    clf_other_paths=clf_other_paths,
                    clf_death_path=clf_death_path,
                    only_segment=bool(only_segment),
                    overwrite_existing=bool(self.overwrite_existing.value),
                    n_workers=int(self.n_workers.value),
                )
                if new_md is not None:
                    self.metadata_loader.metadata = new_md
                    new_md.to_csv(self.metadata_loader.metadata_csv_path, index=False)
                print("✅ Apply finished.")
            except Exception: traceback.print_exc()
            finally:
                self.spinner_apply.layout.display = "none"
                self._lock(False)
class CellposeChannelConfigPanel:
    """Panel for configuring channel labels for Cellpose segmentation.
    
    This panel calculates the number of channels from metadata cell types
    and allows users to configure channel labels either:
    - Same for all samples (global configuration)
    - Per-sample configuration
    
    Channel labels are saved to the config file under the 'cellpose' section.
    """
    
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self._detect_cell_types()
        self._calculate_number_of_channels()
        
        # Load saved config
        cellpose_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {})
        
        # Display calculated number of channels
        self.n_channels_display = widgets.HTML(
            value=f"<b>Number of channels:</b> {self.n_channels} "
                  f"(Organoids: {len(self.organoid_types)}, "
                  f"Immune: {len(self.immune_types)}, "
                  f"Other: {len(self.other_types)}"
                  f"{', Dead: 1' if self.has_death else ''})"
        )
        
        # Mode selection: same for all vs per-sample
        saved_mode = cellpose_cfg.get("labels_mode", "same_for_all")
        self.labels_mode = widgets.RadioButtons(
            options=[
                ('Same labels for all samples', 'same_for_all'),
                ('Configure labels per sample', 'per_sample')
            ],
            value=saved_mode,
            description='Configuration mode:',
            style={'description_width': '140px'}
        )
        
        # Container for dynamic channel label UI
        self.channel_config_container = widgets.VBox()
        
        # Save button
        self.btn_save = widgets.Button(
            description='Save Channel Configuration',
            button_style='success',
            icon='save'
        )
        self.btn_save.on_click(self._on_save)
        
        # Output for status messages
        self.out = widgets.Output()
        
        # Setup observers
        self.labels_mode.observe(self._on_mode_change, names='value')
        
        # Build initial UI
        self._build_channel_config_ui()
        
        # Main UI
        self.ui = widgets.VBox([
            widgets.HTML('<h3>Cellpose Channel Configuration</h3>'),
            self.n_channels_display,
            widgets.HTML('<hr>'),
            self.labels_mode,
            self.channel_config_container,
            self.btn_save,
            self.out
        ])
    
    def _detect_cell_types(self):
        """Detect cell types from metadata"""
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel
        )
        metadata = self.metadata_loader.metadata
        self.organoid_types = detect_organoid_types_from_metadata(metadata)
        self.immune_types = detect_immune_cell_types_from_metadata(metadata)
        self.other_types = detect_other_cell_types_from_metadata(metadata)
        self.has_death = has_dead_channel(metadata)
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types
    
    def _calculate_number_of_channels(self):
        """Calculate total number of channels from cell types"""
        self.n_channels = len(self.organoid_types) + len(self.immune_types) + len(self.other_types)
        if self.has_death:
            self.n_channels += 1
    
    def _get_channel_options(self):
        """Get available channel label options"""
        options = list(self.organoid_types) + list(self.immune_types) + list(self.other_types)
        if self.has_death:
            options.append('dead')
        return options
    
    def _get_sample_names(self):
        """Get list of sample names from metadata"""
        return list(self.metadata_loader.metadata['sample_name'].unique())
    
    def _on_mode_change(self, change):
        """Handle mode change between same_for_all and per_sample"""
        self._build_channel_config_ui()
    
    def _build_channel_config_ui(self):
        """Build the channel configuration UI based on current mode"""
        cellpose_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {})
        mode = self.labels_mode.value
        
        self.channel_dropdowns = {}
        children = []
        
        if self.n_channels == 0:
            children.append(widgets.HTML('<p style="color: orange;">No channels detected. Please check your metadata configuration.</p>'))
            self.channel_config_container.children = children
            return
        
        channel_options = self._get_channel_options()
        
        if mode == 'same_for_all':
            children.append(widgets.HTML('<h4>Set channel labels (applies to all samples)</h4>'))
            
            # Load saved global channel labels
            saved_labels = cellpose_cfg.get("channel_labels", {})
            
            self.channel_dropdowns['global'] = {}
            for i in range(self.n_channels):
                # Default: try to use saved value, or match index to option
                saved_val = saved_labels.get(str(i), saved_labels.get(i, None))
                if saved_val and saved_val in channel_options:
                    default_val = saved_val
                elif i < len(channel_options):
                    default_val = channel_options[i]
                else:
                    default_val = channel_options[0] if channel_options else ''
                
                dropdown = widgets.Dropdown(
                    options=channel_options,
                    value=default_val,
                    description=f'Channel {i}:',
                    style={'description_width': '100px'},
                    layout=widgets.Layout(width='300px')
                )
                self.channel_dropdowns['global'][i] = dropdown
                children.append(dropdown)
        
        else:  # per_sample
            children.append(widgets.HTML('<h4>Set channel labels per sample</h4>'))
            
            # Load saved per-sample labels
            saved_per_sample = cellpose_cfg.get("per_sample_channel_labels", {})
            
            sample_names = self._get_sample_names()
            
            for sample_name in sample_names:
                children.append(widgets.HTML(f'<b>{sample_name}</b>'))
                self.channel_dropdowns[sample_name] = {}
                
                # Get saved labels for this sample
                sample_saved = saved_per_sample.get(sample_name, {})
                
                sample_row = []
                for i in range(self.n_channels):
                    saved_val = sample_saved.get(str(i), sample_saved.get(i, None))
                    if saved_val and saved_val in channel_options:
                        default_val = saved_val
                    elif i < len(channel_options):
                        default_val = channel_options[i]
                    else:
                        default_val = channel_options[0] if channel_options else ''
                    
                    dropdown = widgets.Dropdown(
                        options=channel_options,
                        value=default_val,
                        description=f'Ch {i}:',
                        style={'description_width': '50px'},
                        layout=widgets.Layout(width='180px')
                    )
                    self.channel_dropdowns[sample_name][i] = dropdown
                    sample_row.append(dropdown)
                
                children.append(widgets.HBox(sample_row))
        
        self.channel_config_container.children = children
    
    def _persist_params(self):
        """Save channel configuration to config file"""
        cellpose_cfg = self.metadata_loader.behav3d_parameters.setdefault("cellpose", {})
        
        cellpose_cfg["number_of_channels"] = self.n_channels
        cellpose_cfg["labels_mode"] = self.labels_mode.value
        
        if self.labels_mode.value == 'same_for_all':
            # Save global channel labels
            channel_labels = {}
            if 'global' in self.channel_dropdowns:
                for i, dropdown in self.channel_dropdowns['global'].items():
                    channel_labels[i] = dropdown.value
            cellpose_cfg["channel_labels"] = channel_labels
            cellpose_cfg["per_sample_channel_labels"] = {}
        else:
            # Save per-sample channel labels
            per_sample_labels = {}
            for sample_name, dropdowns in self.channel_dropdowns.items():
                if sample_name != 'global':
                    per_sample_labels[sample_name] = {i: dd.value for i, dd in dropdowns.items()}
            cellpose_cfg["per_sample_channel_labels"] = per_sample_labels
            cellpose_cfg["channel_labels"] = {}
        
        # Save to file
        if hasattr(self.metadata_loader, "behav3d_parameters_path"):
            yaml.safe_dump(
                self.metadata_loader.behav3d_parameters,
                self.metadata_loader.behav3d_parameters_path.open("w"),
                sort_keys=False
            )
    
    def _on_save(self, btn):
        """Save channel configuration"""
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                print("Channel configuration saved successfully!")
                
                # Show summary
                mode = self.labels_mode.value
                if mode == 'same_for_all' and 'global' in self.channel_dropdowns:
                    labels = {i: dd.value for i, dd in self.channel_dropdowns['global'].items()}
                    print(f"   Mode: Same for all samples")
                    print(f"   Labels: {labels}")
                else:
                    print(f"   Mode: Per-sample configuration")
                    for sample_name, dropdowns in self.channel_dropdowns.items():
                        if sample_name != 'global':
                            labels = {i: dd.value for i, dd in dropdowns.items()}
                            print(f"   {sample_name}: {labels}")
            except Exception:
                traceback.print_exc()
    
    def get_channel_labels(self, sample_name=None):
        """Get channel labels for a specific sample or global labels.
        
        This method can be called by other widgets/functions to retrieve
        the configured channel labels.
        
        Args:
            sample_name: If provided and mode is per_sample, returns labels for that sample.
                        If None or mode is same_for_all, returns global labels.
        
        Returns:
            dict: Mapping of channel index to label name
        """
        cellpose_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {})
        mode = cellpose_cfg.get("labels_mode", "same_for_all")
        
        if mode == 'same_for_all' or sample_name is None:
            return cellpose_cfg.get("channel_labels", {})
        else:
            per_sample = cellpose_cfg.get("per_sample_channel_labels", {})
            return per_sample.get(sample_name, cellpose_cfg.get("channel_labels", {}))
