import ipywidgets as widgets
from IPython.display import display, HTML
from pathlib import Path
import os
import pandas as pd

# Configure PyOpenCL (cache off, build chatter silenced).
# MUST RUN BEFORE pyopencl or apoc ARE IMPORTED
from behav3d.core.opencl_env import configure_pyopencl
configure_pyopencl()
import yaml
import traceback
from functools import partial
from .utils import (
    _cfg_get,
    _mk_timepoint_range,
    detect_cell_type_category,
    ExternalImageImporter,
    make_metadata_callback,
    PathPicker,
    spinning_loader,
)
from behav3d.preprocessing.segmentation.napari_pixelclassifier import (
    train_pixel_classifier,
    run_pixel_classifier_segmentation,
    SIGMA_UM,
    format_sigma_table,
)
import numpy as np
import napari


def _filter_merge_types(cell_types):
    """Remove derived output cell types (*_merged, *_grouped) from training lists.
    Keeps *_N_multicolor types since those are real per-channel inputs that
    need segmentation and tracking before being merged."""
    from behav3d.core.metadata import is_combined_multicolor_celltype
    return [ct for ct in cell_types
            if not is_combined_multicolor_celltype(ct)]
from ipyfilechooser import FileChooser
from behav3d.preprocessing.segmentation.cellpose_prediction import (
    run_cellpose_and_sync_metadata,
    run_otsu_threshold_segmentation_from_zarr,
    run_otsu_and_sync_metadata,
)
from behav3d.io.images import load_zarr
from behav3d.core.metadata import (
    has_dead_channel,
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    is_multicolor_celltype,
    multicolor_base_name,
)

# Cellpose channel config: extra raw C indices with no analyzed label
CELLPOSE_CHANNEL_UNUSED = "(unused)"


def _normalize_segmentation_labels_for_view(label_arr, raw_img_shape, layer_name="labels"):
    """
    Normalize a label array to the raw-image viewer shape (T, Z, Y, X).

    Raw images in the viewer are displayed per channel as ``img_arr[:, ch]``,
    so labels must align to ``TZYX`` even when the stored array still carries a
    singleton channel axis.
    """
    raw_shape = tuple(int(v) for v in raw_img_shape)
    label_shape = tuple(int(v) for v in label_arr.shape)

    if len(raw_shape) != 5:
        raise ValueError(f"{layer_name}: expected raw image shape (T,C,Z,Y,X), got {raw_shape}.")

    target_shape = (raw_shape[0],) + tuple(raw_shape[-3:])

    if label_shape == target_shape:
        return label_arr

    if len(label_shape) == 5:
        if label_shape[0] != target_shape[0] or tuple(label_shape[-3:]) != tuple(target_shape[1:]):
            raise ValueError(
                f"{layer_name}: expected label shape {target_shape} or {target_shape[:1] + (1,) + target_shape[1:]}, "
                f"got {label_shape}."
            )
        if int(label_shape[1]) != 1:
            raise ValueError(
                f"{layer_name}: cannot align multi-channel labels with shape {label_shape} to viewer shape {target_shape}."
            )
        return label_arr[:, 0, ...]

    if len(label_shape) != 4:
        raise ValueError(
            f"{layer_name}: expected 4D/5D label array compatible with viewer shape {target_shape}, got {label_shape}."
        )

    raise ValueError(
        f"{layer_name}: expected label shape {target_shape}, got {label_shape}."
    )


def _build_notebook_annotation_legend(label_map, has_death=False):
    """Generate an HTML table displaying the ConvPaint annotation legend."""
    if not label_map or "celltype_to_label" not in label_map:
        return
        
    BACKGROUND_LABEL = 0
    UNKNOWN_LABEL = 1
    CELL_TYPE_COLORS = {
        "tcells": "#00FFFF",  # Cyan
        "nkcells": "#FF00FF", # Magenta
        "macrophages": "#00FF00", # Green
        "tumor": "#FF0000", # Red
        "organoid": "#FF0000", # Red
    }
    
    html = [
        "<div style='margin-top: 10px; margin-bottom: 10px; font-family: sans-serif;'>",
        "<h4>ConvPaint Annotation Legend</h4>",
        "<table style='border-collapse: collapse; width: 300px; text-align: left;'>",
        "<tr style='border-bottom: 1px solid #ddd;'><th>Index</th><th>Color</th><th>Cell Type</th></tr>"
    ]
    
    bg_label = int(label_map.get("background_label", BACKGROUND_LABEL))
    html.append(
        f"<tr><td>{bg_label}</td>"
        f"<td style='background-color: #000000; width: 20px;'></td>"
        f"<td>Background</td></tr>"
    )
    
    ordered_cts = sorted(list(label_map["celltype_to_label"].keys()))
    for ct in ordered_cts:
        idx = int(label_map["celltype_to_label"][ct])
        # Find base name for coloring (e.g. organoid1 -> organoid)
        base_ct = ''.join([i for i in ct if not i.isdigit()]).lower()
        color = CELL_TYPE_COLORS.get(base_ct, "#888888")
        html.append(
            f"<tr><td>{idx}</td>"
            f"<td style='background-color: {color}; width: 20px;'></td>"
            f"<td>{ct}</td></tr>"
        )
        
    if has_death:
        death_label = int(label_map.get("death_label", UNKNOWN_LABEL))
        html.append(
            f"<tr><td>{death_label}</td>"
            f"<td style='background-color: #FFFFFF; width: 20px; border: 1px solid #000;'></td>"
            f"<td>Death</td></tr>"
        )
        
    html.append("</table></div>")
    display(HTML("".join(html)))


class PixelClassifierPanel:
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        
        self._detect_cell_types()
        
        print(f"Detected organoid types: {self.organoid_types}")
        print(f"Detected immune cell types: {self.immune_types}")
        print(f"Detected other cell types: {self.other_types}")
        print(f"Dead channel present: {self.has_death}")

        # Check for APOC availability
        try:
            import apoc
            import pyopencl as cl
            import pyclesperanto_prototype as cle
            APOC_AVAILABLE = len(cl.get_platforms()) > 0
            _gpu_devices = cle.available_device_names()
        except Exception:
            APOC_AVAILABLE = False
            _gpu_devices = []

        # Check for ConvPaint availability
        try:
            from napari_convpaint import ConvpaintModel as _CP  # noqa: F401
            CONVPAINT_AVAILABLE = True
        except Exception:
            CONVPAINT_AVAILABLE = False

        _engine_options = ["Scikit-Learn (CPU)"]
        if APOC_AVAILABLE:
            _engine_options.append("APOC (GPU-Accelerated)")
        if CONVPAINT_AVAILABLE:
            _engine_options.append("ConvPaint")

        _saved_engine = pc.get("classifier_engine", "Scikit-Learn (CPU)")
        # Guard: if the saved engine is no longer available, fall back.
        if _saved_engine not in _engine_options:
            _saved_engine = _engine_options[0]

        self.classifier_engine = widgets.Dropdown(
            options=_engine_options,
            value=_saved_engine,
            description="Engine:"
        )
        if not APOC_AVAILABLE:
            self.classifier_engine.tooltip = "APOC (GPU) disabled: 'apoc' and 'pyopencl' required."
        if not CONVPAINT_AVAILABLE:
            self.classifier_engine.tooltip = (self.classifier_engine.tooltip or "") + " ConvPaint disabled: 'napari-convpaint' required."

        _saved_gpu = pc.get("gpu_device_name", _gpu_devices[0] if _gpu_devices else "")
        if _gpu_devices and _saved_gpu not in _gpu_devices:
            _saved_gpu = _gpu_devices[0]

        self.gpu_device = widgets.Dropdown(
            options=_gpu_devices,
            value=_saved_gpu,
            description="GPU Device:",
            layout=widgets.Layout(width="auto", display="none")  # Hidden by default
        )
        self.gpu_device.observe(self._on_gpu_changed, names="value")

        _apoc_strategy_options = [
            "APOC (Direct Instance Segmentation)",
            "APOC Mask + EDT/Watershed Resegmentation",
            "APOC Mask + Peak EDT/Watershed Resegmentation",
            "APOC Probability Map + Watershed",
        ]
        import os
        if os.environ.get("BEHAV3D_DEV_MODE") != "1":
            _apoc_strategy_options.remove("APOC Mask + Peak EDT/Watershed Resegmentation")
        _saved_apoc_strategy = pc.get("apoc_strategy", "APOC (Direct Instance Segmentation)")
        if _saved_apoc_strategy not in _apoc_strategy_options:
            _saved_apoc_strategy = "APOC (Direct Instance Segmentation)"
        self.apoc_strategy = widgets.Dropdown(
            options=_apoc_strategy_options,
            value=_saved_apoc_strategy,
            description="Strategy:",
            style={'description_width': 'initial'},
            layout=widgets.Layout(width="auto", display="none")
        )
        self.apoc_strategy.observe(self._toggle_engine_widgets, names="value")

        # Migrate legacy Direct value (no longer offered) to the EDT default.
        _convpaint_strategy_options = [
            "ConvPaint Mask + EDT/Watershed",
            "ConvPaint Mask + Peak EDT/Watershed",
            "ConvPaint Probability + Watershed",
        ]
        import os
        if os.environ.get("BEHAV3D_DEV_MODE") != "1":
            _convpaint_strategy_options.remove("ConvPaint Mask + Peak EDT/Watershed")
        _saved_cp_strategy = pc.get(
            "convpaint_strategy", "ConvPaint Mask + EDT/Watershed"
        )
        if _saved_cp_strategy not in _convpaint_strategy_options:
            _saved_cp_strategy = "ConvPaint Mask + EDT/Watershed"
        self.convpaint_strategy = widgets.Dropdown(
            options=_convpaint_strategy_options,
            value=_saved_cp_strategy,
            description="Strategy:",
            style={'description_width': 'initial'},
            layout=widgets.Layout(width="auto", display="none")
        )
        self.convpaint_strategy.observe(self._toggle_engine_widgets, names="value")


        self.examples_per_sample = widgets.IntText(
            description="Examples per sample",
            value=int(pc.get("examples_per_sample", 3))
        )
        self.overwrite_sample_images = widgets.Checkbox(
            description="Overwrite sample images",
            value=False,
            indent=False,
            style={"description_width": "auto"},
        )
        '''self.sample_specific_classifier = widgets.Checkbox(
            description="Sample-specific classifier",
            value=bool(pc.get("sample_specific_classifier", False))
        )'''
        
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

        # Unified ConvPaint model picker (single .pkl for all cell types).
        # Only shown when the engine is ConvPaint.
        self.clf_unified_path = PathPicker(
            mode='file',
            description='Unified ConvPaint clf:',
            default="",
            filter_pattern='*.pkl',
            description_width='160px',
            width='100%',
        )
        # ConvPaint death model is also a .pkl (separate binary classifier).
        self.clf_death_convpaint_path = PathPicker(
            mode='file',
            description='Death ConvPaint clf:',
            default="",
            filter_pattern='*.pkl',
            description_width='160px',
            width='100%',
        )

        if self.manual_clf_paths.value:
            self.clf_dir.value = str(pc.get("clf_dir", "") or "")
            self.clf_death_path.value = str(pc.get("clf_death_path", "") or "")
            self.clf_unified_path.value = str(pc.get("clf_unified_path", "") or "")
            self.clf_death_convpaint_path.value = str(
                pc.get("clf_death_convpaint_path", "") or ""
            )
            for cell_type in self.all_cell_types:
                if cell_type in self.clf_paths:
                    saved_path = pc.get(f"clf_{cell_type}_path", "")
                    if saved_path:
                        self.clf_paths[cell_type].value = str(saved_path)

        # Per-group physical(µm)/pixel unit toggles. "Segmentation
        # parameters" (EDT threshold, min size, peak distance) share one
        # toggle; the separate "post-processing size filtering" section
        # (apoc_min_size_inputs) gets its own, matching the napari plugin's
        # per-group design (behav3d.napari._units.UnitGroupManager).
        from .utils import UnitGroupManager
        self._seg_unit_mgr = UnitGroupManager(self.metadata_loader.metadata, default_physical=True)
        self._apoc_size_unit_mgr = UnitGroupManager(self.metadata_loader.metadata, default_physical=True)

        self.edt_thresholds = {}
        self._create_edt_threshold_inputs()
        self._create_postprocessing_inputs()
        self._create_peak_inputs()
        self._create_probability_threshold_inputs()
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
            description="Run batch segmentation",
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

        ct_options = self.all_cell_types
        if self.has_death:
            ct_options = ct_options + ["death"]

        self.segment_ct_dropdown = widgets.Dropdown(
            options=ct_options,
            description="Cell type:",
            layout=widgets.Layout(width="250px")
        )
        self.btn_run_single = widgets.Button(
            description="Run for selected cell type",
            button_style="info",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )

        self.single_apply_row = widgets.HBox(
            [self.segment_ct_dropdown, self.btn_run_single],
            layout=widgets.Layout(align_items="center", gap="8px", margin="0 0 10px 0")
        )

        self.apply_row = widgets.HBox(
            [self.btn_run, self.btn_resegment, self.spinner_apply],
            layout=widgets.Layout(align_items="center", gap="8px")
        )
        
        self.apply_container = widgets.VBox(
            [self.single_apply_row, self.apply_row],
            layout=widgets.Layout(margin="10px 0")
        )

        self.only_resegment_warning = widgets.HTML(
            "<span style='color:#b26a00;'>"
            "Warning: Only resegment must use the same timepoint range as the masks/segmentation it comes from. "
            "If you want a different range, run full segmentation again."
            "</span>"
        )

        self.btn_train.on_click(self._on_train_clicked)
        self.close_button.on_click(self._on_close_clicked)
        self.btn_run.on_click(partial(self._on_apply_clicked, only_segment=False, target_cell_type=None))
        self.btn_resegment.on_click(partial(self._on_apply_clicked, only_segment=True, target_cell_type=None))
        self.btn_run_single.on_click(lambda _: self._on_apply_clicked(None, only_segment=False, target_cell_type=self.segment_ct_dropdown.value))

        # --- Post-processing: size filtering ---
        # We create a list of min size inputs, one per cell type
        self.apoc_min_size_inputs = {}
        for ct in self.all_cell_types:
            # Re-use values from CPU params if available, else standard defaults
            if ct in self.organoid_types: default_size = 1000
            else: default_size = 10
            
            saved_size = pc.get(f"apoc_{ct}_min_size_voxels", default_size)
            self.apoc_min_size_inputs[ct] = widgets.BoundedIntText(
                description=f"{ct.capitalize()}:",
                value=int(saved_size),
                min=0,
                max=1000000000,
                step=10,
                style={'description_width': '100px'},
                layout=widgets.Layout(width='200px'),
            )
            self._apoc_size_unit_mgr.register(
                self.apoc_min_size_inputs[ct], "volume", saved_size,
            )

        self.apoc_size_inputs_box = widgets.VBox(
            [self._apoc_size_unit_mgr.header_row("Units:")]
            + list(self.apoc_min_size_inputs.values())
        )

        self.btn_size_filter = widgets.Button(
            description="Run size filtering",
            button_style="warning",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto"),
        )
        self.spinner_filter = widgets.HTML(value=spinning_loader)
        self.spinner_filter.layout.display = "none"
        
        # New row for size filtering: inputs on the left, button on the right
        self.filter_row = widgets.HBox(
            [self.apoc_size_inputs_box, self.btn_size_filter, self.spinner_filter],
            layout=widgets.Layout(align_items="flex-end", gap="12px"),
        )
        self.btn_size_filter.on_click(self._on_size_filter_clicked)

        self.post_processing_box = widgets.VBox([
            widgets.HTML("<b>Post-processing</b>"),
            widgets.HTML(
                "<span style='color:#666;font-size:0.9em;'>"
                "Remove segments smaller than threshold (voxels). Holes left by removed embedded objects are filled."
                "</span>"
            ),
            self.filter_row,
        ])

        self.out = widgets.Output()

        self.clf_paths_box = widgets.VBox()
        self._build_clf_paths_box()
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
        segmentation_widgets = [
            widgets.HTML("<b>Segmentation parameters per cell type</b>"),
            self._seg_unit_mgr.header_row("Units:"),
        ]

        def _add_cell_type_rows(ct_list):
            for ct in ct_list:
                if ct not in self.edt_thresholds:
                    continue
                segmentation_widgets.append(widgets.HBox([
                    self.edt_thresholds[ct],
                    self.segment_size_mins[ct],
                ]))
                segmentation_widgets.append(widgets.HBox([
                    self.opening_nr_pixels[ct],
                    self.fill_holes[ct],
                ]))
                # Peak EDT params (shown/hidden via _toggle_engine_widgets)
                peak_row = widgets.HBox([
                    self.peak_min_distances[ct],
                    self.peak_min_ratios[ct],
                ])
                peak_row._behav3d_peak_row = True
                segmentation_widgets.append(peak_row)

        if self.organoid_types:
            segmentation_widgets.append(widgets.HTML("<i>Organoids:</i>"))
            _add_cell_type_rows(self.organoid_types)
        if self.immune_types:
            segmentation_widgets.append(widgets.HTML("<i>Immune cells:</i>"))
            _add_cell_type_rows(self.immune_types)
        if self.other_types:
            segmentation_widgets.append(widgets.HTML("<i>Other cells:</i>"))
            _add_cell_type_rows(self.other_types)

        # Store CPU/EDT widget group so it can be reused for CPU and APOC EDT workflows
        self._cpu_seg_params_box = widgets.VBox(segmentation_widgets)

        # Build probability map parameter widgets
        prob_widgets = [widgets.HTML("<b>Probability Map + Watershed parameters per cell type</b>")]
        for ct in self.all_cell_types:
            if ct in self.prob_mask_thresholds:
                prob_widgets.append(widgets.HBox([
                    self.prob_mask_thresholds[ct],
                    self.prob_seed_thresholds[ct],
                ]))
                prob_widgets.append(widgets.HBox([
                    self.segment_size_mins[ct],
                    self.opening_nr_pixels[ct],
                ]))
        self._prob_params_box = widgets.VBox(prob_widgets)

        self.ui = widgets.VBox([
            widgets.VBox([
                widgets.HTML("<b>Train pixel classifier</b>"),
                widgets.HBox([self.classifier_engine, self.gpu_device, self.apoc_strategy, self.convpaint_strategy]),
                widgets.HBox([self.examples_per_sample, self.overwrite_sample_images, self.n_workers]),
                #self.sample_specific_classifier,
                self.train_row,
            ]),
            widgets.HTML("<hr>"),
            widgets.VBox([
                widgets.HTML("<b>Apply segmentation</b>"),
                self.manual_clf_paths,
                self.clf_paths_box,
                self._cpu_seg_params_box,
                self._prob_params_box,
                self.overwrite_existing,
                widgets.HBox([self.use_all_timepoints, self.tp_start, self.tp_end]),
                self.only_resegment_warning,
                self.apply_container,
            ]),
            widgets.HTML("<hr>"),
            self.post_processing_box,
            widgets.HTML("<hr>"),
            self.out
        ])

        # Engine toggle: hide CPU-only widgets when APOC is selected
        self.classifier_engine.observe(self._toggle_engine_widgets, names='value')
        self._toggle_engine_widgets()  # apply on init

        self._viewer = None

    def _detect_cell_types(self):
        import re
        from behav3d.io.images import load_image, load_zarr, save_as_zarr
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel,
        )

        metadata = self.metadata_loader.metadata

        # Prefer explicit cell-type columns found in metadata (this preserves
        # per-channel multicolor names such as 'tcells_1_multicolor').
        def _extract_types_for_prefix(prefix):
            types = set()
            if metadata is None:
                return []
            suffixes = (
                "_line_condition",
                "_segments_image_path",
                "_tracks_image_path",
                "_tracks_csv_path",
            )
            for col in metadata.columns:
                if not col.startswith(f"{prefix}_"):
                    continue
                for suf in suffixes:
                    if col.endswith(suf):
                        # strip prefix + trailing suffix to get the full cell-type token
                        ct = col[len(prefix) + 1 : -len(suf)]
                        types.add(ct)
                        break
            return sorted(types)

        # Extract explicit names (will include per-channel multicolor names)
        self.organoid_types = _extract_types_for_prefix('or') or detect_organoid_types_from_metadata(metadata)
        self.immune_types = _extract_types_for_prefix('im') or detect_immune_cell_types_from_metadata(metadata)
        self.other_types = _extract_types_for_prefix('ot') or detect_other_cell_types_from_metadata(metadata)

        # Remove derived output types (*_merged, *_grouped) — they are generated
        # downstream and should not appear in training or segmentation parameter UI.
        self.organoid_types = _filter_merge_types(self.organoid_types)
        self.immune_types = _filter_merge_types(self.immune_types)
        self.other_types = _filter_merge_types(self.other_types)

        # Backwards-compatible: if explicit columns not present, fall back to detectors
        # Determine death channel presence
        self.has_death = has_dead_channel(metadata)

        # Combine all detected types; keep order organoid -> immune -> other
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types

    def _apply_multicolor_segment_cleanup(self, metadata, output_dir):
        """Run multicolor redundancy cleanup after segmentation.

        This keeps per-channel segmentations separate while removing shared
        voxels that would otherwise be duplicated across multicolor channels.

        Unified-mode ConvPaint already guarantees mutual exclusivity between
        cell-type masks (a voxel can carry only one class), so the cleanup is
        a no-op there and we skip it for performance.
        """
        if str(self.classifier_engine.value) == "ConvPaint":
            print(
                "Skipping multicolor cleanup: unified ConvPaint already "
                "produces mutually-exclusive per-cell-type masks."
            )
            return metadata

        families = {}
        for cell_type in self.organoid_types + self.immune_types + self.other_types:
            if not is_multicolor_celltype(cell_type):
                continue
            base_name = multicolor_base_name(cell_type)
            families.setdefault(base_name, []).append(cell_type)

        if not families:
            return metadata

        from behav3d.preprocessing.segmentation.multicolor_segment_processing import (
            apply_multicolor_segment_correction_for_base,
        )

        cleaned_metadata = metadata
        for base_name, family_cell_types in sorted(families.items()):
            if len(family_cell_types) < 2:
                continue
            print(
                f"Running multicolor cleanup for '{base_name}' after batch segmentation: "
                f"{sorted(family_cell_types)}"
            )
            cleaned_metadata = apply_multicolor_segment_correction_for_base(
                metadata=cleaned_metadata,
                output_dir=str(output_dir),
                base_cell_type=base_name,
                n_channels=len(family_cell_types),
                overwrite=bool(self.overwrite_existing.value),
                n_workers=max(1, int(self.n_workers.value)),
            )

        return cleaned_metadata
    
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
            self.edt_thresholds[cell_type] = widgets.BoundedFloatText(
                description=f"{cell_type.capitalize()} EDT:",
                value=float(saved_threshold),
                min=0,
                max=1000.0,
                step=0.5,
                style={'description_width': '160px'},
            )
            self._seg_unit_mgr.register(
                self.edt_thresholds[cell_type], "distance", saved_threshold,
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
            self.segment_size_mins[cell_type] = widgets.BoundedIntText(
                description=f"{cell_type.capitalize()} min size:",
                value=int(saved_size),
                min=0,
                max=1000000000,
                step=10,
                style={'description_width': '160px'},
            )
            self._seg_unit_mgr.register(
                self.segment_size_mins[cell_type], "volume", saved_size,
            )
            
            # Opening nr pixels
            saved_opening = pc.get(f"{cell_type}_opening_nr_pixels", default_opening)
            self.opening_nr_pixels[cell_type] = widgets.BoundedIntText(
                description=f"{cell_type.capitalize()} opening px:",
                value=int(saved_opening),
                min=0,
                max=10,
                step=1,
                style={'description_width': '160px'},
            )
            
            # Fill holes
            saved_fill = pc.get(f"{cell_type}_fill_holes", default_fill)
            self.fill_holes[cell_type] = widgets.Checkbox(
                description=f"{cell_type.capitalize()} fill holes",
                value=bool(saved_fill),
                indent=False,
            )

    def _create_peak_inputs(self):
        """Create per-cell-type peak_min_distance and peak_min_ratio inputs (Peak EDT strategy)."""
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})

        self.peak_min_distances = {}
        self.peak_min_ratios = {}

        for cell_type in self.all_cell_types:
            saved_dist = pc.get(f"{cell_type}_peak_min_distance", 0.0)
            self.peak_min_distances[cell_type] = widgets.BoundedFloatText(
                description=f"{cell_type.capitalize()} peak dist:",
                value=float(saved_dist),
                min=0.0,
                max=1000.0,
                step=0.5,
                style={'description_width': '160px'},
            )
            self._seg_unit_mgr.register(
                self.peak_min_distances[cell_type], "distance", saved_dist,
            )
            self.peak_min_distances[cell_type].layout.width = 'auto'

            saved_ratio = pc.get(f"{cell_type}_peak_min_ratio", 0.35)
            self.peak_min_ratios[cell_type] = widgets.BoundedFloatText(
                description=f"{cell_type.capitalize()} peak ratio:",
                value=float(saved_ratio),
                min=0.0,
                max=1.0,
                step=0.05,
                style={'description_width': '160px'},
            )
            self.peak_min_ratios[cell_type].layout.width = 'auto'

    def _create_probability_threshold_inputs(self):
        """Create per-cell-type mask and seed threshold inputs for probability map strategy."""
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})

        self.prob_mask_thresholds = {}
        self.prob_seed_thresholds = {}

        for cell_type in self.all_cell_types:
            saved_mask = pc.get(f"{cell_type}_prob_mask_threshold", 0.5)
            self.prob_mask_thresholds[cell_type] = widgets.BoundedFloatText(
                description=f"{cell_type.capitalize()} mask thr:",
                value=float(saved_mask),
                min=0.0,
                max=1.0,
                step=0.05,
                style={'description_width': '160px'},
                layout=widgets.Layout(width='auto'),
            )

            saved_seed = pc.get(f"{cell_type}_prob_seed_threshold", 0.8)
            self.prob_seed_thresholds[cell_type] = widgets.BoundedFloatText(
                description=f"{cell_type.capitalize()} seed thr:",
                value=float(saved_seed),
                min=0.0,
                max=1.0,
                step=0.05,
                style={'description_width': '160px'},
                layout=widgets.Layout(width='auto'),
            )
    
    def _persist_params(self):
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc["classifier_engine"] = str(self.classifier_engine.value)
        pc["examples_per_sample"] = int(self.examples_per_sample.value)
        pc["workers"] = int(self.n_workers.value)
        pc["use_all_timepoints"] = bool(self.use_all_timepoints.value)
        pc["tp_start"] = int(self.tp_start.value)
        pc["tp_end"]   = int(self.tp_end.value)
        pc["manual_clf_paths"] = bool(self.manual_clf_paths.value)
        pc["overwrite_existing"] = bool(self.overwrite_existing.value)
        pc["gpu_device_name"] = str(self.gpu_device.value)
        pc["apoc_strategy"] = str(self.apoc_strategy.value)
        pc["convpaint_strategy"] = str(self.convpaint_strategy.value)
        
        for ct, w in self.apoc_min_size_inputs.items():
            pc[f"apoc_{ct}_min_size_voxels"] = int(round(self._apoc_size_unit_mgr.get_native(w)))

        for cell_type, threshold_widget in self.edt_thresholds.items():
            pc[f"{cell_type}_edt_threshold"] = float(self._seg_unit_mgr.get_native(threshold_widget))

        for cell_type in self.all_cell_types:
            if cell_type in self.segment_size_mins:
                pc[f"{cell_type}_segment_size_min"] = int(round(self._seg_unit_mgr.get_native(self.segment_size_mins[cell_type])))
            if cell_type in self.opening_nr_pixels:
                pc[f"{cell_type}_opening_nr_pixels"] = int(self.opening_nr_pixels[cell_type].value)
            if cell_type in self.fill_holes:
                pc[f"{cell_type}_fill_holes"] = bool(self.fill_holes[cell_type].value)
            if cell_type in self.prob_mask_thresholds:
                pc[f"{cell_type}_prob_mask_threshold"] = float(self.prob_mask_thresholds[cell_type].value)
            if cell_type in self.prob_seed_thresholds:
                pc[f"{cell_type}_prob_seed_threshold"] = float(self.prob_seed_thresholds[cell_type].value)
            if cell_type in self.peak_min_distances:
                pc[f"{cell_type}_peak_min_distance"] = float(self._seg_unit_mgr.get_native(self.peak_min_distances[cell_type]))
            if cell_type in self.peak_min_ratios:
                pc[f"{cell_type}_peak_min_ratio"] = float(self.peak_min_ratios[cell_type].value)
        
        if self.manual_clf_paths.value:
            pc["clf_dir"] = str(self.clf_dir.value or "")
            pc["clf_death_path"] = str(self.clf_death_path.value or "")
            pc["clf_unified_path"] = str(self.clf_unified_path.value or "")
            pc["clf_death_convpaint_path"] = str(
                self.clf_death_convpaint_path.value or ""
            )
            for cell_type, path_picker in self.clf_paths.items():
                pc[f"clf_{cell_type}_path"] = str(path_picker.value or "")

        if hasattr(self.metadata_loader, "behav3d_parameters_path"):
            yaml.safe_dump(
                self.metadata_loader.behav3d_parameters,
                self.metadata_loader.behav3d_parameters_path.open("w"),
                sort_keys=False
            )

    def _toggle_engine_widgets(self, change=None):
        """Show/hide widgets depending on engine and strategy."""
        engine = str(self.classifier_engine.value)
        is_apoc = engine == "APOC (GPU-Accelerated)"
        is_convpaint = engine == "ConvPaint"

        # APOC strategy
        apoc_strategy = str(self.apoc_strategy.value)
        is_apoc_edt = is_apoc and apoc_strategy == "APOC Mask + EDT/Watershed Resegmentation"
        is_apoc_peak = is_apoc and apoc_strategy == "APOC Mask + Peak EDT/Watershed Resegmentation"
        is_apoc_prob = is_apoc and apoc_strategy == "APOC Probability Map + Watershed"

        # ConvPaint strategy
        cp_strategy = str(self.convpaint_strategy.value)
        is_cp_edt = is_convpaint and cp_strategy == "ConvPaint Mask + EDT/Watershed"
        is_cp_peak = is_convpaint and cp_strategy == "ConvPaint Mask + Peak EDT/Watershed"
        is_cp_prob = is_convpaint and cp_strategy == "ConvPaint Probability + Watershed"

        # CPU-like = needs EDT/postprocessing params
        is_cpu_like = (
            (not is_apoc and not is_convpaint)
            or is_apoc_edt or is_apoc_peak
            or is_cp_edt or is_cp_peak
        )
        # Probability params needed?
        needs_prob = is_apoc_prob or is_cp_prob

        hide_cpu = None if is_cpu_like else 'none'

        # CPU segmentation parameters (EDT, min size, opening, fill holes)
        self._cpu_seg_params_box.layout.display = hide_cpu
        self.n_workers.layout.display = (None if not is_apoc and not is_convpaint else 'none')

        # Peak EDT rows: visible only for APOC Peak or ConvPaint Peak
        show_peak = is_apoc_peak or is_cp_peak
        peak_display = None if show_peak else 'none'
        for w in self._cpu_seg_params_box.children:
            if getattr(w, '_behav3d_peak_row', False):
                w.layout.display = peak_display

        # Probability Map parameters
        self._prob_params_box.layout.display = (None if needs_prob else 'none')

        # GPU Device and APOC Strategy selection
        self.gpu_device.layout.display = (None if is_apoc else 'none')
        self.apoc_strategy.layout.display = (None if is_apoc else 'none')
        self.convpaint_strategy.layout.display = (None if is_convpaint else 'none')

        if is_apoc:
            self._apply_gpu_selection()

        # 'Only resegment' button (available for EDT and Probability strategies)
        show_reseg = is_cpu_like or needs_prob
        self.btn_resegment.layout.display = (None if show_reseg else 'none')
        self.only_resegment_warning.layout.display = (None if show_reseg else 'none')

        # Manual classifier paths: hidden for APOC (no manual override supported);
        # visible for ConvPaint (single unified .pkl) and CPU (per-cell-type).
        self.manual_clf_paths.layout.display = ('none' if is_apoc else None)
        if is_apoc:
            self.clf_paths_box.layout.display = 'none'
        else:
            # Rebuild so the picker set matches the active engine.
            self._build_clf_paths_box()
            self.clf_paths_box.layout.display = (
                None if bool(self.manual_clf_paths.value) else 'none'
            )

        # Post-processing (Size Filtering via standalone button)
        # Shown only for Direct APOC; ConvPaint strategies (EDT / Probability)
        # already include size filtering inside their resegmentation pipeline.
        show_post = is_apoc and not is_apoc_edt and not is_apoc_peak and not is_apoc_prob
        self.post_processing_box.layout.display = (None if show_post else 'none')

    def _on_gpu_changed(self, change):
        if change['new']:
            self._apply_gpu_selection()

    def _apply_gpu_selection(self):
        device = str(self.gpu_device.value)
        if device:
            try:
                import pyclesperanto_prototype as cle
                print(f"APOC: Explicitly selecting GPU device: {device}")
                cle.select_device(device)
            except Exception as e:
                print(f"⚠️ Error selecting GPU device: {e}")

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
            self.clf_unified_path.value = ""
            self.clf_death_convpaint_path.value = ""
            for picker in self.clf_paths.values(): picker.value = ""
        else:
            # Refresh the visible pickers in case the engine changed while hidden.
            self._build_clf_paths_box()

    def _build_clf_paths_box(self):
        engine = str(self.classifier_engine.value)
        is_convpaint = engine == "ConvPaint"

        children = [widgets.HTML("<b>Classifier paths</b>")]

        if is_convpaint:
            children.append(widgets.HTML(
                "<span style='color:#666;font-size:0.9em;'>"
                "ConvPaint uses a single unified multi-class classifier. "
                "Death is a separate binary ConvPaint model.</span>"
            ))
            children.append(self.clf_unified_path)
            if self.has_death:
                children.append(widgets.HTML("<i>Death mask:</i>"))
                children.append(self.clf_death_convpaint_path)
        else:
            children.append(self.clf_dir)
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
            self.btn_train, self.apply_container,
            #self.examples_per_sample, self.sample_specific_classifier, self.n_workers,
            self.examples_per_sample, self.n_workers,
            self.use_all_timepoints, self.tp_start, self.tp_end,
            self.manual_clf_paths, self.overwrite_existing,
        ]
        to_lock.extend(self.edt_thresholds.values())
        to_lock.extend(self.segment_size_mins.values())
        to_lock.extend(self.opening_nr_pixels.values())
        to_lock.extend(self.fill_holes.values())
        to_lock.extend(self.prob_mask_thresholds.values())
        to_lock.extend(self.prob_seed_thresholds.values())
        to_lock.extend(self.peak_min_distances.values())
        to_lock.extend(self.peak_min_ratios.values())
        to_lock.extend([self.gpu_device, self.apoc_strategy, self.convpaint_strategy])
        for w in self.apoc_min_size_inputs.values(): w.disabled = state
        for w in [self.btn_size_filter]: w.disabled = state
        for w in to_lock: w.disabled = state
        try:
            self.clf_dir.text.disabled = state
            self.clf_dir.btn.disabled = state
            if self.has_death:
                self.clf_death_path.text.disabled = state
                self.clf_death_path.btn.disabled = state
                self.clf_death_convpaint_path.text.disabled = state
                self.clf_death_convpaint_path.btn.disabled = state
            self.clf_unified_path.text.disabled = state
            self.clf_unified_path.btn.disabled = state
            for picker in self.clf_paths.values():
                picker.text.disabled = state
                picker.btn.disabled = state
        except Exception: pass

    def _get_current_params(self):
        """Collect current parameter values and widget limits into a dict.
        Also includes saved APOC per-cell-type params from the config YAML
        so the APOC widget can restore previous settings on open.
        """
        params = {}
        for cell_type in self.all_cell_types:
            if cell_type in self.edt_thresholds:
                w = self.edt_thresholds[cell_type]
                params[f"{cell_type}_edt_threshold"] = float(self._seg_unit_mgr.get_native(w))
                params[f"{cell_type}_edt_min"] = float(w.min)
                params[f"{cell_type}_edt_max"] = float(w.max)
                params[f"{cell_type}_edt_step"] = float(w.step)
            if cell_type in self.segment_size_mins:
                w = self.segment_size_mins[cell_type]
                params[f"{cell_type}_segment_size_min"] = int(round(self._seg_unit_mgr.get_native(w)))
                params[f"{cell_type}_segment_size_min_limit"] = int(w.min)
                params[f"{cell_type}_segment_size_max_limit"] = int(w.max)
                params[f"{cell_type}_segment_size_step"] = int(w.step)
            if cell_type in self.opening_nr_pixels:
                w = self.opening_nr_pixels[cell_type]
                params[f"{cell_type}_opening_nr_pixels"] = int(w.value)
                params[f"{cell_type}_opening_min"] = int(w.min)
                params[f"{cell_type}_opening_max"] = int(w.max)
                params[f"{cell_type}_opening_step"] = int(w.step)
            if cell_type in self.fill_holes:
                params[f"{cell_type}_fill_holes"] = bool(self.fill_holes[cell_type].value)
            if cell_type in self.prob_mask_thresholds:
                params[f"{cell_type}_prob_mask_threshold"] = float(self.prob_mask_thresholds[cell_type].value)
            if cell_type in self.prob_seed_thresholds:
                params[f"{cell_type}_prob_seed_threshold"] = float(self.prob_seed_thresholds[cell_type].value)

        # Include global GPU selection
        params["gpu_device_name"] = str(self.gpu_device.value)
        params["apoc_strategy"] = str(self.apoc_strategy.value)
        params["convpaint_strategy"] = str(self.convpaint_strategy.value)

        for cell_type in self.all_cell_types:
            if cell_type in self.peak_min_distances:
                params[f"{cell_type}_peak_min_distance"] = float(self._seg_unit_mgr.get_native(self.peak_min_distances[cell_type]))
            if cell_type in self.peak_min_ratios:
                params[f"{cell_type}_peak_min_ratio"] = float(self.peak_min_ratios[cell_type].value)

        # Include saved APOC per-cell-type params from config YAML
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        all_apoc_types = list(self.all_cell_types)
        if self.has_death:
            all_apoc_types.append("dead")
        for ct in all_apoc_types:
            for key in [
                "feature_preset",
                "feature_string",
                "sigmas",
                "grid_sigmas",
                "custom_feature_string",
                "checked_features",
                "consider_original",
                "max_depth",
                "num_ensembles",
                "channels",
            ]:
                saved = pc.get(f"apoc_{ct}_{key}")
                if saved is not None:
                    params[f"apoc_{ct}_{key}"] = saved
        return params
    
    def _update_widgets_from_params(self, params):
        """Update widget values from a params dict and persist to config.
        Accepts both CPU segmentation params and APOC per-cell-type params.
        """
        for cell_type in self.all_cell_types:
            if f"{cell_type}_edt_threshold" in params and cell_type in self.edt_thresholds:
                self._seg_unit_mgr.set_native(self.edt_thresholds[cell_type], float(params[f"{cell_type}_edt_threshold"]))
            if f"{cell_type}_segment_size_min" in params and cell_type in self.segment_size_mins:
                self._seg_unit_mgr.set_native(self.segment_size_mins[cell_type], int(params[f"{cell_type}_segment_size_min"]))
            if f"{cell_type}_opening_nr_pixels" in params and cell_type in self.opening_nr_pixels:
                self.opening_nr_pixels[cell_type].value = int(params[f"{cell_type}_opening_nr_pixels"])
            if f"{cell_type}_fill_holes" in params and cell_type in self.fill_holes:
                self.fill_holes[cell_type].value = bool(params[f"{cell_type}_fill_holes"])
            if f"{cell_type}_prob_mask_threshold" in params and cell_type in self.prob_mask_thresholds:
                self.prob_mask_thresholds[cell_type].value = float(params[f"{cell_type}_prob_mask_threshold"])
            if f"{cell_type}_prob_seed_threshold" in params and cell_type in self.prob_seed_thresholds:
                self.prob_seed_thresholds[cell_type].value = float(params[f"{cell_type}_prob_seed_threshold"])
        if "apoc_strategy" in params:
            apoc_val = str(params["apoc_strategy"])
            if apoc_val in self.apoc_strategy.options:
                self.apoc_strategy.value = apoc_val
        if "convpaint_strategy" in params:
            cp_val = str(params["convpaint_strategy"])
            if cp_val in self.convpaint_strategy.options:
                self.convpaint_strategy.value = cp_val
        for cell_type in self.all_cell_types:
            if f"{cell_type}_peak_min_distance" in params and cell_type in self.peak_min_distances:
                self._seg_unit_mgr.set_native(self.peak_min_distances[cell_type], float(params[f"{cell_type}_peak_min_distance"]))
            if f"{cell_type}_peak_min_ratio" in params and cell_type in self.peak_min_ratios:
                self.peak_min_ratios[cell_type].value = float(params[f"{cell_type}_peak_min_ratio"])

        # Cache APOC and ConvPaint per-cell-type params into the YAML config dict
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        for key, value in params.items():
            if key.startswith("apoc_") or key.startswith("convpaint_"):
                pc[key] = value

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
                
                engine = str(self.classifier_engine.value)
                
                if engine == "APOC (GPU-Accelerated)":
                    # Use the custom APOC widget
                    from behav3d.preprocessing.segmentation.apoc_train import train_pixel_classifier_apoc
                    print("Starting APOC pixel classifier training...")
                    print("Napari viewer will open with the APOC widget.")
                    ret = train_pixel_classifier_apoc(
                        output_dir=str(odir),
                        metadata=self.metadata_loader.metadata,
                        examples_per_sample=int(self.examples_per_sample.value),
                        overwrite_images=bool(self.overwrite_sample_images.value),
                        organoid_types=_filter_merge_types(self.organoid_types),
                        immune_types=_filter_merge_types(self.immune_types),
                        other_types=_filter_merge_types(self.other_types),
                        initial_params=initial_params,
                        on_params_changed=self._update_widgets_from_params,
                    )
                elif engine == "ConvPaint":
                    from behav3d.preprocessing.segmentation.convpaint_train import train_pixel_classifier_convpaint
                    print("Starting ConvPaint pixel classifier training...")
                    print("Napari viewer will open with the ConvPaint widget.")
                    ret = train_pixel_classifier_convpaint(
                        output_dir=str(odir),
                        metadata=self.metadata_loader.metadata,
                        examples_per_sample=int(self.examples_per_sample.value),
                        overwrite_images=bool(self.overwrite_sample_images.value),
                        organoid_types=_filter_merge_types(self.organoid_types),
                        immune_types=_filter_merge_types(self.immune_types),
                        other_types=_filter_merge_types(self.other_types),
                        initial_params=initial_params,
                        on_params_changed=self._update_widgets_from_params,
                    )
                else:
                    # Use the original scikit-learn pixel classifier
                    print("Starting pixel classifier training...")
                    print("Napari viewer will open. Close it when done labeling.")
                    
                    # Print sigma table for the pixel sizes in metadata
                    try:
                        md = self.metadata_loader.metadata
                        px_xy = float(md['pixel_distance_xy'].iloc[0])
                        px_z  = float(md['pixel_distance_z'].iloc[0])
                        print(format_sigma_table(px_xy, px_z))
                    except Exception:
                        pass
                    
                    ret = train_pixel_classifier(
                        output_dir=str(odir),
                        metadata=self.metadata_loader.metadata,
                        examples_per_sample=int(self.examples_per_sample.value),
                        overwrite_images=bool(self.overwrite_sample_images.value),
                        #sample_specific_classifier=bool(self.sample_specific_classifier.value),
                        n_workers=int(self.n_workers.value),
                        organoid_types=_filter_merge_types(self.organoid_types),
                        immune_types=_filter_merge_types(self.immune_types),
                        other_types=_filter_merge_types(self.other_types),
                        initial_params=initial_params,
                        on_params_changed=self._update_widgets_from_params,
                    )
                self._viewer = ret if (ret is not None and hasattr(ret, "close")) else None
                self.spinner_train.layout.display = "none"
                self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"
                if self._viewer is None:
                    self.out.clear_output()
                    print("Training finished.")
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
                print("Training finished.")

    def _print_output_completeness(self, output_dir, only_cell_types=None):
        """Report how complete each existing segmentation output is, before running.

        The napari plugin asks this question in its overwrite dialog; the notebook
        widget has only a checkbox, so the same information is printed instead.
        "Already exists" is not the same as "already finished" - an interrupted run
        leaves a real but partial array, and only the progress journal knows which
        is which.
        """
        from behav3d.preprocessing.segmentation.segment_journal import (
            describe_state, journal_state,
        )

        md = self.metadata_loader.metadata
        if md is None:
            return
        cts = list(only_cell_types) if only_cell_types else list(
            self.organoid_types + self.immune_types + self.other_types
        )
        rows = []
        for sn in md["sample_name"].unique():
            for ct in cts:
                path = Path(output_dir) / "images" / sn / f"{sn}_{ct}_segments.zarr"
                state, n_done, t_total = journal_state(path)
                if state == "missing":
                    continue
                rows.append(f"  - {ct} segments for {sn}: {describe_state(state, n_done, t_total)}")
        if not rows:
            return
        print("Existing outputs:")
        for row in rows:
            print(row)
        if not self.overwrite_existing.value:
            print(
                "  ('Overwrite existing' is off: complete outputs are skipped and "
                "incomplete ones are resumed.)"
            )

    def _on_apply_clicked(self, _, only_segment=False, target_cell_type=None):
        self._lock(True)
        with self.out:
            self.out.clear_output()
            self._persist_params()
            try:
                odir = Path(self.metadata_loader.output_dir).expanduser()
                tpr = _mk_timepoint_range(bool(self.use_all_timepoints.value), int(self.tp_start.value), int(self.tp_end.value))

                # Build only_cell_types filter: None means run all
                only_cell_types = None if target_cell_type is None else [target_cell_type]

                clf_organoid_paths = {ct: str(self.clf_paths[ct].value) for ct in self.organoid_types if ct in self.clf_paths and self.clf_paths[ct].value} if self.manual_clf_paths.value else None
                clf_immune_paths = {ct: str(self.clf_paths[ct].value) for ct in self.immune_types if ct in self.clf_paths and self.clf_paths[ct].value} if self.manual_clf_paths.value else None
                clf_other_paths = {ct: str(self.clf_paths[ct].value) for ct in self.other_types if ct in self.clf_paths and self.clf_paths[ct].value} if self.manual_clf_paths.value else None
                clf_death_path = str(self.clf_death_path.value) if (self.manual_clf_paths.value and self.has_death and self.clf_death_path.value) else None

                self._print_output_completeness(odir, only_cell_types)

                self.spinner_apply.layout.display = None
                
                if str(self.classifier_engine.value) == "APOC (GPU-Accelerated)":
                    from behav3d.preprocessing.segmentation.apoc_segment import run_apoc_segmentation
                    pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
                    new_md = run_apoc_segmentation(
                        output_dir=str(odir),
                        metadata=self.metadata_loader.metadata,
                        metadata_csv_path=str(self.metadata_loader.metadata_csv_path),
                        apoc_config=pc,
                        timepoint_range=tpr,
                        clf_organoid_paths=clf_organoid_paths,
                        clf_immune_paths=clf_immune_paths,
                        clf_other_paths=clf_other_paths,
                        clf_death_path=clf_death_path,
                        only_segment=bool(only_segment),
                        overwrite_existing=bool(self.overwrite_existing.value),
                        n_workers=int(self.n_workers.value),
                        gpu_device=str(self.gpu_device.value),
                        apoc_strategy=str(self.apoc_strategy.value),
                        only_cell_types=only_cell_types,
                    )
                elif str(self.classifier_engine.value) == "ConvPaint":
                    from behav3d.preprocessing.segmentation.convpaint_segment import run_convpaint_segmentation
                    pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})

                    # Unified ConvPaint mode: a single .pkl drives every cell type.
                    clf_unified_path = (
                        str(self.clf_unified_path.value)
                        if self.manual_clf_paths.value and self.clf_unified_path.value
                        else None
                    )
                    clf_death_convpaint_path = (
                        str(self.clf_death_convpaint_path.value)
                        if (self.manual_clf_paths.value and self.has_death
                            and self.clf_death_convpaint_path.value)
                        else None
                    )

                    new_md = run_convpaint_segmentation(
                        output_dir=str(odir),
                        metadata=self.metadata_loader.metadata,
                        metadata_csv_path=str(self.metadata_loader.metadata_csv_path),
                        convpaint_config=pc,
                        timepoint_range=tpr,
                        clf_unified_path=clf_unified_path,
                        clf_death_path=clf_death_convpaint_path,
                        only_segment=bool(only_segment),
                        overwrite_existing=bool(self.overwrite_existing.value),
                        n_workers=int(self.n_workers.value),
                        convpaint_strategy=str(self.convpaint_strategy.value),
                        only_cell_types=only_cell_types,
                    )
                else:
                    # Original scikit-learn pixel classifier
                    try:
                        md = self.metadata_loader.metadata
                        px_xy = float(md['pixel_distance_xy'].iloc[0])
                        px_z  = float(md['pixel_distance_z'].iloc[0])
                        print(format_sigma_table(px_xy, px_z))
                    except Exception:
                        pass
                    
                    new_md = run_pixel_classifier_segmentation(
                        output_dir=str(odir),
                        metadata=self.metadata_loader.metadata,
                        metadata_csv_path=str(self.metadata_loader.metadata_csv_path),
                        organoid_edt_thresholds={ct: float(self._seg_unit_mgr.get_native(self.edt_thresholds[ct])) for ct in self.organoid_types if ct in self.edt_thresholds},
                        immune_edt_thresholds={ct: float(self._seg_unit_mgr.get_native(self.edt_thresholds[ct])) for ct in self.immune_types if ct in self.edt_thresholds},
                        other_edt_thresholds={ct: float(self._seg_unit_mgr.get_native(self.edt_thresholds[ct])) for ct in self.other_types if ct in self.edt_thresholds},
                        organoid_segment_size_mins={ct: int(round(self._seg_unit_mgr.get_native(self.segment_size_mins[ct]))) for ct in self.organoid_types if ct in self.segment_size_mins},
                        immune_segment_size_mins={ct: int(round(self._seg_unit_mgr.get_native(self.segment_size_mins[ct]))) for ct in self.immune_types if ct in self.segment_size_mins},
                        other_segment_size_mins={ct: int(round(self._seg_unit_mgr.get_native(self.segment_size_mins[ct]))) for ct in self.other_types if ct in self.segment_size_mins},
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
                        only_cell_types=only_cell_types,
                    )
                if new_md is not None:
                    new_md = self._apply_multicolor_segment_cleanup(new_md, odir)
                    if isinstance(new_md, dict):
                        # Multicolor cleanup helpers can return path maps; keep metadata as DataFrame.
                        new_md = self.metadata_loader.metadata
                    self.metadata_loader.metadata = new_md
                    new_md.to_csv(self.metadata_loader.metadata_csv_path, index=False)
                print("Apply finished.")
            except Exception: traceback.print_exc()
            finally:
                self.spinner_apply.layout.display = "none"
                self._lock(False)

    def _on_size_filter_clicked(self, _):
        """Remove segments smaller than the user-specified voxel threshold."""
        self._lock(True)
        with self.out:
            self.out.clear_output()
            self.spinner_filter.layout.display = None
            try:
                from behav3d.preprocessing.segmentation.size_filter import filter_segments_by_size

                odir = Path(self.metadata_loader.output_dir).expanduser()
                metadata = self.metadata_loader.metadata
                all_cell_types = self.organoid_types + self.immune_types + self.other_types
                sample_names = metadata['sample_name'].unique()

                # Persist the setting
                pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
                for ct, w in self.apoc_min_size_inputs.items():
                    pc[f"apoc_{ct}_min_size_voxels"] = int(round(self._apoc_size_unit_mgr.get_native(w)))

                if hasattr(self.metadata_loader, "behav3d_parameters_path"):
                    yaml.safe_dump(
                        self.metadata_loader.behav3d_parameters,
                        self.metadata_loader.behav3d_parameters_path.open("w"),
                        sort_keys=False,
                    )

                print(f"🔍 Size filtering: removing small segments and filling holes")

                for sample_name in sample_names:
                    img_dir = odir / "images" / sample_name
                    found_any = False

                    for ct in all_cell_types:
                        seg_path = img_dir / f"{sample_name}_{ct}_segments.zarr"
                        if not seg_path.exists():
                            continue

                        min_size = int(round(self._apoc_size_unit_mgr.get_native(self.apoc_min_size_inputs[ct])))
                        if min_size <= 0:
                            continue

                        found_any = True
                        print(f"  {sample_name} / {ct} (min_size={min_size}):")
                        removed = filter_segments_by_size(
                            seg_path, 
                            min_size, 
                            fill_holes=True
                        )
                        print(f"    → removed {removed} segments/holes")

                    if not found_any:
                        print(f"  ⚠️ No segmentation files found for {sample_name}")

                print("\n✅ Size filtering complete!")

            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_filter.layout.display = "none"
                self._lock(False)


def build_external_seg_import_ui(metadata_loader):
    """Build a standalone widget for importing external segmentation masks (+ dead mask) to Zarr."""
    metadata = metadata_loader.metadata
    sample_names = list(metadata["sample_name"].astype(str))
    out_dir = Path(metadata_loader.output_dir) / "images"

    all_cell_types = (
        detect_organoid_types_from_metadata(metadata)
        + detect_immune_cell_types_from_metadata(metadata)
        + detect_other_cell_types_from_metadata(metadata)
    )

    children = []
    for ct in all_cell_types:
        cat = detect_cell_type_category(ct, metadata)
        prefix = {"organoid": "or", "immune": "im", "other": "ot"}[cat]
        col = f"{prefix}_{ct}_segments_image_path"
        children.append(ExternalImageImporter(
            label=f"{ct} segmentation",
            output_dir=str(out_dir),
            on_converted=make_metadata_callback(metadata_loader, col),
            sample_names=sample_names,
            output_filename="{sample_name}_" + ct + "_segments.zarr",
            is_label_image=True,
        ))

    if has_dead_channel(metadata):
        children.append(ExternalImageImporter(
            label="Dead mask",
            output_dir=str(out_dir),
            on_converted=make_metadata_callback(metadata_loader, "dead_mask_path"),
            sample_names=sample_names,
            output_filename="{sample_name}_mask_dead.zarr",
            is_label_image=True,
        ))

    accordion = widgets.Accordion(children=[widgets.VBox(children)])
    accordion.set_title(0, "Import external segmentation")
    accordion.selected_index = None

    return widgets.VBox([
        widgets.HTML("<b>Import external segmentation</b>"),
        accordion,
    ])


class CellposeChannelConfigPanel:

    """Panel for configuring channel labels for Cellpose segmentation.

    Channel count starts from metadata (organoid / immune / other + optional dead).
    Use **Extra slots** for additional raw C indices (e.g. not listed in metadata).
    Assign **(unused)** when that index should not map to a Cellpose label.
    """

    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self._detect_cell_types()
        self._calculate_base_channels()

        cellpose_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {})
        self.n_extra_slots = int(cellpose_cfg.get("extra_channel_slots", 0) or 0)
        self.n_channels = self.n_base_channels + self.n_extra_slots

        self.n_channels_display = widgets.HTML(value=self._channels_summary_html())

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

        self.channel_config_container = widgets.VBox()

        self.extra_slots_input = widgets.BoundedIntText(
            value=self.n_extra_slots,
            min=0,
            max=64,
            description='Extra slots:',
            style={'description_width': '100px'},
            layout=widgets.Layout(width='260px'),
        )
        self.btn_apply_extra = widgets.Button(
            description='Apply extra channels',
            button_style='info',
            tooltip='Rebuild dropdowns; click Save to write behav3d_parameters.yml',
            layout=widgets.Layout(width='200px'),
        )
        self.btn_apply_extra.on_click(self._on_apply_extra_slots)
        self.extra_help = widgets.HTML(
            '<p style="font-size:12px;color:#555;margin:4px 0 0 0;">'
            'Extra slots add Channel 5, 6, … when the raw image has more C planes than '
            'metadata types. Use <b>(unused)</b> for indices that exist only for ordering '
            'or reference, not for segmentation.</p>'
        )

        self.btn_save = widgets.Button(
            description='Save Channel Configuration',
            button_style='success',
            icon='save'
        )
        self.btn_save.on_click(self._on_save)

        self.out = widgets.Output()

        self.labels_mode.observe(self._on_mode_change, names='value')

        self._build_channel_config_ui()

        self.ui = widgets.VBox([
            widgets.HTML('<h3>Cellpose Channel Configuration</h3>'),
            self.n_channels_display,
            widgets.HTML('<hr>'),
            widgets.HBox([self.extra_slots_input, self.btn_apply_extra]),
            self.extra_help,
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
    
    def _calculate_base_channels(self):
        """Metadata-derived channel count (organoid + immune + other + optional dead)."""
        self.n_base_channels = (
            len(self.organoid_types) + len(self.immune_types) + len(self.other_types)
        )
        if self.has_death:
            self.n_base_channels += 1

    def _channels_summary_html(self):
        extra_txt = f" + {self.n_extra_slots} extra" if self.n_extra_slots else ""
        return (
            f"<p><b>Channel slots:</b> {self.n_channels} total "
            f"({self.n_base_channels} from metadata{extra_txt}) — "
            f"Organoids: {len(self.organoid_types)}, Immune: {len(self.immune_types)}, "
            f"Other: {len(self.other_types)}"
            f"{', Dead: 1' if self.has_death else ''}</p>"
        )

    def _on_apply_extra_slots(self, _):
        self.n_extra_slots = int(self.extra_slots_input.value)
        self.n_channels = self.n_base_channels + self.n_extra_slots
        self.n_channels_display.value = self._channels_summary_html()
        self._build_channel_config_ui()

    def _analytic_channel_options(self):
        """Real cell-type labels only (no (unused) sentinel)."""
        options = list(self.organoid_types) + list(self.immune_types) + list(self.other_types)
        if self.has_death:
            options.append('dead')
        return options

    def _get_channel_options(self):
        """Dropdown options: (unused) first, then analytic labels."""
        return [CELLPOSE_CHANNEL_UNUSED] + self._analytic_channel_options()
    
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
            children.append(widgets.HTML(
                '<p style="color: orange;">No channel slots (metadata types empty and extra slots 0). '
                'Configure populations or add extra slots.</p>'
            ))
            self.channel_config_container.children = children
            return

        channel_options = self._get_channel_options()
        analytic = self._analytic_channel_options()
        n_analytic = len(analytic)

        if mode == 'same_for_all':
            children.append(widgets.HTML('<h4>Set channel labels (applies to all samples)</h4>'))

            saved_labels = cellpose_cfg.get("channel_labels", {})

            self.channel_dropdowns['global'] = {}
            for i in range(self.n_channels):
                saved_val = saved_labels.get(str(i), saved_labels.get(i, None))
                if saved_val and saved_val in channel_options:
                    default_val = saved_val
                elif i < n_analytic:
                    default_val = analytic[i]
                else:
                    default_val = CELLPOSE_CHANNEL_UNUSED

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

            saved_per_sample = cellpose_cfg.get("per_sample_channel_labels", {})

            sample_names = self._get_sample_names()

            for sample_name in sample_names:
                children.append(widgets.HTML(f'<b>{sample_name}</b>'))
                self.channel_dropdowns[sample_name] = {}

                sample_saved = saved_per_sample.get(sample_name, {})

                sample_row = []
                for i in range(self.n_channels):
                    saved_val = sample_saved.get(str(i), sample_saved.get(i, None))
                    if saved_val and saved_val in channel_options:
                        default_val = saved_val
                    elif i < n_analytic:
                        default_val = analytic[i]
                    else:
                        default_val = CELLPOSE_CHANNEL_UNUSED

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

        self.n_extra_slots = int(self.extra_slots_input.value)
        self.n_channels = self.n_base_channels + self.n_extra_slots
        cellpose_cfg["extra_channel_slots"] = self.n_extra_slots
        cellpose_cfg["number_of_channels"] = self.n_channels
        self.n_channels_display.value = self._channels_summary_html()

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


class CellposeInferencePanel:
    """Widget panel for Cellpose model loading and inference."""
    
    def __init__(self, metadata_loader, category='organoid'):
        self.metadata_loader = metadata_loader
        self.category = category  # 'organoid' or 'immune'
        self._detect_cell_types()
        
        # Model Selection
        self.model_path_text = widgets.Text(
            value='',
            placeholder='Type path or click Browse...',
            description=f'{category.capitalize()} model:',
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='70%')
        )

        self.browse_btn = widgets.Button(
            description='Browse…',
            icon='folder-open',
            layout=widgets.Layout(width='100px')
        )
        self.browse_btn.on_click(self._on_browse_clicked)
        self.browse_output = widgets.Output()

        self.load_btn = widgets.Button(
            description=f"Load {category.capitalize()} Model", 
            button_style='info', 
            layout={'width': 'max-content'}
        )
        self.load_btn.on_click(self._on_load_clicked)
        
        # Inference
        self.label_selector = widgets.Dropdown(
            description='Cell type to segment:',
            style={'description_width': 'initial'},
            layout={'width': 'max-content'},
        )
        self._update_label_options()

        self.apply_btn = widgets.Button(
            description=f"Apply {category.capitalize()} Segmentation", 
            button_style='success', 
            layout={'width': 'max-content'}
        )
        self.apply_btn.on_click(self._on_apply_clicked)
        
        self.out = widgets.Output()
        
        # State variables
        self.pretrained_model_dir = None
        self.channel_order = []
        self.available_labels = []

        # Build UI
        parts = [
            widgets.HBox([self.model_path_text, self.browse_btn]),
            self.browse_output,
            self.load_btn,
            widgets.HTML("<hr>"),
            self.label_selector,
            self.apply_btn,
            self.out
        ]
        
        self.ui = widgets.VBox(parts, layout=widgets.Layout(grid_gap='10px'))

    def _detect_cell_types(self):
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
        
        if self.category == 'organoid':
            self.cell_types = self.organoid_types
        else:
            self.cell_types = self.immune_types + self.other_types

    def _update_label_options(self):
        if self.cell_types:
            self.label_selector.options = self.cell_types
            self.label_selector.value = self.cell_types[0]
        else:
            self.label_selector.options = ["None detected"]

    def _on_browse_clicked(self, b):
        with self.browse_output:
            self.browse_output.clear_output()
            fc = FileChooser('.', select_default=True, show_only_dirs=False)
            def on_select(chooser):
                if chooser.selected:
                    self.model_path_text.value = str(chooser.selected)
                self.browse_output.clear_output()
            fc.register_callback(on_select)
            widgets.display(fc)

    def _on_load_clicked(self, b):
        self.pretrained_model_dir = self.model_path_text.value
        if not self.pretrained_model_dir:
            with self.out:
                self.out.clear_output()
                print("Error: Please select a model path first.")
            return

        # Extract channel order from model name
        try:
            self.channel_order = self.pretrained_model_dir.split(os.path.sep)[-1].split('__')[-1].split('_')
            self.channel_order = [ch.split('-')[-1] for ch in self.channel_order]
        except:
            self.channel_order = ["Unknown"]

        # Get labels from config
        self.available_labels = []
        config_path = getattr(self.metadata_loader, 'behav3d_parameters_path', None)
        if config_path and config_path.exists():
            try:
                with open(config_path, 'r') as f:
                    config = yaml.safe_load(f) or {}
                cellpose_cfg = config.get("cellpose", {})
                mode = cellpose_cfg.get("labels_mode", "same_for_all")
                if mode == "same_for_all":
                    self.available_labels = list(cellpose_cfg.get("channel_labels", {}).values())
                else:
                    labels_set = set()
                    for sl in cellpose_cfg.get("per_sample_channel_labels", {}).values():
                        labels_set.update(sl.values())
                    self.available_labels = sorted(list(labels_set))
            except: pass
        
        self.available_labels = [
            l for l in self.available_labels
            if l and l != 'dead' and l != CELLPOSE_CHANNEL_UNUSED
        ]

        with self.out:
            self.out.clear_output()
            print(f"Loaded {self.category.capitalize()} model")
            print(f"   Channels detected: {self.channel_order}")
            print(f"   Cell types: {self.cell_types}")

    def _on_apply_clicked(self, b):
        if not self.pretrained_model_dir:
            with self.out:
                print("Error: Load a model first.")
            return
            
        label = self.label_selector.value
        with self.out:
            self.out.clear_output()
            try:
                # updates memory and saves CSV
                _, summary = run_cellpose_and_sync_metadata(
                    output_dir=self.metadata_loader.output_dir,
                    metadata_loader=self.metadata_loader,
                    pretrained_model_dir=self.pretrained_model_dir,
                    input_channels=[label], # assign to ch0
                    label_name=label,
                    timepoint_range= None
                )
                
                # summary
                n_proc = len(summary['processed'])
                n_skip = len(summary['skipped'])
                if n_skip > 0:
                    print(f"Cellpose segmentation finished: {n_proc} processed, {n_skip} skipped.")
                    print(f"Skipped samples: {', '.join(summary['skipped'])}")
                else:
                    print(f"Cellpose segmentation finished for all {n_proc} samples.")
                print("Metadata updated.")
            except Exception:
                traceback.print_exc()

    def display(self):
        display(self.ui)


class SAMChannelSelectorBox(widgets.VBox):
    """Per-cell-type channel picker for Cellpose-SAM.

    Notebook twin of ``behav3d.napari._segmentation.SAMChannelSelector``: one
    checkbox per raw channel, for whichever cell type is currently selected in
    ``cell_type_dropdown`` - "which channels feed this cell type", separate
    from classic Cellpose's per-channel :class:`CellposeChannelConfigPanel`.
    Stores to its own config block, ``cellpose_sam.channel_selection:
    {cell_type: [channel_idx, ...]}``, uniform across samples (like the
    existing ``cellpose_sam.size_filter``).
    """

    def __init__(self, metadata_loader, cell_type_dropdown, get_cell_types):
        self.metadata_loader = metadata_loader
        self._cell_type_dropdown = cell_type_dropdown
        self._get_cell_types = get_cell_types
        self._checkboxes = []

        self.hint = widgets.HTML(
            "<i>Select which raw channels feed Cellpose-SAM for the cell type "
            "chosen above. Multiple channels are stacked as separate planes "
            "(up to 3).</i>"
        )
        self.checkbox_box = widgets.VBox()
        super().__init__([self.hint, self.checkbox_box])

        self._detect_channels()
        self._build_checkboxes()
        cell_type_dropdown.observe(lambda _c: self._load_for_current_cell_type(), names="value")

    # ── channel count + label hints (read-only reference to classic Cellpose config) ──
    def _detect_channels(self):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel,
        )
        metadata = self.metadata_loader.metadata
        if metadata is None:
            self._n_channels = 0
            return
        n = (len(detect_organoid_types_from_metadata(metadata))
             + len(detect_immune_cell_types_from_metadata(metadata))
             + len(detect_other_cell_types_from_metadata(metadata)))
        if has_dead_channel(metadata):
            n += 1
        classic_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {}) or {}
        manual_n = classic_cfg.get("manual_n_channels")
        self._n_channels = manual_n if manual_n is not None else n

    # ── UI ──────────────────────────────────────────────────────────────
    def _build_checkboxes(self):
        self._checkboxes = []
        if self._n_channels == 0:
            self.checkbox_box.children = [widgets.HTML(
                '<span style="color: orange;">No channels detected. Load metadata first.</span>'
            )]
            return

        self._ensure_migrated()

        boxes = []
        for i in range(self._n_channels):
            cb = widgets.Checkbox(value=False, description=f"Channel {i}", indent=False)
            cb.observe(lambda _c: self._sync_current_cell_type_to_memory(), names="value")
            boxes.append(cb)
            self._checkboxes.append(cb)
        self.checkbox_box.children = boxes

        self._load_for_current_cell_type()

    # ── config: cellpose_sam.channel_selection ────────────────────────────
    def _current_cell_type(self):
        v = self._cell_type_dropdown.value
        return v if v and v != "None detected" else None

    def _ensure_migrated(self):
        """Seed a default for any cell type that has no channel_selection
        entry yet, from classic Cellpose's single-channel assignment (if
        any) - so a working config exists immediately for users who already
        configured channel labels there, without visiting every cell type."""
        from behav3d.preprocessing.segmentation.cellpose_prediction import (
            _label_to_channel_from_stored_map,
        )
        cfg = self.metadata_loader.behav3d_parameters.setdefault("cellpose_sam", {})
        selection = cfg.setdefault("channel_selection", {})
        classic_cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose", {}) or {}
        label_to_channel = _label_to_channel_from_stored_map(classic_cfg.get("channel_labels", {}))
        for ct in self._get_cell_types():
            if ct not in selection:
                ch = label_to_channel.get(ct)
                selection[ct] = [ch] if ch is not None else []

    def _load_for_current_cell_type(self):
        cell_type = self._current_cell_type()
        if cell_type is None or not self._checkboxes:
            return
        selected = set(self.channels_for(cell_type))
        for i, cb in enumerate(self._checkboxes):
            cb.value = i in selected

    def channels_for(self, cell_type) -> list:
        """Return the persisted channel indices for *cell_type*."""
        cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose_sam", {}) or {}
        return [int(c) for c in ((cfg.get("channel_selection") or {}).get(cell_type) or [])]

    def _sync_current_cell_type_to_memory(self):
        cell_type = self._current_cell_type()
        if cell_type is None:
            return
        cfg = self.metadata_loader.behav3d_parameters.setdefault("cellpose_sam", {})
        selection = cfg.setdefault("channel_selection", {})
        selection[cell_type] = [i for i, cb in enumerate(self._checkboxes) if cb.value]

    # ── refresh ──────────────────────────────────────────────────────────
    def refresh(self):
        """Re-detect channel count/hints and rebuild - e.g. when metadata
        changes or classic Cellpose's channel config is edited elsewhere."""
        self._detect_channels()
        self._build_checkboxes()


class CellposeSAMInferencePanel:
    """Notebook panel for Cellpose-SAM (cellpose v4) zero-shot segmentation.

    Notebook twin of ``behav3d.napari._segmentation.CellposeSAMWidget``; both call
    ``run_cellpose_sam_and_sync_metadata``, so they stay behaviourally identical.

    There is no model to load or train. Channel selection is its own,
    Cellpose-SAM-specific per-cell-type checkbox picker
    (:class:`SAMChannelSelectorBox`), separate from classic Cellpose's
    per-channel ``CellposeChannelConfigPanel``.
    """

    def __init__(self, metadata_loader, category=None):
        """*category* optionally restricts the cell-type dropdown to
        ``'organoid'`` or ``'immune/other cell types'``; ``None`` offers all.
        """
        self.metadata_loader = metadata_loader
        self.category = category
        self._preview_labels = None
        self._preview_sizes = None
        self._detect_cell_types()

        cfg = _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose_sam", {}) or {}

        # ── Environment ──────────────────────────────────────────────
        self.env_status = widgets.HTML("Checking Cellpose-SAM environment…")
        self.setup_btn = widgets.Button(
            description="Set up Cellpose-SAM",
            icon="wrench",
            layout={"width": "max-content"},
        )
        self.setup_btn.on_click(self._on_setup_clicked)

        # ── What to segment ──────────────────────────────────────────
        self.label_selector = widgets.Dropdown(
            description="Cell type to segment:",
            options=self.cell_types or ["None detected"],
            style={"description_width": "initial"},
            layout={"width": "max-content"},
        )
        self.label_selector.observe(self._on_cell_type_changed, names="value")

        self.all_cell_types_cb = widgets.Checkbox(
            value=bool(cfg.get("all_cell_types", False)),
            description="Run all cell types in one batch",
            indent=False,
        )

        self.channel_panel = SAMChannelSelectorBox(
            self.metadata_loader, self.label_selector,
            get_cell_types=lambda: self.cell_types,
        )

        self.use_all_tp = widgets.Checkbox(
            value=bool(cfg.get("use_all_timepoints", True)),
            description="Process all timepoints",
            indent=False,
        )
        self.tp_start = widgets.IntText(value=int(cfg.get("tp_start", 0)),
                                        description="From:", layout={"width": "150px"})
        self.tp_end = widgets.IntText(value=int(cfg.get("tp_end", 0)),
                                      description="To:", layout={"width": "150px"})

        # ── Device ───────────────────────────────────────────────────
        from behav3d.preprocessing.segmentation.convpaint_train import _detect_torch_devices
        self._devices = _detect_torch_devices()
        self.device_selector = widgets.Dropdown(
            description="GPU device:",
            options=[(label, value) for label, value in self._devices],
            value=cfg.get("device", "auto") if any(
                v == cfg.get("device", "auto") for _l, v in self._devices
            ) else "auto",
            style={"description_width": "initial"},
            layout={"width": "max-content"},
        )
        self.force_cpu_cb = widgets.Checkbox(
            value=bool(cfg.get("force_cpu", False)),
            description="Force CPU-only processing",
            indent=False,
        )

        # ── Segmentation settings ────────────────────────────────────
        from behav3d.preprocessing.segmentation.cellpose_sam_prediction import CPSAM_MODELS
        self.model_selector = widgets.Dropdown(
            description="Model:", options=list(CPSAM_MODELS),
            value=cfg.get("model_name", "cpsam") if cfg.get("model_name", "cpsam") in CPSAM_MODELS else "cpsam",
            style={"description_width": "initial"}, layout={"width": "max-content"},
        )
        self.diameter = widgets.FloatText(value=float(cfg.get("diameter", 0.0)),
                                          description="Diameter (0=auto):",
                                          style={"description_width": "initial"})
        self.do_3d = widgets.Checkbox(value=bool(cfg.get("do_3D", True)),
                                      description="3D segmentation (do_3D)", indent=False)
        self.do_3d.observe(self._on_do_3d_changed, names="value")
        self.stitch_threshold = widgets.FloatSlider(
            value=float(cfg.get("stitch_threshold", 0.0)), min=0.0, max=1.0, step=0.05,
            description="Stitch threshold:", style={"description_width": "initial"},
        )
        self.flow_threshold = widgets.FloatSlider(
            value=float(cfg.get("flow_threshold", 0.4)), min=0.0, max=3.0, step=0.1,
            description="Flow threshold:", style={"description_width": "initial"},
        )
        self.cellprob_threshold = widgets.FloatSlider(
            value=float(cfg.get("cellprob_threshold", 0.0)), min=-6.0, max=6.0, step=0.5,
            description="Cell prob threshold:", style={"description_width": "initial"},
        )
        self.batch_size = widgets.IntText(value=int(cfg.get("batch_size", 8)),
                                          description="Batch size:",
                                          style={"description_width": "initial"})
        self.drop_2d_cb = widgets.Checkbox(
            value=bool(cfg.get("drop_2d_segments", True)),
            description="Remove flat (single-slice) segments", indent=False,
        )

        from behav3d.preprocessing.segmentation.cellpose_sam_prediction import (
            power_safe_settings,
        )
        self._power_safe_preset = power_safe_settings()

        # Advanced / normalization
        self.norm_low = widgets.FloatText(value=float(cfg.get("norm_percentile_low", 1.0)),
                                          description="Norm pct low:",
                                          style={"description_width": "initial"})
        self.norm_high = widgets.FloatText(value=float(cfg.get("norm_percentile_high", 99.0)),
                                           description="Norm pct high:",
                                           style={"description_width": "initial"})
        self.norm3d = widgets.Checkbox(value=bool(cfg.get("norm3D", True)),
                                       description="Normalize across whole Z stack (norm3D)",
                                       indent=False)
        self.niter = widgets.IntText(value=int(cfg.get("niter", 0)), description="niter:",
                                     style={"description_width": "initial"})
        self.flow3d_smooth = widgets.FloatText(value=float(cfg.get("flow3D_smooth", 0.0)),
                                               description="flow3D smooth:",
                                               style={"description_width": "initial"})
        self.max_size_fraction = widgets.FloatText(
            value=float(cfg.get("max_size_fraction", 1.0)),
            description="Max size fraction:", style={"description_width": "initial"},
        )
        # Cellpose GUI's "custom filter settings" group.
        self.sharpen_radius = widgets.FloatText(value=float(cfg.get("sharpen_radius", 0.0)),
                                                description="Sharpen radius:",
                                                style={"description_width": "initial"})
        self.smooth_radius = widgets.FloatText(value=float(cfg.get("smooth_radius", 0.0)),
                                               description="Smooth radius:",
                                               style={"description_width": "initial"})
        self.tile_norm_blocksize = widgets.FloatText(
            value=float(cfg.get("tile_norm_blocksize", 0.0)),
            description="Tile norm blocksize:", style={"description_width": "initial"})
        self.tile_norm_smooth3D = widgets.FloatText(
            value=float(cfg.get("tile_norm_smooth3D", 0.0)),
            description="Tile norm smooth3D:", style={"description_width": "initial"})
        self.n_threads = widgets.IntText(
            value=int(cfg.get("n_threads", 0)), description="CPU threads (0=all):",
            style={"description_width": "initial"},
        )
        self.cooldown_s = widgets.FloatText(
            value=float(cfg.get("cooldown_s", 0.0)), description="Cooldown between frames (s):",
            style={"description_width": "initial"},
        )
        advanced = widgets.Accordion(children=[widgets.VBox([
            self.norm_low, self.norm_high, self.norm3d, self.niter,
            self.flow3d_smooth, self.max_size_fraction,
            self.sharpen_radius, self.smooth_radius,
            self.tile_norm_blocksize, self.tile_norm_smooth3D,
            widgets.HTML(
                "<i>max_size_fraction: cellpose's own default is 0.4, which silently "
                "deletes any object covering &gt;40% of the frame (large organoids). "
                "BEHAV3D defaults to 1.0 and filters by size instead.</i>"
            ),
            widgets.HBox([self.n_threads, self.cooldown_s]),
            widgets.HTML(
                "<i>CPU threads caps mask post-processing/compression threads; "
                "cooldown pauses between timepoints to let voltage regulators and "
                "battery recover. A full Cellpose-SAM run saturates the GPU and every "
                "CPU core at once, which on some laptops draws more power than the "
                "machine can sustain and causes an abrupt shutdown mid-run. If that "
                f"happens, try around {self._power_safe_preset['n_threads']} threads "
                f"and a {self._power_safe_preset['cooldown_s']:g}s cooldown.</i>"
            ),
        ])])
        advanced.set_title(0, "Advanced / normalization")
        advanced.selected_index = None

        # ── Size filter + preview ────────────────────────────────────
        self.preview_sample = widgets.Dropdown(
            description="Preview sample:", options=self._sample_names() or ["(no samples)"],
            style={"description_width": "initial"}, layout={"width": "max-content"},
        )
        self.preview_tp = widgets.IntText(value=int(cfg.get("preview_timepoint", 0)),
                                          description="t:", layout={"width": "120px"})
        self.preview_btn = widgets.Button(description="Preview segmentation", icon="search",
                                          button_style="info", layout={"width": "max-content"})
        self.preview_btn.on_click(self._on_preview_clicked)

        self.size_min = widgets.IntText(value=0, description="Min size (voxels):",
                                        style={"description_width": "initial"})
        self.size_max = widgets.IntText(value=0, description="Max size (voxels):",
                                        style={"description_width": "initial"})
        self.size_min.observe(self._on_size_changed, names="value")
        self.size_max.observe(self._on_size_changed, names="value")
        self.filter_stats = widgets.HTML("<i>Run a preview to see object sizes.</i>")
        self.hist_out = widgets.Output()
        self.show_napari_btn = widgets.Button(
            description="Show preview in napari", icon="eye",
            layout={"width": "max-content"},
        )
        self.show_napari_btn.on_click(self._on_show_napari_clicked)
        self._load_size_filter_for_cell_type()

        # ── Apply ────────────────────────────────────────────────────
        self.apply_btn = widgets.Button(description="Apply Cellpose-SAM Segmentation",
                                        button_style="success", layout={"width": "max-content"})
        self.apply_btn.on_click(self._on_apply_clicked)
        self.spinner = widgets.HTML(value=spinning_loader)
        self.spinner.layout.display = "none"
        self.out = widgets.Output()

        self.ui = widgets.VBox([
            widgets.HTML("<b>1. Environment</b>"),
            self.env_status, self.setup_btn,
            widgets.HTML("<hr><b>2. What to segment</b>"),
            self.label_selector, self.channel_panel, self.all_cell_types_cb,
            self.use_all_tp, widgets.HBox([self.tp_start, self.tp_end]),
            widgets.HTML("<hr><b>3. Device</b>"),
            self.device_selector, self.force_cpu_cb,
            widgets.HTML("<hr><b>4. Segmentation settings</b>"),
            self.model_selector, self.diameter, self.do_3d, self.stitch_threshold,
            self.flow_threshold, self.cellprob_threshold, self.batch_size, self.drop_2d_cb,
            advanced,
            widgets.HTML("<hr><b>5. Size filter (preview before batch)</b>"),
            widgets.HBox([self.preview_sample, self.preview_tp]),
            self.preview_btn,
            widgets.HBox([self.size_min, self.size_max]),
            self.filter_stats, self.hist_out, self.show_napari_btn,
            widgets.HTML("<hr>"),
            widgets.HBox([self.apply_btn, self.spinner]),
            self.out,
        ], layout=widgets.Layout(grid_gap="6px"))

        self._on_do_3d_changed({"new": self.do_3d.value})
        self._refresh_env_status()

    # ── helpers ─────────────────────────────────────────────────────
    def _detect_cell_types(self):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
        )
        metadata = self.metadata_loader.metadata
        self.organoid_types = _filter_merge_types(detect_organoid_types_from_metadata(metadata))
        self.immune_types = _filter_merge_types(detect_immune_cell_types_from_metadata(metadata))
        self.other_types = _filter_merge_types(detect_other_cell_types_from_metadata(metadata))
        if self.category == "organoid":
            self.cell_types = self.organoid_types
        elif self.category:
            self.cell_types = self.immune_types + self.other_types
        else:
            self.cell_types = self.organoid_types + self.immune_types + self.other_types

    def _sample_names(self):
        md = self.metadata_loader.metadata
        return list(md["sample_name"].unique()) if md is not None else []

    def _cfg(self):
        return _cfg_get(self.metadata_loader.behav3d_parameters, "cellpose_sam", {}) or {}

    def _refresh_env_status(self):
        from behav3d.preprocessing.segmentation.cpsam_env import cpsam_env_status
        status = cpsam_env_status(self.metadata_loader.behav3d_parameters)
        if status["available"]:
            gpus = ", ".join(status.get("gpus") or []) or "no CUDA GPU"
            self.env_status.value = (
                f"<span style='color:green'>&#10003; Ready — cellpose "
                f"{status['cellpose_version']}, torch {status.get('torch_version')} ({gpus})</span>"
            )
            self.setup_btn.disabled = True
            self.apply_btn.disabled = False
            self.preview_btn.disabled = False
        else:
            self.env_status.value = f"<span style='color:#b8860b'>&#9888; {status['error']}</span>"
            self.setup_btn.disabled = False
            self.apply_btn.disabled = True
            self.preview_btn.disabled = True

    def _on_setup_clicked(self, _b):
        from behav3d.preprocessing.segmentation.cpsam_env import create_cpsam_env
        self.setup_btn.disabled = True
        self.spinner.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                create_cpsam_env(progress_cb=print)
                print("Cellpose-SAM environment ready.")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner.layout.display = "none"
                self._refresh_env_status()

    def _on_do_3d_changed(self, change):
        # cellpose ignores flow_threshold entirely when do_3D=True.
        is_3d = bool(change.get("new"))
        self.flow_threshold.disabled = is_3d
        self.flow_threshold.description = (
            "Flow threshold (ignored in 3D):" if is_3d else "Flow threshold:"
        )
        self.stitch_threshold.disabled = is_3d

    def _on_cell_type_changed(self, _change):
        self._load_size_filter_for_cell_type()

    def _load_size_filter_for_cell_type(self):
        entry = (self._cfg().get("size_filter") or {}).get(self.label_selector.value, {})
        self.size_min.unobserve(self._on_size_changed, names="value")
        self.size_max.unobserve(self._on_size_changed, names="value")
        self.size_min.value = int(entry.get("size_min", 0) or 0)
        self.size_max.value = int(entry.get("size_max", 0) or 0)
        self.size_min.observe(self._on_size_changed, names="value")
        self.size_max.observe(self._on_size_changed, names="value")

    def _collect_sam_params(self):
        diameter = float(self.diameter.value)
        return {
            "diameter": diameter if diameter > 0 else None,
            "flow_threshold": float(self.flow_threshold.value),
            "cellprob_threshold": float(self.cellprob_threshold.value),
            "niter": int(self.niter.value),
            "do_3D": bool(self.do_3d.value),
            "stitch_threshold": float(self.stitch_threshold.value),
            "flow3D_smooth": float(self.flow3d_smooth.value),
            "batch_size": int(self.batch_size.value),
            "max_size_fraction": float(self.max_size_fraction.value),
            "norm_percentile_low": float(self.norm_low.value),
            "norm_percentile_high": float(self.norm_high.value),
            "norm3D": bool(self.norm3d.value),
            "sharpen_radius": float(self.sharpen_radius.value),
            "smooth_radius": float(self.smooth_radius.value),
            "tile_norm_blocksize": float(self.tile_norm_blocksize.value),
            "tile_norm_smooth3D": float(self.tile_norm_smooth3D.value),
            "drop_2d_segments": bool(self.drop_2d_cb.value),
        }

    def _selected_device(self):
        return "cpu" if self.force_cpu_cb.value else str(self.device_selector.value or "auto")

    def _persist_params(self):
        cfg = self.metadata_loader.behav3d_parameters.setdefault("cellpose_sam", {})
        cfg.update(self._collect_sam_params())
        cfg["diameter"] = float(self.diameter.value)
        cfg["model_name"] = self.model_selector.value
        cfg["device"] = str(self.device_selector.value or "auto")
        cfg["force_cpu"] = bool(self.force_cpu_cb.value)
        cfg["all_cell_types"] = bool(self.all_cell_types_cb.value)
        cfg["use_all_timepoints"] = bool(self.use_all_tp.value)
        cfg["tp_start"] = int(self.tp_start.value)
        cfg["tp_end"] = int(self.tp_end.value)
        cfg["preview_timepoint"] = int(self.preview_tp.value)
        cfg["n_threads"] = int(self.n_threads.value)
        cfg["cooldown_s"] = float(self.cooldown_s.value)
        cfg.setdefault("size_filter", {})[self.label_selector.value] = {
            "size_min": int(self.size_min.value),
            "size_max": int(self.size_max.value),
        }
        if hasattr(self.metadata_loader, "behav3d_parameters_path"):
            yaml.safe_dump(
                self.metadata_loader.behav3d_parameters,
                self.metadata_loader.behav3d_parameters_path.open("w"),
                sort_keys=False,
            )

    # ── preview ─────────────────────────────────────────────────────
    def _on_preview_clicked(self, _b):
        from behav3d.preprocessing.segmentation.cellpose_sam_prediction import preview_cellpose_sam
        from behav3d.preprocessing.segmentation.size_filter import compute_label_sizes

        self.preview_btn.disabled = True
        self.spinner.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                labels = preview_cellpose_sam(
                    output_dir=self.metadata_loader.output_dir,
                    metadata=self.metadata_loader.metadata,
                    sample_name=self.preview_sample.value,
                    cell_type=self.label_selector.value,
                    timepoint=int(self.preview_tp.value),
                    model_name=self.model_selector.value,
                    device=self._selected_device(),
                    force_cpu=bool(self.force_cpu_cb.value),
                    sam_params=self._collect_sam_params(),
                    config=self.metadata_loader.behav3d_parameters,
                )
                self._preview_labels = labels
                _lbl, sizes = compute_label_sizes(labels)
                self._preview_sizes = sizes
                print(f"Preview complete: {sizes.size} objects.")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner.layout.display = "none"
                self.preview_btn.disabled = False
        self._update_preview_filter()

    def _on_size_changed(self, _change):
        self._update_preview_filter()

    def _update_preview_filter(self):
        """Re-filter the cached preview. No inference re-run, so this is instant."""
        if self._preview_sizes is None:
            return
        import matplotlib.pyplot as plt

        sizes = self._preview_sizes
        size_min = int(self.size_min.value) or None
        size_max = int(self.size_max.value) or None
        keep = np.ones(sizes.size, dtype=bool)
        if size_min:
            keep &= sizes >= size_min
        if size_max:
            keep &= sizes <= size_max
        n_keep = int(keep.sum())

        if sizes.size:
            self.filter_stats.value = (
                f"<b>{sizes.size}</b> objects &rarr; <b>{n_keep}</b> kept "
                f"({sizes.size - n_keep} removed). Volumes: min={sizes.min()}, "
                f"median={int(np.median(sizes))}, max={sizes.max()} voxels."
            )
        else:
            self.filter_stats.value = "<i>No objects in the preview.</i>"
            return

        with self.hist_out:
            self.hist_out.clear_output(wait=True)
            fig, ax = plt.subplots(figsize=(6, 2.4))
            ax.hist(sizes, bins=min(40, max(5, sizes.size)), color="#6699cc", edgecolor="white")
            if size_min:
                ax.axvline(size_min, color="crimson", ls="--", label=f"min={size_min}")
            if size_max:
                ax.axvline(size_max, color="darkorange", ls="--", label=f"max={size_max}")
            ax.set_xlabel("object volume (voxels)")
            ax.set_ylabel("count")
            ax.set_title(f"{sizes.size} objects → {n_keep} kept")
            if size_min or size_max:
                ax.legend(fontsize=8)
            fig.tight_layout()
            plt.show()

    def _on_show_napari_clicked(self, _b):
        from behav3d.preprocessing.segmentation.size_filter import filter_labels_by_size
        if self._preview_labels is None:
            with self.out:
                print("Run a preview first.")
            return
        filtered = filter_labels_by_size(
            self._preview_labels,
            size_min=int(self.size_min.value) or None,
            size_max=int(self.size_max.value) or None,
        )
        viewer = napari.current_viewer() or napari.Viewer()
        name = "Cellpose-SAM preview"
        if name in viewer.layers:
            viewer.layers[name].data = filtered
        else:
            viewer.add_labels(filtered, name=name)

    # ── apply ───────────────────────────────────────────────────────
    def _on_apply_clicked(self, _b):
        from behav3d.preprocessing.segmentation.cellpose_sam_prediction import (
            run_cellpose_sam_and_sync_metadata,
        )
        cell_types = (list(self.cell_types) if self.all_cell_types_cb.value
                      else [self.label_selector.value])

        self.apply_btn.disabled = True
        self.spinner.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                timepoint_range = _mk_timepoint_range(
                    self.use_all_tp.value, self.tp_start.value, self.tp_end.value
                )
                size_filters = self._cfg().get("size_filter") or {}
                total = {"processed": [], "skipped": []}
                for ct in cell_types:
                    _, summary = run_cellpose_sam_and_sync_metadata(
                        output_dir=self.metadata_loader.output_dir,
                        metadata_loader=self.metadata_loader,
                        label_name=ct,
                        model_name=self.model_selector.value,
                        timepoint_range=timepoint_range,
                        device=self._selected_device(),
                        force_cpu=bool(self.force_cpu_cb.value),
                        sam_params=self._collect_sam_params(),
                        size_filter=size_filters.get(ct),
                        resume=True,
                        n_threads=int(self.n_threads.value) or None,
                        cooldown_s=float(self.cooldown_s.value),
                    )
                    total["processed"].extend(summary["processed"])
                    total["skipped"].extend(summary["skipped"])

                n_proc, n_skip = len(total["processed"]), len(total["skipped"])
                if n_skip:
                    print(f"Cellpose-SAM finished: {n_proc} processed, {n_skip} skipped.")
                    print(f"Skipped: {', '.join(total['skipped'])}")
                else:
                    print(f"Cellpose-SAM finished: {n_proc} processed.")
                print("Metadata updated.")
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner.layout.display = "none"
                self.apply_btn.disabled = False

    def display(self):
        display(self.ui)


class DeadMaskPanel:
    """Widget panel for Run Dead Mask (Otsu) segmentation."""

    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self.has_death = has_dead_channel(metadata_loader.metadata)
        
        self.btn = widgets.Button(
            description="Run Dead Mask (Otsu)", 
            button_style='warning', 
            layout={'width': 'max-content'},
            disabled=not self.has_death
        )
        self.btn.on_click(self._on_clicked)
        self.out = widgets.Output()
        
        if not self.has_death:
            with self.out:
                print("No dead channel detected in metadata.")

        self.ui = widgets.VBox([self.btn, self.out])

    def _on_clicked(self, b):
        with self.out:
            self.out.clear_output()
            print("Starting Dead Mask (Otsu) segmentation...")
            try:
                _, summary = run_otsu_and_sync_metadata(
                    output_dir=self.metadata_loader.output_dir,
                    metadata_loader=self.metadata_loader,
                    mask_suffix="_mask_dead",
                    timepoint_range=None
                )
                # summary
                n_proc = len(summary['processed'])
                n_skip = len(summary['skipped'])
                if n_skip > 0:
                    print(f"Dead mask segmentation finished: {n_proc} processed, {n_skip} skipped.")
                    print(f"Skipped samples: {', '.join(summary['skipped'])}")
                else:
                    print(f"Dead mask segmentation finished for all {n_proc} samples.")
                print("Metadata updated.")
            except Exception:
                traceback.print_exc()

    def display(self):
        display(self.ui)


class SegmentationVisualizationPanel:
    """Panel for visualizing segmentation results (pixel classifier, Cellpose, Otsu, …) in Napari."""

    _CHANNEL_COLORS = (
        "cyan", "yellow", "red", "green", "magenta", "blue",
        "gray", "turbo", "viridis", "plasma", "inferno", "twilight",
    )

    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self._viewer = None

        self._status = widgets.HTML("<b>Waiting for metadata…</b>")
        self._refresh_btn = widgets.Button(
            description="Refresh",
            tooltip="Reload sample list from metadata",
        )
        self._refresh_btn.on_click(self._on_refresh_clicked)

        self.sample_dropdown = widgets.Dropdown(
            options=[],
            value=None,
            description="Sample:",
            layout=widgets.Layout(width="350px"),
            disabled=True,
        )

        self.use_range = widgets.Checkbox(
            description="Use custom time range",
            value=False,
            indent=False,
        )
        self.start_t = widgets.IntText(value=0, description="Start T:", layout=widgets.Layout(width="180px"))
        self.end_t   = widgets.IntText(value=2, description="End T:",   layout=widgets.Layout(width="180px"))
        self.range_box = widgets.HBox([self.start_t, self.end_t])
        self.range_box.layout.display = "none"
        self.use_range.observe(
            lambda ch: setattr(self.range_box.layout, "display", "flex" if ch["new"] else "none"),
            names="value",
        )

        self.open_button = widgets.Button(
            description="Visualize segmentation results",
            button_style="primary",
            icon="eye",
            layout=widgets.Layout(width="max-content"),
            disabled=True,
        )
        self.close_button = widgets.Button(
            description="Close viewer",
            button_style="danger",
            icon="stop",
            layout=widgets.Layout(width="max-content", display="none"),
        )
        self.msg = widgets.Output()

        self.open_button.on_click(self._on_open_clicked)
        self.close_button.on_click(self._on_close_clicked)

        self._panel = widgets.VBox([
            widgets.HBox([self._status, self._refresh_btn]),
            self.sample_dropdown,
            self.use_range,
            self.range_box,
            widgets.HBox([self.open_button, self.close_button]),
            self.msg,
        ])
        self._maybe_build_from_loader()

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------
    @property
    def ui(self):
        return self._panel

    def display(self):
        display(self._panel)

    # ------------------------------------------------------------------
    # Core viewer logic
    # ------------------------------------------------------------------
    def open_viewer(self):
        """Build mask list from metadata and launch Napari."""
        sample_name = self.sample_dropdown.value
        time_range  = (int(self.start_t.value), int(self.end_t.value)) if self.use_range.value else None
        metadata    = self.metadata_loader.metadata
        output_dir  = Path(self.metadata_loader.output_dir)

        with self.msg:
            self.msg.clear_output()
            try:
                print(f"Sample: {sample_name}")
                print(f"Timepoint range: {time_range}")

                raw_image_zarr = output_dir / "images" / sample_name / f"{sample_name}.zarr"

                # ---- collect mask paths from metadata ----
                mask_files = {}
                organoid_types = detect_organoid_types_from_metadata(metadata)
                immune_types   = detect_immune_cell_types_from_metadata(metadata)
                other_types    = detect_other_cell_types_from_metadata(metadata)
                print(f"Detected cell types – Organoids: {organoid_types}, Immune: {immune_types}, Other: {other_types}")

                sample_row = metadata[metadata["sample_name"] == sample_name]
                if sample_row.empty:
                    print(f"[Warning] '{sample_name}' not found in metadata")
                else:
                    row = sample_row.iloc[0]
                    for cell_type, prefix in (
                        [(ct, "or") for ct in organoid_types]
                        + [(ct, "im") for ct in immune_types]
                        + [(ct, "ot") for ct in other_types]
                    ):
                        col = f"{prefix}_{cell_type}_segments_image_path"
                        if col in metadata.columns:
                            p = row.get(col)
                            if p and not pd.isna(p):
                                mask_files[f"{cell_type}_segments"] = Path(p)

                    if has_dead_channel(metadata) and "dead_mask_path" in metadata.columns:
                        p = row.get("dead_mask_path")
                        if p and not pd.isna(p):
                            mask_files["mask_dead"] = Path(p)

                if not mask_files:
                    # Fallback: scan directory for any *_segments.zarr
                    sample_dir = output_dir / "images" / sample_name
                    if sample_dir.exists():
                        for zarr_path in sample_dir.glob("*_segments.zarr"):
                            name = zarr_path.stem.replace(f"{sample_name}_", "").replace("_segments", "")
                            mask_files[f"{name}_segments"] = zarr_path
                        dead = sample_dir / f"{sample_name}_mask_dead.zarr"
                        if dead.exists():
                            mask_files["mask_dead"] = dead

                for name, path in mask_files.items():
                    print(f"Segments/Mask '{name}': {path}")

                # ---- load arrays ----
                loaded_masks = {}
                for name, path in mask_files.items():
                    if path.exists():
                        loaded_masks[name] = load_zarr(path)
                    else:
                        print(f"[Warning] Segments/Mask '{name}' not found: {path}")

                img_arr = load_zarr(raw_image_zarr)
                print(f"Image shape: {img_arr.shape}")

                # ---- time slicing ----
                if time_range is not None:
                    s, e = time_range
                    s = max(0, min(s, img_arr.shape[0] - 1))
                    e = max(s,  min(e, img_arr.shape[0] - 1))
                    print(f"Slicing T={s}…{e}")
                    img_arr = img_arr[s:e + 1]
                    loaded_masks = {k: v[s:e + 1] for k, v in loaded_masks.items()}

                # ---- build viewer ----
                if self._viewer is not None:
                    try:
                        self._viewer.close()
                    except Exception:
                        pass

                viewer = napari.Viewer()

                if img_arr.ndim != 5:
                    raise ValueError(f"Expected (T,C,Z,Y,X), got {img_arr.shape}")

                raw_shape = tuple(int(v) for v in img_arr.shape)
                viewer_label_shape = (raw_shape[0],) + tuple(raw_shape[-3:])
                print(f"Raw viewer shape for labels: {viewer_label_shape}")

                T, C, Z, Y, X = raw_shape
                for ch in range(C):
                    viewer.add_image(
                        img_arr[:, ch],
                        name=f"channel_{ch}",
                        colormap=self._CHANNEL_COLORS[ch] if ch < len(self._CHANNEL_COLORS) else "gray",
                        scale=(1, 1, 1, 1),
                        blending="additive",
                        opacity=0.8,
                        channel_axis=None,
                    )

                added_label_layers = 0
                for name, mask_arr in loaded_masks.items():
                    source_path = mask_files.get(name)
                    source_text = str(source_path) if source_path is not None else "<unknown>"
                    original_shape = tuple(int(v) for v in mask_arr.shape)
                    try:
                        mask_view = _normalize_segmentation_labels_for_view(
                            mask_arr,
                            raw_shape,
                            layer_name=name,
                        )
                        normalized_shape = tuple(int(v) for v in mask_view.shape)
                        print(
                            f"Adding layer '{name}' from {source_text} "
                            f"(original shape {original_shape} -> viewer shape {normalized_shape})"
                        )
                        viewer.add_labels(
                            mask_view,
                            name=name,
                            scale=(1, 1, 1, 1),
                            opacity=0.6,
                        )
                        added_label_layers += 1
                    except Exception as e:
                        print(
                            f"Skipping layer '{name}' from {source_text} "
                            f"(original shape {original_shape}, expected viewer shape {viewer_label_shape}): {e}"
                        )

                if not loaded_masks:
                    print("No segments/mask layers found.")
                elif added_label_layers == 0:
                    print("No segments/mask layers added after shape normalization.")

                # Try to load and display the ConvPaint label legend if present
                try:
                    from behav3d.preprocessing.segmentation.convpaint_label_map import load_label_map, unified_label_map_path
                    pixel_class_outdir = output_dir / "images" / "PixelClassification"
                    map_path = unified_label_map_path(pixel_class_outdir)
                    if map_path.exists():
                        label_map = load_label_map(map_path)
                        has_death = has_dead_channel(metadata)
                        _build_notebook_annotation_legend(label_map, has_death=has_death)
                except Exception as e:
                    print(f"Could not load ConvPaint legend: {e}")

                napari.run()
                self._viewer = viewer
                self.close_button.layout.display = "inline-block"

            except Exception:
                traceback.print_exc()

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _maybe_build_from_loader(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if isinstance(df, pd.DataFrame) and not df.empty and "sample_name" in df.columns:
            names = df["sample_name"].astype(str).unique().tolist()
            self.sample_dropdown.options = names
            self.sample_dropdown.value   = names[0]
            self.sample_dropdown.disabled = False
            self.open_button.disabled     = False
            self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
        else:
            self._status.value = "<b>Waiting for metadata…</b>"
            self.sample_dropdown.disabled = True
            self.open_button.disabled     = True

    def _on_open_clicked(self, _):    self.open_viewer()
    def _on_refresh_clicked(self, _): self._maybe_build_from_loader()
    def _on_close_clicked(self, _):
        with self.msg:
            try:
                if self._viewer:
                    self._viewer.close()
                    self._viewer = None
            finally:
                self.close_button.layout.display = "none"
