import os
import re
import traceback
from copy import deepcopy
from pathlib import Path

import ipywidgets as widgets
import numpy as np
import pandas as pd
import yaml
from IPython.display import display
from ipyfilechooser import FileChooser
from traitlets import Any, Unicode, Bool

from behav3d.core.metadata import (
    check_behav3d_metadata,
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    load_behav3d_metadata,
)
from behav3d.io.formats.zarr import save_as_zarr
from behav3d.io.images import (
    convert_label_file_to_zarr,
    convert_raw_file_to_zarr,
    get_image_dimension_order,
    get_image_shape,
    load_image,
)
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv

# ===============================
# CONSTANTS
# ===============================
CONFIG_PATH = Path("behav3d_config.yml")

_DEFAULT_CONFIG = {
    "seed": 42,
    "paths": {
        "metadata_csv": "",
        "output_dir": ""
    },
    "dim_order": {"default_apply_all": "TCZYX"},
    "signal_unmixing": {
        "sink_channel": 2,
        "source_channels": "0, 1",
        "train_z": "4,5,6",
        "train_t": "1,3,6",
        "bg_percentile": 1,
        "median_size": 3,
        "gaussian_sigma": 4,
        "visualize_napari": True,
        "training_log": False,
        "timepoints": 0,
        "use_all_timepoints": True,
        "use_range": False,
        "start_t": 0,
        "end_t": 0,
        "use_unmix_path": False,
        "channel_colors": ["cyan", "yellow", "red", "green", "magenta", "blue",
                           "gray", "turbo", "viridis", "plasma", "inferno", "twilight"]
    },
    "pixel_classifier": {
        "examples_per_sample": 3,
        "sample_specific_classifier": False,
        "workers": 8,
        "use_all_timepoints": True,
        "tp_start": 0,
        "tp_end": 0,
    },
    "cellpose": {
        "number_of_channels": 0,
        "labels_mode": "same_for_all",  # "same_for_all" or "per_sample"
        "channel_labels": {},  # {0: "organoid1", 1: "tcell", ...}
        "per_sample_channel_labels": {},  # {sample_name: {0: "organoid1", ...}, ...}
    },
    "tracking": {
        "track_organoids_together": False,
        "all_organoids": {
            "method": "propagation_all_organoids",
            "overwrite": False,
        },
        "immune": {
            "method": "trackpy",
            "overwrite": False,
            "lap": {
                "track_cost_px": 60,
                "gap_close_cost_px": 45,
                "gap_close_max_frames": 5,
                "merging_cost_px": 0,
                "splitting_cost_px": 0
            },
            "trackpy": {
                "search_range_px": 31,
                "memory_frames": 5,
                "adaptive_stop": 5.0,
                "adaptive_step": 0.95
            },
            "btrack": {
                "config_preset": "cell",
                "config_path": "",
                "max_search_radius": 100,
                "update_method": "EXACT",
                "step_size": 100,
                "n_workers": 16,
                "use_optimize": False,
                "hypotheses": ["P_FP", "P_init", "P_term", "P_link"],
                "dist_thresh": 60,
                "time_thresh": 3
            }
        },
        "other": {
            "method": "trackpy",
            "overwrite": False,
            "lap": {
                "track_cost_px": 60,
                "gap_close_cost_px": 45,
                "gap_close_max_frames": 5,
                "merging_cost_px": 0,
                "splitting_cost_px": 0
            },
            "trackpy": {
                "search_range_px": 31,
                "memory_frames": 5,
                "adaptive_stop": 5.0,
                "adaptive_step": 0.95
            },
            "btrack": {
                "config_preset": "cell",
                "config_path": "",
                "max_search_radius": 100,
                "update_method": "EXACT",
                "step_size": 100,
                "n_workers": 8,
                "use_optimize": False,
                "hypotheses": ["P_FP", "P_init", "P_term", "P_link"],
                "dist_thresh": 60,
                "time_thresh": 3
            }
        },
        "organoid": {
            "method": "propagation",  # also supports "propagation_all_organoids"
            "overwrite": False,
            "lap": {
                "track_cost_px": 60,
                "gap_close_cost_px": 80,
                "gap_close_max_frames": 3,
                "merging_cost_px": 0,
                "splitting_cost_px": 0
            },
            "trackpy": {
                "search_range_px": 35,
                "memory_frames": 2,
                "adaptive_stop": 10.0,
                "adaptive_step": 0.95
            },
            "btrack": {
                "config_preset": "cell",
                "config_path": "",
                "max_search_radius": 100,
                "update_method": "EXACT",
                "step_size": 100,
                "n_workers": 8,
                "use_optimize": False,
                "hypotheses": ["P_FP", "P_init", "P_term", "P_link"],
                "dist_thresh": 60,
                "time_thresh": 3
            }
        }
    },
    "tracking_visualization": {
        "sample_name": None,
        "use_range": False,
        "start_t": 0,
        "end_t": 0,
        "channel_colors": ["cyan", "yellow", "red", "green", "magenta", "blue",
                           "gray", "turbo", "viridis", "plasma", "inferno", "twilight"]
    },
    "features": {
        "immune": {
            "features_choice": ["movement", "intensity", "contact", "death"],
            "contact_threshold": 0,
            "n_workers": 16,
            "overwrite": False
        },
        "organoid": {
            "features_choice": ["intensity", "death", "morphology"],
            "contact_threshold": 0,
            "n_workers": 8,
            "overwrite": False
        },
        "other": {
            "features_choice": ["movement", "intensity", "contact"],
            "contact_threshold": 0,
            "n_workers": 8,
            "overwrite": False
        }
    },
    "track_filtering": {
        "immune": {
            "exp_duration": 24.0,
            "exp_duration_enabled": False,
            "min_track_length": 0,
            "min_track_length_enabled": True,
            "max_track_length": 999999,
            "max_track_length_enabled": True,
        },
        "organoid": {
            "exp_duration_enabled": False,
            "exp_duration": 24.0,
            "min_track_length_enabled": False,
            "min_track_length": 50,
            "max_track_length_enabled": False,
            "max_track_length": 999999,
        },
        "other": {
            "exp_duration": 24.0,
            "exp_duration_enabled": False,
            "min_track_length": 0,
            "min_track_length_enabled": True,
            "max_track_length": 999999,
            "max_track_length_enabled": True,
        }
    },
    "analysis": {
        "immune": {
            "seed": 42,
            "umap_min_dist": 0.1,
            "umap_n_neighbors": 15,
            "nr_of_clusters": 5,
            "dtw_feature_groups_enabled": {
                "morphology": False,
                "movement": False,
                "intensity": False,
                "death": False,
                "contact": False,
            },
            'dtw_features_input': [],
            "dtw_features_resolved": [],
            'z_normalize': {},
        },
        "organoid": {
            "seed": 42,
            "umap_min_dist": 0.1,
            "umap_n_neighbors": 15,
            "nr_of_clusters": 3,
            "dtw_feature_groups_enabled": {
                "morphology": False,
                "movement": False,
                "intensity": False,
                "death": False,
                "contact": False,
            },
            'dtw_features_input': [],
            "dtw_features_resolved": [],
            'z_normalize': {},
        },
        "other": {
            "seed": 42,
            "umap_min_dist": 0.1,
            "umap_n_neighbors": 15,
            "nr_of_clusters": 5,
            "dtw_feature_groups_enabled": {
                "morphology": False,
                "movement": False,
                "intensity": False,
                "death": False,
                "contact": False,
            },
            'dtw_features_input': [],
            "dtw_features_resolved": [],
            'z_normalize': {},
        }
    },
    "backprojection": {
        "mode": "mean",
        "save": False,
        "last_sample": None,
        "feature_groups_enabled": {
            "morphology": False,
            "movement": False,
            "intensity": False,
            "death": False,
            "contact": False,
        },
        "columns_input": [],
        "columns_resolved": [],
    },
    "active_killing": {
        "observation_window": 5,
        "death_signal_column": "percentage_dead_mask",
        "killing_threshold_multiplier": 1.5,
        "min_contact_duration": 1,
        "use_absolute_threshold": False,
        "absolute_killing_threshold": None,
        "save_results": True,
    },
    "death_dynamics": {
        "organoid": {
            "dead_perc_threshold": 0.02
        },
        "immune": {
            "dead_perc_threshold": 0.25
        },
        "other": {
            "dead_perc_threshold": 0.10
        }
    }
}

behav3d_calculated_features = {
    "morphology": [
        "nr_pixels", "volume", "bbox_volume", "extent", "solidity",
        "equivalent_diameter", "major_axis_length", "minor_axis_length",
        "elongation", "surface_area", "sphericity", "convex_volume", "orientation_vector",
    ],
    "movement": [
        "displacement", "cumulative_displacement", "displacement_from_origin",
        "mean_square_displacement", "speed", "mean_speed",
    ],
    "intensity": [
        "mean_intensity_*", "mean_dead_dye",
    ],
    "death": [
        "percentage_dead_mask", "nr_dead_mask_pixels", "increase_dead_mask", "dead",
    ],
    "contact": [
        "*_contact", "*_contact_pixels", "touching_*", "active_*_contact",
    ],
    "active_killing": [
        "is_active_killing", "killing_efficiency",
    ],
}

spinning_loader = (
    '<div style="display:flex;align-items:center;gap:6px;">'
    '<svg width="18" height="18" viewBox="0 0 44 44" stroke="#1f77b4">'
    '  <g fill="none" fill-rule="evenodd" stroke-width="4">'
    '    <circle cx="22" cy="22" r="18" stroke-opacity=".2"></circle>'
    '    <path d="M40 22c0-9.94-8.06-18-18-18">'
    '      <animateTransform attributeName="transform" type="rotate"'
    '        from="0 22 22" to="360 22 22" dur="0.8s" repeatCount="indefinite"/>'
    '    </path>'
    '  </g>'
    '</svg>'
    '<span style="font-size:12px;color:#555;">Running…</span>'
    '</div>'
)

# ===============================
# HELPER FUNCTIONS
# ===============================
def _deep_merge(a: dict, b: dict) -> dict:
    out = deepcopy(a)
    for k, v in (b or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _deep_merge(out[k], v)
        else:
            out[k] = v
    return out

def _load_config(path: Path = CONFIG_PATH) -> dict:
    cfg = deepcopy(_DEFAULT_CONFIG)
    if path.exists():
        try:
            user = yaml.safe_load(path.read_text(encoding="utf-8"))
            if isinstance(user, dict):
                cfg = _deep_merge(cfg, user)
        except Exception as e:
            print(f"⚠️ Failed to load config {path}: {e}")
    return cfg

def _cfg_get(cfg: dict, dotted_key: str, default=None):
    cur = cfg
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur

CONFIG = _load_config()

def _mk_timepoint_range(use_all: bool, start: int, end: int):
    if use_all:
        return None
    return (int(start), int(end))

def detect_cell_type_category(cell_type, metadata):
    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)
    
    if cell_type in organoid_types:
        return 'organoid'
    elif cell_type in immune_types:
        return 'immune'
    elif cell_type in other_types:
        return 'other'
    else:
        if 'organoid' in cell_type.lower():
            return 'organoid'
        else:
            raise ValueError(
                f"Cell type '{cell_type}' not found in metadata. "
                f"Detected: organoid={organoid_types}, immune={immune_types}, other={other_types}."
            )

# ===============================
# BASE WIDGETS
# ===============================
class PathPicker(widgets.HBox):
    def __init__(
        self,
        mode='file',
        start_dir='.',
        default='',
        description='Path:',
        button_description='Browse…',
        placeholder='Type a path or click Browse…',
        filter_pattern=None,
        description_width='90px',
        width='100%',
    ):
        super().__init__()
        self._mode = mode
        self._start_dir = start_dir
        
        self.text = widgets.Text(
            value=default,
            placeholder=placeholder,
            description=description,
            style={'description_width': description_width},
            layout=widgets.Layout(flex='1 1 auto', width='auto')
        )
        
        self.btn = widgets.Button(
            description=button_description,
            icon='folder-open',
            layout=widgets.Layout(width='auto', flex='0 0 auto')
        )
        
        self.chooser_out = widgets.Output()
        self.chooser = None
        
        self.btn.on_click(self._toggle_chooser)
        self.text.observe(self._on_text_change, names='value')
        
        self.children = [self.text, self.btn, self.chooser_out]
        self.layout = widgets.Layout(width=width, display='flex', flex_flow='row wrap')

    @property
    def value(self):
        return self.text.value

    @value.setter
    def value(self, path: str):
        self.text.value = path

    @property
    def mode(self):
        return self._mode

    def _toggle_chooser(self, _=None):
        if self.chooser is None:
            with self.chooser_out:
                self.chooser = FileChooser(
                    self._start_dir,
                    select_default=True,
                    show_only_dirs=(self._mode == 'dir'),
                    filter_pattern=getattr(self, 'filter_pattern', None)
                )
                self.chooser.register_callback(self._on_select)
                display(self.chooser)
        else:
            self.chooser_out.clear_output()
            self.chooser = None

    def _on_select(self, chooser):
        path = chooser.selected_path if self._mode == 'dir' else chooser.selected
        if path:
            self.text.value = str(path)
        self._toggle_chooser()

    def _on_text_change(self, change):
        pass


DIM_ORDER_OPTIONS = [
    # 5D
    "TCZYX", "TZCYX", "CTZYX", "CZTYX", "ZCTYX", "ZTCYX",
    # 4D (missing one of T / C / Z)
    "CZYX", "ZCYX", "TZYX", "ZTYX", "TCYX", "CTYX",
    # 3D / 2D
    "ZYX", "CYX", "TYX",
    "YX",
]


def make_metadata_callback(metadata_loader, col):
    """Return an on_converted callback that writes *zarr_path* into *col* for the given sample."""
    def _cb(sample_name, src, zarr_path):
        md = metadata_loader.metadata
        if col not in md.columns:
            md[col] = pd.NA
        md[col] = md[col].astype("object")
        mask = md["sample_name"].astype(str) == str(sample_name)
        md.loc[mask, col] = str(zarr_path)
        csv_path = getattr(metadata_loader, "metadata_csv_path", None)
        if csv_path:
            md.to_csv(csv_path, index=False)
        print(f"  Metadata '{col}' for sample '{sample_name}' -> {zarr_path}")
    return _cb


def make_tracking_metadata_callback(metadata_loader, img_col, csv_col, trackdata_dir):
    """
    Return an on_converted callback for tracking imports.

    After the zarr is created, automatically generates the tracks CSV
    (same format as laptracking, trackpy, etc.) using the existing
    convert_tracked_image_to_csv(), and stores both paths in metadata.
    """
    def _cb(sample_name, src, zarr_path):
        md = metadata_loader.metadata

        # --- 1. Store zarr path ---
        for col in [img_col, csv_col]:
            if col not in md.columns:
                md[col] = pd.NA
            md[col] = md[col].astype("object")

        mask = md["sample_name"].astype(str) == str(sample_name)
        md.loc[mask, img_col] = str(zarr_path)

        # --- 2. Generate CSV (same as all tracking methods) ---
        row = md.loc[mask].iloc[0]
        el_xy = float(row.get("pixel_distance_xy", 1))
        el_z = float(row.get("pixel_distance_z", 1))

        csv_dir = Path(trackdata_dir)
        csv_dir.mkdir(parents=True, exist_ok=True)
        csv_path = csv_dir / f"{Path(zarr_path).stem.replace('_tracked', '_tracks')}.csv"

        print(f"  Generating tracks CSV: {csv_path}")
        convert_tracked_image_to_csv(
            img_path=zarr_path,
            outpath=csv_path,
            element_size_x=el_xy,
            element_size_y=el_xy,
            element_size_z=el_z,
        )
        md.loc[mask, csv_col] = str(csv_path)

        # --- 3. Save metadata ---
        csv_meta_path = getattr(metadata_loader, "metadata_csv_path", None)
        if csv_meta_path:
            md.to_csv(csv_meta_path, index=False)

        print(f"  Metadata '{img_col}' -> {zarr_path}")
        print(f"  Metadata '{csv_col}' -> {csv_path}")
    return _cb


class ExternalImageImporter(widgets.VBox):
    """
    Widget: pick a sample, browse/paste an image path, probe its shape,
    select axis order, and convert to Zarr.  Reuses PathPicker for the
    file input.  *on_converted(sample_name, src_path, zarr_path)* is
    called after a successful conversion so callers can update metadata.

    Parameters
    ----------
    label : str
        Display label (e.g. "tcell segmentation").
    output_dir : str
        Base output directory (zarr stored inside {output_dir}/{sample_name}/).
    on_converted : callable
        Callback(sample_name, src_path, zarr_path) after successful conversion.
    sample_names : list[str]
        Sample names for the dropdown.
    output_filename : str or None
        If set, a Python format string for the zarr filename using {sample_name}.
        E.g. "{sample_name}_tcell_segments.zarr".  If None, uses "{src_stem}.zarr".
    is_label_image : bool
        If True, cast output to uint16 (same as all internal segmentation/tracking
        methods). Raises ValueError if any value exceeds 65535.
    """

    def __init__(self, label="Image", output_dir=None, on_converted=None,
                 sample_names=None, output_filename=None, is_label_image=False):
        super().__init__()
        self._get_shape, self._get_dim_order = get_image_shape, get_image_dimension_order
        self._output_dir, self._on_converted = output_dir, on_converted
        self._output_filename = output_filename
        self._is_label_image = is_label_image
        self.last_zarr_path = None

        self.sample_dd = widgets.Dropdown(
            options=sample_names or [], description="Sample:",
            style={"description_width": "80px"}, layout=widgets.Layout(width="300px"))
        self.path_picker = PathPicker(
            description="File:", placeholder="Path to image (TIFF, H5, CZI, …)",
            description_width="50px", width="100%")
        self.btn_probe = widgets.Button(
            description="Probe", icon="search", button_style="info",
            layout=widgets.Layout(width="90px"))
        self.info = widgets.HTML(value="")
        self.axis_dd = widgets.Dropdown(
            options=["-- probe first --"], description="Axis order:",
            style={"description_width": "90px"}, layout=widgets.Layout(width="260px"))
        n_cpu = os.cpu_count() or 4
        self.n_workers = widgets.IntText(
            value=max(1, n_cpu // 2), description="Workers:",
            style={"description_width": "70px"}, layout=widgets.Layout(width="160px"))
        self.btn_convert = widgets.Button(
            description=f"Convert to Zarr", button_style="success",
            icon="exchange", disabled=True, layout=widgets.Layout(width="180px"))
        self.out = widgets.Output()

        self.btn_probe.on_click(self._on_probe)
        self.btn_convert.on_click(self._on_convert)

        self.children = [
            widgets.HTML(f"<b>Import external {label}</b>"),
            widgets.HBox([self.sample_dd, self.path_picker, self.btn_probe]),
            self.info,
            widgets.HBox([self.axis_dd, self.n_workers, self.btn_convert]),
            self.out,
        ]

    def _on_probe(self, _=None):
        self.out.clear_output()
        p = Path(self.path_picker.value.strip())
        if not p.exists():
            self.info.value = "<span style='color:red'>File not found.</span>"
            self.btn_convert.disabled = True
            return
        try:
            shape, detected = self._get_shape(p), self._get_dim_order(p)
        except Exception as exc:
            self.info.value = f"<span style='color:red'>Error: {exc}</span>"
            self.btn_convert.disabled = True
            return

        shape_str = " x ".join(str(s) for s in shape)
        det_str = detected or "not detected"
        self.info.value = (
            f"<b>Shape:</b> {shape_str} ({len(shape)}D) &nbsp;|&nbsp; "
            f"<b>Detected order:</b> {det_str}")

        matching = [o for o in DIM_ORDER_OPTIONS if len(o) == len(shape)] or list(DIM_ORDER_OPTIONS)
        self.axis_dd.options = matching
        self.axis_dd.value = detected if (detected and detected in matching) else matching[0]
        self.btn_convert.disabled = False

    def _on_convert(self, _=None):
        self.out.clear_output()
        self.btn_convert.disabled = True
        src = Path(self.path_picker.value.strip())
        if not src.exists():
            with self.out: print(f"File not found: {src}")
            self.btn_convert.disabled = False
            return

        sample_name = str(self.sample_dd.value)

        # Output path: {output_dir}/{sample_name}/{filename}.zarr
        base_dir = Path(self._output_dir) if self._output_dir else src.parent
        sample_dir = base_dir / sample_name
        sample_dir.mkdir(parents=True, exist_ok=True)

        if self._output_filename:
            zarr_name = self._output_filename.format(sample_name=sample_name)
        else:
            zarr_name = f"{src.stem}.zarr"
        zarr_path = sample_dir / zarr_name

        with self.out:
            try:
                if self._is_label_image:
                    convert_label_file_to_zarr(
                        path=src,
                        outpath=zarr_path,
                        axis_order=self.axis_dd.value,
                        overwrite=True,
                    )
                else:
                    convert_raw_file_to_zarr(
                        path=src,
                        outpath=zarr_path,
                        axis_order=self.axis_dd.value,
                        overwrite=True,
                        n_workers=max(1, int(self.n_workers.value)),
                    )

                # dtype check for label images (segmentation/tracking)
                # Prefer uint16 (standard for pipeline), fall back to int32
                # for very large label sets (>65535 objects).
                if self._is_label_image:
                    arr = load_image(zarr_path)
                    max_val = int(np.max(np.asarray(arr)))
                    if arr.dtype == np.uint16:
                        pass  # Already correct
                    elif max_val <= np.iinfo(np.uint16).max:
                        print(f"  Casting from {arr.dtype} to uint16 (max value {max_val} fits)")
                        save_as_zarr(np.asarray(arr).astype(np.uint16), zarr_path)
                    elif arr.dtype != np.int32:
                        print(f"  Casting from {arr.dtype} to int32 (max value {max_val} exceeds uint16)")
                        save_as_zarr(np.asarray(arr).astype(np.int32), zarr_path)

                self.last_zarr_path = str(zarr_path)
                if self._on_converted:
                    self._on_converted(sample_name, str(src), str(zarr_path))
                print(f"Zarr saved: {zarr_path}")
            except Exception:
                traceback.print_exc()
        self.btn_convert.disabled = False



class MetadataLoader(widgets.VBox):
    def __init__(
        self,
        metadata_path_picker,
        output_dir_picker,
        func=False,
        button_description: str = "Load metadata",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.load_behav3d_metadata = load_behav3d_metadata
        self.check_behav3d_metadata = check_behav3d_metadata
        
        self.metadata_path_picker = metadata_path_picker
        self.output_dir_picker = output_dir_picker
        self.func = func
        
        self.btn_load = widgets.Button(
            description=button_description,
            button_style="primary",
            icon="upload",
            layout=widgets.Layout(width="200px")
        )
        self.out = widgets.Output()
        self.btn_load.on_click(self._on_click)
        
        self.value = None
        self.path = None
        self.behav3d_parameters = _load_config()
        self.output_dir = self.output_dir_picker.value if self.output_dir_picker else ""
        
        self.children = [self.btn_load, self.out]

    def set_file_picker(self, file_picker: widgets.Widget):
        self.metadata_path_picker = file_picker

    def load(self, path=None):
        if path is None:
            path = self.metadata_path_picker.value
        
        if not path or not Path(path).exists():
            with self.out:
                print(f"❌ Error: Metadata file does not exist: {path}")
            return
        
        try:
            df = self.load_behav3d_metadata(path)
            with self.out:
                self.out.clear_output()
                self.check_behav3d_metadata(df, func=self.func)
            
            self.value = df
            self.path = path
            self.output_dir = self.output_dir_picker.value if self.output_dir_picker else str(Path(path).parent)
            
        except Exception as e:
            with self.out:
                print(f"❌ Error loading metadata: {e}")
                traceback.print_exc()

    def _on_click(self, _):
        self.load()
