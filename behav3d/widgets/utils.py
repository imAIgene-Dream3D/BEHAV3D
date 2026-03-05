import ipywidgets as widgets
from ipyfilechooser import FileChooser
from pathlib import Path
import yaml
from copy import deepcopy
import pandas as pd
import re
from traitlets import Any, Unicode, Bool

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
            }
        },
        "organoid": {
            "method": "propagation",
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
        "death_signal_column": "mean_dead_dye",
        "killing_threshold_multiplier": 1.5,
        "min_contact_duration": 1,
        "use_absolute_threshold": False,
        "absolute_killing_threshold": None,
        "save_results": True,
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
    from behav3d.core.metadata import (
        detect_organoid_types_from_metadata,
        detect_immune_cell_types_from_metadata,
        detect_other_cell_types_from_metadata
    )
    
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
        from behav3d.core.metadata import load_behav3d_metadata, check_behav3d_metadata
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
                import traceback
                traceback.print_exc()

    def _on_click(self, _):
        self.load()
