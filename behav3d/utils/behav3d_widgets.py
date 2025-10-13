# behav3d_widgets.py
import random

import ipywidgets as widgets
from ipyfilechooser import FileChooser
from IPython.display import display, clear_output
from behav3d.utils import load_behav3d_metadata, check_behav3d_metadata, expand_column_patterns
from behav3d.preprocessing import convert_input_files_to_zarr
import builtins
from pathlib import Path
from traitlets import Any, Unicode, Bool
import os
from behav3d .utils.fileio import load_image, get_image_shape, get_image_dimension_order
import asyncio
import pandas as pd
from behav3d.preprocessing.unmixing import visualize_unmix
from behav3d.preprocessing.unmixing.signal_unmixing import signal_unmixing
from behav3d.preprocessing.segmentation.napari_pixelclassifier import train_pixel_classifier, run_pixel_classifier_segmentation
import traceback
from behav3d.preprocessing.tracking import visualize_tracks
from behav3d.preprocessing.tracking.laptracking import run_tcell_laptracking
from behav3d.preprocessing.tracking.trackpy_tracking import run_tcell_trackpy_tracking
from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking

import json
from copy import deepcopy
import yaml
import fnmatch
from behav3d.analysis import summarize_track_features
from behav3d.analysis.feature_extraction import run_feature_extraction
from behav3d.analysis.tcell_analysis import filter_tcell_tracks, run_tcell_analysis 
from behav3d.analysis.organoid_analysis import filter_organoid_tracks, run_organoid_analysis 

from behav3d.analysis.backprojection import backproject_mean_features_behav3d, backproject_time_features_behav3d
import napari

behav3d_calculated_features = {
    "morphology": [
        "nr_pixels",
        "volume",
        "bbox_volume",
        "extent",
        "solidity",
        "equivalent_diameter",
        "major_axis_length",
        "minor_axis_length",
        "elongation",
        "surface_area",
        "sphericity",
        "convex_volume",
        "orientation_vector",
    ],
    "movement": [
        "displacement",
        "cumulative_displacement",
        "displacement_from_origin",
        "mean_square_displacement",
        "speed",
        "mean_speed",
    ],
    "intensity": [
        "mean_intensity_*",
        "mean_dead_dye",
    ],
    "death": [
        "percentage_dead_mask",
        "nr_dead_mask_pixels",
        "increase_dead_mask",
        "dead",
    ],
    "contact": [
        "organoid_contact",
        "organoid_contact_pixels",
        "tcell_contact",
        "tcell_contact_pixels",
        "active_tcell_contact",
    ],
}
# ===============================
# JSON CONFIG (loader + defaults)
# ===============================
CONFIG_PATH = Path("behav3d_config.yml")  # adjust if you prefer a different location

_DEFAULT_CONFIG = {
    "seed": 42,
    "paths": {
        "metadata_csv": r"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv",
        "output_dir": r"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE"
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
        "channel_colors": ["cyan", "yellow", "red", "green", "magenta", "blue"]
    },
    "pixel_classifier": {
        "examples_per_sample": 3,
        "sample_specific_classifier": False,
        "workers": 8,
        "organoid_edt_threshold": 12.0,
        "use_all_timepoints": True,
        "tp_start": 0,
        "tp_end": 0
    },
    "tracking": {
        "tcell": {
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
                "adaptive_stop": 5,
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
        },
        "organoid_2": {
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
        "channel_colors": ["cyan", "yellow", "red", "green", "magenta", "blue"]
    },
    "features": {
        "tcell": {
            "dead_mask_percentage_threshold": 0.25,
            "features_choice": ["movement", "intensity", "contact", "death"],
            "contact_threshold": 0,
            "n_workers": 16,
            "overwrite": False
        },
        "organoid": {
            "dead_mask_percentage_threshold": 0.02,
            "features_choice": ["intensity", "death", "morphology"],
            "contact_threshold": 0,
            "n_workers": 8,
            "overwrite": False
        },
        "organoid_2": {
            "dead_mask_percentage_threshold": 0.02,
            "features_choice": ["intensity", "death", "morphology"],
            "contact_threshold": 0,
            "n_workers": 8,
            "overwrite": False
        }
    },
    "track_filtering": {
        "tcell": {
            "exp_duration": 24.0,
            "exp_duration_enabled": False,
            "min_track_length": 0, # frames
            "min_track_length_enabled": True,
            "max_track_length": 999999, # frames
            "max_track_length_enabled": True,
            "filter_t0_dead": True,
            "filter_t0_dead_enabled": False
        },
        "organoid": {
            "exp_duration_enabled": False,
            "exp_duration": 999999,      # timepoints
            "min_track_length_enabled": False,
            "min_track_length": 100,       # frames
            "max_track_length_enabled": False,
            "max_track_length": 100,  # frames
            "min_size_enabled": True,
            "min_size": 1000,
        },
        "organoid_2": {
            "exp_duration_enabled": False,
            "exp_duration": 999999,      # timepoints
            "min_track_length_enabled": False,
            "min_track_length": 100,       # frames
            "max_track_length_enabled": False,
            "max_track_length": 100,  # frames
            "min_size_enabled": True,
            "min_size": 1000,
        }
    },
    "analysis": {
        "tcell": {
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
            'dtw_features_input': [
                'mean_square_displacement',
                'speed',
                'mean_dead_dye',
                'organoid_contact',
                'tcell_contact'
            ],      # patterns allowed (e.g. "mean_intensity_*")
            "dtw_features_resolved": [],   # expanded at run
            'z_normalize': {
                'mean_square_displacement': True,
                'speed': True,
                'mean_dead_dye': True,
                'organoid_contact': False,
                'tcell_contact': False
            },             # {feature_name: bool}
        },
        "organoid": {
            "dead_perc_threshold": 0.02
        },
        "organoid_2": {
            "dead_perc_threshold": 0.02
        }
    },
    "backprojection": {
        "mode": "mean",              # "mean" or "time"
        "save": False,               # if False, widget prevents writing .zarr files
        "last_sample": None,         # remember last picked sample
        "feature_groups_enabled": {
            "morphology": False,
            "movement": False,
            "intensity": False,
            "death": False,
            "contact": False,
        },
        "columns_input": [],         # patterns selected in the UI (e.g., "mean_intensity_*")
        "columns_resolved": [],      # expanded exact column names (filled at run)
    }
}

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
            # ignore malformed files; keep defaults
            print(f"⚠️ Failed to load config {path}: {e}")
    return cfg

def _cfg_get(cfg: dict, dotted_key: str, default=None):
    cur = cfg
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur

CFG = _load_config()

spinning_loader = (
                '<div style="display:flex;align-items:center;gap:6px;">'
                # circular spinner (animated via SVG, no CSS needed)
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
def _mk_timepoint_range(use_all: bool, start: int, end: int):
    if use_all:
        return None
    return (int(start), int(end))
class PathPicker(widgets.VBox):
    """
    Displayable widget for picking a file or directory.
    - mode: 'file' or 'dir'
    - start_dir: starting directory for the chooser
    - default: default filename (file mode) or starting path (dir mode)
    - filter_pattern: e.g. '*.csv' (file mode only)
    """
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
        assert mode in ('file', 'dir')
        self._mode = mode

        # Prefer JSON defaults if caller didn't provide one
        if not default:
            if mode == 'file':
                default = _cfg_get(CFG, "paths.metadata_csv", "") or ""
            else:  # dir
                default = _cfg_get(CFG, "paths.output_dir", "") or ""

        # Text box (the only persistent UI)
        self.text = widgets.Textarea(
            value=str(default or ''),
            description=description,
            placeholder=placeholder,
            style={'description_width': description_width},
            layout=widgets.Layout(width=width, height='32px'),
        )

        # Browse button
        self.button = widgets.Button(
            description=button_description, icon='folder-open',
            tooltip=f"Choose a {'folder' if mode=='dir' else 'file'}"
        )

        # Resolve starting directory & filename
        start_dir = start_dir or '.'
        filename = ''
        if mode == 'file' and default:
            p = Path(default)
            start_dir = str(p.parent if p.parent.exists() else start_dir)
            filename = p.name
            # if p.is_absolute():
            #     start_dir = str(p.parent if p.parent.exists() else start_dir)
            #     filename = p.name
            # else:
            #     filename = str(default)
        elif mode == 'dir' and default and os.path.isdir(default):
            start_dir = str(default)

        # FileChooser (hidden until button click)
        self.fc = FileChooser(path=start_dir, filename=filename)
        self.fc.title = f"<b>Select {'directory' if mode=='dir' else 'file'}</b>"
        self.fc.show_only_dirs = (mode == 'dir')
        self.fc.use_dir_icons = True
        if filter_pattern and mode == 'file':
            self.fc.filter_pattern = filter_pattern

        # Container: row with Text + Button (chooser injected below when needed)
        self._row = widgets.HBox([self.text, self.button])
        super().__init__([self._row])

        # Wire events
        self.button.on_click(self._toggle_chooser)
        self.fc.register_callback(self._on_select)
        self.text.observe(self._on_text_change, names='value')

    # ---- public API ----
    @property
    def value(self) -> str:
        """Current selected path (string)."""
        return self.text.value

    @value.setter
    def value(self, path: str):
        """Set the path programmatically (updates chooser location if possible)."""
        self.text.value = str(path or '')

    @property
    def mode(self) -> str:
        return self._mode

    # ---- internal helpers ----
    def _toggle_chooser(self, _=None):
        if self.fc in self.children:
            # hide
            self.children = tuple(c for c in self.children if c is not self.fc)
        else:
            # show
            self.children = tuple(list(self.children) + [self.fc])

    def _on_select(self, chooser):
        # Fill text and hide chooser
        if self._mode == 'dir':
            self.text.value = chooser.selected_path
        else:
            self.text.value = chooser.selected
        if self.fc in self.children:
            self.children = tuple(c for c in self.children if c is not self.fc)

    def _on_text_change(self, change):
        # Keep chooser in sync when user types/pastes a path
        newv = (change.get('new') or '').strip()
        if not newv:
            return
        try:
            if self._mode == 'dir':
                if os.path.isdir(newv):
                    self.fc.reset(path=newv)
            else:
                # file mode: split into directory + filename
                parent = os.path.dirname(newv) or '.'
                fname = os.path.basename(newv)
                if os.path.isdir(parent):
                    if fname:
                        self.fc.reset(path=parent, filename=fname)
                    else:
                        self.fc.reset(path=parent)
        except Exception:
            # Don't crash UI if reset fails; ignore silently
            pass
class MetadataLoader(widgets.VBox):
    """
    A VBox widget that wraps a file PathPicker and a 'Load' button.
    After clicking the button, `value` holds the loaded pandas DataFrame,
    and `path` holds the CSV path. Works without globals; read .value later.
    """
    _busy = Bool(default_value=False).tag(sync=False)
    _handler_bound = Bool(default_value=False).tag(sync=False)
    
    def __init__(
        self,
        metadata_path_picker,
        output_dir_picker,
        func=False,
        button_description: str = "Load metadata",
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.metadata = None
        self.output_dir = None
        self.metadata_csv_path = None

        self.func = func
        self.file_picker = metadata_path_picker
        self.output_dir_picker = output_dir_picker

        # initialize pickers with JSON defaults if empty
        try:
            if hasattr(self.file_picker, "value") and not self.file_picker.value:
                self.file_picker.value = _cfg_get(CFG, "paths.metadata_csv", "") or ""
        except Exception:
            pass
        try:
            if hasattr(self.output_dir_picker, "value") and not self.output_dir_picker.value:
                self.output_dir_picker.value = _cfg_get(CFG, "paths.output_dir", "") or ""
        except Exception:
            pass
        
        self.button = widgets.Button(description=button_description, button_style='success', layout=widgets.Layout(
        width="fit-content",   # size to content
        flex="0 0 auto"        # don't stretch in HBox/VBox
    ))
        self.out = widgets.Output()

        self._click_handler = self._on_click
        self.button.on_click(self._click_handler)
        # Build UI (file_picker may be None initially)
        children = []
        if self.output_dir_picker is not None:
            children.append(self.output_dir_picker)
        if self.file_picker is not None:
            children.append(self.file_picker)  
        children += [self.button, self.out]
        self.children = tuple(children)

    def set_file_picker(self, file_picker: widgets.Widget):
        """Attach/replace the file picker widget (must have a `.value` path)."""
        self.file_picker = file_picker

    def load(self, path = None):
        """Programmatic load (same as clicking the button)."""
        if path is None:
            if self.file_picker is None:
                raise ValueError("No file_picker attached and no path provided.")
            path = self.file_picker.value

        with self.out:
            clear_output(wait=True)
            if not path or not str(path).lower().endswith(".csv"):
                print("Please choose a .csv file.")
                return

            self.metadata = load_behav3d_metadata(path)
            self.output_dir = self.output_dir_picker.value if self.output_dir_picker is not None else ""
            self.metadata_csv_path = str(path)
            
            self.behav3d_parameters_path = Path(self.output_dir, "behav3d_parameters.yml")
            
            if self.behav3d_parameters_path.exists():
                self.behav3d_parameters = _load_config(path = self.behav3d_parameters_path)
            else:
                self.behav3d_parameters = deepcopy(_DEFAULT_CONFIG)
                
            self.behav3d_parameters["paths"]["metadata_csv"] = str(self.metadata_csv_path)
            self.behav3d_parameters["paths"]["output_dir"]   = str(self.output_dir)
            yaml.safe_dump(self.behav3d_parameters, self.behav3d_parameters_path.open("w"), sort_keys=False)
            
            check_behav3d_metadata(self.metadata, self.func)
            print("✅ Checks passed!")
            display(self.metadata)

            
    # Button handler
    def _on_click(self, _):
        if self._busy:
            return  # re-entrancy guard prevents double execution
        self._busy = True
        self.load()
        # try:
        #     self.load()
        # except Exception as e:
        #     with self.out:
        #         print(f"❌ Error: {e}")
        # finally:
        self._busy = False

def convert_zarr_button(metadata_loader, dim_order_widget):
    btn = widgets.Button(
        description="Convert to Zarr",
        button_style="success",
        icon="cogs",
        layout=widgets.Layout(
            width="fit-content",   # size to content
            flex="0 0 auto"        # don't stretch in HBox/VBox
        )
    )
    spinner = widgets.HTML(value=spinning_loader)  # uses the same spinner HTML used elsewhere
    spinner.layout.display = "none"                # hidden by default
    out = widgets.Output()

    # Button row: button + spinner
    row = widgets.HBox([btn, spinner], layout=widgets.Layout(align_items="center", gap="8px"))

    def _run_conversion(_):
        with out:
            out.clear_output()
            print("Starting conversion…")
            btn.disabled = True
            spinner.layout.display = None  # show spinner

            # Optionally persist current dim-order choices before converting
            try:
                dim_order_widget.write_dimorder_to_metadata()
            except Exception:
                pass

            # ensure `result` is always defined even if conversion errors
            result = metadata_loader.metadata
            try:
                result = convert_input_files_to_zarr(
                    metadata=metadata_loader.metadata,
                    output_dir=metadata_loader.output_dir
                )
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                metadata_loader.metadata = result
                try:
                    metadata_loader.metadata.to_csv(metadata_loader.metadata_csv_path, index=False)
                except Exception:
                    import traceback; traceback.print_exc()
                print("Done ✅")
                btn.disabled = False
                spinner.layout.display = "none"  # hide spinner

    btn.on_click(_run_conversion)
    # Return the row (button + spinner) and the log output area
    return widgets.VBox([row, out])

class DimOrderTable:
    DIM_ORDER_OPTIONS = [
        "TCZYX",
        "TZCYX",
        "ZCTYX",
        "ZTCYX",
        "CZTYX",
        "CTZYX",
    ]
    DEFAULT_ORDER = "TCZYX"

    def __init__(
        self,
        metadata_loader,
        sample_col="sample_name",
        path_col="raw_image_path",
        widths=("40%","30%","30%"),
        dim_col="dimension_order",
        auto_write=False
    ):
        self.metadata_loader = metadata_loader
        self.sample_col = sample_col
        self.path_col = path_col
        self.dim_col = dim_col
        self.auto_write = auto_write

        self._width_sample, self._width_shape, self._width_order = widths

        # State
        self._rows = []
        self._selections = {}

        # Precompute allowed tuples
        self._allowed_orders = self.DIM_ORDER_OPTIONS

        # UI skeleton
        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(description="Refresh", tooltip="Build/Update table from metadata")
        self._refresh_btn.on_click(self._on_refresh)

        _apply_default = str(_cfg_get(self.metadata_loader.behav3d_parameters, "dim_order.default_apply_all", self.DEFAULT_ORDER))
        if _apply_default not in self.DIM_ORDER_OPTIONS:
            _apply_default = self.DEFAULT_ORDER

        self._apply_all_dd = widgets.Dropdown(
            options=self.DIM_ORDER_OPTIONS,
            value=_apply_default,
            description="Apply to all:",
            style={'description_width': 'initial'},
            layout=widgets.Layout(width="350px", margin="6px 8px 6px 0")
        )
        self._apply_all_btn = widgets.Button(description="Apply", layout=widgets.Layout(width="100px"))
        self._apply_all_btn.on_click(lambda _: self.set_all(self._apply_all_dd.value))

        self._header = widgets.HBox([
            widgets.Label("Sample", layout=widgets.Layout(width=self._width_sample, overflow="hidden")),
            widgets.Label("Image shape", layout=widgets.Layout(width=self._width_shape, overflow="hidden")),
            widgets.Label("Dim order", layout=widgets.Layout(width=self._width_order, overflow="hidden")),
        ])
        self._table_body = widgets.VBox()
        self._out = widgets.Output()

        self.widget = widgets.VBox([
            self._status,
            widgets.HBox([self._refresh_btn]),
            self._header,
            self._table_body,
            widgets.HBox([self._apply_all_dd, self._apply_all_btn]),
            self._out
        ])

        # NEW: auto-load immediately if metadata already exists
        self._maybe_autoload()

    # NEW: helper that only shows “waiting” if metadata is None; otherwise loads
    def _maybe_autoload(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if df is None:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            return
        if not isinstance(df, pd.DataFrame):
            self._status.value = "<b>metadata_loader.metadata must be a pandas DataFrame.</b>"
            return
        missing = [c for c in (self.sample_col, self.path_col) if c not in df.columns]
        if missing:
            self._status.value = f"<b>Metadata missing columns: {missing}</b>"
            return
        self._build_rows_from(df)
        self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"

    # UPDATED: only show “waiting” when metadata is None; otherwise load/refresh
    def _on_refresh(self, _btn):
        df = getattr(self.metadata_loader, "metadata", None)
        if df is None:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            return
        if not isinstance(df, pd.DataFrame):
            self._status.value = "<b>metadata_loader.metadata must be a pandas DataFrame.</b>"
            return

        missing = [c for c in (self.sample_col, self.path_col) if c not in df.columns]
        if missing:
            self._status.value = f"<b>Metadata missing columns: {missing}</b>"
            return

        self._build_rows_from(df)
        self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
        
    def display(self):
        display(self.widget)

    def get_selections(self) -> dict:
        """Return dict: sample_name -> tuple like ('T','C','Z','Y','X')"""
        return dict(self._selections)

    def get_selections_str(self) -> dict:
        """Return dict: sample_name -> 'TCZYX'"""
        return {s: self._order_to_str(o) for s, o in self._selections.items()}

    def write_dimorder_to_metadata(self, col=None):
        col = col or self.dim_col
        df = getattr(self.metadata_loader, "metadata", None)
        if not isinstance(df, pd.DataFrame):
            raise ValueError("metadata_loader.metadata is not a DataFrame yet.")
        if df.empty:
            raise ValueError("metadata_loader.metadata is empty.")
        if self.sample_col not in df.columns:
            raise ValueError(f"metadata is missing '{self.sample_col}' column.")
        # Map selections by sample name
        str_map = self.get_selections_str()
        self.metadata_loader.metadata[col] = df[self.sample_col].map(str_map)
        return self.metadata_loader.metadata

    def set_all(self, order_tuple):
        """Programmatically set all rows to the given order tuple or string."""
        if isinstance(order_tuple, str):
            order_tuple = tuple(order_tuple)
        for row in self._rows:
            row['dd'].value = order_tuple  # triggers observers

    def refresh_shapes(self):
        """Recompute shape labels (e.g., after mounting a drive)."""
        for row in self._rows:
            row['shape_label'].value = self._probe_image_shape(row['path'])

    def to_dataframe(self):
        data = []
        for row in self._rows:
            tup = row["dd"].value
            data.append({
                "sample_name": row["sample"],
                "raw_image_path": row["path"],
                "image_shape": row["shape_label"].value,
                "dim_order": tup,
                self.dim_col: self._order_to_str(tup),
            })
        return pd.DataFrame(data)

    def _on_refresh(self, _btn):
        df = getattr(self.metadata_loader, "metadata", None)
        if not (isinstance(df, pd.DataFrame) and not df.empty):
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            return
        missing = [c for c in (self.sample_col, self.path_col) if c not in df.columns]
        if missing:
            self._status.value = f"<b>Metadata missing columns: {missing}</b>"
            return

        self._build_rows_from(df)
        self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"

    def _build_rows_from(self, df: pd.DataFrame):
        # Reset state
        self._rows.clear()
        self._selections.clear()
        row_boxes = []

        # Optional: prefill from existing string column if present
        prefill = None
        if self.dim_col in df.columns:
            prefill = df[self.dim_col].astype(str).to_dict()  # index -> string

        for idx, r in df.iterrows():
            sample = str(r[self.sample_col])
            path = r[self.path_col]

            probed_shape = self._probe_image_shape(path)
            sample_lbl = widgets.Label(sample, layout=widgets.Layout(width=self._width_sample, overflow="hidden"))
            shape_lbl  = widgets.Label(probed_shape, layout=widgets.Layout(width=self._width_shape, overflow="hidden"))

            # Determine default dropdown value
            default_val = get_image_dimension_order(path)
            if default_val is None:
                default_val = self.DEFAULT_ORDER
            if prefill is not None:
                s = prefill.get(idx, None)
                if isinstance(s, str):
                    tu = tuple(s.upper())
                    if tu in self._allowed_orders:
                        default_val = tu

            dd = widgets.Dropdown(
                options=self.DIM_ORDER_OPTIONS,
                value=default_val,
                layout=widgets.Layout(width=self._width_order),
            )

            # keep live selection
            self._selections[sample] = dd.value
            dd.observe(lambda ch, s=sample: self._on_dd_changed(s, ch), names='value')

            self._rows.append({"sample": sample, "path": path, "shape_label": shape_lbl, "dd": dd})
            row_boxes.append(widgets.HBox([sample_lbl, shape_lbl, dd]))

        # Swap UI body
        self._table_body.children = tuple(row_boxes)

    def _on_dd_changed(self, sample_name, change):
        if change.get('name') == 'value':
            new_tuple = change['new']
            self._selections[sample_name] = new_tuple
            if self.auto_write:
                df = getattr(self.metadata_loader, "metadata", None)
                if isinstance(df, pd.DataFrame) and self.sample_col in df.columns:
                    mask = df[self.sample_col].astype(str) == str(sample_name)
                    df.loc[mask, self.dim_col] = self._order_to_str(new_tuple)

    def _probe_image_shape(self, path):
        path = Path(path)
        if not path.exists():
            return "⛔ missing"
        try:
            shape = get_image_shape(path)
            return "×".join(map(str, shape)) if shape else "unknown"
        except Exception as e:
            return f"⚠️ {type(e).__name__}"

    @staticmethod
    def _order_to_str(order_tuple):
        return "".join(order_tuple)
    


class SignalUnmixingPanel:
    def __init__(
            self, 
            metadata_loader,
            default_time_range=None,
            channel_colors=None):
        
        self.metadata_loader = metadata_loader
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing", {})

        # -------- Run params --------
        self.sink_channel = widgets.IntText(
            description="Sink channel",
            value=int(pc.get("sink_channel", 2))
        )
        self.source_channels = widgets.Text(
            value=pc.get("source_channels", "0,1"), 
            description="Source channels",
            style={"description_width": "initial"}, 
        )
        self.train_z = widgets.Text(
            value=pc.get("train_z", "4,5,6"), 
            description="Train Z", 
        )
        self.train_t = widgets.Text(
            value=pc.get("train_t", "1,3,6"), 
            description="Train T", 
        )
        self.bg_percentile = widgets.IntText(
            value=int(pc.get("bg_percentile", 1)), 
            description="BG percentile", 
        )
        self.median_size = widgets.IntText(
            value=int(pc.get("median_size", 3)), 
            description="Median size",
        )
        self.gaussian_sigma = widgets.IntText(
            value=int(pc.get("gaussian_sigma", 4)), 
            description="Gaussian sigma",
            style={"description_width": "initial"}, 
        )  
        self.use_all_timepoints = widgets.Checkbox(
            description="Process ALL timepoints",
            value=bool(pc.get("use_all_timepoints", True))
        )
        self.tp_n = widgets.IntText(
            description="Number of Timepoints", 
            value=int(pc.get("timepoints", 0)),
            style={"description_width": "initial"}  # Show full description
        )

        # Start/End visibility
        self.use_all_timepoints.observe(self._toggle_timepoint_inputs, names='value')
        self._toggle_timepoint_inputs()

        self.training_log = widgets.Checkbox(
            description="Show Training Log",
            value=bool(pc.get("training_log", True))
        )

        # -------- sum_scale table --------
        sum_scale_widgets = []
        rows = []
        file_list = self.metadata_loader.metadata.sample_name.to_list()

        # Header (first row)
        header = widgets.HBox([
            widgets.Label(value="File Name", layout=widgets.Layout(width='200px')),
            widgets.Label(value="sum_scale (source1, source2)", layout=widgets.Layout(width='200px')),
            widgets.Label(value="")  # Placeholder for the button column
        ])
        rows.append(header) 

        # Shared function to apply all
        def apply_to_all_callback(btn):
            first_value = sum_scale_widgets[0]['input'].value
            for i, row in enumerate(sum_scale_widgets):
                if i == 0:
                    continue  # Skip first
                row['input'].value = first_value
                # row['input'].style = {'description_width': 'initial', 'background': '#f0f0f0'}  # grey

        sum_scale_config = pc.get("sum_scale", "20,30")

        for i, file in enumerate(file_list):
            file_label = widgets.Label(value=file, layout=widgets.Layout(width='200px'))

            # Determine the value for this file
            if isinstance(sum_scale_config, dict):
                file_sum = sum_scale_config.get(file, [20, 30])  # fallback default if missing
                value_str = ", ".join(str(x) for x in file_sum)
            elif isinstance(sum_scale_config, str):
                value_str = sum_scale_config  # global default
            else:
                value_str = "20,30"  # final fallback

            self.sum_scale_input = widgets.Text(
                value=value_str,
                description='',
                layout=widgets.Layout(width='150px')
            )

            self.sum_scale_input.style = {'description_width': 'initial'}

            if i == 0:
                # Add "Apply to all" button only on the first row
                self.apply_button = widgets.Button(description="Apply to all", layout=widgets.Layout(width='120px'))
                self.apply_button.on_click(apply_to_all_callback)
            else:
                self.apply_button = widgets.Label(value="")  # Empty placeholder

            # Track widgets
            sum_scale_widgets.append({
                'label': file_label,
                'input': self.sum_scale_input,
                'button': self.apply_button
            })

            # Add the row as an HBox
            row = widgets.HBox([file_label, self.sum_scale_input, self.apply_button])
            rows.append(row)

        # Store for later use
        self.sum_scale_widgets = sum_scale_widgets


        self.btn_unmix = widgets.Button(
            description="Run signal unmixing",
            button_style="primary",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )

        # Uses the same global spinner HTML you use elsewhere
        self.spinner_unmix = widgets.HTML(value=spinning_loader)
        self.spinner_unmix.layout.display = "none"

        # Row: Unmix | spinner
        self.unmix_row = widgets.HBox(
            [self.btn_unmix, self.spinner_unmix],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        self.channel_colors = tuple(channel_colors or _cfg_get(
            self.metadata_loader.behav3d_parameters, "signal_unmixing.channel_colors",
            ["cyan", "yellow", "red", "green", "magenta", "blue"]
        ))

        self._viewer = None
        
        # --- UI: status + refresh
        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(
            description="Refresh",
            tooltip="Build/Update selector from metadata_loader.metadata",
        )
        self._refresh_btn.on_click(self._on_refresh_clicked)

        # --- Main controls (disabled until metadata present)
        self.sample_dropdown = widgets.Dropdown(
            options=[],
            value=None,
            description="Sample:",
            layout=widgets.Layout(width="350px"),
            disabled=True,
        )

        # Tickbox to enable/disable time range (from JSON)
        self.use_range = widgets.Checkbox(
            description="Use custom time range",
            value=bool(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.use_range", False)),
            indent=False,
        )

        # Start/End boxes (defaults from constructor OR JSON)
        if isinstance(default_time_range, (tuple, list)) and len(default_time_range) == 2:
            _start_default, _end_default = map(int, default_time_range)
        else:
            _start_default = int(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.start_t", 0))
            _end_default   = int(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.end_t", 0))

        self.start_t = widgets.IntText(
            description="Start T:",
            value=_start_default,
            layout=widgets.Layout(width="180px")
        )
        self.end_t   = widgets.IntText(
            description="End T:",
            value=_end_default,
            layout=widgets.Layout(width="180px")
        )

        self.range_box = widgets.HBox([self.start_t, self.end_t])
        self.range_box.layout.display = "flex" if self.use_range.value else "none"

        def _toggle_range_visibility(change):
            self.range_box.layout.display = "flex" if change["new"] else "none"
        self.use_range.observe(_toggle_range_visibility, names="value")

        self.open_button = widgets.Button(
            description="Visualize Signal Unmixing results",
            button_style="primary",
            tooltip="Launch Napari for the selected sample",
            icon="eye",
            layout=widgets.Layout(width="300px"),
            disabled=True,
        )

        self.close_button = widgets.Button(
            description="Close viewer",
            button_style="danger",
            icon="stop",
            tooltip="Close the active Napari viewer",
            layout=widgets.Layout(width="200px", display="none"),  # hidden by default
        )
        

        # Try immediate build if metadata already present
        self._maybe_build_from_loader()

        self.btn_save = widgets.Button(
            description="Save Signal Unmixing",
            button_style="success",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        self.use_unmix_path = widgets.Checkbox(
            description="Use Unmix as raw path",
            value=bool(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.use_unmix_path", False)),
            indent=False,
        )
        self.spinner_save = widgets.HTML(value=spinning_loader)
        self.spinner_save.layout.display = "none"
        self.save_row = widgets.HBox(
            [self.btn_save, self.spinner_save],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        # Wire handlers
        self.btn_unmix.on_click(self._on_btn_unmix_clicked)
        self.open_button.on_click(self._on_open_clicked)
        self.close_button.on_click(self._on_close_clicked)
        self.btn_save.on_click(self._on_save_clicked)

        self.out = widgets.Output()

        # Layout
        unmix_box = widgets.VBox([
            widgets.HTML("<b>Run Signal Unmixing</b>"),
            widgets.VBox(rows),
            widgets.HBox([self.sink_channel, self.source_channels]),
            widgets.HBox([self.train_z, self.train_t]),
            widgets.HBox([self.bg_percentile, self.median_size, self.gaussian_sigma]),
            widgets.HBox([self.training_log]),
            widgets.HBox([self.use_all_timepoints, self.tp_n]),
            self.unmix_row,
        ])

        visualize_box = widgets.VBox([
                widgets.HTML("<b>Visualize Signal Unmixing Result</b>"),
                widgets.HBox([self._status, self._refresh_btn]),
                self.sample_dropdown,
                self.use_range,
                self.range_box,
                widgets.HBox([self.open_button, self.close_button]),
        ])
        
        save_box = widgets.VBox([
            widgets.HTML("<b>Save Signal Unmixing Result</b>"),
            widgets.HTML("Modifies the tcell channel to  the new unmixed channel and saves a new metadata file."),
            self.use_unmix_path,
            self.save_row,
        ])


        self.ui = widgets.VBox([unmix_box, widgets.HTML("<hr>"), visualize_box, widgets.HTML("<hr>"), save_box, widgets.HTML("<hr>"), self.out])

        ####### Viewer section
        
        # viewer handle (if training opens napari and returns a viewer)
        self._viewer = None

 #######################BUTTON FUNCTIONS####################
    def display(self):
        display(self.ui)


    def _toggle_timepoint_inputs(self, change=None):
        show = not self.use_all_timepoints.value
        disp = None if show else 'none'
        self.tp_n.layout.display = disp
        self.tp_n.disabled = not show
        self.tp_n.value = 0 if self.use_all_timepoints.value else self.tp_n.value

    def _persist_params(self):
        self.metadata_loader.behav3d_parameters.setdefault("signal_unmixing", {})
        pc = self.metadata_loader.behav3d_parameters["signal_unmixing"]
        pc["sink_channel"] = int(self.sink_channel.value)
        pc["source_channels"] = str(self.source_channels.value)
        pc["train_z"] = str(self.train_z.value)
        pc["train_t"] = str(self.train_t.value)
        pc["bg_percentile"] = int(self.bg_percentile.value)
        pc["median_size"] = int(self.median_size.value)
        pc["gaussian_sigma"]   = int(self.gaussian_sigma.value)
        pc["training_log"] = bool(self.training_log.value)
        pc["tp_n"] = int(self.tp_n.value)
        pc["use_all_timepoints"] = bool(self.use_all_timepoints.value)
        pc["timepoints"] = int(self.tp_n.value)
        # Save per-file sum_scale values
        sum_scale = {}
        for row in self.sum_scale_widgets:
            filename = row['label'].value
            input_str = row['input'].value
            try:
                # Parse as a pair of integers
                parts = [int(s.strip()) for s in input_str.split(",") if s.strip().isdigit()]
                if len(parts) == 2:
                    sum_scale[filename] = parts
                else:
                    print(f"⚠️ Skipping file '{filename}' — invalid sum_scale: '{input_str}'")
            except Exception as e:
                print(f"⚠️ Error parsing sum_scale for '{filename}': {e}")

        pc["sum_scale"] = sum_scale

        yaml.safe_dump(
            self.metadata_loader.behav3d_parameters,
            self.metadata_loader.behav3d_parameters_path.open("w"),
            sort_keys=False
        )

    def _on_btn_unmix_clicked(self, _):
        with self.out:
            self.out.clear_output()
            self._lock(True)
            # Save previous parameters to compare
            previous_pc = deepcopy(_cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing", {}))
            self._persist_params()
            try:
                odir = Path(self.metadata_loader.output_dir).expanduser()
                odir.mkdir(parents=True, exist_ok=True)
                sum_scale_per_file = {}
                for row in self.sum_scale_widgets:
                    filename = row['label'].value
                    val = row['input'].value
                    try:
                        parts = [int(s.strip()) for s in val.split(",") if s.strip().isdigit()]
                        if len(parts) == 2:
                            sum_scale_per_file[filename] = parts
                        else:
                            print(f"⚠️ Skipping file '{filename}' due to invalid sum_scale: {val}")
                    except Exception as e:
                        print(f"⚠️ Error parsing sum_scale for file '{filename}': {e}")


                print(f"▶️ Applying signal unmixing to channel {self.sink_channel.value}…", flush=True)
                print(f"  output_dir={odir}\\images\\SignalUnmixing", flush=True)
                print(f"  sum_scales=\n{sum_scale_per_file}", flush=True)
                print(f"  sink_channel={self.sink_channel.value}")
                print(f"  source_channels={self.source_channels.value}")
                print(f"  train_z={self.train_z.value}")
                print(f"  train_t={self.train_t.value}")
                print(f"  bg_percentile={self.bg_percentile.value}")
                print(f"  median_size={self.median_size.value}")
                print(f"  gaussian_sigma={self.gaussian_sigma.value}")
                print(f"  tp_n={self.tp_n.value}")
                print(f"  training_log={self.training_log.value}", flush=True)

                self.spinner_unmix.layout.display = None

                new_md = signal_unmixing(
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(odir),
                    sink_channel=int(self.sink_channel.value),
                    source_channels=[int(s.strip()) for s in str(self.source_channels.value).split(",") if s.strip().isdigit()],
                    sum_scale=sum_scale_per_file, # dict: filename -> [int, int]
                    train_z=[int(s.strip()) for s in str(self.train_z.value).split(",") if s.strip().isdigit()],
                    train_t=[int(s.strip()) for s in str(self.train_t.value).split(",") if s.strip().isdigit()],
                    bg_percentile=int(self.bg_percentile.value),    
                    median_size=int(self.median_size.value),
                    gaussian_sigma=int(self.gaussian_sigma.value),
                    training_log=bool(self.training_log.value),
                    timepoints=int(self.tp_n.value),
                    previous_pc=previous_pc,
                )
                try:
                    if new_md is not None:
                        self.metadata_loader.metadata = new_md
                        new_md.to_csv(self.metadata_loader.metadata_csv_path, index=False)
                except Exception:
                    import traceback; traceback.print_exc()

                print("✅ Unmixing finished.", flush=True)
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self.spinner_unmix.layout.display = "none"
                self._lock(False)


    ######### CHECKKKKKKK##########################``

    # ---------------- Public API ----------------
    def display(self):
        display(self._panel)

    def get_selected_row(self) -> pd.Series:
        self._ensure_metadata_ready()
        name = self.sample_dropdown.value
        df = self.metadata_loader.metadata
        row = df[df["sample_name"].astype(str) == str(name)]
        if row.empty:
            raise KeyError(f"sample_name '{name}' not found in metadata.")
        return row.iloc[0]
    
    
    def open_viewer(self):
        row = self.get_selected_row()

        # only use range if tickbox checked
        if self.use_range.value:
            start_t, end_t = int(self.start_t.value), int(self.end_t.value)
            if end_t < start_t:
                start_t, end_t = end_t, start_t
            time_range = (start_t, end_t)
        else:
            time_range = None

        with self.out:
            self.out.clear_output()
            try:
                print(f"Opening viewer for: {row['sample_name']}")
                self._viewer = visualize_unmix(
                    metadata_row=row,
                    timepoint_range=time_range,
                    channel_colors=self.channel_colors,
                )
            except Exception as e:
                print(f"Error: {e}")
            self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"

    # ---------------- Internal helpers ----------------
    def _ensure_metadata_ready(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if not isinstance(df, pd.DataFrame) or df.empty:
            raise RuntimeError("Metadata not loaded yet. Click 'Refresh' once metadata_loader.metadata is set.")
   
    def _maybe_build_from_loader(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if isinstance(df, pd.DataFrame) and not df.empty and ("sample_name" in df.columns):
            self._build_from_metadata(df)
            self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
        else:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            self.sample_dropdown.disabled = True
            self.open_button.disabled = True

    def _build_from_metadata(self, df: pd.DataFrame):
        sample_names = df["sample_name"].astype(str).unique().tolist()
        if not sample_names:
            self._status.value = "<b>No sample_name values found in metadata.</b>"
            self.sample_dropdown.options = []
            self.sample_dropdown.value = None
            self.sample_dropdown.disabled = True
            self.open_button.disabled = True
            return

        desired = _cfg_get(self.metadata_loader.behav3d_parameters, "signal_unmixing.sample_name", None)
        self.sample_dropdown.options = sample_names
        self.sample_dropdown.value = desired if (desired in sample_names) else sample_names[0]
        self.sample_dropdown.disabled = False
        self.open_button.disabled = False

    
    def _lock(self, state: bool):
        # keep close_button enabled so user can close at any time
        for w in [
            self.sample_dropdown, self.sink_channel, self.sum_scale_input, self.apply_button,
            self.source_channels, self.gaussian_sigma, self.median_size, self.bg_percentile,
            self.btn_save, self.btn_unmix, self.use_all_timepoints, self.tp_n, self.use_range, self.train_t, self.end_t,
            self.start_t, self.train_z, self.training_log, self.open_button, self.close_button, self.range_box, self.save_row, self.use_unmix_path
        ]:
            w.disabled = state

        # # also lock/unlock the path pickers
        # try:
        #     for p in [self.clf_dir, self.clf_org_path, self.clf_tcell_path, self.clf_death_path, self.clf_org2_path]:
        #         p.text.disabled = state
        #         p.button.disabled = state
        # except Exception:
        #     pass

    # ---------------- Callbacks ----------------
    def _on_open_clicked(self, _):
        self._lock(True)
        self.open_viewer()
        self._lock(False)

    def _on_close_clicked(self, _):
        # Only relevant when we’re in non-blocking mode and hold a viewer
        with self.out:
            try:
                if self._viewer is not None:
                    print("Closing viewer…")
                    self._viewer.close()
                    self._viewer = None
                else:
                    print("No active viewer to close (viewer was likely opened in blocking mode).")
            finally:
                # Hide the button regardless
                self.close_button.layout.display = "none"
                
    def _on_refresh_clicked(self, _):
        with self.out:
            self.out.clear_output()
            try:
                self._maybe_build_from_loader()
            except Exception as e:
                print(f"Refresh failed: {e}")

    def _on_save_clicked(self, _):
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                metadata = self.metadata_loader.metadata
                for idx, sample in metadata.iterrows():
                    unmix_path = Path(sample["signal_unmixing_image_path"]).expanduser()
                    image = load_image(unmix_path)
                    # Unmixed channel is always the last one
                    new_tcell_channel = image.shape[1] - 1  # assuming shape is (T, C, Z, Y, X)
                    print(f"Sample '{sample['sample_name']}': setting tcell_channel to {new_tcell_channel}") 
                    metadata.at[idx, 'tcell_channel'] = new_tcell_channel
                    if self.use_unmix_path:
                        if sample["signal_unmixing_image_path"] != metadata.at[idx, 'raw_image_path']:
                            metadata.at[idx, 'original_raw_image_path'] = sample["raw_image_path"]
                            metadata.at[idx, 'raw_image_path'] = sample["signal_unmixing_image_path"]
                    try:
                        if metadata is not None:
                            self.metadata_loader.metadata = metadata
                            metadata.to_csv(self.metadata_loader.metadata_csv_path, index=False)
                    except Exception:
                        import traceback; traceback.print_exc()

                print(f"✅ Saved updated metadata")
            except Exception as e:
                print(f"Error saving metadata: {e}")


class PixelClassifierPanel:
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})

        # -------- Train controls --------
        self.examples_per_sample = widgets.IntText(
            description="Examples per sample",
            value=int(pc.get("examples_per_sample", 3))
        )
        self.sample_specific_classifier = widgets.Checkbox(
            description="Sample-specific classifier",
            value=bool(pc.get("sample_specific_classifier", False))
        )        
        
        # Checkbox for two organoid types (e.g., WT vs. KO)
        self.two_org_types = widgets.Checkbox(
            description="Segment 2 organoid types",
              value=bool(pc.get("two_org_types", False))  
        )
        self.n_workers = widgets.IntText(
            description="Workers",
            value=int(pc.get("workers", (os.cpu_count() or 8))),
            max=max(8, (os.cpu_count() or 8))
        )

        # ---- Manual classifier toggle ----
        self.manual_clf_paths = widgets.Checkbox(
            description="Manually supply classifiers",
            value=bool(pc.get("manual_clf_paths", False)),
            indent=False
        )

        # -------- Classifier path pickers (default EMPTY) --------
        self.clf_dir = PathPicker(
            mode='dir',
            description='Classifier dir:',
            default="", 
            description_width='160px',
            width='100%',
        )
        self.clf_org_path = PathPicker(
            mode='file',
            description='Organoid clf:',
            default="", 
            filter_pattern='*.joblib',
            description_width='160px',
            width='100%',
        )
        self.clf_org2_path = PathPicker(
            mode='file',
            description='Organoid 2 clf:',
            default="", 
            filter_pattern='*.joblib',
            description_width='160px',
            width='100%',
        )
        self.clf_tcell_path = PathPicker(
            mode='file',
            description='T-cell clf:',
            default="",   
            filter_pattern='*.joblib',
            description_width='160px',
            width='100%',
        )
        self.clf_death_path = PathPicker(
            mode='file',
            description='Death clf:',
            default="",         
            filter_pattern='*.joblib',
            description_width='160px',
            width='100%',
        )

        # If user previously saved manual paths AND toggle is True, restore them now.
        if self.manual_clf_paths.value:
            self.clf_dir.value        = str(pc.get("clf_dir", "") or "")
            self.clf_org_path.value   = str(pc.get("clf_org_path", "") or "")
            self.clf_tcell_path.value = str(pc.get("clf_tcell_path", "") or "")
            self.clf_death_path.value = str(pc.get("clf_death_path", "") or "")
            if self.two_org_types.value:
                self.clf_org2_path.value   = str(pc.get("clf_org2_path", "") or "")

        # -------- Apply controls --------
        self.organoid_edt_threshold = widgets.FloatText(
            description="Organoid EDT thr",
            value=float(pc.get("organoid_edt_threshold", 12.0)),
            style={'description_width': '160px'}
        )
        self.use_all_timepoints = widgets.Checkbox(
            description="Process ALL timepoints",
            value=bool(pc.get("use_all_timepoints", True))
        )
        self.tp_start = widgets.IntText(description="Start t", value=int(pc.get("tp_start", 0)))
        self.tp_end   = widgets.IntText(description="End t",   value=int(pc.get("tp_end", 0)))

        # Start/End visibility
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
            layout=widgets.Layout(width="200px", display="none")  # hidden until a viewer is open
        )

        self.spinner_train = widgets.HTML(value=spinning_loader)
        self.spinner_train.layout.display = "none"

        # Row: Train | Close viewer | spinner
        self.train_row = widgets.HBox(
            [self.btn_train, self.close_button, self.spinner_train],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        # Apply button + spinner
        self.btn_apply = widgets.Button(
            description="Apply pixel classifier",
            button_style="success",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        self.spinner_apply = widgets.HTML(value=spinning_loader)
        self.spinner_apply.layout.display = "none"
        self.apply_row = widgets.HBox(
            [self.btn_apply, self.spinner_apply],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        # Wire handlers
        self.btn_train.on_click(self._on_train_clicked)
        self.close_button.on_click(self._on_close_clicked)
        self.btn_apply.on_click(self._on_apply_clicked)

        self.out = widgets.Output()

        clf_path_widgets = [
            widgets.HTML("<b>Classifier paths</b>"),
            self.clf_dir,
            self.clf_org_path,
            self.clf_tcell_path,
            self.clf_death_path,
        ]

        # Conditionally add the second organization path
        if self.two_org_types.value:
            clf_path_widgets.append(self.clf_org2_path)

        self.clf_paths_box = widgets.VBox(clf_path_widgets)
        
        # Show/hide according to the checkbox
        self.two_org_types.observe(self._on_two_org_types_changed, names='value')
        self.manual_clf_paths.observe(self._toggle_clf_path_section, names='value')

        self._build_clf_paths_box()
        self._toggle_clf_path_section()

        # When directory changes, auto-fill file pickers (only when manual mode is enabled)
        def _dir_changed(change):
            if not self.manual_clf_paths.value:
                return
            newd = (change.get('new') or '').strip()
            if newd:
                self._apply_dir_to_clf_paths(newd)
        self.clf_dir.text.observe(_dir_changed, names='value')

        # Layout
        train_box = widgets.VBox([
            widgets.HTML("<b>Train pixel classifier</b>"),
            widgets.HBox([self.examples_per_sample, self.n_workers]),
            self.sample_specific_classifier,
            self.two_org_types,
            self.train_row,
        ])

        self.tp_row = widgets.HBox([self.use_all_timepoints, self.tp_start, self.tp_end])

        apply_box = widgets.VBox([
            widgets.HTML("<b>Apply segmentation</b>"),
            self.manual_clf_paths,               
            self.clf_paths_box,   
            widgets.HTML("<b>Segmentation setting</b>"),  
            self.organoid_edt_threshold,
            self.tp_row,
            self.apply_row,
        ])

        self.ui = widgets.VBox([train_box, widgets.HTML("<hr>"), apply_box, widgets.HTML("<hr>"), self.out])

        # viewer handle (if training opens napari and returns a viewer)
        self._viewer = None

    # ---------- config ----------
    def _persist_params(self):
        self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc = self.metadata_loader.behav3d_parameters["pixel_classifier"]
        pc["examples_per_sample"] = int(self.examples_per_sample.value)
        pc["sample_specific_classifier"] = bool(self.sample_specific_classifier.value)
        pc["workers"] = int(self.n_workers.value)
        pc["organoid_edt_threshold"] = float(self.organoid_edt_threshold.value)
        pc["use_all_timepoints"] = bool(self.use_all_timepoints.value)
        pc["tp_start"] = int(self.tp_start.value)
        pc["tp_end"]   = int(self.tp_end.value)
        pc["two_org_types"] = bool(self.two_org_types.value)

        pc["manual_clf_paths"] = bool(self.manual_clf_paths.value)
        if self.manual_clf_paths.value:
            pc["clf_dir"]        = str(self.clf_dir.value or "")
            pc["clf_org_path"]   = str(self.clf_org_path.value or "")
            pc["clf_tcell_path"] = str(self.clf_tcell_path.value or "")
            pc["clf_death_path"] = str(self.clf_death_path.value or "")
            if self.two_org_types.value:
                pc["clf_org2_path"]   = str(self.clf_org2_path.value or "")


        yaml.safe_dump(
            self.metadata_loader.behav3d_parameters,
            self.metadata_loader.behav3d_parameters_path.open("w"),
            sort_keys=False
        )

    # ---------- UI helpers ----------
    def _toggle_timepoint_inputs(self, change=None):
        show = not self.use_all_timepoints.value
        disp = None if show else 'none'
        self.tp_start.layout.display = disp
        self.tp_end.layout.display   = disp
        self.tp_start.disabled = not show
        self.tp_end.disabled   = not show

    # >>> helper: compute default file paths for a directory
    def _default_clf_paths_for_dir(self, d: str):
        if not d:
            if self.two_org_types.value:
                return ("", "", "", "")
            else:
                return ("", "", "")
            
        p = Path(d).expanduser()
        org_path = str(p / 'PixelClassifier_Organoid.joblib')
        tcell_path = str(p / 'PixelClassifier_TCell.joblib')
        death_path = str(p / 'PixelClassifier_Death.joblib')

        if self.two_org_types.value:
            org2_path = str(p / 'PixelClassifier_Organoid2.joblib')
            return (org_path, tcell_path, death_path, org2_path)
        else:
            return (org_path, tcell_path, death_path)

    def _apply_dir_to_clf_paths(self, d: str):
        paths = self._default_clf_paths_for_dir(d)

        self.clf_org_path.value = paths[0]
        self.clf_tcell_path.value = paths[1]
        self.clf_death_path.value = paths[2]

        if self.two_org_types.value and len(paths) > 3:
            self.clf_org2_path.value = paths[3]

    def _toggle_clf_path_section(self, change=None):
        manual = bool(self.manual_clf_paths.value)
        self.clf_paths_box.layout.display = (None if manual else 'none')
        
        if not manual:
            for p in [self.clf_dir, self.clf_org_path, self.clf_tcell_path, 
                    self.clf_death_path, self.clf_org2_path]:
                p.value = ""

    def _build_clf_paths_box(self):
        children = [
            widgets.HTML("<b>Classifier paths</b>"),
            self.clf_dir,
            self.clf_org_path,
        ]

        if self.two_org_types.value:
            children.append(self.clf_org2_path)

        children.extend([
            self.clf_tcell_path,
            self.clf_death_path,
        ])

        self.clf_paths_box.children = children

    def _on_two_org_types_changed(self, change=None):
        self._build_clf_paths_box()

        if self.manual_clf_paths.value and self.clf_dir.value:
            self._apply_dir_to_clf_paths(self.clf_dir.value)




    def display(self):
        display(self.ui)

    def _lock(self, state: bool):
        # keep close_button enabled so user can close at any time
        for w in [
            self.btn_train, self.btn_apply,
            self.examples_per_sample, self.sample_specific_classifier, self.n_workers, self.two_org_types,
            self.organoid_edt_threshold, self.use_all_timepoints, self.tp_start, self.tp_end, self.manual_clf_paths,
        ]:
            w.disabled = state

        # also lock/unlock the path pickers
        try:
            for p in [self.clf_dir, self.clf_org_path, self.clf_tcell_path, self.clf_death_path, self.clf_org2_path]:
                p.text.disabled = state
                p.button.disabled = state
        except Exception:
            pass

    # ---------- callbacks ----------
    def _on_train_clicked(self, _):
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()
                odir = Path(self.metadata_loader.output_dir).expanduser()
                odir.mkdir(parents=True, exist_ok=True)

                print("▶️ Training pixel classifier…", flush=True)
                print(f"  output_dir={odir}")
                print(f"  examples_per_sample={self.examples_per_sample.value}")
                print(f"  sample_specific_classifier={self.sample_specific_classifier.value}")
                print(f"  two_org_types={self.two_org_types.value}")
                print(f"  n_workers={self.n_workers.value}")

                # UI state
                self._lock(True)
                self.spinner_train.layout.display = None

                # Run training (synchronously). If it opens napari non-blocking and
                # returns a viewer, store it so the Close button can close it.
                try:
                    from behav3d.preprocessing.segmentation.napari_pixelclassifier import train_pixel_classifier
                except Exception:
                    # If your train function is imported elsewhere, rely on that name in globals
                    pass

                ret = train_pixel_classifier(
                    output_dir=str(odir),
                    metadata=self.metadata_loader.metadata,
                    examples_per_sample=int(self.examples_per_sample.value),
                    sample_specific_classifier=bool(self.sample_specific_classifier.value),
                    n_workers=int(self.n_workers.value),
                    two_org_types=bool(self.two_org_types.value),
                )

                # Try to capture a viewer handle
                self._viewer = None
                try:
                    if ret is not None and hasattr(ret, "close"):
                        self._viewer = ret
                    else:
                        # best-effort: some environments expose current viewer getter
                        get_curr = getattr(napari, "current_viewer", None)
                        if callable(get_curr):
                            self._viewer = get_curr()
                except Exception:
                    pass

                # UI finalize for train
                self.spinner_train.layout.display = "none"
                # Show Close button only if we have a viewer to close
                self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"
                if self._viewer is None:
                    print("✅ Training finished.", flush=True)
                else:
                    print("✅ Training UI opened in Napari (use 'Close viewer' to close).", flush=True)

            except Exception:
                import traceback; traceback.print_exc()
                self.spinner_train.layout.display = "none"
                self.close_button.layout.display = "none"
            finally:
                # Re-enable inputs regardless; close button stays visible if viewer exists
                self._lock(False)

    def _on_close_clicked(self, _):
        with self.out:
            try:
                if self._viewer is not None:
                    print("Closing viewer…")
                    try:
                        self._viewer.close()
                    except Exception as e:
                        print(f"Error closing viewer: {e}")
                    self._viewer = None
                else:
                    print("No active viewer to close (viewer may have been opened in blocking mode).")
            finally:
                self.close_button.layout.display = "none"

    def _on_apply_clicked(self, _):
        self._lock(True)
        with self.out:
            self.out.clear_output()
            self._persist_params()
            try:
                odir = Path(self.metadata_loader.output_dir).expanduser()
                odir.mkdir(parents=True, exist_ok=True)

                tpr = _mk_timepoint_range(
                    use_all=bool(self.use_all_timepoints.value),
                    start=int(self.tp_start.value),
                    end=int(self.tp_end.value),
                )

                print("▶️ Applying pixel classifier segmentation…", flush=True)
                print(f"  output_dir={odir}")
                print(f"  organoid_edt_threshold={self.organoid_edt_threshold.value}")
                print(f"  timepoint_range={tpr}", flush=True)
                # Echo chosen classifier paths if manual mode is on
                if self.manual_clf_paths.value:
                    print(f"  clf_dir={self.clf_dir.value}")
                    print(f"  clf_org_path={self.clf_org_path.value}")
                    if self.two_org_types.value:
                        print(f"  clf_org2_path={self.clf_org2_path.value}")
                    print(f"  clf_tcell_path={self.clf_tcell_path.value}")
                    print(f"  clf_death_path={self.clf_death_path.value}")

                self.spinner_apply.layout.display = None

                new_md = run_pixel_classifier_segmentation(
                    output_dir=str(odir),
                    metadata=self.metadata_loader.metadata,
                    organoid_edt_threshold=float(self.organoid_edt_threshold.value),
                    timepoint_range=tpr,
                    two_org_types=bool(self.two_org_types.value),
                    clf_org_path=self.clf_org_path.value,
                    clf_tcell_path=self.clf_tcell_path.value,
                    clf_death_path=self.clf_death_path.value,
                    clf_org2_path=self.clf_org2_path.value if self.two_org_types.value else None,

                )

                try:
                    if new_md is not None:
                        self.metadata_loader.metadata = new_md
                        new_md.to_csv(self.metadata_loader.metadata_csv_path, index=False)
                except Exception:
                    import traceback; traceback.print_exc()

                print("✅ Apply finished.", flush=True)
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self.spinner_apply.layout.display = "none"
                self._lock(False)

class TrackingPanel:
    """
    Minimal UI for tracking (LAP, TrackPy, or Propagation).
    - Uses metadata_loader.output_dir and metadata_loader.metadata_csv_path directly
    - Updates metadata_loader.metadata and always saves CSV
    - cell_type is supplied in __init__
    """
    def __init__(self, metadata_loader, cell_type="tcell"):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()

        # Ensure config skeleton exists, then read profile
        self._ensure_cfg_skeleton()
        params = dict(self.metadata_loader.behav3d_parameters or {})
        tcfg = params["tracking"][self.cell_type]

        # Method selector
        self.method = widgets.Dropdown(
            description="Tracking method",
            options=[("LAP (laptrack)", "lap"),
                     ("TrackPy", "trackpy"),
                     ("Propagation", "propagation")],
            value=str(tcfg.get("method", "lap")),
            style={'description_width': '160px'}
        )

        # Shared controls
        self.overwrite = widgets.Checkbox(
            description="Overwrite existing",
            value=bool(tcfg.get("overwrite", False))
        )

        # LAP params (all as text boxes)
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

        # TrackPy params
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

        # Propagation params (none needed; just an info box)
        self.prop_params = widgets.VBox([
            widgets.HTML("<b>Propagation tracking</b>"),
            widgets.HTML("<i>No tunable parameters.</i>")
        ])

        # Swap param groups on method change
        self.param_container = widgets.VBox()
        self.method.observe(self._toggle_param_groups, names='value')
        self._toggle_param_groups()

        # Run + spinner + log
        self.btn_run = widgets.Button(
            description=f"Run {cell_type} tracking",
            button_style="success",
            layout=widgets.Layout(
                width="fit-content",  # shrink-to-fit text
                flex="0 0 auto"       # don't stretch
            )
        )
        self.spinner_run = widgets.HTML(value=spinning_loader)
        self.spinner_run.layout.display = "none"
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_run],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        self.btn_run.on_click(self._on_run_clicked)
        self.out = widgets.Output()

        # Layout
        self.ui = widgets.VBox([
            widgets.HTML(f"<b>{cell_type} Tracking</b>"),
            self.method,
            widgets.HBox([self.overwrite]),
            self.param_container,
            self.run_row,
            widgets.HTML("<hr>"),
            self.out
        ])

    # ---------------- UI helpers ----------------
    def display(self):
        display(self.ui)

    def _toggle_param_groups(self, change=None):
        if self.method.value == "lap":
            self.param_container.children = (self.lap_params,)
        elif self.method.value == "trackpy":
            self.param_container.children = (self.tp_params,)
        else:  # "propagation"
            self.param_container.children = (self.prop_params,)

    def _lock(self, state: bool):
        for w in [
            self.method, self.overwrite,
            self.track_cost_dist, self.gap_cost_dist, self.gap_max_frames,
            self.merging_cost_dist, self.splitting_cost_dist,
            self.tp_search_range, self.tp_memory, self.tp_adaptive_stop, self.tp_adaptive_step,
            self.btn_run
        ]:
            if hasattr(w, "disabled"):
                w.disabled = state

    def _force_commit_pending_changes(self):
        """Ensure IntText/FloatText edits are committed (they only send on blur/Enter)."""
        try:
            from IPython.display import Javascript, display
            display(Javascript("document.activeElement && document.activeElement.blur();"))
            import time as _t; _t.sleep(0.08)  # tiny pause for the comm to deliver
        except Exception:
            pass

    # ---------------- config I/O ----------------
    def _ensure_cfg_skeleton(self):
        p = dict(self.metadata_loader.behav3d_parameters or {})
        p.setdefault("tracking", {})
        p["tracking"].setdefault(self.cell_type, {})
        prof = p["tracking"][self.cell_type]
        prof.setdefault("method", "lap")
        prof.setdefault("overwrite", False)
        prof.setdefault("lap", {
            "track_cost_px": 45,
            "gap_close_cost_px": 60,
            "gap_close_max_frames": 3,
            "merging_cost_px": 0,
            "splitting_cost_px": 0
        })
        prof.setdefault("trackpy", {
            "search_range_px": 31,
            "memory_frames": 2,
            "adaptive_stop": 10.0,
            "adaptive_step": 0.95
        })
        self.metadata_loader.behav3d_parameters = p  # keep loader in sync

    def _persist_params(self):
        # Make sure the skeleton exists
        self._ensure_cfg_skeleton()
        params = dict(self.metadata_loader.behav3d_parameters or {})
        params.setdefault("tracking", {})
        params["tracking"].setdefault(self.cell_type, {"lap": {}, "trackpy": {}})
        prof = params["tracking"][self.cell_type]

        prof["method"] = str(self.method.value)
        prof["overwrite"] = bool(self.overwrite.value)

        prof.setdefault("lap", {})
        prof["lap"]["track_cost_px"]        = int(self.track_cost_dist.value)
        prof["lap"]["gap_close_cost_px"]    = int(self.gap_cost_dist.value)
        prof["lap"]["gap_close_max_frames"] = int(self.gap_max_frames.value)
        prof["lap"]["merging_cost_px"]      = int(self.merging_cost_dist.value)
        prof["lap"]["splitting_cost_px"]    = int(self.splitting_cost_dist.value)

        prof.setdefault("trackpy", {})
        prof["trackpy"]["search_range_px"]  = int(self.tp_search_range.value)
        prof["trackpy"]["memory_frames"]    = int(self.tp_memory.value)
        prof["trackpy"]["adaptive_stop"]    = float(self.tp_adaptive_stop.value)
        prof["trackpy"]["adaptive_step"]    = float(self.tp_adaptive_step.value)

        # Write back to loader + disk
        self.metadata_loader.behav3d_parameters = params
        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(params, f, sort_keys=False)

    # ---------------- run ----------------
    def _on_run_clicked(self, _):
        self._lock(True)
        self.spinner_run.layout.display = None  # show spinner
        with self.out:
            self.out.clear_output()
            try:
                # Commit any in-progress text edits BEFORE reading values / persisting
                self._force_commit_pending_changes()

                # Persist current UI to YAML
                self._persist_params()

                out_dir = Path(self.metadata_loader.output_dir).expanduser()
                out_dir.mkdir(parents=True, exist_ok=True)
                csv_path = Path(self.metadata_loader.metadata_csv_path).expanduser()

                method = self.method.value
                if method == "lap":
                    # Square distances for laptracking (function expects px^2)
                    tc = int(self.track_cost_dist.value) ** 2
                    gc = int(self.gap_cost_dist.value) ** 2
                    mc = int(self.merging_cost_dist.value)
                    sc = int(self.splitting_cost_dist.value)
                    merging  = (mc ** 2) if mc > 0 else False
                    splitting = (sc ** 2) if sc > 0 else False

                    print("▶️ LAP tracking…", flush=True)
                    print(f"  track_cost_cutoff={tc}, gap_closing_cost_cutoff={gc}")
                    print(f"  gap_closing_max_frame_count={self.gap_max_frames.value}")
                    print(f"  merging_cost_cutoff={merging}, splitting_cost_cutoff={splitting}", flush=True)

                    new_md = run_tcell_laptracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        track_cost_cutoff=tc,
                        gap_closing_cost_cutoff=gc,
                        gap_closing_max_frame_count=int(self.gap_max_frames.value),
                        merging_cost_cutoff=merging,
                        splitting_cost_cutoff=splitting,
                        cell_type=self.cell_type,
                        overwrite=bool(self.overwrite.value),
                    )

                elif method == "trackpy":
                    print("▶️ TrackPy tracking…", flush=True)
                    print(f"  search_range={self.tp_search_range.value}, memory={self.tp_memory.value}")
                    print(f"  adaptive_stop={self.tp_adaptive_stop.value}, adaptive_step={self.tp_adaptive_step.value}", flush=True)

                    new_md = run_tcell_trackpy_tracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        overwrite=bool(self.overwrite.value),
                        cell_type=self.cell_type,
                        search_range=int(self.tp_search_range.value),
                        memory=int(self.tp_memory.value),
                        adaptive_stop=float(self.tp_adaptive_stop.value),
                        adaptive_step=float(self.tp_adaptive_step.value),
                    )

                else:  # "propagation"
                    print("▶️ Propagation tracking…", flush=True)
                    new_md = run_propagation_tracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        cell_type=self.cell_type,
                        overwrite=bool(self.overwrite.value),
                    )

                # Update loader + ALWAYS save CSV
                self.metadata_loader.metadata = new_md
                print(f"💾 Saving metadata to: {csv_path}", flush=True)
                new_md.to_csv(csv_path, sep=",", index=False)
                print("✅ Tracking finished.", flush=True)

            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self.spinner_run.layout.display = "none"  # hide spinner
                self._lock(False)

class TrackingVisualizationPanel:
    def __init__(
        self,
        metadata_loader,
        default_time_range=None,
        channel_colors=None,
    ):
        self.metadata_loader = metadata_loader
        # prefer JSON channel colors if not passed
        self.channel_colors = tuple(channel_colors or _cfg_get(
            self.metadata_loader.behav3d_parameters, "tracking_visualization.channel_colors",
            ["cyan", "yellow", "red", "green", "magenta", "blue"]
        ))

        self._viewer = None
        
        # --- UI: status + refresh
        self._status = widgets.HTML("<b>Waiting for user to load metadata…</b>")
        self._refresh_btn = widgets.Button(
            description="Refresh",
            tooltip="Build/Update selector from metadata_loader.metadata",
        )
        self._refresh_btn.on_click(self._on_refresh_clicked)

        # --- Main controls (disabled until metadata present)
        self.sample_dropdown = widgets.Dropdown(
            options=[],
            value=None,
            description="Sample:",
            layout=widgets.Layout(width="350px"),
            disabled=True,
        )

        # Tickbox to enable/disable time range (from JSON)
        self.use_range = widgets.Checkbox(
            description="Use custom time range",
            value=bool(_cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.use_range", False)),
            indent=False,
        )

        # Start/End boxes (defaults from constructor OR JSON)
        if isinstance(default_time_range, (tuple, list)) and len(default_time_range) == 2:
            _start_default, _end_default = map(int, default_time_range)
        else:
            _start_default = int(_cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.start_t", 0))
            _end_default   = int(_cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.end_t", 0))

        self.start_t = widgets.IntText(
            description="Start T:",
            value=_start_default,
            layout=widgets.Layout(width="180px")
        )
        self.end_t   = widgets.IntText(
            description="End T:",
            value=_end_default,
            layout=widgets.Layout(width="180px")
        )

        self.range_box = widgets.HBox([self.start_t, self.end_t])
        self.range_box.layout.display = "flex" if self.use_range.value else "none"

        def _toggle_range_visibility(change):
            self.range_box.layout.display = "flex" if change["new"] else "none"
        self.use_range.observe(_toggle_range_visibility, names="value")

        self.open_button = widgets.Button(
            description="Visualize segmentation and tracking results",
            button_style="primary",
            tooltip="Launch Napari for the selected sample",
            icon="eye",
            layout=widgets.Layout(width="300px"),
            disabled=True,
        )

        self.close_button = widgets.Button(
            description="Close viewer",
            button_style="danger",
            icon="stop",
            tooltip="Close the active Napari viewer",
            layout=widgets.Layout(width="200px", display="none"),  # hidden by default
        )
        
        self.msg = widgets.Output()

        # Events
        self.open_button.on_click(self._on_open_clicked)
        self.close_button.on_click(self._on_close_clicked)

        # Layout
        self._panel = widgets.VBox(
            [
                widgets.HBox([self._status, self._refresh_btn]),
                self.sample_dropdown,
                self.use_range,
                self.range_box,
                widgets.HBox([self.open_button, self.close_button]),
                self.msg,
            ]
        )

        # Try immediate build if metadata already present
        self._maybe_build_from_loader()

    # ---------------- Public API ----------------
    def display(self):
        display(self._panel)

    def get_selected_row(self) -> pd.Series:
        self._ensure_metadata_ready()
        name = self.sample_dropdown.value
        df = self.metadata_loader.metadata
        row = df[df["sample_name"].astype(str) == str(name)]
        if row.empty:
            raise KeyError(f"sample_name '{name}' not found in metadata.")
        return row.iloc[0]

    def open_viewer(self):
        row = self.get_selected_row()

        # only use range if tickbox checked
        if self.use_range.value:
            start_t, end_t = int(self.start_t.value), int(self.end_t.value)
            if end_t < start_t:
                start_t, end_t = end_t, start_t
            time_range = (start_t, end_t)
        else:
            time_range = None

        with self.msg:
            self.msg.clear_output()
            try:
                print(f"Opening viewer for: {row['sample_name']}")
                self._viewer = visualize_tracks(
                    metadata_row=row,
                    timepoint_range=time_range,
                    channel_colors=self.channel_colors,
                )
            except Exception as e:
                print(f"Error: {e}")
            self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"

    # ---------------- Internal helpers ----------------
    def _ensure_metadata_ready(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if not isinstance(df, pd.DataFrame) or df.empty:
            raise RuntimeError("Metadata not loaded yet. Click 'Refresh' once metadata_loader.metadata is set.")
   
    def _maybe_build_from_loader(self):
        df = getattr(self.metadata_loader, "metadata", None)
        if isinstance(df, pd.DataFrame) and not df.empty and ("sample_name" in df.columns):
            self._build_from_metadata(df)
            self._status.value = "<span style='color:green'>Metadata loaded ✅</span>"
        else:
            self._status.value = "<b>Waiting for user to load metadata…</b>"
            self.sample_dropdown.disabled = True
            self.open_button.disabled = True

    def _build_from_metadata(self, df: pd.DataFrame):
        sample_names = df["sample_name"].astype(str).unique().tolist()
        if not sample_names:
            self._status.value = "<b>No sample_name values found in metadata.</b>"
            self.sample_dropdown.options = []
            self.sample_dropdown.value = None
            self.sample_dropdown.disabled = True
            self.open_button.disabled = True
            return

        desired = _cfg_get(self.metadata_loader.behav3d_parameters, "tracking_visualization.sample_name", None)
        self.sample_dropdown.options = sample_names
        self.sample_dropdown.value = desired if (desired in sample_names) else sample_names[0]
        self.sample_dropdown.disabled = False
        self.open_button.disabled = False

    # ---------------- Callbacks ----------------
    def _on_open_clicked(self, _):
        self.open_viewer()

    def _on_close_clicked(self, _):
        # Only relevant when we’re in non-blocking mode and hold a viewer
        with self.msg:
            try:
                if self._viewer is not None:
                    print("Closing viewer…")
                    self._viewer.close()
                    self._viewer = None
                else:
                    print("No active viewer to close (viewer was likely opened in blocking mode).")
            finally:
                # Hide the button regardless
                self.close_button.layout.display = "none"
                
    def _on_refresh_clicked(self, _):
        with self.msg:
            self.msg.clear_output()
            try:
                self._maybe_build_from_loader()
            except Exception as e:
                print(f"Refresh failed: {e}")
       
class FeatureExtractionPanel:
    """
    Minimal UI for feature extraction (per-cell-type), mirroring TrackingPanel style.
    - Uses metadata_loader.output_dir and metadata_loader.metadata_csv_path directly
    - Updates metadata_loader.metadata and always saves CSV
    - Persists used settings to behav3d_parameters.yml under features.<cell_type>
    """
    def __init__(self, metadata_loader, cell_type="tcell"):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()

        fcfg = _cfg_get(self.metadata_loader.behav3d_parameters,
                        f"features.{self.cell_type}", {}) or {}

        # --- Controls ---
        self.dead_mask_threshold = widgets.FloatText(
            description="Dead mask % thr",
            value=float(fcfg.get("dead_mask_percentage_threshold", 0.05)),
            style={'description_width':'160px'}
        )

        self._all_features = ["movement", "intensity", "morphology", "contact", "death"]
        default_feats = fcfg.get("features_choice", self._all_features)
        if not isinstance(default_feats, (list, tuple)):
            default_feats = self._all_features
        default_feats = [f for f in default_feats if f in self._all_features] or self._all_features

        self.feature_checks = {
            f: widgets.Checkbox(description=f.capitalize(), value=(f in default_feats), indent=False)
            for f in self._all_features
        }

        contact_thr_default = float(fcfg.get("contact_threshold", 0.0))
        self.contact_threshold = widgets.FloatText(
            description="Contact Threshold (µm)",
            value=contact_thr_default,
            layout=widgets.Layout(width="240px"),
            style={'description_width': '180px'}
        )

        # Build "Features" box; put contact_threshold to the right of the Contact checkbox
        contact_row = widgets.HBox([
            self.feature_checks["contact"],
            self.contact_threshold
        ], layout=widgets.Layout(align_items="center", gap="12px"))

        def _toggle_contact_threshold(change=None):
            show = bool(self.feature_checks["contact"].value)
            self.contact_threshold.layout.display = None if show else "none"

        _toggle_contact_threshold()  # set initial visibility
        self.feature_checks["contact"].observe(_toggle_contact_threshold, names="value")

        feat_rows = []
        for f in self._all_features:
            if f == "contact":
                feat_rows.append(contact_row)
            else:
                feat_rows.append(self.feature_checks[f])
                
        self.features_box = widgets.VBox(
            [widgets.HTML("<b>Features</b>")] + feat_rows
        )

        self.n_workers = widgets.IntText(
            description="Workers",
            value=int(fcfg.get("n_workers", max(8, (os.cpu_count() or 8)))),
            max=max(8, (os.cpu_count() or 8)),
            style={'description_width':'160px'}
        )

        # NEW: overwrite checkbox
        self.overwrite = widgets.Checkbox(
            description="Overwrite existing",
            value=bool(fcfg.get("overwrite", False))
        )

        # Run + log
        self.btn_run = widgets.Button(description=f"Run {cell_type} feature extraction", button_style="success", layout=widgets.Layout(
        width="fit-content",   # size to content
        flex="0 0 auto"        # don't stretch in HBox/VBox
    ))
        self.btn_run.on_click(self._on_run_clicked)
        
        self.spinner_html = widgets.HTML(
            value=spinning_loader
        )
        self.spinner_html.layout.display = "none"  # hidden until run starts

        
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html],
            layout=widgets.Layout(align_items="center", gap="10px")
        )
        
        self.out = widgets.Output()

        # Layout
        self.ui = widgets.VBox([
            widgets.HTML(f"<b> {self.cell_type} Feature Extraction</b>"),
            widgets.HBox([self.dead_mask_threshold, self.n_workers]),
            self.features_box,
            self.overwrite,             # <- new row
            self.run_row,
            widgets.HTML("<hr>"),
            self.out
        ])

    def _selected_features(self):
        return [f for f in self._all_features if self.feature_checks[f].value]

    def _persist_params(self):
        params = self.metadata_loader.behav3d_parameters
        params.setdefault("features", {})
        params["features"].setdefault(self.cell_type, {})
        prof = params["features"][self.cell_type]

        prof["dead_mask_percentage_threshold"] = float(self.dead_mask_threshold.value)
        prof["features_choice"] = self._selected_features()
        prof["n_workers"] = int(self.n_workers.value)
        prof["overwrite"] = bool(self.overwrite.value)   # <- save overwrite
        prof["contact_threshold"] = float(self.contact_threshold.value)


        params.setdefault("paths", {})
        if getattr(self.metadata_loader, "metadata_csv_path", None):
            params["paths"]["metadata_csv"] = str(Path(self.metadata_loader.metadata_csv_path).expanduser())
        if getattr(self.metadata_loader, "output_dir", None):
            params["paths"]["output_dir"] = str(Path(self.metadata_loader.output_dir).expanduser())

        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(params, f, sort_keys=False)

    def _lock(self, state: bool):
        for w in [self.dead_mask_threshold, self.n_workers, self.overwrite,
                  self.btn_run, *self.feature_checks.values()]:
            if hasattr(w, "disabled"):
                w.disabled = state

    def _on_run_clicked(self, _):
        self._lock(True)
        self.spinner_html.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                self._persist_params()

                out_dir = Path(self.metadata_loader.output_dir).expanduser()
                out_dir.mkdir(parents=True, exist_ok=True)
                csv_path = Path(self.metadata_loader.metadata_csv_path).expanduser()

                thr = float(self.dead_mask_threshold.value)
                feats = self._selected_features()
                workers = int(self.n_workers.value)
                ow = bool(self.overwrite.value)

                print("▶️ Running feature extraction…", flush=True)
                print(f"  dead_mask_percentage_threshold={thr}")
                print(f"  features_choice={feats}")
                print(f"  n_workers={workers}, overwrite={ow}")
                print(f"  output_dir={out_dir}", flush=True)

                new_md = run_feature_extraction(
                    dead_mask_percentage_threshold=thr,
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(out_dir),
                    features_choice=feats,
                    cell_type=self.cell_type,
                    n_workers=workers,
                    overwrite=ow                # <- pass overwrite here
                )
                print("✅ Feature extraction finished.", flush=True)
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)

class TcellFilterPanel:
    def __init__(self, metadata_loader, cell_type="tcell"):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()

        # seed config
        params = self.metadata_loader.behav3d_parameters
        params.setdefault("track_filtering", {})
        if self.cell_type not in params["track_filtering"]:
            from copy import deepcopy
            params["track_filtering"][self.cell_type] = deepcopy(
                _DEFAULT_CONFIG["track_filtering"][self.cell_type]
            )
            with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                yaml.safe_dump(params, f, sort_keys=False)

        cfg = params["track_filtering"][self.cell_type]

        # ---- exp_duration ----
        self.en_exp_duration = widgets.Checkbox(
            description="Trim down full time series to supplied duration",
            value=bool(cfg.get("exp_duration_enabled", False)),
            indent=False,
            style={'description_width': '240px'}
        )
        self.exp_duration = widgets.IntText(
            description="Max timepoints",
            value=int(cfg.get("exp_duration", 999999)),
            style={'description_width': '160px'},
            continuous_update=False
        )

        self.row_exp = widgets.HBox([self.exp_duration], layout=widgets.Layout(display="none"))
        self.en_exp_duration.observe(
            lambda c: self.row_exp.layout.__setattr__("display", None if self.en_exp_duration.value else "none"),
            names="value"
        )
        self.en_exp_duration.observe(lambda c: _update_unit_toggle_visibility(), names="value")

        if self.en_exp_duration.value:
            self.row_exp.layout.display = None

        # ---- min_track_length ----
        self.en_min_len = widgets.Checkbox(
            description="Select only tracks with minimal length",
            value=bool(cfg.get("min_track_length_enabled", False)),
            indent=False
        )
        # REPLACED IntSlider -> IntText
        self.min_track_length = widgets.IntText(
            description="Minimal length (frames)",
            value=int(cfg.get("min_track_length", 0)),
            style={'description_width': '160px'},
            continuous_update=False
        )
        self.row_min = widgets.HBox([self.min_track_length], layout=widgets.Layout(display="none"))
        self.en_min_len.observe(
            lambda c: self.row_min.layout.__setattr__("display", None if self.en_min_len.value else "none"),
            names="value"
        )
        self.en_min_len.observe(lambda c: _update_unit_toggle_visibility(), names="value")

        if self.en_min_len.value:
            self.row_min.layout.display = None

        # ---- max_track_length ----
        self.en_max_len = widgets.Checkbox(
            description="Trim down tracks to supplied length",
            value=bool(cfg.get("max_track_length_enabled", False)),
            indent=False
        )
        # REPLACED IntSlider -> IntText
        self.max_track_length = widgets.IntText(
            description="Maximal length (frames)",
            value=int(cfg.get("max_track_length", 999999)),
            style={'description_width': '160px'},
            continuous_update=False
        )
        self.row_max = widgets.HBox([self.max_track_length], layout=widgets.Layout(display="none"))
        self.en_max_len.observe(
            lambda c: self.row_max.layout.__setattr__("display", None if self.en_max_len.value else "none"),
            names="value"
        )
        self.en_max_len.observe(lambda c: _update_unit_toggle_visibility(), names="value")

        if self.en_max_len.value:
            self.row_max.layout.display = None

        # ---- filter_t0_dead ----
        self.filter_t0_dead = widgets.Checkbox(
            description="Filter tracks that are dead at t=0",
            value=bool(cfg.get("filter_t0_dead", False)),
            indent=False
        )
        
        # Run button + output
        self.btn_run = widgets.Button(description=f"Filter {cell_type} tracks & summarize", button_style="success", layout=widgets.Layout(
        width="fit-content",   # size to content
        flex="0 0 auto"        # don't stretch in HBox/VBox
    ))
        self.btn_run.on_click(self._on_run_clicked)
        
        self.spinner_html = widgets.HTML(
            value=spinning_loader
        )
        self.spinner_html.layout.display = "none"  # hidden until run starts

                # ----- time unit toggle (Frames / Hours) -----
        _TIME_TYPE_OPTIONS = [("frames", "Frames"), ("hours", "real_time")]
        default_time_type = str(cfg.get("time_type", "Frames"))
        default_label = next((lbl for (lbl, val) in _TIME_TYPE_OPTIONS if val == default_time_type), "frames")

        self.time_type_toggle = widgets.ToggleButtons(
            options=_TIME_TYPE_OPTIONS,  # [('frames','Frames'), ('hours','real_time')]
            value=next(val for (lbl, val) in _TIME_TYPE_OPTIONS if lbl == default_label),
            description="Unit to use for time-based filters:",
            button_style="info",
            layout=widgets.Layout(width="auto")
        )

        # Update field labels to match units
        def _refresh_labels(time_type_value: str):
            if time_type_value == "real_time":
                self.exp_duration.description = "Max duration (hours)"
                self.min_track_length.description = "Minimal length (hours)"
                self.max_track_length.description = "Maximal length (hours)"
            else:
                self.exp_duration.description = "Max timepoints (frames)"
                self.min_track_length.description = "Minimal length (frames)"
                self.max_track_length.description = "Maximal length (frames)"

        _refresh_labels(self.time_type_toggle.value)
        self.time_type_toggle.observe(lambda c: _refresh_labels(self.time_type_toggle.value), names="value")

        # Wrap the toggle in a row that we can show/hide
        self.row_unit = widgets.HBox([self.time_type_toggle], layout=widgets.Layout(display="none"))

        # Show the row only if any filter is enabled
        def _update_unit_toggle_visibility():
            visible = self.en_exp_duration.value or self.en_min_len.value or self.en_max_len.value
            self.row_unit.layout.display = None if visible else "none"

        # Call once now, and whenever any enabling checkbox changes
        _update_unit_toggle_visibility()
        
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html],
            layout=widgets.Layout(align_items="center", gap="10px")
        )
        
        self.out = widgets.Output()

        # Full layout
        self.ui = widgets.VBox([
            widgets.HTML(f"<b>{self.cell_type} Track Filtering</b>"),
            self.en_exp_duration, self.row_exp,
            self.en_min_len, self.row_min,
            self.en_max_len, self.row_max,
            self.filter_t0_dead,
            self.row_unit,
            self.run_row,
            widgets.HTML("<hr>"),
            self.out
        ])

        # results
        self.df_tracks_filt = None
        self.df_tracks_summ = None

    def display(self):
        display(self.ui)

    def _effective_values(self):
        exp_duration = float(self.exp_duration.value) if self.en_exp_duration.value else 999999.0
        min_len      = int(self.min_track_length.value) if self.en_min_len.value else 0
        max_len      = int(self.max_track_length.value) if self.en_max_len.value else 999999
        filt_dead    = bool(self.filter_t0_dead.value)
        time_type    = str(self.time_type_toggle.value)  # "frames" or "real_time"
        return exp_duration, min_len, max_len, filt_dead, time_type

    def _persist_params(self):
        params = self.metadata_loader.behav3d_parameters
        prof = params["track_filtering"][self.cell_type]

        # Save enable flags
        prof["exp_duration_enabled"]      = bool(self.en_exp_duration.value)
        prof["min_track_length_enabled"]  = bool(self.en_min_len.value)
        prof["max_track_length_enabled"]  = bool(self.en_max_len.value)

        # Always save the text values exactly as the user set them
        prof["exp_duration"]      = float(self.exp_duration.value)
        prof["min_track_length"]  = int(self.min_track_length.value)
        prof["max_track_length"]  = int(self.max_track_length.value)

        # Persist checkbox directly
        prof["filter_t0_dead"] = bool(self.filter_t0_dead.value)
        prof["time_type"] = str(self.time_type_toggle.value)

        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(params, f, sort_keys=False)

    def _lock(self, state: bool):
        for w in [
            self.en_exp_duration, self.exp_duration,
            self.en_min_len, self.min_track_length,
            self.en_max_len, self.max_track_length,
            self.filter_t0_dead, self.btn_run
        ]:
            if hasattr(w, "disabled"):
                w.disabled = state

    def _on_run_clicked(self, _):
        self._lock(True)
        self.spinner_html.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                out_dir = Path(self.metadata_loader.output_dir).expanduser()
                out_dir.mkdir(parents=True, exist_ok=True)

                exp_duration, min_len, max_len, filt_t0, time_type = self._effective_values()
                print("▶️ Filtering tracks…")
                print(f"  exp_duration={exp_duration}, min_len={min_len}, max_len={max_len}, "
                      f"filter_t0_dead={filt_t0}, unit={'frames' if time_type=='Frames' else 'hours'}")
                self._persist_params()
                
                self.df_tracks_filt = filter_tcell_tracks(
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(out_dir),
                    exp_duration=exp_duration,
                    min_track_length=min_len,
                    max_track_length=max_len,
                    filter_t0_dead=filt_t0,
                    cell_type=self.cell_type,
                    time_type=time_type,
                )
                print("✅ Filtering done. Summarizing…")

                self.df_tracks_summ = summarize_track_features(
                    output_dir=str(out_dir),
                    cell_type=self.cell_type
                )

                from IPython.display import display as _disp
                print("— Filtered tracks (head) —")
                try: _disp(self.df_tracks_filt.head(10))
                except Exception: print(f"[df_tracks_filt] shape={getattr(self.df_tracks_filt,'shape',None)}")
                print("\n— Summary (head) —")
                try: _disp(self.df_tracks_summ.head(10))
                except Exception: print(f"[df_tracks_summ] shape={getattr(self.df_tracks_summ,'shape',None)}")
                print("✅ Finished.")

            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)

class TCellAnalysisPanel:
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        # ---- config bootstrap (load + ensure defaults) ----
        try:
            groups = behav3d_calculated_features
            md = self.metadata_loader.metadata
            # Detect if organoid_2 features should be added
            self.two_org_types = False
            if md is not None and 'organoid_2_tracks_csv_path' in md.columns:
                s = md['organoid_2_tracks_csv_path'].dropna()
                if not s.empty and Path(str(s.iloc[0])).exists():
                    self.two_org_types = True 
                    # Add organoid_2 contact features
                    groups["contact"].extend([
                        "organoid_2_contact",
                        "organoid_2_contact_pixels",
            ])
        except NameError:
            raise RuntimeError("Define behav3d_calculated_features before creating TCellAnalysisPanel.")

        params = dict(self.metadata_loader.behav3d_parameters or {})
        params.setdefault("analysis", {})
        params["analysis"].setdefault("tcell", self._default_panel_config(groups))
        self._params = params
        self._panel_cfg = self._params["analysis"]["tcell"]

        # ---- headings ----
        self.section_title = widgets.HTML('<div style="font-size:22px;font-weight:700;">T cell analysis</div>')
        self.sel_title     = widgets.HTML('<div style="font-size:20px;font-weight:700;">Select features to use for Dynamic Time Warping (DTW):</div>')
        self.norm_title    = widgets.HTML('<div style="font-size:20px;font-weight:700;">Normalize</div>')
        self.umap_title    = widgets.HTML('<div style="font-size:20px;font-weight:700;">UMAP settings</div>')
        # ---- Seed ----
        self.seed_widget = widgets.IntText(
            description="Seed",
            value=int(self._panel_cfg.get("seed", self._params.get("seed", 42))),
            style={"description_width": "80px"},
        )

        # ---------- layout helpers (shared by both sections) ----------
        def _grid_for(children, indent_px="24px", columns=3):
            # Force 3 columns for every group
            return widgets.GridBox(
                children,
                layout=widgets.Layout(
                    grid_template_columns=" ".join(["max-content"] * columns),
                    grid_gap="6px 18px",
                    margin=f"0 0 0 {indent_px}",
                )
            )

        # ---- SELECTION (grouped; always visible) ----
        self._group_rows = {}  # group -> {"group_cb","child_cbs","grid","container","group_handler","child_handler"}
        sel_group_boxes = []

        preset = list(self._panel_cfg.get("dtw_features_input", []))
        preset_set = set(preset)
        groups_enabled = dict(self._panel_cfg.get("dtw_feature_groups_enabled", {}))

        for group_name, feats in groups.items():
            # Children init from preset; group init from saved flag or "all children on"
            child_cbs = [
                widgets.Checkbox(value=(f in preset_set), description=f, indent=True)
                for f in feats
            ]
            g_init = bool(groups_enabled.get(group_name, all(cb.value for cb in child_cbs)))

            gcb = widgets.Checkbox(value=g_init, indent=False)
            glabel = widgets.HTML(f"<b>{group_name}</b>")
            header = widgets.HBox([gcb, glabel])

            grid = _grid_for(child_cbs, columns=3)   # ALWAYS 3 columns
            container = widgets.VBox([header, grid])

            self._group_rows[group_name] = {
                "group_cb": gcb, "child_cbs": child_cbs, "grid": grid, "container": container
            }
            sel_group_boxes.append(container)

        self.checkboxes_box = widgets.VBox(sel_group_boxes)
        # Wire selection events (batched)
        def make_group_handler(grp_name):
            def _on_group(change):
                if change["name"] != "value":
                    return
                val = bool(change["new"])
                row = self._group_rows[grp_name]
                # batch set children values only (no hiding/disable)
                for cb in row["child_cbs"]:
                    cb.unobserve(row["child_handler"], names="value")
                for cb in row["child_cbs"]:
                    cb.value = val
                for cb in row["child_cbs"]:
                    cb.observe(row["child_handler"], names="value")
                # sync normalize for this group (once)
                self._sync_normalize_visibility_batched([grp_name])
            return _on_group

        def make_child_handler(grp_name):
            def _on_child(change):
                if change["name"] != "value":
                    return
                row = self._group_rows[grp_name]
                all_on = all(cb.value for cb in row["child_cbs"])
                row["group_cb"].unobserve(row["group_handler"], names="value")
                row["group_cb"].value = all_on
                row["group_cb"].observe(row["group_handler"], names="value")
                self._sync_normalize_visibility_batched([grp_name])
            return _on_child

        for grp_name, row in self._group_rows.items():
            row["group_handler"] = make_group_handler(grp_name)
            row["child_handler"] = make_child_handler(grp_name)
            row["group_cb"].observe(row["group_handler"], names="value")
            for cb in row["child_cbs"]:
                cb.observe(row["child_handler"], names="value")

        # ---- top-level selection actions ----
        self.btn_select_all = widgets.Button(description="Select all")
        self.btn_clear_all  = widgets.Button(description="Clear")
        self.btn_select_all.on_click(lambda *_: self._set_all(True))
        self.btn_clear_all.on_click(lambda *_: self._set_all(False))

        self.output_area_features = widgets.Output()

        # ---- NORMALIZE (3 columns like selection; only selected visible) ----
        self._norm_group_rows = {}  # group -> {"group_cb","child_cbs","grid","container","group_handler"}
        norm_group_boxes = []
        for group_name, feats in groups.items():
            gcb = widgets.Checkbox(value=False, indent=False)  # master for visible children
            glabel = widgets.HTML(f"<b>{group_name}</b>")
            header = widgets.HBox([gcb, glabel])

            # Same child order & widgets; visibility controls position
            child_cbs = []
            for f in feats:
                cb = widgets.Checkbox(
                    value=bool(self._panel_cfg.get("z_normalize", {}).get(f, False)),
                    description=f,
                    indent=True
                )
                cb.layout.visibility = 'hidden'  # start hidden; will show if selected
                child_cbs.append(cb)

            grid = _grid_for(child_cbs, columns=3)  # force 3 columns
            container = widgets.VBox([header, grid])

            self._norm_group_rows[group_name] = {
                "group_cb": gcb, "child_cbs": child_cbs, "grid": grid, "container": container
            }
            norm_group_boxes.append(container)

        # normalize group master toggles only currently visible children
        def make_norm_group_handler(grp_name):
            def _on_group(change):
                if change["name"] != "value":
                    return
                val = bool(change["new"])
                row = self._norm_group_rows[grp_name]
                for cb in row["child_cbs"]:
                    if cb.layout.visibility != 'hidden':
                        cb.value = val
            return _on_group
        for grp_name, row in self._norm_group_rows.items():
            row["group_handler"] = make_norm_group_handler(grp_name)
            row["group_cb"].observe(row["group_handler"], names="value")

        # Global normalize controls (like selection)
        self.norm_select_all_btn = widgets.Button(description="Select all")
        self.norm_clear_all_btn  = widgets.Button(description="Clear")
        self.norm_select_all_btn.on_click(self._on_norm_select_all)
        self.norm_clear_all_btn.on_click(self._on_norm_clear_all)

        self.normalize_section = widgets.VBox([
            self.norm_title,
            widgets.HBox([self.norm_select_all_btn, self.norm_clear_all_btn],
                         layout=widgets.Layout(gap="8px")),
            widgets.VBox(norm_group_boxes)
        ])

        
        # ---- UMAP / clustering ----
        self.umap_distance_widget  = widgets.FloatText(description="UMAP min_dist",
                                                       value=float(self._panel_cfg.get("umap_min_dist", 0.1)),
                                                       style={"description_width": "140px"})
        self.umap_neighbors_widget = widgets.IntText(description="UMAP n_neighbors",
                                                     value=int(self._panel_cfg.get("umap_n_neighbors", 15)),
                                                     style={"description_width": "140px"})
        self.num_clusters_widget   = widgets.IntText(description="# clusters",
                                                     value=int(self._panel_cfg.get("nr_of_clusters", 5)),
                                                     style={"description_width": "140px"})
        
        self.umap_box = widgets.VBox([
            self.umap_title,
            widgets.HBox([self.umap_distance_widget, self.umap_neighbors_widget, self.num_clusters_widget],
                         layout=widgets.Layout(flex_flow="row wrap", gap="12px"))
        ])
        
        self.output_area_params = widgets.Output()

        # ---- Run ----
        self.btn_run = widgets.Button(
            description=f"Run T cell analysis",
            button_style="success",
            layout=widgets.Layout(width="300px")
        )
        self.btn_run.on_click(self._on_run_clicked)
        
        
        # Tiny inline SVG spinner + label (hidden by default)
        self.spinner_html = widgets.HTML(
            value=spinning_loader
        )
        self.spinner_html.layout.display = "none"  # hidden until run starts

        
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html],
            layout=widgets.Layout(align_items="center", gap="10px")
        )
        
        self.out_run = widgets.Output()

        # ---- UI ----
        self.ui = widgets.VBox([
            self.section_title,
            self.seed_widget,
            self.sel_title,
            widgets.HBox([self.btn_select_all, self.btn_clear_all], layout=widgets.Layout(gap="8px")),
            self.checkboxes_box,
            self.output_area_features,
            self.normalize_section,
            widgets.HTML("<hr>"),
            self.umap_box,
            self.output_area_params,
            widgets.HTML("<hr>"),
            self.run_row,
            self.out_run,
        ])

        # internal state
        self.dtw_features = self._selected_features()
        self.umap_minimal_distance = float(self.umap_distance_widget.value)
        self.umap_n_neighbors = int(self.umap_neighbors_widget.value)
        self.nr_of_clusters = int(self.num_clusters_widget.value)
        self.df_tracks_clustered = None

        # initial normalize visibility sync (keep positions; only selected visible)
        self._sync_normalize_visibility_batched(list(self._group_rows.keys()))

    # ---------- defaults & persistence ----------
    def _selected_normalize_patterns(self):
        """
        Return the list of pattern names whose normalize checkbox is ON (visible items only).
        """
        pats = []
        for grp_name, row in self._norm_group_rows.items():
            for cb in row["child_cbs"]:
                # visible in Normalize panel means it's part of the DTW selection
                if cb.layout.visibility != 'hidden' and cb.value:
                    pats.append(cb.description)
        # de-dup preserve order
        seen, out = set(), []
        for p in pats:
            if p not in seen:
                out.append(p); seen.add(p)
        return out


    def _default_panel_config(self, groups_dict):
        return {
            "seed": 42,
            "umap_min_dist": 0.1,
            "umap_n_neighbors": 15,
            "nr_of_clusters": 5,
            "dtw_feature_groups_enabled": {g: True for g in groups_dict.keys()},
            "dtw_features_input": [],
            "dtw_features_resolved": [],
            "z_normalize": {},
        }

    def _persist_params(self, *, resolved_use=None, resolved_norm=None):
        # snapshot current UI state to parameter dict
        self._panel_cfg["seed"] = int(self.seed_widget.value)
        self._panel_cfg["umap_min_dist"] = float(self.umap_distance_widget.value)
        self._panel_cfg["umap_n_neighbors"] = int(self.umap_neighbors_widget.value)
        self._panel_cfg["nr_of_clusters"] = int(self.num_clusters_widget.value)

        # DTW selection (patterns as shown in the UI)
        self._panel_cfg["dtw_features_input"] = list(self._selected_features())
        self._panel_cfg["dtw_feature_groups_enabled"] = {
            g: row["group_cb"].value for g, row in self._group_rows.items()
        }

        # Normalize selections (pattern-level switches as shown in the UI)
        self._panel_cfg["z_normalize"] = self._collect_znorm_map()  # {feature(pattern): bool}
        self._panel_cfg["columns_to_normalize_input"] = self._selected_normalize_patterns()

        # Optionally store resolved (expanded) columns
        if resolved_use is not None:
            self._panel_cfg["dtw_features_resolved"] = list(resolved_use)
            self._panel_cfg["columns_to_use_resolved"] = list(resolved_use)
        if resolved_norm is not None:
            self._panel_cfg["columns_to_normalize_resolved"] = list(resolved_norm)

        # write back to loader + disk
        self._params["analysis"]["tcell"] = self._panel_cfg
        self.metadata_loader.behav3d_parameters = self._params
        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(self._params, f, sort_keys=False)

    # ---------- helpers ----------
    def display(self):
        display(self.ui)

    def _set_all(self, val: bool):
        # Bulk toggle: set all group masters and all children
        for grp_name, row in self._group_rows.items():
            row["group_cb"].unobserve(row["group_handler"], names="value")
            for cb in row["child_cbs"]:
                cb.unobserve(row["child_handler"], names="value")
            row["group_cb"].value = val
            for cb in row["child_cbs"]:
                cb.value = val
            for cb in row["child_cbs"]:
                cb.observe(row["child_handler"], names="value")
            row["group_cb"].observe(row["group_handler"], names="value")
        # One batched normalize refresh
        self._sync_normalize_visibility_batched(list(self._group_rows.keys()))

    def _all_child_checkboxes(self):
        out = []
        for row in self._group_rows.values():
            out.extend(row["child_cbs"])
        return out

    def _selected_features(self):
        seen, selected = set(), []
        for cb in self._all_child_checkboxes():
            if cb.value and cb.description not in seen:
                selected.append(cb.description); seen.add(cb.description)
        return selected

    def _infer_available_columns(self):
        for attr in ("available_columns", "feature_columns"):
            cols = getattr(self.metadata_loader, attr, None)
            if cols:
                return list(cols)
        meta = getattr(self.metadata_loader, "metadata", None)
        if meta is not None and hasattr(meta, "columns"):
            try:
                return list(meta.columns)
            except Exception:
                pass
        return None

    def _expand_patterns(self, selected):
        cols = self._infer_available_columns()
        if not cols:
            return selected
        out = []
        for name in selected:
            if any(ch in name for ch in "*?["):
                matches = [c for c in cols if fnmatch.fnmatchcase(c, name)]
                out.extend(matches if matches else [name])
            else:
                out.append(name)
        seen, uniq = set(), []
        for f in out:
            if f not in seen:
                uniq.append(f); seen.add(f)
        return uniq

    # ---- normalize layout sync ----
    def _visible_norm_children(self, grp_name):
        row = self._norm_group_rows[grp_name]
        return [cb for cb in row["child_cbs"] if cb.layout.visibility != 'hidden']

    def _sync_normalize_visibility_batched(self, groups_to_update):
        """
        Mirror DTW selection to normalize:
        - keep grid structure identical using visibility
        - show only selected features in their original slots
        - hide a normalize group container if none of its features are selected
        """
        selected = set(self._selected_features())

        for grp_name in groups_to_update:
            sel_row = self._group_rows[grp_name]
            nrow = self._norm_group_rows[grp_name]

            any_visible = False
            for cb_sel in sel_row["child_cbs"]:
                feat = cb_sel.description
                cb_norm = next(c for c in nrow["child_cbs"] if c.description == feat)
                is_selected = cb_sel.value
                cb_norm.layout.visibility = None if is_selected else 'hidden'
                any_visible = any_visible or is_selected

            # normalize group master reflects only visible children
            vis_children = self._visible_norm_children(grp_name)
            gh = nrow.get("group_handler")
            nrow["group_cb"].unobserve(gh, names="value")
            nrow["group_cb"].value = bool(vis_children) and all(cb.value for cb in vis_children)
            nrow["group_cb"].observe(gh, names="value")

            # show container only if something visible
            nrow["container"].layout.display = None if any_visible else "none"

        # hide whole section if nothing visible
        any_group_visible = any(self._norm_group_rows[g]["container"].layout.display != "none"
                                for g in self._norm_group_rows)
        self.normalize_section.layout.display = None if any_group_visible else "none"

    def _collect_znorm_map(self):
        # Persist only visible (selected) normalize choices
        out = {}
        for grp_name, row in self._norm_group_rows.items():
            for cb in row["child_cbs"]:
                if cb.layout.visibility != 'hidden':
                    out[cb.description] = bool(cb.value)
        return out

    # ---- normalize global actions ----
    def _on_norm_select_all(self, *_):
        # Set all visible normalize checkboxes to True, then update each group's master
        for grp_name, row in self._norm_group_rows.items():
            for cb in row["child_cbs"]:
                if cb.layout.visibility != 'hidden':
                    cb.value = True
            gh = row.get("group_handler")
            vis_children = self._visible_norm_children(grp_name)
            row["group_cb"].unobserve(gh, names="value")
            row["group_cb"].value = bool(vis_children) and all(cb.value for cb in vis_children)
            row["group_cb"].observe(gh, names="value")

    def _on_norm_clear_all(self, *_):
        # Set all visible normalize checkboxes to False, then update each group's master
        for grp_name, row in self._norm_group_rows.items():
            for cb in row["child_cbs"]:
                if cb.layout.visibility != 'hidden':
                    cb.value = False
            gh = row.get("group_handler")
            vis_children = self._visible_norm_children(grp_name)
            row["group_cb"].unobserve(gh, names="value")
            row["group_cb"].value = False if vis_children else False
            row["group_cb"].observe(gh, names="value")

    def _lock(self, locked):
        # Selection widgets
        for row in self._group_rows.values():
            row["group_cb"].disabled = locked
            for cb in row["child_cbs"]:
                cb.disabled = locked
        # Normalize widgets
        for row in self._norm_group_rows.values():
            row["group_cb"].disabled = locked
            for cb in row["child_cbs"]:
                cb.disabled = locked
        for w in [
            self.seed_widget, self.btn_select_all, self.btn_clear_all,
            self.umap_distance_widget, self.umap_neighbors_widget, self.num_clusters_widget,
            self.norm_select_all_btn, self.norm_clear_all_btn, self.btn_run
        ]:
            w.disabled = locked

    # ---------- run ----------
    def _on_run_clicked(self, *_):
        self._lock(True)
        self.out_run.clear_output()
        self.spinner_html.layout.display = None
        with self.out_run:
            try:
                print("▶️ Running T-cell behavioral analysis…")

                # Read current UI state
                seed = int(self.seed_widget.value)
                self.umap_minimal_distance = float(self.umap_distance_widget.value)
                self.umap_n_neighbors      = int(self.umap_neighbors_widget.value)
                self.nr_of_clusters        = int(self.num_clusters_widget.value)

                # Selection patterns (always visible)
                dtw_patterns = self._selected_features()
                if not dtw_patterns:
                    print("⚠️ Please select at least one DTW feature.")
                    return

                # Normalize patterns (only those you ticked in Normalize)
                norm_patterns = self._selected_normalize_patterns()

                # Expand any globs (e.g. mean_intensity_*)
                columns_to_use       = self._expand_patterns(dtw_patterns)
                columns_to_normalize = self._expand_patterns(norm_patterns)

                # Persist both pattern-level and resolved lists
                self._persist_params(resolved_use=columns_to_use, resolved_norm=columns_to_normalize)

                # Sanity printout
                print(f"  output_dir                = {self.output_dir}")
                print(f"  seed                     = {seed}")
                print(f"  UMAP: n_neighbors={self.umap_n_neighbors}, min_dist={self.umap_minimal_distance}")
                print(f"  clusters                 = {self.nr_of_clusters}")
                print(f"  columns_to_use (resolved)       [{len(columns_to_use)}]: {columns_to_use}")
                print(f"  columns_to_normalize (resolved) [{len(columns_to_normalize)}]: {columns_to_normalize}")

                # Ensure output dir
                Path(self.output_dir).mkdir(parents=True, exist_ok=True)
                random.seed(seed)

                # Call with the new signature
                self.df_tracks_clustered = run_tcell_analysis(
                    output_dir=self.output_dir,
                    umap_minimal_distance=self.umap_minimal_distance,
                    umap_n_neighbors=self.umap_n_neighbors,
                    nr_of_clusters=self.nr_of_clusters,
                    columns_to_use=columns_to_use,
                    columns_to_normalize=columns_to_normalize,
                    seed=seed
                )

                # Show preview if it quacks like a DataFrame
                try:
                    from IPython.display import display as _disp
                    # print("✅ Done. Preview (head):")
                    _disp(self.df_tracks_clustered.head(10))
                except Exception:
                    pass
                # print("✅ Done. Result object:", type(self.df_tracks_clustered))

            except Exception:
                import traceback
                print("❌ Error while running analysis:")
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)

class OrganoidFilterPanel:
    """
    Styled like your TrackFilterPanel example, but wired to filter_organoid_tracks(...).
    """
    def __init__(self, metadata_loader, cell_type="organoid"):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip() or "organoid"

        # ---- bootstrap config ----
        params = self.metadata_loader.behav3d_parameters or {}
        params.setdefault("track_filtering", {})
        if self.cell_type not in params["track_filtering"]:
            # try to seed from global _DEFAULT_CONFIG if present
            try:
                base = _DEFAULT_CONFIG["track_filtering"][self.cell_type]  # type: ignore[name-defined]
            except Exception:
                base = _DEFAULT_TRACK_FILTERING[self.cell_type]
            params["track_filtering"][self.cell_type] = dict(base)
            # persist immediately
            with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                yaml.safe_dump(params, f, sort_keys=False)
        self._params = params
        cfg = params["track_filtering"][self.cell_type]

        # ---- exp_duration ----
        self.en_exp_duration = widgets.Checkbox(
            description="Trim down full time series to supplied duration",
            value=bool(cfg.get("exp_duration_enabled", False)),
            indent=False,
            style={'description_width': '240px'}
        )
        self.exp_duration = widgets.IntText(
            description="Max timepoints",
            value=int(cfg.get("exp_duration", 999999)),
            style={'description_width': '160px'},
            continuous_update=False
        )
        self.row_exp = widgets.HBox([self.exp_duration], layout=widgets.Layout(display="none"))
        self.en_exp_duration.observe(
            lambda c: self.row_exp.layout.__setattr__("display", None if self.en_exp_duration.value else "none"),
            names="value"
        )
        if self.en_exp_duration.value:
            self.row_exp.layout.display = None

        # ---- min_track_length ----
        self.en_min_len = widgets.Checkbox(
            description="Select only tracks with minimal length",
            value=bool(cfg.get("min_track_length_enabled", False)),
            indent=False
        )
        self.min_track_length = widgets.IntText(
            description="Minimal length (frames)",
            value=int(cfg.get("min_track_length", 0)),
            style={'description_width': '160px'},
            continuous_update=False
        )
        self.row_min = widgets.HBox([self.min_track_length], layout=widgets.Layout(display="none"))
        self.en_min_len.observe(
            lambda c: self.row_min.layout.__setattr__("display", None if self.en_min_len.value else "none"),
            names="value"
        )
        if self.en_min_len.value:
            self.row_min.layout.display = None

        # ---- max_track_length ----
        self.en_max_len = widgets.Checkbox(
            description="Trim down tracks to supplied length",
            value=bool(cfg.get("max_track_length_enabled", False)),
            indent=False
        )
        self.max_track_length = widgets.IntText(
            description="Maximal length (frames)",
            value=int(cfg.get("max_track_length", 999999)),
            style={'description_width': '160px'},
            continuous_update=False
        )
        self.row_max = widgets.HBox([self.max_track_length], layout=widgets.Layout(display="none"))
        self.en_max_len.observe(
            lambda c: self.row_max.layout.__setattr__("display", None if self.en_max_len.value else "none"),
            names="value"
        )
        if self.en_max_len.value:
            self.row_max.layout.display = None

        # ---- Frames / Hours toggle (shown only if any filter is enabled) ----
        _TIME_TYPE_OPTIONS = [("frames", "Frames"), ("hours", "real_time")]
        default_time_type = str(cfg.get("time_type", "Frames"))
        self.time_type_toggle = widgets.ToggleButtons(
            options=_TIME_TYPE_OPTIONS,
            value=default_time_type if default_time_type in ("Frames", "real_time") else "Frames",
            description="Unit",
            button_style="info",
            layout=widgets.Layout(width="auto")
        )
        def _refresh_labels(val: str):
            if val == "real_time":
                self.exp_duration.description = "Max duration (hours)"
                self.min_track_length.description = "Minimal length (hours)"
                self.max_track_length.description = "Maximal length (hours)"
            else:
                self.exp_duration.description = "Max timepoints (frames)"
                self.min_track_length.description = "Minimal length (frames)"
                self.max_track_length.description = "Maximal length (frames)"
        _refresh_labels(self.time_type_toggle.value)
        self.time_type_toggle.observe(lambda c: _refresh_labels(self.time_type_toggle.value), names="value")

        self.row_unit = widgets.HBox([self.time_type_toggle], layout=widgets.Layout(display="none"))
        def _update_unit_toggle_visibility(_=None):
            show = self.en_exp_duration.value or self.en_min_len.value or self.en_max_len.value
            self.row_unit.layout.display = None if show else "none"
        for cb in (self.en_exp_duration, self.en_min_len, self.en_max_len):
            cb.observe(_update_unit_toggle_visibility, names="value")
        _update_unit_toggle_visibility()
        
        self.en_min_size = widgets.Checkbox(
            description="Filter by minimal size at t=1",
            value=bool(cfg.get("min_size_enabled", True)),
            indent=False
        )
        self.min_size = widgets.IntText(
            description="Minimal size (px @ t=1)",
            value=int(cfg.get("min_size", 1000)),
            style={'description_width': '160px'},
            continuous_update=False
        )
        self.row_min_size = widgets.HBox([self.min_size], layout=widgets.Layout(display="none"))
        self.en_min_size.observe(
            lambda c: self.row_min_size.layout.__setattr__("display", None if self.en_min_size.value else "none"),
            names="value"
        )
        if self.en_min_size.value:
            self.row_min_size.layout.display = None

        # ---- run row ----
        self.btn_run = widgets.Button(description=f"Filter {cell_type} tracks & summarize", button_style="success", layout=widgets.Layout(
            width="fit-content", flex="0 0 auto"
    ))
        self.btn_run.on_click(self._on_run_clicked)
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"

        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html],
            layout=widgets.Layout(align_items="center", gap="10px")
        )

        self.out = widgets.Output()

        # ---- full UI ----
        self.ui = widgets.VBox([
            widgets.HTML(f"<b>{self.cell_type.title()} Track Filtering</b>"),
            self.en_exp_duration, self.row_exp,
            self.en_min_len, self.row_min,
            self.en_max_len, self.row_max,
            self.en_min_size, self.row_min_size,
            self.row_unit,
            self.run_row,
            widgets.HTML("<hr>"),
            self.out
        ])

        # results
        self.df_tracks_filt = None

    # ---------- public ----------
    def display(self):
        display(self.ui)

    # ---------- internals ----------
    def _effective_values(self):
        exp_duration = int(self.exp_duration.value) if self.en_exp_duration.value else 999999
        min_len      = int(self.min_track_length.value) if self.en_min_len.value else 0
        max_len      = int(self.max_track_length.value) if self.en_max_len.value else 999999
        min_size     = int(self.min_size.value) if self.en_min_size.value else 0
        time_type    = str(self.time_type_toggle.value)
        return exp_duration, min_len, max_len, min_size, time_type

    def _persist_params(self):
        params = self._params
        prof = params["track_filtering"][self.cell_type]

        prof["exp_duration_enabled"]      = bool(self.en_exp_duration.value)
        prof["min_track_length_enabled"]  = bool(self.en_min_len.value)
        prof["max_track_length_enabled"]  = bool(self.en_max_len.value)
        prof["min_size_enabled"]          = bool(self.en_min_size.value)

        prof["exp_duration"]     = int(self.exp_duration.value)
        prof["min_track_length"] = int(self.min_track_length.value)
        prof["max_track_length"] = int(self.max_track_length.value)
        prof["min_size"]         = int(self.min_size.value)
        prof["time_type"]        = str(self.time_type_toggle.value)

        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(params, f, sort_keys=False)

    def _lock(self, state: bool):
        for w in [
            self.en_exp_duration, self.exp_duration,
            self.en_min_len, self.min_track_length,
            self.en_max_len, self.max_track_length,
            self.time_type_toggle,
            self.en_min_size, self.min_size,
            self.btn_run
        ]:
            if hasattr(w, "disabled"):
                w.disabled = state

    def _on_run_clicked(self, _):
        self._lock(True)
        self.spinner_html.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                out_dir = Path(self.metadata_loader.output_dir).expanduser()
                out_dir.mkdir(parents=True, exist_ok=True)

                exp_duration, min_len, max_len, min_size, time_type = self._effective_values()
                print("▶️ Filtering tracks…")
                print(f"  exp_duration={exp_duration}, min_len={min_len}, max_len={max_len}, min_size={min_size}, unit={'hours' if time_type=='real_time' else 'frames'}")

                # persist UI state
                self._persist_params()

                # call your function
                self.df_tracks_filt = filter_organoid_tracks(
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(out_dir),
                    exp_duration=exp_duration,
                    min_track_length=min_len,
                    max_track_length=max_len,
                    min_size=min_size,
                    cell_type=self.cell_type,
                    time_type=time_type
                )
                print("✅ Filtering complete.")

                # locations
                analysis_outdir = out_dir / "analysis" / self.cell_type
                feature_outdir  = analysis_outdir / "track_features"
                qc_outdir       = analysis_outdir / "quality_control"

                combined_csv = feature_outdir / f"BEHAV3D_{self.cell_type}_combined_track_features.csv"
                filtered_csv = feature_outdir / f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv"
                filter_plot  = qc_outdir      / "BEHAV3D_filter_counts.pdf"

                print("\nOutputs:")
                print(f"  Combined tracks: {combined_csv}")
                print(f"  Filtered tracks: {filtered_csv}")
                print(f"  Filter counts:   {filter_plot}")

                # preview
                try:
                    from IPython.display import display as _disp
                    print("\n— Filtered tracks (head) —")
                    _disp(self.df_tracks_filt.head(10))
                except Exception:
                    print(f"(preview unavailable) shape={getattr(self.df_tracks_filt,'shape',None)}")

                # quick per-sample summary before/after
                try:
                    group_cols = ['sample_name', 'organoid_line', 'tcell_line', 'exp_nr', 'well']

                    def _count(df):
                        return (
                            df.groupby(group_cols, dropna=False)["TrackID"]
                              .nunique()
                              .reset_index(name="nr_tracks")
                              .groupby("sample_name", dropna=False)["nr_tracks"]
                              .sum()
                              .reset_index()
                        )

                    before_df = pd.read_csv(combined_csv)
                    after_df  = pd.read_csv(filtered_csv)
                    before_counts = _count(before_df).rename(columns={"nr_tracks": "tracks_before"})
                    after_counts  = _count(after_df).rename(columns={"nr_tracks": "tracks_after"})
                    summary = before_counts.merge(after_counts, on="sample_name", how="outer").fillna(0)
                    summary["removed"] = summary["tracks_before"] - summary["tracks_after"]

                    print("\n— Track counts by sample —")
                    _disp(summary)
                except Exception as e:
                    print(f"(Could not compute summary counts: {e})")

                print("\n✅ Finished.")
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)
                
class OrganoidAnalysisPanel:
    """Dead-%-only panel. Persists to behav3d_parameters['analysis']['organoid'] and writes YAML on run."""
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = cell_type

        # --- load + seed config from _DEFAULT_CONFIG (preserve user values) ---
        params_on_disk = self.metadata_loader.behav3d_parameters or {}
        params = deepcopy(params_on_disk)
        _deep_merge(params, deepcopy(_DEFAULT_CONFIG))
        if params != params_on_disk:
            with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                yaml.safe_dump(params, f, sort_keys=False)
        self._params = params
        self._org = params["analysis"][self.cell_type]

        # --- UI ---
        self.title = widgets.HTML("<b>Organoid Death Analysis</b>")

        self.dead_perc_threshold = widgets.FloatText(
            description="Dead % threshold",
            value=float(self._org.get("dead_perc_threshold", 1e-7)),
            style={'description_width': '160px'},
            continuous_update=False,
        )

        self.btn_run = widgets.Button(
            description="Run organoid analysis",
            button_style="success",
            layout=widgets.Layout(width="360px")   # wider button
        )
        self.btn_run.on_click(self._on_run_clicked)

        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"

        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html],
            layout=widgets.Layout(align_items="center", gap="10px")
        )

        self.out = widgets.Output()

        self.ui = widgets.VBox([
            self.title,
            self.dead_perc_threshold,
            self.run_row,
            widgets.HTML("<hr>"),
            self.out
        ])

    # ---------- public ----------
    def display(self):
        display(self.ui)

    # ---------- internals ----------
    def _resolve_output_dir(self) -> Path:
        ml_out = getattr(self.metadata_loader, "output_dir", None)
        if ml_out:
            return Path(ml_out).expanduser()
        return Path(self._params["paths"]["output_dir"]).expanduser()

    def _persist_params(self):
        self._org["dead_perc_threshold"] = float(self.dead_perc_threshold.value)
        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(self._params, f, sort_keys=False)

    def _lock(self, locked: bool):
        self.dead_perc_threshold.disabled = locked
        self.btn_run.disabled = locked

    # ---------- run ----------
    def _on_run_clicked(self, _):
        self._lock(True)
        self.spinner_html.layout.display = None
        with self.out:
            self.out.clear_output()
            try:
                out_dir = self._resolve_output_dir()
                out_dir.mkdir(parents=True, exist_ok=True)

                # overwrite config values with current UI selection
                self._persist_params()

                dth = float(self.dead_perc_threshold.value)
                print(f"▶️ Running {self.cell_type} analysis…")
                print(f"  dead_perc_threshold = {dth}")

                run_organoid_analysis(
                    dead_perc_threshold=dth,
                    output_dir=str(out_dir),
                    df_tracks_path=None,
                    org_type=self.cell_type
                )

                # show outputs (optional preview)
                results_outdir = out_dir / "analysis" / self.cell_type / "results"
                general_csv = results_outdir / f"combined_general_{self.cell_type}_dynamics_analysis.csv"
                general_pdf = results_outdir / f"combined_general_{self.cell_type}_dynamics_analysis.pdf"

                print("\nOutputs:")
                print(f"  Results CSV: {general_csv}")
                print(f"  Results PDF: {general_pdf}")

                if general_csv.exists():
                    try:
                        df_general = pd.read_csv(general_csv)
                        print("\n— General results (head) —")
                        display(df_general.head(10))
                    except Exception as e:
                        print(f"(Could not read {general_csv}: {e})")

                print("\n✅ Finished.")
            except Exception:
                import traceback; traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)
                
                
class BackprojectionPanel:
    def __init__(self, metadata_loader, cell_type="tcell"):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        # ---- load groups ----
        try:
            self._base_groups = behav3d_calculated_features  # dict[group] -> list[str/pattern]
            md = self.metadata_loader.metadata
            self.two_org_types = False
            # Conditional addition of organoid_2
            if md is not None and 'organoid_2_tracks_csv_path' in md.columns:
                s = md['organoid_2_tracks_csv_path'].dropna()
                if not s.empty and Path(str(s.iloc[0])).exists():
                    self.two_org_types = True
                    # Add organoid_2 contact features
                    self._base_groups["contact"].extend([
                        "organoid_2_contact",
                        "organoid_2_contact_pixels",
                    ])
        except NameError:
            raise RuntimeError("Define behav3d_calculated_features before creating BackprojectionPanel.")

        # ---- config bootstrap ----
        params = dict(self.metadata_loader.behav3d_parameters or {})
        params.setdefault("backprojection", self._default_backproj_config(self._base_groups, cell_type))
        self._params = params
        self._cfg = self._params["backprojection"]

        # ---- headings ----
        self.title     = widgets.HTML('<div style="font-size:22px;font-weight:700;">Backproject features to images (napari)</div>')
        self.sel_title = widgets.HTML('<div style="font-size:20px;font-weight:700;">Select feature columns</div>')

        # ---- sample selector ----
        md = self.metadata_loader.metadata
        if md is None or "sample_name" not in md.columns:
            raise RuntimeError("metadata_loader.metadata must have a 'sample_name' column.")
        sample_list = sorted(map(str, md["sample_name"].unique().tolist()))
        last_sample = self._cfg.get("last_sample") or (sample_list[0] if sample_list else None)
        self.sample_dd = widgets.Dropdown(description="Sample", options=sample_list, value=last_sample,
                                          layout=widgets.Layout(width="360px"))
        self.sample_dd.style.description_width = "80px"

        # ---- cell type toggle ----
        # Safely check for the column and for at least one existing, non-empty path before using it.
        if self.two_org_types:
           self._celltype_map = {"T cell": "tcell", "Organoid": "organoid", "Organoid 2": "organoid_2"}
        else:
            self._celltype_map = {"T cell": "tcell", "Organoid": "organoid"}

        inv_ct_map = {v: k for k, v in self._celltype_map.items()}
        ct_default = inv_ct_map.get(str(self._cfg.get("cell_type", cell_type)).lower(), "T cell")
        self.celltype_tb = widgets.ToggleButtons(
            options=list(self._celltype_map.keys()),
            value=ct_default,
            description="Cell type",
            layout=widgets.Layout(width="360px")
        )
        self.celltype_tb.style.button_width = "160px"
        self.celltype_tb.style.description_width = "80px"

        # ---- mode (mean/time) ----
        self._mode_map = {"Mean features": "mean", "Time features": "time"}
        inv_mode_map = {v: k for k, v in self._mode_map.items()}
        self.mode_tb = widgets.ToggleButtons(
            options=list(self._mode_map.keys()),
            value=inv_mode_map.get(self._cfg.get("mode", "mean"), "Mean features"),
            description="Mode",
            layout=widgets.Layout(width="360px")
        )
        self.mode_tb.style.button_width = "160px"
        self.mode_tb.style.description_width = "80px"

        # ---- save toggle ----
        self.save_cb = widgets.Checkbox(
            description="Save .zarr backprojection to disk",
            value=bool(self._cfg.get("save", False)),
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )

        # ---- status + refresh ----
        self.status_html = widgets.HTML("")
        self.refresh_btn = widgets.Button(description="Refresh", icon="refresh", layout=widgets.Layout(width="120px"))
        self.refresh_btn.on_click(self._on_refresh_clicked)

        # ---- TIME MODE selection (grouped patterns; 3 columns) ----
        self._time_group_rows = {}  # group -> {group_cb, child_cbs, grid, container}
        self.time_selection_box = self._build_time_selection_box()

        # ---- MEAN MODE selection (grouped actual columns from UMAP_clusters.csv; 3 columns) ----
        self._mean_group_rows = {}  # group -> {group_cb, child_cbs, grid, container}
        self.mean_selection_box = widgets.VBox([])
        self._rebuild_mean_selection_from_file()  # initial

        # ---- selection container (swaps by mode) ----
        self.selection_container = widgets.VBox([])
        self._swap_selection_ui()

        # react to changes
        self.mode_tb.observe(lambda ch: (self._swap_selection_ui(), self._update_status_for_mode()), names='value')
        self.celltype_tb.observe(lambda ch: (self._on_celltype_changed(), self._update_status_for_mode()), names='value')

        # ---- Select/Clear (affects current UI) ----
        self.btn_select_all = widgets.Button(description="Select all")
        self.btn_clear_all  = widgets.Button(description="Clear")
        self.btn_select_all.on_click(lambda *_: self._set_all(True))
        self.btn_clear_all.on_click(lambda *_: self._set_all(False))

        # ---- Run + spinner + Stop viewer ----
        self.btn_run = widgets.Button(
            description="Run backprojection",
            button_style="success",
            layout=widgets.Layout(width="fit-content", flex="0 0 auto")
        )
        self.btn_run.on_click(self._on_run_clicked)

        # uses a global `spinning_loader` HTML snippet you've defined elsewhere
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"

        self.btn_stop = widgets.Button(
            description="Close viewer",
            button_style="danger",
            icon="stop",
            tooltip="Close the active Napari viewer",
            layout=widgets.Layout(width="fit-content", display="none")
        )
        self.btn_stop.on_click(self._on_stop_clicked)

        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html, self.btn_stop],
            layout=widgets.Layout(align_items="center", gap="10px")
        )

        self.out = widgets.Output()

        # ---- UI ----
        self.ui = widgets.VBox([
            self.title,
            self.sample_dd,
            self.celltype_tb,   # under sample
            self.mode_tb,       # under cell type
            self.save_cb,
            widgets.HBox([self.status_html, self.refresh_btn], layout=widgets.Layout(align_items="center", gap="12px")),
            self.sel_title,
            widgets.HBox([self.btn_select_all, self.btn_clear_all], layout=widgets.Layout(gap="8px")),
            self.selection_container,
            widgets.HTML("<hr>"),
            self.run_row,
            self.out
        ])

        self._viewer = None
        self._update_status_for_mode()

    # ---------- defaults ----------
    def _default_backproj_config(self, groups_dict, cell_type_default):
        return {
            "cell_type": str(cell_type_default),
            "mode": "mean",                # "mean" or "time"
            "save": False,
            "last_sample": None,
            "columns_input_time": [],
            "columns_input_mean": [],
        }

    # ---------- helpers for paths & status ----------
    def _current_results_csv_path(self):
        cell_type = self._celltype_map[self.celltype_tb.value]
        mode = self._mode_map[self.mode_tb.value]
        results_dir = Path(self.output_dir, "analysis", cell_type, "results")
        if mode == "mean":
            return Path(results_dir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv")
        else:
            # not strictly required for selection, but useful for status
            return Path(results_dir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv")

    def _update_status_for_mode(self):
        p = self._current_results_csv_path()
        if self._mode_map[self.mode_tb.value] == "mean":
            if not p.exists():
                self.status_html.value = f"<span style='color:#b00'>Cannot find results yet at:<br><code>{p}</code></span>"
            else:
                self.status_html.value = f"<span style='color:#070'>Results detected:<br><code>{p.name}</code></span>"
        else:
            # time mode does not need a file to show patterns
            if p.exists():
                self.status_html.value = f"<span style='color:#070'>Data present (header available):<br><code>{p.name}</code></span>"
            else:
                self.status_html.value = "<span>Ready (no file required for selection)</span>"

    def _on_refresh_clicked(self, *_):
        # Rebuild mean selection (since it depends on the file), then update status.
        if self._mode_map[self.mode_tb.value] == "mean":
            self._rebuild_mean_selection_from_file()
            self._swap_selection_ui()
        self._update_status_for_mode()

    # ---------- UI builders ----------
    def _grid_for(self, children, indent_px="24px", columns=3):
        return widgets.GridBox(
            children,
            layout=widgets.Layout(
                grid_template_columns=" ".join(["max-content"] * columns),
                grid_gap="6px 18px",
                margin=f"0 0 0 {indent_px}",
            )
        )

    def _build_time_selection_box(self):
        preset_time = list(self._cfg.get("columns_input_time", self._cfg.get("columns_input", [])))
        preset_set = set(preset_time)

        rows = {}
        boxes = []
        for group_name, feats in self._base_groups.items():
            child_cbs = [widgets.Checkbox(value=(f in preset_set), description=f, indent=True) for f in feats]
            gcb = widgets.Checkbox(value=all(cb.value for cb in child_cbs), indent=False)
            glabel = widgets.HTML(f"<b>{group_name}</b>")
            header = widgets.HBox([gcb, glabel])

            grid = self._grid_for(child_cbs, columns=3)
            container = widgets.VBox([header, grid])

            # handlers
            def make_group_handler(gc, cc_list):
                def _on(change):
                    if change["name"] != "value": return
                    val = bool(change["new"])
                    for cb in cc_list: cb.unobserve(child_handler, names="value")
                    for cb in cc_list: cb.value = val
                    for cb in cc_list: cb.observe(child_handler, names="value")
                return _on

            def make_child_handler(gc, cc_list):
                def _on(change):
                    if change["name"] != "value": return
                    all_on = all(cb.value for cb in cc_list)
                    gc.unobserve(group_handler, names="value")
                    gc.value = all_on
                    gc.observe(group_handler, names="value")
                return _on

            group_handler = make_group_handler(gcb, child_cbs)
            child_handler = make_child_handler(gcb, child_cbs)
            gcb.observe(group_handler, names="value")
            for cb in child_cbs:
                cb.observe(child_handler, names="value")

            rows[group_name] = {"group_cb": gcb, "child_cbs": child_cbs, "grid": grid, "container": container}
            boxes.append(container)

        self._time_group_rows = rows
        return widgets.VBox(boxes)

    def _rebuild_mean_selection_from_file(self):
        """Build grouped checkboxes from the UMAP_clusters.csv header, using mean_ prefixed versions of base names."""
        # File for current cell type
        clusters_csv = self._current_results_csv_path()
        cols = []
        if clusters_csv.exists():
            try:
                cols = pd.read_csv(clusters_csv, nrows=1).columns.tolist()
            except Exception:
                cols = []

        # Prepare target patterns per group: add mean_ if not already present
        def to_mean_pattern(name: str) -> str:
            return name if str(name).startswith("mean_") else ("mean_" + str(name))

        group_to_columns = {}
        for gname, base_feats in self._base_groups.items():
            present = []
            for raw in base_feats:
                patt = to_mean_pattern(raw)
                if any(ch in patt for ch in "*?["):
                    present.extend([c for c in cols if fnmatch.fnmatchcase(str(c), patt)])
                else:
                    if patt in cols:
                        present.append(patt)
            # de-dup keep order
            seen = set(); present = [x for x in present if not (x in seen or seen.add(x))]
            group_to_columns[gname] = present

        # Build rows or placeholder
        preset_mean = set(self._cfg.get("columns_input_mean", []))
        rows = {}
        boxes = []

        if not cols:  # file missing or unreadable
            msg = widgets.HTML("<i>Cannot find results yet. Run analysis, then click <b>Refresh</b>.</i>")
            self._mean_group_rows = {}
            self.mean_selection_box.children = (msg,)
            return

        for gname, gcols in group_to_columns.items():
            if not gcols:
                # Keep empty groups visible for consistent layout (can hide if you prefer)
                child_cbs = []
            else:
                child_cbs = [widgets.Checkbox(value=(c in preset_mean), description=c, indent=True) for c in gcols]
            gcb = widgets.Checkbox(value=(bool(child_cbs) and all(cb.value for cb in child_cbs)), indent=False)
            glabel = widgets.HTML(f"<b>{gname}</b>")
            header = widgets.HBox([gcb, glabel])

            grid = self._grid_for(child_cbs, columns=3)
            container = widgets.VBox([header, grid])

            # handlers
            def make_group_handler(gc, cc_list):
                def _on(change):
                    if change["name"] != "value": return
                    val = bool(change["new"])
                    for cb in cc_list: cb.value = val
                return _on

            def make_child_handler(gc, cc_list):
                def _on(change):
                    if change["name"] != "value": return
                    gc.value = (bool(cc_list) and all(cb.value for cb in cc_list))
                return _on

            gh = make_group_handler(gcb, child_cbs)
            ch = make_child_handler(gcb, child_cbs)
            gcb.observe(gh, names="value")
            for cb in child_cbs:
                cb.observe(ch, names="value")

            rows[gname] = {"group_cb": gcb, "child_cbs": child_cbs, "grid": grid, "container": container}
            boxes.append(container)

        self._mean_group_rows = rows
        self.mean_selection_box.children = boxes

    def _swap_selection_ui(self):
        mode = self._mode_map[self.mode_tb.value]
        if mode == "mean":
            # make sure mean UI is up-to-date with current cell type's file
            self._rebuild_mean_selection_from_file()
            self.selection_container.children = (self.mean_selection_box,)
        else:
            self.selection_container.children = (self.time_selection_box,)

    def _on_celltype_changed(self):
        # Rebuild mean-selection UI because the source CSV depends on cell type
        if self._mode_map[self.mode_tb.value] == "mean":
            self._rebuild_mean_selection_from_file()

    # ---------- general helpers ----------
    def display(self): 
        display(self.ui)

    def _set_all(self, val: bool):
        if self._mode_map[self.mode_tb.value] == "mean":
            rows = self._mean_group_rows
        else:
            rows = self._time_group_rows
        for row in rows.values():
            row["group_cb"].value = val
            for cb in row["child_cbs"]:
                cb.value = val

    def _selected_time_patterns(self):
        seen, sel = set(), []
        for row in self._time_group_rows.values():
            for cb in row["child_cbs"]:
                if cb.value and cb.description not in seen:
                    sel.append(cb.description); seen.add(cb.description)
        return sel

    def _selected_mean_columns(self):
        seen, sel = set(), []
        for row in self._mean_group_rows.values():
            for cb in row["child_cbs"]:
                if cb.value and cb.description not in seen:
                    sel.append(cb.description); seen.add(cb.description)
        return sel

    def _expand_patterns(self, selected, available_columns):
        if not available_columns:
            return list(dict.fromkeys(selected))
        avail = list(available_columns.tolist() if hasattr(available_columns, "tolist") else list(available_columns))
        out = []
        for name in selected:
            if any(ch in name for ch in "*?["):
                matches = [c for c in avail if fnmatch.fnmatchcase(str(c), name)]
                out.extend(matches if matches else [name])
            else:
                out.append(name)
        # de-dup preserve order
        return list(dict.fromkeys(out))

    def _lock(self, state: bool):
        for w in [self.sample_dd, self.celltype_tb, self.mode_tb, self.save_cb,
                  self.btn_select_all, self.btn_clear_all, self.btn_run, self.refresh_btn]:
            w.disabled = state
        for rows in (self._time_group_rows, self._mean_group_rows):
            for row in rows.values():
                row["group_cb"].disabled = state
                for cb in row["child_cbs"]:
                    cb.disabled = state

    def _persist(self, *, resolved=None):
        mode = self._mode_map[self.mode_tb.value]
        self._cfg["cell_type"] = self._celltype_map[self.celltype_tb.value]
        self._cfg["mode"] = mode
        self._cfg["save"] = bool(self.save_cb.value)
        self._cfg["last_sample"] = self.sample_dd.value

        if mode == "time":
            self._cfg["columns_input_time"] = list(self._selected_time_patterns())
        else:
            self._cfg["columns_input_mean"] = list(self._selected_mean_columns())

        if resolved is not None:
            self._cfg["columns_resolved"] = list(resolved)

        self._params["backprojection"] = self._cfg
        self.metadata_loader.behav3d_parameters = self._params
        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(self._params, f, sort_keys=False)

    # ---------- run ----------
    def _on_run_clicked(self, *_):
        self._lock(True)
        self.spinner_html.layout.display = None
        self.btn_stop.layout.display = "inline-block"  # allow user to try closing viewer
        self.out.clear_output()

        with self.out:
            patched = False
            try:
                sample = self.sample_dd.value
                cell_type = self._celltype_map[self.celltype_tb.value]
                mode = self._mode_map[self.mode_tb.value]
                save = bool(self.save_cb.value)

                print(f"▶️ Backprojecting… sample={sample} | cell_type={cell_type} | mode={mode} | save={save}")

                # Resolve columns to use
                if mode == "mean":
                    cols_to_use = self._selected_mean_columns()
                    resolved = cols_to_use
                else:
                    results_dir = Path(self.output_dir, "analysis", cell_type, "results")
                    clustered_csv = Path(results_dir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv")
                    try:
                        header_cols = pd.read_csv(clustered_csv, nrows=1).columns
                    except Exception:
                        header_cols = None
                    patterns = self._selected_time_patterns()
                    resolved = self._expand_patterns(patterns, header_cols)
                    cols_to_use = resolved

                self._persist(resolved=resolved)

                # Respect Save checkbox by temporarily disabling save_as_zarr when save=False
                if not save:
                    try:
                        import behav3d.utils.fileio as fio
                        self._orig_save_as_zarr = fio.save_as_zarr
                        def _noop(*args, **kwargs): return None
                        fio.save_as_zarr = _noop
                        patched = True
                    except Exception:
                        print("⚠️ Could not patch save_as_zarr; results may be written.")

                # Choose function
                fn = backproject_mean_features_behav3d if mode == "mean" else backproject_time_features_behav3d

                # Run (these functions launch napari internally)
                _ = fn(
                    metadata=self.metadata_loader.metadata,
                    sample_name=sample,
                    config=None,
                    output_dir=self.output_dir,
                    cell_type=cell_type,
                    columns=cols_to_use if cols_to_use else [],  # [] lets function choose default
                    save=save
                )

                print("✅ Backprojection finished (napari was launched inside the function).")

            except Exception:
                import traceback; traceback.print_exc()
            finally:
                try:
                    if patched:
                        import behav3d.utils.fileio as fio
                        fio.save_as_zarr = self._orig_save_as_zarr
                except Exception:
                    pass
                self.spinner_html.layout.display = "none"
                self._lock(False)

    def _on_stop_clicked(self, *_):
        """Best-effort close of any open napari window."""
        try:
            import napari
            from qtpy.QtWidgets import QApplication
            app = QApplication.instance()
            if app is not None:
                for w in list(app.topLevelWidgets()):
                    name = w.__class__.__name__.lower()
                    if "napari" in name or "mainwindow" in name:
                        try: w.close()
                        except Exception: pass
                app.processEvents()
            else:
                try:
                    v = napari.current_viewer()
                    if v is not None:
                        v.close()
                except Exception:
                    pass
        except Exception as e:
            with self.out:
                print(f"Close viewer attempt: {e}")