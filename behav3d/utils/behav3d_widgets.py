# behav3d_widgets.py
import random
import re

import ipywidgets as widgets
from ipyfilechooser import FileChooser
from IPython.display import display, clear_output
from behav3d.utils import detect_organoid_types_from_metadata, detect_immune_cell_types_from_metadata, detect_other_cell_types_from_metadata
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
from behav3d.preprocessing.tracking.laptracking import run_laptracking
from behav3d.preprocessing.tracking.trackpy_tracking import run_trackpy_tracking_generic
from behav3d.preprocessing.tracking.propagation_tracking import run_propagation_tracking
import json
from copy import deepcopy
import yaml
import fnmatch
from behav3d.analysis import summarize_track_features
from behav3d.analysis.feature_extraction import run_feature_extraction
# Analysis functions imported from adapted tcell_analysis and organoid_analysis
from behav3d.analysis.tcell_analysis import filter_cell_tracks, run_tcell_analysis
from behav3d.analysis.organoid_analysis import filter_organoid_tracks, run_organoid_analysis

from behav3d.analysis.backprojection import backproject_mean_features_behav3d, backproject_time_features_behav3d
import napari

from functools import partial


# ===============================
# UNIFIED CATEGORY DETECTION
# ===============================
def detect_cell_type_category(cell_type, metadata):
    """
    Unified category detection for ANY panel in the pipeline.
    
    Determines if a cell_type belongs to 'organoid', 'immune', or 'other' category
    by checking metadata column prefixes (or_, im_, ot_).
    
    Parameters
    ----------
    cell_type : str
        The cell type name (e.g., 'organoid1', 'tcell', 'macro', 'tum')
    metadata : pd.DataFrame
        The loaded metadata with prefixed columns
        
    Returns
    -------
    str
        One of: 'organoid', 'immune', 'other'
        
    Raises
    ------
    ValueError
        If cell_type not found in any detected category
        
    Examples
    --------
    >>> detect_cell_type_category('organoid1', metadata)
    'organoid'
    >>> detect_cell_type_category('tcell', metadata)
    'immune'
    >>> detect_cell_type_category('tumor', metadata)
    'other'
    """
    from behav3d.utils import (
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
        # Fallback: check if "organoid" substring exists
        if 'organoid' in cell_type.lower():
            return 'organoid'
        else:
            raise ValueError(
                f"Cell type '{cell_type}' not found in metadata. "
                f"Detected: organoid={organoid_types}, immune={immune_types}, other={other_types}. "
                f"Ensure metadata has columns with prefixes: or_{cell_type}_*, im_{cell_type}_*, or ot_{cell_type}_*"
            )


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
        # Contact features are DYNAMIC - generated for ALL cell types in metadata:
        # - {cell_type}_contact           (bool)
        # - {cell_type}_contact_pixels    (int)  
        # - touching_{cell_type}s         (str)
        # - active_{cell_type}_contact    (float) - works for ANY cell type
        "*_contact",
        "*_contact_pixels",
        "touching_*",
        "active_*_contact",
    ],
    "active_killing": [
        # Active killing features from advanced_feature_extraction
        # Only available if Active Killing Analysis has been run
        # Global features (across all target types):
        "is_active_killing",           # Boolean: True if this timepoint is active killing
        "killing_efficiency",          # Ratio of actual vs expected background death
        # Note: death_signal_increase_*tp excluded from DTW - not suitable for time-series comparison
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
        # Dynamic EDT thresholds are stored as "{celltype}_edt_threshold"

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
        # Specific cell types will be added dynamically
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
        # Category templates - specific cell types inherit these settings
        # Cell types are auto-detected from metadata (im_*, or_*, ot_* prefixes)
        "immune": {
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
        "other": {
            "dead_mask_percentage_threshold": 0.10,
            "features_choice": ["movement", "intensity", "contact"],
            "contact_threshold": 0,
            "n_workers": 8,
            "overwrite": False
        }
    },
    "track_filtering": {
        # Category templates - specific cell types inherit these settings
        "immune": {
            "exp_duration": 24.0,
            "exp_duration_enabled": False,
            "min_track_length": 0, # frames
            "min_track_length_enabled": True,
            "max_track_length": 999999, # frames
            "max_track_length_enabled": True,
            "filter_t0_dead": True,
            "filter_t0_dead_enabled": True,  # Enabled - remove dead cells added at start
        },
        "organoid": {
            "exp_duration_enabled": False,
            "exp_duration": 24.0,        # hours (disabled by default - organoids persist)
            "min_track_length_enabled": False,
            "min_track_length": 50,      # frames (disabled - adjust based on experiment length)
            "max_track_length_enabled": False,
            "max_track_length": 999999,  # frames (disabled - organoids should persist)
            "filter_t0_dead": True,
            "filter_t0_dead_enabled": True,
        },
        "other": {
            "exp_duration": 24.0,
            "exp_duration_enabled": False,
            "min_track_length": 0, # frames
            "min_track_length_enabled": True,
            "max_track_length": 999999, # frames
            "max_track_length_enabled": True,
            "filter_t0_dead": False,
            "filter_t0_dead_enabled": False
        }
    },
    "analysis": {
        # Category templates - specific cell types inherit these settings
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
            'dtw_features_input': [],  # Empty list = use groups
            "dtw_features_resolved": [],   # expanded at run
            'z_normalize': {},  # Auto-determined at runtime
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
            'dtw_features_input': [],  # Empty list = use groups
            "dtw_features_resolved": [],
            'z_normalize': {},  # Auto-determined at runtime
            "dead_perc_threshold": 0.02
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
            'dtw_features_input': [],  # Empty list = use groups
            "dtw_features_resolved": [],
            'z_normalize': {},  # Auto-determined at runtime
            "dead_perc_threshold": 0.10
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
    },
    "active_killing": {
        "observation_window": 5,
        "death_signal_column": "mean_dead_dye",
        "killing_threshold_multiplier": 1.5,
        "min_contact_duration": 1,
        "save_results": True,
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

class MetadataBuilder(widgets.VBox):
    """
    Interactive widget for creating/editing BEHAV3D metadata CSV files.
    Supports:
    - Dynamic number of samples
    - Organoid, immune, and other cell type populations
    - Dead channel only (yes/no toggle)
    - Fill-down from sample 1
    - Load existing CSV for editing
    """
    def __init__(self, output_dir_picker=None, **kwargs):
        super().__init__(**kwargs)
        self.output_dir_picker = output_dir_picker
        
        # Population configuration
        self.n_organoid_types = 0
        self.n_immune_types = 0
        self.n_other_types = 0
        self.organoid_names = []
        self.immune_names = []
        self.other_names = []
        
        # Channel names configuration
        self.include_channels = False
        self.n_channels = 0
        
        # Sample data storage
        self.sample_forms = []
        
        # Store loaded dataframe for re-population
        self.loaded_df = None
        
        self._build_ui()
    
    def _build_ui(self):
        """Build the initial configuration UI"""
        # Header
        header = widgets.HTML('<h3>Metadata Builder</h3>')
        instructions = widgets.HTML(
            '<p>Define the number of samples and configure cell type populations. '
            'Fields marked with * are required.</p>'
        )
        
        # Number of samples
        self.n_samples_input = widgets.IntText(
            value=1,
            description='Number of samples*:',
            style={'description_width': '150px'}
        )
        
        # Load CSV button (NEW)
        self.btn_load_csv = widgets.Button(
            description='Load Existing CSV',
            button_style='info',
            icon='upload',
            tooltip='Load an existing metadata CSV to edit'
        )
        self.btn_load_csv.on_click(self._on_load_csv)
        
        # Configure populations button
        self.btn_configure = widgets.Button(
            description='Configure Cell Types',
            button_style='primary',
            icon='cog'
        )
        self.btn_configure.on_click(self._on_configure_populations)
        
        # Container for dynamic content
        self.main_container = widgets.VBox()
        
        # Assembly
        self.children = [
            header,
            instructions,
            widgets.HBox([self.n_samples_input, self.btn_load_csv]),
            self.btn_configure,
            self.main_container
        ]
    
    def _on_load_csv(self, btn):
        """Load existing CSV and populate the form"""
        if self.output_dir_picker is None:
            self.main_container.children = [widgets.HTML('<p style="color:red">⚠️ No output directory picker provided</p>')]
            return
        
        output_dir = Path(self.output_dir_picker.value)
        if not output_dir.exists():
            self.main_container.children = [widgets.HTML('<p style="color:red">⚠️ Output directory does not exist</p>')]
            return
        
        # Look for metadata.csv in the output directory ***pending: create a path to load the name.cvs you want
        csv_path = output_dir / 'metadata.csv'
        if not csv_path.exists():
            self.main_container.children = [widgets.HTML(f'<p style="color:red">⚠️ No metadata.csv found in {output_dir}</p>')] #
            return
        
        try:
            df = pd.read_csv(csv_path)
            self._populate_from_dataframe(df)
            self.main_container.children = [widgets.HTML('<p style="color:green">✅ CSV loaded successfully! Scroll down to edit.</p>')]
        except Exception as e:
            self.main_container.children = [widgets.HTML(f'<p style="color:red">⚠️ Error loading CSV: {e}</p>')]
    
    def _populate_from_dataframe(self, df):
        """Populate the form from an existing DataFrame"""
        
        # Store loaded dataframe for re-population
        self.loaded_df = df.copy()
        
        # Detect cell types from columns (with new prefix system: or_, im_, ot_)
        organoid_types = detect_organoid_types_from_metadata(df)
        immune_types = detect_immune_cell_types_from_metadata(df)
        other_types = detect_other_cell_types_from_metadata(df)
               
        self.n_organoid_types = len(organoid_types)
        self.n_immune_types = len(immune_types)
        self.n_other_types = len(other_types)
        self.organoid_names = organoid_types
        self.immune_names = immune_types
        self.other_names = other_types
        
        # Detect channel columns (pattern: channel_N_label)
        channel_cols = [col for col in df.columns if re.match(r'^channel_\d+_label$', col)]
        if channel_cols:
            self.include_channels = True
            self.n_channels = len(channel_cols)
        else:
            self.include_channels = False
            self.n_channels = 0
        
        # Set number of samples
        self.n_samples_input.value = len(df)
        
        # Build the form
        self._build_data_entry_form()
    
    def _populate_form_fields(self, df):
        """Populate form fields from DataFrame"""
        # Populate fields from DataFrame
        for idx, row in df.iterrows():
            if idx >= len(self.sample_forms):
                break
            
            form = self.sample_forms[idx]
            
            # Channel name fields
            for field_name, widget in form['channels'].items():
                if field_name in row and pd.notna(row[field_name]):
                    widget.value = str(row[field_name]).strip()
            
            # Basic fields
            for field_name, widget in form['basic'].items():
                if field_name in row and pd.notna(row[field_name]):
                    try:
                        if isinstance(widget, widgets.IntText):
                            widget.value = int(row[field_name])
                        elif isinstance(widget, widgets.FloatText):
                            widget.value = float(row[field_name])
                        else:
                            widget.value = str(row[field_name])
                    except (ValueError, TypeError):
                        # If conversion fails, try as string
                        widget.value = str(row[field_name]) if row[field_name] != '' else ''
            
            # Cell type fields with line_condition splitting
            # Need to determine prefix for each cell type
            for cell_type, fields_dict in form['cell_types'].items():
                # Determine prefix based on which category this cell type belongs to
                if cell_type in self.organoid_names:
                    prefix = 'or'
                elif cell_type in self.immune_names:
                    prefix = 'im'
                elif cell_type in self.other_names:
                    prefix = 'ot'
                else:
                    continue  # Skip unknown types
                
                # Handle line_condition splitting with prefix
                line_cond_col = f'{prefix}_{cell_type}_line_condition'
                if line_cond_col in row and pd.notna(row[line_cond_col]) and str(row[line_cond_col]).strip() != '':
                    # Split "32_ko" -> line="32", condition="ko"
                    value = str(row[line_cond_col]).strip()
                    parts = value.split('_', 1)
                    if 'line' in fields_dict:
                        fields_dict['line'].value = parts[0] if parts[0] else ''
                    if 'condition' in fields_dict and len(parts) > 1:
                        fields_dict['condition'].value = parts[1] if parts[1] else ''
                    elif 'condition' in fields_dict:
                        fields_dict['condition'].value = ''
                
                # Load optional path fields if they exist
                for path_field in ['segments_image_path', 'tracks_image_path', 'tracks_csv_path']:
                    col_name = f'{prefix}_{cell_type}_{path_field}'
                    if col_name in row and pd.notna(row[col_name]) and path_field in fields_dict:
                        fields_dict[path_field].value = str(row[col_name]).strip()
            
            # Dead channel (including mask_path)
            if 'dead_channel' in row and pd.notna(row['dead_channel']) and str(row['dead_channel']).strip() != '':
                try:
                    form['dead_channel']['enabled'].value = True
                    form['dead_channel']['number'].value = int(row['dead_channel'])
                    # Load dead_mask_path if present
                    if 'dead_mask_path' in row and pd.notna(row['dead_mask_path']) and 'mask_path' in form['dead_channel']:
                        form['dead_channel']['mask_path'].value = str(row['dead_mask_path']).strip()
                except (ValueError, TypeError):
                    # If conversion fails, leave dead channel disabled
                    form['dead_channel']['enabled'].value = False
    
    def _on_configure_populations(self, btn):
        """Show population configuration UI"""
        # Organoid types
        org_label = widgets.HTML('<h4>Organoid Populations</h4>')
        self.n_organoid_input = widgets.IntText(
            value=self.n_organoid_types,
            description='Number of types:',
            style={'description_width': '120px'}
        )
        
        # Immune cell types
        immune_label = widgets.HTML('<h4>Immune Cell Populations</h4>')
        self.n_immune_input = widgets.IntText(
            value=self.n_immune_types,
            description='Number of types:',
            style={'description_width': '120px'}
        )
        
        # Other cell types  
        other_label = widgets.HTML('<h4>Other Cell Types</h4>')
        self.n_other_input = widgets.IntText(
            value=self.n_other_types,
            description='Number of types:',
            style={'description_width': '120px'}
        )
        
        # Include channels section
        channels_label = widgets.HTML('<h4>Channel Names</h4>')
        self.include_channels_checkbox = widgets.Checkbox(
            description='Include channels?',
            value=self.include_channels,
            style={'description_width': '120px'}
        )
        self.n_channels_input = widgets.IntText(
            value=self.n_channels if self.n_channels > 0 else 1,
            description='Number of channels:',
            style={'description_width': '140px'}
        )
        
        # Toggle visibility of channel count based on checkbox
        def _toggle_channels_input(change):
            self.n_channels_input.layout.display = None if change['new'] else 'none'
        self.include_channels_checkbox.observe(_toggle_channels_input, names='value')
        _toggle_channels_input({'new': self.include_channels_checkbox.value})
        
        btn_confirm = widgets.Button(description='Next: Name Cell Types', button_style='success')
        btn_confirm.on_click(self._show_cell_type_naming)
        
        self.main_container.children = [
            org_label,
            self.n_organoid_input,
            immune_label,
            self.n_immune_input,
            other_label,
            self.n_other_input,
            channels_label,
            self.include_channels_checkbox,
            self.n_channels_input,
            btn_confirm
        ]
    
    def _show_cell_type_naming(self, btn):
        """Show UI to name each cell type"""
        self.n_organoid_types = self.n_organoid_input.value
        self.n_immune_types = self.n_immune_input.value
        self.n_other_types = self.n_other_input.value
        
        # Capture channel configuration
        self.include_channels = self.include_channels_checkbox.value
        self.n_channels = self.n_channels_input.value if self.include_channels else 0
        
        naming_widgets = []
        self.organoid_name_inputs = []
        self.immune_name_inputs = []
        self.other_name_inputs = []
        
        # Organoid naming
        if self.n_organoid_types > 0:
            naming_widgets.append(widgets.HTML('<h4>Name Your Organoid Types</h4>'))
            for i in range(self.n_organoid_types):
                default = self.organoid_names[i] if i < len(self.organoid_names) else f'organoid{i+1}'
                w = widgets.Text(
                    value=default,
                    description=f'Organoid {i+1}:',
                    placeholder='e.g., organoid1, organoidWT',
                    style={'description_width': '100px'}
                )
                self.organoid_name_inputs.append(w)
                naming_widgets.append(w)
        
        # Immune cell naming
        if self.n_immune_types > 0:
            naming_widgets.append(widgets.HTML('<h4>Name Your Immune Cell Types</h4>'))
            for i in range(self.n_immune_types):
                default = self.immune_names[i] if i < len(self.immune_names) else f'immune{i+1}'
                w = widgets.Text(
                    value=default,
                    description=f'Immune {i+1}:',
                    placeholder='e.g., tcells, nk',
                    style={'description_width': '100px'}
                )
                self.immune_name_inputs.append(w)
                naming_widgets.append(w)
        
        # Other cell naming
        if self.n_other_types > 0:
            naming_widgets.append(widgets.HTML('<h4>Name Your Other Cell Types</h4>'))
            for i in range(self.n_other_types):
                default = self.other_names[i] if i < len(self.other_names) else f'other{i+1}'
                w = widgets.Text(
                    value=default,
                    description=f'Other {i+1}:',
                    placeholder='e.g., tumorcells, fibroblast',
                    style={'description_width': '100px'}
                )
                self.other_name_inputs.append(w)
                naming_widgets.append(w)
        
        btn_next = widgets.Button(description='Create Data Entry Form', button_style='success')
        btn_next.on_click(self._on_names_confirmed)
        naming_widgets.append(btn_next)
        
        self.main_container.children = naming_widgets
    
    def _on_names_confirmed(self, btn):
        """Collect names and build data entry form"""
        self.organoid_names = [w.value.strip() for w in self.organoid_name_inputs]
        self.immune_names = [w.value.strip() for w in self.immune_name_inputs]
        self.other_names = [w.value.strip() for w in self.other_name_inputs]
        
        self._build_data_entry_form()
    
    def _build_data_entry_form(self):
        """Build the main data entry form for all samples"""
        n_samples = self.n_samples_input.value
        self.sample_forms = []
        
        form_widgets = [widgets.HTML('<h3> Sample Data Entry</h3>')]
        
        # Fill-down button (NEW)
        btn_fill_down = widgets.Button(
            description='Fill All from Sample 1',
            button_style='warning',
            icon='arrow-down',
            tooltip='Copy all values from Sample 1 to all other samples (except sample names)'
        )
        btn_fill_down.on_click(self._on_fill_down)
        form_widgets.append(btn_fill_down)
        
        # Create form for each sample
        for i in range(n_samples):
            sample_form = self._create_sample_form(i)
            self.sample_forms.append(sample_form)
            form_widgets.append(sample_form['widget'])
        
        # Save button
        btn_save = widgets.Button(
            description='Save Metadata to CSV',
            button_style='success',
            icon='save'
        )
        btn_save.on_click(self._on_save_metadata)
        
        self.save_output = widgets.Output()
        
        form_widgets.extend([btn_save, self.save_output])
        self.main_container.children = form_widgets
        
        # If we have a loaded dataframe, re-populate the forms
        if self.loaded_df is not None:
            self._populate_form_fields(self.loaded_df)
    
    def _on_fill_down(self, btn):
        """Copy all values from Sample 1 to all other samples (except sample_name)"""
        if len(self.sample_forms) < 2:
            return
        
        source_form = self.sample_forms[0]
        
        for i in range(1, len(self.sample_forms)):
            target_form = self.sample_forms[i]
            
            # Copy channel fields
            for field_name, src_widget in source_form['channels'].items():
                if field_name in target_form['channels']:
                    target_form['channels'][field_name].value = src_widget.value
            
            # Copy basic fields (except sample_name)
            for field_name, src_widget in source_form['basic'].items():
                if field_name != 'sample_name':
                    target_form['basic'][field_name].value = src_widget.value
            
            # Copy cell type fields
            for cell_type in source_form['cell_types']:
                for field_name, src_widget in source_form['cell_types'][cell_type].items():
                    target_form['cell_types'][cell_type][field_name].value = src_widget.value
            
            # Copy dead channel (including mask_path)
            target_form['dead_channel']['enabled'].value = source_form['dead_channel']['enabled'].value
            if 'mask_path' in source_form['dead_channel'] and 'mask_path' in target_form['dead_channel']:
                target_form['dead_channel']['mask_path'].value = source_form['dead_channel']['mask_path'].value
            target_form['dead_channel']['number'].value = source_form['dead_channel']['number'].value
        
        with self.save_output:
            self.save_output.clear_output()
            print('✅ Filled all samples from Sample 1!')
    
    def _create_sample_form(self, sample_idx):
        """Create form widgets for a single sample"""
        form_data = {
            'basic': {},
            'cell_types': {},
            'dead_channel': {},
            'channels': {},  # For channel name fields
            'widget': None
        }
        
        # Sample header
        header = widgets.HTML(f'<h4 style="background:#e1f5fe;padding:8px;">Sample {sample_idx + 1}</h4>')
        
        # Basic fields
        form_data['basic']['sample_name'] = widgets.Text(
            description='Sample name*:',
            placeholder='e.g., Sample001',
            style={'description_width': '150px'},
            layout=widgets.Layout(width='400px')
        )
        
        form_data['basic']['exp_nr'] = widgets.IntText(
            description='Exp number*:',
            value=1,
            style={'description_width': '150px'}
        )
        
        form_data['basic']['well'] = widgets.Text(
            description='Well*:',
            placeholder='e.g., A1',
            style={'description_width': '150px'}
        )
        
        form_data['basic']['raw_image_path'] = widgets.Text(
            description='Raw image path*:',
            placeholder='/path/to/image.czi',
            style={'description_width': '150px'},
            layout=widgets.Layout(width='600px')
        )
        
        form_data['basic']['dimension_order'] = widgets.Text(
            description='Dimension order:',
            placeholder='Optional - e.g., TCZYX',
            style={'description_width': '150px'}
        )
        
        # Pixel/time metadata
        form_data['basic']['pixel_distance_xy'] = widgets.FloatText(
            description='Pixel xy (μm)*:',
            value=0.5,
            style={'description_width': '150px'}
        )
        
        form_data['basic']['pixel_distance_z'] = widgets.FloatText(
            description='Pixel z (μm)*:',
            value=2.0,
            style={'description_width': '150px'}
        )
        
        form_data['basic']['distance_unit'] = widgets.Text(
            description='Distance unit*:',
            value='μm',
            style={'description_width': '150px'}
        )
        
        form_data['basic']['time_interval'] = widgets.FloatText(
            description='Time interval*:',
            value=1.0,
            style={'description_width': '150px'}
        )
        
        form_data['basic']['time_unit'] = widgets.Text(
            description='Time unit*:',
            value='s',
            style={'description_width': '150px'}
        )
        
        # Dead channel
        dead_channel_label = widgets.HTML('<h5>Dead Channel</h5>')
        dead_enabled = widgets.Checkbox(
            description='Include dead channel',
            value=True
        )
        dead_mask_path = widgets.Text(
            description='Dead mask path:',
            placeholder='Optional - path to dead cell mask',
            style={'description_width': '130px'},
            layout=widgets.Layout(width='600px')
        )
        dead_number = widgets.IntText(
            description='Dead channel #:',
            value=0,
            style={'description_width': '120px'}
        )
        
        # Toggle visibility - show/hide both mask path and channel number
        def _toggle_dead(change):
            display_val = None if change['new'] else 'none'
            dead_mask_path.layout.display = display_val
            dead_number.layout.display = display_val
        dead_enabled.observe(_toggle_dead, names='value')
        _toggle_dead({'new': dead_enabled.value})
        
        form_data['dead_channel']['enabled'] = dead_enabled
        form_data['dead_channel']['mask_path'] = dead_mask_path
        form_data['dead_channel']['number'] = dead_number
        
        # Cell type sections
        cell_type_widgets = []
        
        # Organoids
        for org_name in self.organoid_names:
            cell_type_widgets.append(widgets.HTML(f'<h5>{org_name.capitalize()}</h5>'))
            fields = {}
            fields['line'] = widgets.Text(
                description='Line*:',
                placeholder='e.g., 32',
                style={'description_width': '100px'}
            )
            fields['condition'] = widgets.Text(
                description='Condition*:',
                placeholder='e.g., ko, wt',
                style={'description_width': '100px'}
            )
            # Optional path fields
            fields['segments_image_path'] = widgets.Text(
                description='Segments path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_image_path'] = widgets.Text(
                description='Tracks img path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_csv_path'] = widgets.Text(
                description='Tracks csv path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            form_data['cell_types'][org_name] = fields
            cell_type_widgets.extend([
                fields['line'], 
                fields['condition'],
                fields['segments_image_path'],
                fields['tracks_image_path'],
                fields['tracks_csv_path']
            ])
        
        # Immune cells
        for immune_name in self.immune_names:
            cell_type_widgets.append(widgets.HTML(f'<h5>{immune_name.capitalize()}</h5>'))
            fields = {}
            fields['line'] = widgets.Text(
                description='Line*:',
                placeholder='e.g., CD8',
                style={'description_width': '100px'}
            )
            fields['condition'] = widgets.Text(
                description='Condition*:',
                placeholder='e.g., activated',
                style={'description_width': '100px'}
            )
            # Optional path fields
            fields['segments_image_path'] = widgets.Text(
                description='Segments path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_image_path'] = widgets.Text(
                description='Tracks img path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_csv_path'] = widgets.Text(
                description='Tracks csv path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            form_data['cell_types'][immune_name] = fields
            cell_type_widgets.extend([
                fields['line'], 
                fields['condition'],
                fields['segments_image_path'],
                fields['tracks_image_path'],
                fields['tracks_csv_path']
            ])
        
        # Other cells
        for other_name in self.other_names:
            cell_type_widgets.append(widgets.HTML(f'<h5>{other_name.capitalize()}</h5>'))
            fields = {}
            fields['line'] = widgets.Text(
                description='Line*:',
                placeholder='e.g., WT',
                style={'description_width': '100px'}
            )
            fields['condition'] = widgets.Text(
                description='Condition*:',
                placeholder='e.g., control',
                style={'description_width': '100px'}
            )
            # Optional path fields
            fields['segments_image_path'] = widgets.Text(
                description='Segments path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_image_path'] = widgets.Text(
                description='Tracks img path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            fields['tracks_csv_path'] = widgets.Text(
                description='Tracks csv path:',
                placeholder='Optional - filled by pipeline',
                style={'description_width': '130px'},
                layout=widgets.Layout(width='500px')
            )
            form_data['cell_types'][other_name] = fields
            cell_type_widgets.extend([
                fields['line'], 
                fields['condition'],
                fields['segments_image_path'],
                fields['tracks_image_path'],
                fields['tracks_csv_path']
            ])
        
        # Channel label fields (if include_channels is enabled)
        channel_widgets = []
        if self.include_channels and self.n_channels > 0:
            for i in range(self.n_channels):
                channel_field = widgets.Text(
                    description=f'Channel {i+1} label:',
                    placeholder=f'e.g., Tcell, Organoid',
                    style={'description_width': '130px'},
                    layout=widgets.Layout(width='400px')
                )
                form_data['channels'][f'channel_{i+1}_label'] = channel_field
                channel_widgets.append(channel_field)
        
        # Assemble form - Basic Information includes channel labels
        form_children = [
            header,
            widgets.HTML('<h5>Basic Information</h5>'),
            form_data['basic']['sample_name'],
            form_data['basic']['exp_nr'],
            form_data['basic']['well'],
            form_data['basic']['raw_image_path'],
            form_data['basic']['dimension_order'],
        ]
        
        # Add channel label fields inside Basic Information section
        if channel_widgets:
            form_children.append(widgets.HTML('<b>Channel Labels:</b>'))
            form_children.extend(channel_widgets)
        
        # Continue with rest of the form
        form_children.extend([
            widgets.HTML('<h5>Imaging Parameters</h5>'),
            form_data['basic']['pixel_distance_xy'],
            form_data['basic']['pixel_distance_z'],
            form_data['basic']['distance_unit'],
            form_data['basic']['time_interval'],
            form_data['basic']['time_unit'],
            dead_channel_label,
            dead_enabled,
            dead_mask_path,
            dead_number,
            widgets.HTML('<h5>Cell Type Configuration</h5>'),
            *cell_type_widgets,
            widgets.HTML('<hr>')
        ])
        
        form_widget = widgets.VBox(form_children)
        
        form_data['widget'] = form_widget
        return form_data
    
    def _on_save_metadata(self, btn):
        """Save all form data to CSV"""
        if self.output_dir_picker is None:
            with self.save_output:
                self.save_output.clear_output()
                print('⚠️ No output directory picker provided')
            return
        
        output_dir = Path(self.output_dir_picker.value)
        if not output_dir.exists():
            with self.save_output:
                self.save_output.clear_output()
                print(f'⚠️ Output directory does not exist: {output_dir}')
            return
        
        # Collect data from all forms
        rows = []
        for form in self.sample_forms:
            row = {}
            
            # Channel names (if any)
            for field_name, widget in form['channels'].items():
                row[field_name] = widget.value.strip()
            
            # Basic fields (includes dimension_order now)
            for field_name, widget in form['basic'].items():
                row[field_name] = widget.value
            
            # Dead channel (including mask_path) - only add columns if enabled
            # If disabled, columns won't exist and pipeline interprets as "no death channel"
            if form['dead_channel']['enabled'].value:
                if 'mask_path' in form['dead_channel']:
                    row['dead_mask_path'] = form['dead_channel']['mask_path'].value.strip()
                row['dead_channel'] = form['dead_channel']['number'].value
            # Note: If disabled, we do NOT create these columns at all
            
            # Cell types with line_condition merging and prefixes
            for cell_type, fields in form['cell_types'].items():
                # Determine prefix based on category
                if cell_type in self.organoid_names:
                    prefix = 'or'
                elif cell_type in self.immune_names:
                    prefix = 'im'
                elif cell_type in self.other_names:
                    prefix = 'ot'
                else:
                    continue  # Skip unknown types
                
                line = fields['line'].value.strip()
                condition = fields['condition'].value.strip()
                
                # Merge line_condition as "{prefix}_{celltype}_line_condition"
                col_name = f'{prefix}_{cell_type}_line_condition'
                if line and condition:
                    row[col_name] = f'{line}_{condition}'
                elif line:
                    row[col_name] = line
                elif condition:
                    row[col_name] = f'_{condition}'
                else:
                    row[col_name] = ''
                
                # Optional path fields - save if filled, otherwise empty
                row[f'{prefix}_{cell_type}_segments_image_path'] = fields.get('segments_image_path', widgets.Text()).value.strip()
                row[f'{prefix}_{cell_type}_tracks_image_path'] = fields.get('tracks_image_path', widgets.Text()).value.strip()
                row[f'{prefix}_{cell_type}_tracks_csv_path'] = fields.get('tracks_csv_path', widgets.Text()).value.strip()
            
            rows.append(row)
        
        # Create DataFrame
        df = pd.DataFrame(rows)
        
        # Save to CSV
        csv_path = output_dir / 'metadata.csv'
        df.to_csv(csv_path, index=False)
        
        with self.save_output:
            self.save_output.clear_output()
            print(f'Metadata saved to: {csv_path}')
            print(f'{len(df)} samples, {len(df.columns)} columns')
            display(df)

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
                traceback.print_exc()
            finally:
                metadata_loader.metadata = result
                try:
                    metadata_loader.metadata.to_csv(metadata_loader.metadata_csv_path, index=False)
                except Exception:
                    traceback.print_exc()
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
                    traceback.print_exc()

                print("✅ Unmixing finished.", flush=True)
            except Exception:
                traceback.print_exc()
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
                        traceback.print_exc()

                print(f"✅ Saved updated metadata")
            except Exception as e:
                print(f"Error saving metadata: {e}")

class PixelClassifierPanel:
    def __init__(self, metadata_loader):
        self.metadata_loader = metadata_loader
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        
        # Detect cell types from metadata
        self._detect_cell_types()
        
        print(f"Detected organoid types: {self.organoid_types}")
        print(f"Detected immune cell types: {self.immune_types}")
        print(f"Detected other cell types: {self.other_types}")
        print(f"Dead channel present: {self.has_death}")

        # -------- Train controls --------
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

        # ---- Manual classifier toggle ----
        self.manual_clf_paths = widgets.Checkbox(
            description="Manually supply classifiers",
            value=bool(pc.get("manual_clf_paths", False)),
            indent=False
        )

        # NEW: Overwrite existing outputs toggle
        self.overwrite_existing = widgets.Checkbox(
            description="Overwrite existing",
            value=bool(pc.get("overwrite_existing", False))
        )

        # -------- Dynamic classifier path pickers --------
        self.clf_dir = PathPicker(
            mode='dir',
            description='Classifier dir:',
            default="",
            description_width='160px',
            width='100%',
        )
        
        # Create path pickers for each cell type dynamically
        self.clf_paths = {}
        self._create_clf_path_pickers()
        
        # Death mask classifier (conditionally displayed)
        self.clf_death_path = PathPicker(
            mode='file',
            description='Death clf:',
            default="",
            filter_pattern='*.joblib',
            description_width='160px',
            width='100%',
        )
        
        # Restore saved paths if manual mode is enabled
        if self.manual_clf_paths.value:
            self.clf_dir.value = str(pc.get("clf_dir", "") or "")
            self.clf_death_path.value = str(pc.get("clf_death_path", "") or "")
            # Restore dynamic cell type paths
            for cell_type in self.all_cell_types:
                if cell_type in self.clf_paths:
                    saved_path = pc.get(f"clf_{cell_type}_path", "")
                    if saved_path:
                        self.clf_paths[cell_type].value = str(saved_path)

        # -------- Dynamic EDT thresholds for each cell type --------
        self.edt_thresholds = {}
        self._create_edt_threshold_inputs()
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

        # Two Apply buttons + spinner
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
        
        self.btn_apply = self.btn_run
        
        self.spinner_apply = widgets.HTML(value=spinning_loader)
        self.spinner_apply.layout.display = "none"
        self.apply_row = widgets.HBox(
            [self.btn_run, self.btn_resegment, self.spinner_apply],
            layout=widgets.Layout(align_items="center", gap="8px")
        )

        # Wire handlers
        self.btn_train.on_click(self._on_train_clicked)
        self.close_button.on_click(self._on_close_clicked)
        self.btn_run.on_click(partial(self._on_apply_clicked, only_segment=False))
        self.btn_resegment.on_click(partial(self._on_apply_clicked, only_segment=True))

        self.out = widgets.Output()

        # Build the classifier paths box
        self.clf_paths_box = widgets.VBox()
        self._build_clf_paths_box()
        
        # Wire observers
        self.manual_clf_paths.observe(self._toggle_clf_path_section, names='value')
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
            self.train_row,
        ])

        self.tp_row = widgets.HBox([self.use_all_timepoints, self.tp_start, self.tp_end])
        
        # Build dynamic threshold widgets layout
        threshold_widgets = [widgets.HTML("<b>EDT Thresholds per cell type</b>")]
        
        if self.organoid_types:
            threshold_widgets.append(widgets.HTML("<i>Organoids:</i>"))
            for org_type in self.organoid_types:
                if org_type in self.edt_thresholds:
                    threshold_widgets.append(self.edt_thresholds[org_type])
        
        if self.immune_types:
            threshold_widgets.append(widgets.HTML("<i>Immune cells:</i>"))
            for immune_type in self.immune_types:
                if immune_type in self.edt_thresholds:
                    threshold_widgets.append(self.edt_thresholds[immune_type])
        
        if self.other_types:
            threshold_widgets.append(widgets.HTML("<i>Other cells:</i>"))
            for other_type in self.other_types:
                if other_type in self.edt_thresholds:
                    threshold_widgets.append(self.edt_thresholds[other_type])

        apply_box = widgets.VBox([
            widgets.HTML("<b>Apply segmentation</b>"),
            self.manual_clf_paths,
            self.clf_paths_box,
            *threshold_widgets,
            self.overwrite_existing,
            self.tp_row,
            self.apply_row,
        ])

        self.ui = widgets.VBox([train_box, widgets.HTML("<hr>"), apply_box, widgets.HTML("<hr>"), self.out])

        # viewer handle (if training opens napari and returns a viewer)
        self._viewer = None

    # ---------- Cell type detection ----------
    def _detect_cell_types(self):
        """Detect all cell types from metadata"""
        from behav3d.utils import (
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
        
        # Combined list for iteration
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types
    
    def _create_clf_path_pickers(self):
        """Create path pickers for each detected cell type"""
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
        """Create EDT threshold inputs for each cell type with appropriate defaults"""
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
        
        for cell_type in self.all_cell_types:
            # Determine default threshold based on cell type category
            if cell_type in self.organoid_types:
                default_threshold = 12.0  # Large objects
            elif cell_type in self.immune_types:
                default_threshold = 2.5   # Small objects
            else:  # other types
                default_threshold = 1.0   # User must adjust
            
            # Check if saved in config
            saved_threshold = pc.get(f"{cell_type}_edt_threshold", default_threshold)
            
            self.edt_thresholds[cell_type] = widgets.FloatText(
                description=f"{cell_type.capitalize()} EDT:",
                value=float(saved_threshold),
                style={'description_width': '160px'},
                tooltip=f"EDT threshold for {cell_type} segmentation"
            )
    
   #def _update_immune_cell_paths(self):
        #"""Update classifier paths when directory changes (for backward compatibility)"""
        # This method helps auto-populate paths based on directory
        #pass
    
    # ---------- config ----------
    def _persist_params(self):
        self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc = self.metadata_loader.behav3d_parameters["pixel_classifier"]
        pc["examples_per_sample"] = int(self.examples_per_sample.value)
        pc["sample_specific_classifier"] = bool(self.sample_specific_classifier.value)
        pc["workers"] = int(self.n_workers.value)
        pc["use_all_timepoints"] = bool(self.use_all_timepoints.value)
        pc["tp_start"] = int(self.tp_start.value)
        pc["tp_end"]   = int(self.tp_end.value)
        pc["manual_clf_paths"] = bool(self.manual_clf_paths.value)
        pc["overwrite_existing"] = bool(self.overwrite_existing.value)
        
        # Save dynamic EDT thresholds
        for cell_type, threshold_widget in self.edt_thresholds.items():
            pc[f"{cell_type}_edt_threshold"] = float(threshold_widget.value)
        
        # Save dynamic classifier paths if manual mode
        if self.manual_clf_paths.value:
            pc["clf_dir"] = str(self.clf_dir.value or "")
            pc["clf_death_path"] = str(self.clf_death_path.value or "")
            for cell_type, path_picker in self.clf_paths.items():
                pc[f"clf_{cell_type}_path"] = str(path_picker.value or "")

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

    # >>> helper: compute default file paths for a directory (dynamic)
    def _default_clf_paths_for_dir(self, d: str):
        """Generate default classifier paths for all cell types"""
        if not d:
            return {}
        
        p = Path(d).expanduser()
        paths = {}
        
        # Generate path for each cell type
        for cell_type in self.all_cell_types:
            # Capitalize first letter for filename
            filename = f'PixelClassifier_{cell_type.capitalize()}.joblib'
            paths[cell_type] = str(p / filename)
        
        # Death classifier
        paths['death'] = str(p / 'PixelClassifier_Death.joblib')
        
        return paths

    def _apply_dir_to_clf_paths(self, d: str):
        """Apply directory to all classifier path pickers"""
        paths = self._default_clf_paths_for_dir(d)
        
        for cell_type, path in paths.items():
            if cell_type == 'death':
                self.clf_death_path.value = path
            elif cell_type in self.clf_paths:
                self.clf_paths[cell_type].value = path

    def _toggle_clf_path_section(self, change=None):
        manual = bool(self.manual_clf_paths.value)
        self.clf_paths_box.layout.display = (None if manual else 'none')
        if not manual:
            # Clear all paths
            self.clf_dir.value = ""
            self.clf_death_path.value = ""
            for picker in self.clf_paths.values():
                picker.value = ""

    def _build_clf_paths_box(self):
        """Build classifier paths box with dynamic cell types"""
        children = [
            widgets.HTML("<b>Classifier paths</b>"),
            self.clf_dir,
        ]
        
        # Add organoid classifiers
        if self.organoid_types:
            children.append(widgets.HTML("<i>Organoids:</i>"))
            for org_type in self.organoid_types:
                if org_type in self.clf_paths:
                    children.append(self.clf_paths[org_type])
        
        # Add immune cell classifiers
        if self.immune_types:
            children.append(widgets.HTML("<i>Immune cells:</i>"))
            for immune_type in self.immune_types:
                if immune_type in self.clf_paths:
                    children.append(self.clf_paths[immune_type])
        
        # Add other cell classifiers
        if self.other_types:
            children.append(widgets.HTML("<i>Other cells:</i>"))
            for other_type in self.other_types:
                if other_type in self.clf_paths:
                    children.append(self.clf_paths[other_type])
        
        # Death classifier - only if dead channel is present
        if self.has_death:
            children.append(widgets.HTML("<i>Death mask:</i>"))
            children.append(self.clf_death_path)
        
        self.clf_paths_box.children = children

    def display(self):
        display(self.ui)

    def _lock(self, state: bool):
        # keep close_button enabled so user can close at any time
        widgets_to_lock = [
            self.btn_train, self.btn_run, self.btn_resegment,
            self.examples_per_sample, self.sample_specific_classifier, self.n_workers,
            self.use_all_timepoints, self.tp_start, self.tp_end,
            self.manual_clf_paths, self.overwrite_existing,
        ]
        
        # Add dynamic threshold widgets
        widgets_to_lock.extend(self.edt_thresholds.values())
        
        for w in widgets_to_lock:
            w.disabled = state

        # Lock/unlock path pickers
        try:
            self.clf_dir.text.disabled = state
            self.clf_dir.button.disabled = state
            
            # Death path - only if present
            if self.has_death:
                self.clf_death_path.text.disabled = state
                self.clf_death_path.button.disabled = state
            
            # Dynamic cell type pickers
            for picker in self.clf_paths.values():
                picker.text.disabled = state
                picker.button.disabled = state
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
                print(f"  organoid_types={self.organoid_types}")
                print(f"  immune_types={self.immune_types}")
                print(f"  other_types={self.other_types}")
                print(f"  n_workers={self.n_workers.value}")

                # UI state
                self._lock(True)
                self.spinner_train.layout.display = None

                # Import where available; allow fallback to global
                try:
                    from behav3d.preprocessing.segmentation.napari_pixelclassifier import train_pixel_classifier
                except Exception:
                    pass

                ret = train_pixel_classifier(
                    output_dir=str(odir),
                    metadata=self.metadata_loader.metadata,
                    examples_per_sample=int(self.examples_per_sample.value),
                    sample_specific_classifier=bool(self.sample_specific_classifier.value),
                    n_workers=int(self.n_workers.value),
                    organoid_types=self.organoid_types,
                    immune_types=self.immune_types,
                    other_types=self.other_types,
                )

                # Try to capture a viewer handle
                self._viewer = None
                try:
                    if ret is not None and hasattr(ret, "close"):
                        self._viewer = ret
                    else:
                        get_curr = getattr(napari, "current_viewer", None)
                        if callable(get_curr):
                            self._viewer = get_curr()
                except Exception:
                    pass

                # UI finalize for train
                self.spinner_train.layout.display = "none"
                self.close_button.layout.display = "inline-block" if self._viewer is not None else "none"
                if self._viewer is None:
                    print("✅ Training finished.", flush=True)
                else:
                    print("✅ Training UI opened in Napari (use 'Close viewer' to close).", flush=True)

            except Exception:
                traceback.print_exc()
                self.spinner_train.layout.display = "none"
                self.close_button.layout.display = "none"
            finally:
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

    def _on_apply_clicked(self, _, only_segment=False):
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

                # Build EDT threshold dictionaries
                organoid_edt_thresholds = {
                    cell_type: float(self.edt_thresholds[cell_type].value)
                    for cell_type in self.organoid_types
                    if cell_type in self.edt_thresholds
                }
                immune_edt_thresholds = {
                    cell_type: float(self.edt_thresholds[cell_type].value)
                    for cell_type in self.immune_types
                    if cell_type in self.edt_thresholds
                }
                other_edt_thresholds = {
                    cell_type: float(self.edt_thresholds[cell_type].value)
                    for cell_type in self.other_types
                    if cell_type in self.edt_thresholds
                }

                print("▶️ Applying pixel classifier segmentation…", flush=True)
                print(f"  output_dir={odir}")
                print(f"  organoid_edt_thresholds={organoid_edt_thresholds}")
                print(f"  immune_edt_thresholds={immune_edt_thresholds}")
                print(f"  other_edt_thresholds={other_edt_thresholds}")
                print(f"  timepoint_range={tpr}", flush=True)
                print(f"  only_segment={only_segment}")
                print(f"  overwrite_existing={bool(self.overwrite_existing.value)}")
                
                # Build classifier path dictionaries if manual mode
                clf_organoid_paths = None
                clf_immune_paths = None
                clf_other_paths = None
                clf_death_path = None
                
                if self.manual_clf_paths.value:
                    clf_organoid_paths = {
                        cell_type: str(self.clf_paths[cell_type].value)
                        for cell_type in self.organoid_types
                        if cell_type in self.clf_paths and self.clf_paths[cell_type].value
                    }
                    clf_immune_paths = {
                        cell_type: str(self.clf_paths[cell_type].value)
                        for cell_type in self.immune_types
                        if cell_type in self.clf_paths and self.clf_paths[cell_type].value
                    }
                    clf_other_paths = {
                        cell_type: str(self.clf_paths[cell_type].value)
                        for cell_type in self.other_types
                        if cell_type in self.clf_paths and self.clf_paths[cell_type].value
                    }
                    # Death path - only if death channel is present
                    if self.has_death:
                        clf_death_path = str(self.clf_death_path.value) if self.clf_death_path.value else None
                    
                    print(f"  Manual classifier paths:")
                    print(f"    Organoids: {clf_organoid_paths}")
                    print(f"    Immune: {clf_immune_paths}")
                    print(f"    Other: {clf_other_paths}")
                    print(f"    Death: {clf_death_path}")

                self.spinner_apply.layout.display = None

                # Build kwargs for new API (with dictionaries)
                call_kwargs = dict(
                    output_dir=str(odir),
                    metadata=self.metadata_loader.metadata,
                    organoid_edt_thresholds=organoid_edt_thresholds,
                    immune_edt_thresholds=immune_edt_thresholds,
                    other_edt_thresholds=other_edt_thresholds,
                    timepoint_range=tpr,
                    clf_organoid_paths=clf_organoid_paths,
                    clf_immune_paths=clf_immune_paths,
                    clf_other_paths=clf_other_paths,
                    clf_death_path=clf_death_path,
                    only_segment=bool(only_segment),
                    overwrite_existing=bool(self.overwrite_existing.value),
                    n_workers=int(self.n_workers.value),
                )

                # Call the segmentation function
                try:
                    new_md = run_pixel_classifier_segmentation(**call_kwargs)
                except TypeError as e:
                    # Try backward compatibility fallback (if needed)
                    print(f"⚠️ TypeError: {e}")
                    print("Attempting backward compatibility mode...")
                    call_kwargs.pop("overwrite_existing", None)
                    new_md = run_pixel_classifier_segmentation(**call_kwargs)

                try:
                    if new_md is not None:
                        self.metadata_loader.metadata = new_md
                        new_md.to_csv(self.metadata_loader.metadata_csv_path, index=False)
                except Exception:
                    traceback.print_exc()

                print("✅ Apply finished.", flush=True)
            except Exception:
                traceback.print_exc()
            finally:
                self.spinner_apply.layout.display = "none"
                self._lock(False)
                


class TrackingPanel:
    """
    Minimal UI for tracking (LAP, TrackPy, or Propagation) for any cell type.
    - Automatically detects cell type category (organoid/immune/other)
    - Loads appropriate default configuration based on category
    - Uses metadata_loader.output_dir and metadata_loader.metadata_csv_path directly
    - Updates metadata_loader.metadata and always saves CSV
    - cell_type is supplied in __init__
    """
    def __init__(self, metadata_loader, cell_type):
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()

        # Detect ALL cell types from metadata
        self._detect_all_cell_types()
        
        # Detect category for THIS specific cell_type
        self._detect_category()

        # Ensure config skeleton exists, then read profile
        self._ensure_cfg_skeleton()
        params = dict(self.metadata_loader.behav3d_parameters or {})
        
        # Load config based on category
        tcfg = self._get_config_for_cell_type(params)

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

    def _detect_all_cell_types(self):
        """Detect all cell types from metadata"""
        from behav3d.utils import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata
        )
        
        metadata = self.metadata_loader.metadata
        if metadata is None:
            raise ValueError(
                f"Metadata must be loaded before creating TrackingPanel for '{self.cell_type}'. "
                "Please load metadata first using MetadataLoader."
            )
        
        self.organoid_types = detect_organoid_types_from_metadata(metadata)
        self.immune_types = detect_immune_cell_types_from_metadata(metadata)
        self.other_types = detect_other_cell_types_from_metadata(metadata)
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types

    def _detect_category(self):
        """Detect if THIS cell_type is organoid, immune, or other category."""
        # Use unified category detection function
        self.category = detect_cell_type_category(self.cell_type, self.metadata_loader.metadata)
    
    def _get_config_for_cell_type(self, params):
        """Get tracking config based on cell type category.
        
        - Organoid types: use 'organoid' config
        - Immune types: use 'immune' config
        - Other types: use 'other' config
        """
        tracking_cfg = params.get("tracking", {})
        
        # Check if cell_type has specific config
        if self.cell_type in tracking_cfg:
            return tracking_cfg[self.cell_type]
        
        # Otherwise use category-specific defaults
        if self.category == 'organoid':
            return tracking_cfg.get('organoid', {})
        elif self.category == 'immune':
            return tracking_cfg.get('immune', {})
        else:  # other
            return tracking_cfg.get('other', {})

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
        """Ensure config skeleton exists with category-appropriate defaults."""
        p = dict(self.metadata_loader.behav3d_parameters or {})
        p.setdefault("tracking", {})
        
        # Ensure base configs exist for all three categories
        p["tracking"].setdefault("immune", {})
        p["tracking"].setdefault("other", {})
        p["tracking"].setdefault("organoid", {})
        
        # Set defaults for immune cell types
        immune_cfg = p["tracking"]["immune"]
        immune_cfg.setdefault("method", "trackpy")
        immune_cfg.setdefault("overwrite", False)
        immune_cfg.setdefault("lap", {
            "track_cost_px": 60,
            "gap_close_cost_px": 45,
            "gap_close_max_frames": 5,
            "merging_cost_px": 0,
            "splitting_cost_px": 0
        })
        immune_cfg.setdefault("trackpy", {
            "search_range_px": 31,
            "memory_frames": 5,
            "adaptive_stop": 5.0,
            "adaptive_step": 0.95
        })
        
        # Set defaults for other cell types (same as immune for now)
        other_cfg = p["tracking"]["other"]
        other_cfg.setdefault("method", "trackpy")
        other_cfg.setdefault("overwrite", False)
        other_cfg.setdefault("lap", {
            "track_cost_px": 60,
            "gap_close_cost_px": 45,
            "gap_close_max_frames": 5,
            "merging_cost_px": 0,
            "splitting_cost_px": 0
        })
        other_cfg.setdefault("trackpy", {
            "search_range_px": 31,
            "memory_frames": 5,
            "adaptive_stop": 5.0,
            "adaptive_step": 0.95
        })
        
        # Set defaults for organoid types
        organoid_cfg = p["tracking"]["organoid"]
        organoid_cfg.setdefault("method", "propagation")
        organoid_cfg.setdefault("overwrite", False)
        organoid_cfg.setdefault("lap", {
            "track_cost_px": 60,
            "gap_close_cost_px": 80,
            "gap_close_max_frames": 3,
            "merging_cost_px": 0,
            "splitting_cost_px": 0
        })
        organoid_cfg.setdefault("trackpy", {
            "search_range_px": 35,
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

                    new_md = run_laptracking(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        cell_type=self.cell_type,
                        track_cost_cutoff=tc,
                        gap_closing_cost_cutoff=gc,
                        gap_closing_max_frame_count=int(self.gap_max_frames.value),
                        merging_cost_cutoff=merging,
                        splitting_cost_cutoff=splitting,
                        overwrite=bool(self.overwrite.value),
                    )

                elif method == "trackpy":
                    print("▶️ TrackPy tracking…", flush=True)
                    print(f"  search_range={self.tp_search_range.value}, memory={self.tp_memory.value}")
                    print(f"  adaptive_stop={self.tp_adaptive_stop.value}, adaptive_step={self.tp_adaptive_step.value}", flush=True)

                    new_md = run_trackpy_tracking_generic(
                        metadata=self.metadata_loader.metadata,
                        output_dir=str(out_dir),
                        cell_type=self.cell_type,
                        overwrite=bool(self.overwrite.value),
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
                traceback.print_exc()
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
        # Check if any sample has a dead_channel defined in metadata
        has_dead_channel = False
        if hasattr(metadata_loader, 'metadata') and metadata_loader.metadata is not None:
            has_dead_channel = 'dead_channel' in metadata_loader.metadata.columns and \
                              metadata_loader.metadata['dead_channel'].notna().any()
        
        self.dead_mask_threshold = widgets.FloatText(
            description="Dead mask % thr",
            value=float(fcfg.get("dead_mask_percentage_threshold", 0.05)),
            style={'description_width':'160px'}
        )
        # Hide dead mask threshold if no dead_channel in metadata
        if not has_dead_channel:
            self.dead_mask_threshold.layout.display = "none"

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
            elif f == "death" and not has_dead_channel:
                # Skip death feature if no dead_channel in metadata
                continue
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
                contact_thr = float(self.contact_threshold.value)
                ow = bool(self.overwrite.value)

                print("▶️ Running feature extraction…", flush=True)
                print(f"  dead_mask_percentage_threshold={thr}")
                print(f"  features_choice={feats}")
                print(f"  n_workers={workers}, overwrite={ow}")
                print(f"  contact_threshold={contact_thr}")
                print(f"  output_dir={out_dir}", flush=True)

                new_md = run_feature_extraction(
                    dead_mask_percentage_threshold=thr,
                    contact_threshold=contact_thr,
                    metadata=self.metadata_loader.metadata,
                    output_dir=str(out_dir),
                    features_choice=feats,
                    cell_type=self.cell_type,
                    n_workers=workers,
                    overwrite=ow                # <- pass overwrite here
                )
                print("✅ Feature extraction finished.", flush=True)
            except Exception:
                traceback.print_exc()
            finally:


                self.spinner_html.layout.display = "none"
                self._lock(False)


class ActiveKillingPanel:
    """
    Advanced feature extraction panel for Active Killing Analysis.
    
    Detects functional immune cell killing events by analyzing death signal 
    changes after cell-cell contact. Analyzes killing against ALL organoid
    types automatically detected from metadata.
    """
    
    def __init__(self, metadata_loader):
        """
        Parameters
        ----------
        metadata_loader : MetadataLoader
            Metadata loader instance with loaded metadata
        """
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        # Detect all cell types from metadata
        from behav3d.utils import (
            detect_immune_cell_types_from_metadata,
            detect_organoid_types_from_metadata,
            detect_other_cell_types_from_metadata
        )
        
        md = self.metadata_loader.metadata
        if md is None:
            raise RuntimeError("metadata_loader.metadata must be loaded before creating ActiveKillingPanel.")
        
        self.immune_types = detect_immune_cell_types_from_metadata(md)
        self.organoid_types = detect_organoid_types_from_metadata(md)
        self.other_types = detect_other_cell_types_from_metadata(md)
        
        # Potential immune cells (attackers): immune + other
        self.potential_immune = self.immune_types + self.other_types
        # All target types (organoids) will be analyzed automatically
        self.target_types = self.organoid_types
        
        # Load config
        params = dict(self.metadata_loader.behav3d_parameters or {})
        params.setdefault("active_killing", deepcopy(_DEFAULT_CONFIG.get("active_killing", {})))
        self._params = params
        self._cfg = self._params["active_killing"]
        
        # ---- Section Title ----
        self.section_title = widgets.HTML(
            '<div style="font-size:22px;font-weight:700;">Active Killing Analysis</div>'
        )
        
        self.description = widgets.HTML(
            '<div style="color:#555;font-size:13px;margin-bottom:10px;">'
            'Detects functional killing events by analyzing death signal changes after immune-target contact.<br>'
            '<b>Targets:</b> All organoid types will be analyzed automatically.'
            '</div>'
        )
        
        # ---- Cell Type Selection ----
        # Immune cell dropdown
        immune_options = self.potential_immune if self.potential_immune else ["(none detected)"]
        self.immune_dd = widgets.Dropdown(
            options=immune_options,
            value=immune_options[0] if immune_options else None,
            description="Immune cell:",
            style={'description_width': '120px'},
            layout=widgets.Layout(width="280px")
        )
        
        # Show detected target types (read-only info)
        target_info = ", ".join(self.target_types) if self.target_types else "(none detected)"
        self.target_info_html = widgets.HTML(
            f'<div style="padding:5px;background:#f0f0f0;border-radius:4px;">'
            f'<b>Target cell types:</b> {target_info}</div>'
        )
        
        self.cell_selection_row = widgets.VBox([
            self.immune_dd,
            self.target_info_html
        ], layout=widgets.Layout(gap="10px"))
        
        # ---- Parameters ----
        self.observation_window = widgets.IntText(
            description="Observation window:",
            value=int(self._cfg.get("observation_window", 5)),
            style={'description_width': '150px'},
            layout=widgets.Layout(width="220px")
        )
        self.observation_window_label = widgets.HTML(
            '<span style="color:#666;font-size:12px;">timepoints after contact</span>'
        )
        
        # Death signal column dropdown
        death_signal_options = ["mean_dead_dye", "percentage_dead_mask", "nr_dead_mask_pixels"]
        self.death_signal_dd = widgets.Dropdown(
            options=death_signal_options,
            value=self._cfg.get("death_signal_column", "mean_dead_dye"),
            description="Death signal:",
            style={'description_width': '150px'},
            layout=widgets.Layout(width="300px")
        )
        
        self.killing_threshold = widgets.FloatText(
            description="Killing threshold:",
            value=float(self._cfg.get("killing_threshold_multiplier", 1.5)),
            style={'description_width': '150px'},
            layout=widgets.Layout(width="220px")
        )
        self.killing_threshold_label = widgets.HTML(
            '<span style="color:#666;font-size:12px;">× background rate</span>'
        )
        
        self.min_contact_duration = widgets.IntText(
            description="Min contact duration:",
            value=int(self._cfg.get("min_contact_duration", 1)),
            style={'description_width': '150px'},
            layout=widgets.Layout(width="220px")
        )
        self.min_contact_duration_label = widgets.HTML(
            '<span style="color:#666;font-size:12px;">timepoints</span>'
        )
        
        self.save_results = widgets.Checkbox(
            description="Save results to CSV",
            value=bool(self._cfg.get("save_results", True)),
            indent=False
        )
        
        # Parameter rows
        self.param_row1 = widgets.HBox([
            self.observation_window, self.observation_window_label,
            widgets.HTML("&nbsp;&nbsp;&nbsp;"),
            self.killing_threshold, self.killing_threshold_label
        ], layout=widgets.Layout(align_items="center"))
        
        self.param_row2 = widgets.HBox([
            self.min_contact_duration, self.min_contact_duration_label,
            widgets.HTML("&nbsp;&nbsp;&nbsp;"),
            self.death_signal_dd
        ], layout=widgets.Layout(align_items="center"))
        
        # ---- Run Button ----
        self.btn_run = widgets.Button(
            description="Run Active Killing Analysis",
            button_style="danger",
            icon="bolt",
            layout=widgets.Layout(width="260px")
        )
        self.btn_run.on_click(self._on_run_clicked)
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html, self.save_results],
            layout=widgets.Layout(align_items="center", gap="15px")
        )
        
        # ---- Output ----
        self.out = widgets.Output()
        
        # ---- Validation message ----
        self.validation_html = widgets.HTML("")
        self._validate_inputs()
        
        # Observe changes to update validation
        self.immune_dd.observe(lambda _: self._validate_inputs(), names="value")
        
        # ---- Build UI ----
        self.ui = widgets.VBox([
            self.section_title,
            self.description,
            widgets.HTML("<b>Cell Type Selection</b>"),
            self.cell_selection_row,
            self.validation_html,
            widgets.HTML("<hr>"),
            widgets.HTML("<b>Analysis Parameters</b>"),
            self.param_row1,
            self.param_row2,
            widgets.HTML("<hr>"),
            self.run_row,
            self.out
        ])
    
    def _validate_inputs(self):
        """Validate that required track files exist"""
        immune = self.immune_dd.value
        
        messages = []
        valid = True
        
        if immune == "(none detected)":
            messages.append("⚠️ No immune cell types detected in metadata")
            valid = False
        
        if not self.target_types:
            messages.append("⚠️ No organoid types detected in metadata")
            valid = False
        
        if valid:
            # Check if immune track feature files exist
            immune_path = Path(self.output_dir, "analysis", immune, "track_features",
                              f"BEHAV3D_{immune}_combined_track_features_filtered.csv")
            
            if not immune_path.exists():
                alt_path = immune_path.with_name(f"BEHAV3D_{immune}_combined_track_features.csv")
                if not alt_path.exists():
                    messages.append(f"⚠️ {immune} track features not found. Run feature extraction first.")
                    valid = False
                else:
                    immune_path = alt_path
            
            # Check at least one target has tracks and contact columns exist
            targets_found = []
            for target in self.target_types:
                target_path = Path(self.output_dir, "analysis", target, "track_features",
                                  f"BEHAV3D_{target}_combined_track_features_filtered.csv")
                if not target_path.exists():
                    target_path = target_path.with_name(f"BEHAV3D_{target}_combined_track_features.csv")
                
                if target_path.exists():
                    targets_found.append(target)
            
            if not targets_found:
                messages.append(f"⚠️ No target track features found. Run feature extraction first.")
                valid = False
            
            # Check contact columns exist for at least one target
            if valid and immune_path.exists():
                try:
                    df_sample = pd.read_csv(immune_path, nrows=1)
                    contacts_found = []
                    for target in targets_found:
                        contact_col = f"{target}_contact"
                        if contact_col in df_sample.columns:
                            contacts_found.append(target)
                    
                    if not contacts_found:
                        messages.append(f"⚠️ No contact columns found in {immune} tracks. "
                                       f"Ensure contact features were calculated.")
                        valid = False
                except Exception:
                    pass
        
        if valid:
            self.validation_html.value = '<span style="color:green;">✓ Ready to run</span>'
        else:
            self.validation_html.value = '<br>'.join([f'<span style="color:#c00;">{m}</span>' for m in messages])
        
        self.btn_run.disabled = not valid
        return valid
    
    def _persist_params(self):
        """Save parameters to config file"""
        self._cfg["observation_window"] = int(self.observation_window.value)
        self._cfg["death_signal_column"] = str(self.death_signal_dd.value)
        self._cfg["killing_threshold_multiplier"] = float(self.killing_threshold.value)
        self._cfg["min_contact_duration"] = int(self.min_contact_duration.value)
        self._cfg["save_results"] = bool(self.save_results.value)
        self._cfg["last_immune_cell"] = str(self.immune_dd.value)
        # Target types are now auto-detected, no need to save
        
        self._params["active_killing"] = self._cfg
        with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
            yaml.safe_dump(self._params, f, sort_keys=False)
    
    def _lock(self, locked):
        """Lock/unlock controls"""
        for w in [self.immune_dd, self.observation_window,
                  self.death_signal_dd, self.killing_threshold, 
                  self.min_contact_duration, self.save_results, self.btn_run]:
            w.disabled = locked
    
    def _on_run_clicked(self, *_):
        """Run active killing analysis for all organoid types"""
        self._lock(True)
        self.spinner_html.layout.display = None
        self.out.clear_output()
        
        with self.out:
            try:
                self._persist_params()
                
                immune_cell = str(self.immune_dd.value)
                observation_window = int(self.observation_window.value)
                death_signal = str(self.death_signal_dd.value)
                threshold_mult = float(self.killing_threshold.value)
                min_contact = int(self.min_contact_duration.value)
                save = bool(self.save_results.value)
                
                print(f"▶️ Running Active Killing Analysis...")
                print(f"  Immune cell type: {immune_cell}")
                print(f"  Target cell types: {', '.join(self.target_types)}")
                print(f"  Observation window: {observation_window} timepoints")
                print(f"  Death signal column: {death_signal}")
                print(f"  Killing threshold: {threshold_mult}× background")
                print(f"  Min contact duration: {min_contact} timepoints")
                print(f"  Save results: {save}")
                print()
                
                # Import and run the analysis
                from behav3d.analysis.advanced_feature_extraction import run_active_killing_analysis
                
                # Run analysis for all target types (None triggers auto-detection)
                df_killing_events, df_summary, stats = run_active_killing_analysis(
                    metadata=self.metadata_loader.metadata,
                    output_dir=self.output_dir,
                    immune_cell_type=immune_cell,
                    target_cell_types=None,  # Auto-detect all organoid types
                    observation_window=observation_window,
                    death_signal_column=death_signal,
                    min_contact_duration=min_contact,
                    killing_threshold_multiplier=threshold_mult,
                    save_results=save
                )
                
                print()
                print("=" * 60)
                print("RESULTS SUMMARY")
                print("=" * 60)
                print(f"Total qualifying contact events: {stats.get('total_contact_events', 0)}")
                print(f"Total contact timepoints analyzed: {stats.get('total_contact_timepoints', 0)}")
                print(f"Active killing timepoints: {stats.get('total_active_killing_timepoints', 0)}")
                print(f"Overall killing rate: {stats.get('overall_killing_rate', 0):.1%}")
                print()
                
                # Show per-sample background rates if available
                if 'background_death_rates' in stats:
                    print("Per-sample background death rates:")
                    for sample_name, rate in stats['background_death_rates'].items():
                        print(f"  {sample_name}: {rate:.6f} per timepoint")
                    print()
                
                if not df_summary.empty:
                    print("Per-sample and target type summary:")
                    display(df_summary)
                
                print()
                print(f"✅ Active Killing Analysis complete!")
                
                if save:
                    results_dir = Path(self.output_dir, "analysis", immune_cell, "active_killing")
                    print(f"   Results saved to: {results_dir}")
                
            except Exception:
                print("❌ Error during Active Killing Analysis:")
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)

class TrackFilterPanel:
    """
    Generic track filtering panel that works for ANY cell type.
    
    If Active Killing Analysis was run, this panel will automatically use the advanced
    features CSV which includes killing data. The filtered output will be saved to the
    active_killing folder to preserve the killing features for downstream analysis.
    """
    
    def __init__(self, metadata_loader, cell_type):
        """
        Parameters
        ----------
        metadata_loader : object
            Metadata loader instance with behav3d_parameters
        cell_type : str
            Name of the cell type to filter (e.g., "tcell", "organoid", "macrophage")
        """
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        # Use unified category detection
        self.category = detect_cell_type_category(self.cell_type, metadata_loader.metadata)
        
        # Check if advanced features exist (Active Killing Analysis was run)
        active_killing_dir = Path(self.output_dir, "analysis", self.cell_type, "active_killing")
        self._advanced_features_path = Path(active_killing_dir, f"BEHAV3D_{self.cell_type}_advanced_track_features.csv")
        self._use_advanced_features = self._advanced_features_path.exists()
        
        if self._use_advanced_features:
            print(f"[TrackFilterPanel] Advanced features FOUND for {self.cell_type} - will include active killing data")
        
        # Load/init config
        params = self.metadata_loader.behav3d_parameters
        params.setdefault("track_filtering", {})
        if self.cell_type not in params["track_filtering"]:
            # Use category defaults from _DEFAULT_CONFIG
            cat_cfg = _DEFAULT_CONFIG["track_filtering"].get(self.category, {})
            params["track_filtering"][self.cell_type] = deepcopy(cat_cfg)
            with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                yaml.safe_dump(params, f, sort_keys=False)
        
        cfg = params["track_filtering"][self.cell_type]
        
        # ---- Title ----
        self.section_title = widgets.HTML(
            f'<div style="font-size:22px;font-weight:700;">{self.cell_type.capitalize()} Track Filtering</div>'
        )
        
        # ---- Experiment duration ----
        self.en_exp_duration = widgets.Checkbox(
            description="Trim down full time series to supplied duration",
            value=bool(cfg.get("exp_duration_enabled", True)),
            indent=False
        )
        self.exp_duration = widgets.IntText(
            description="Max timepoints (frames)",
            value=int(cfg.get("exp_duration", 350)),
            style={'description_width': '180px'}
        )
        self.row_exp = widgets.HBox([self.exp_duration])
        if not self.en_exp_duration.value:
            self.row_exp.layout.display = "none"
        self.en_exp_duration.observe(
            lambda c: self.row_exp.layout.__setattr__("display", None if c["new"] else "none"),
            names="value"
        )
        
        # ---- Minimum track length ----
        self.en_min_length = widgets.Checkbox(
            description="Select only tracks with minimal length",
            value=bool(cfg.get("min_length_enabled", True)),
            indent=False
        )
        self.min_track_length = widgets.IntText(
            description="Minimal length (frames)",
            value=int(cfg.get("min_track_length", 30)),
            style={'description_width': '180px'}
        )
        self.row_min = widgets.HBox([self.min_track_length])
        if not self.en_min_length.value:
            self.row_min.layout.display = "none"
        self.en_min_length.observe(
            lambda c: self.row_min.layout.__setattr__("display", None if c["new"] else "none"),
            names="value"
        )
        
        # ---- Maximum track length (trim) ----
        self.en_max_length = widgets.Checkbox(
            description="Trim down tracks to supplied length",
            value=bool(cfg.get("max_length_enabled", True)),
            indent=False
        )
        self.max_track_length = widgets.IntText(
            description="Maximal length (frames)",
            value=int(cfg.get("max_track_length", 30)),
            style={'description_width': '180px'}
        )
        self.row_max = widgets.HBox([self.max_track_length])
        if not self.en_max_length.value:
            self.row_max.layout.display = "none"
        self.en_max_length.observe(
            lambda c: self.row_max.layout.__setattr__("display", None if c["new"] else "none"),
            names="value"
        )
        
        # ---- Filter t0 dead OR min size at t=1 ----
        from behav3d.utils import has_dead_channel
        self.has_dead = has_dead_channel(self.metadata_loader.metadata)
        
        if self.category == "organoid":
            # Organoids: filter by minimal size at t=1
            self.filter_t0_dead = widgets.Checkbox(
                description="Filter by minimal size at t=1",
                value=bool(cfg.get("filter_min_size_t1", True)),
                indent=False
            )
            self.min_size_t1 = widgets.IntText(
                description="Minimal size (px @ t=1)",
                value=int(cfg.get("min_size_t1", 1000)),
                style={'description_width': '180px'}
            )
            self.row_size_t1 = widgets.HBox([self.min_size_t1])
            if not self.filter_t0_dead.value:
                self.row_size_t1.layout.display = "none"
            self.filter_t0_dead.observe(
                lambda c: self.row_size_t1.layout.__setattr__("display", None if c["new"] else "none"),
                names="value"
            )
        else:
            # Immune/other: filter by dead at t=0
            self.filter_t0_dead = widgets.Checkbox(
                description="Filter tracks that are dead at t=0",
                value=bool(cfg.get("filter_t0_dead", True)),
                indent=False
            )
            self.min_size_t1 = None
            self.row_size_t1 = None
            # Show only if dead channel exists
            if not self.has_dead:
                self.filter_t0_dead.layout.display = "none"
        
        # ---- Time unit toggle buttons (not shown for organoid category) ----
        self.time_unit_label = widgets.HTML(
            value='<div style="margin-top:10px;"><b>Unit to use for time-based filters:</b></div>'
        )
        self.time_type = widgets.ToggleButtons(
            options=["frames", "hours"],
            value=cfg.get("time_type", "frames"),
            button_style="",
            style={'button_width': '100px'}
        )
        
        # Hide time unit toggle for organoid category
        if self.category == "organoid":
            self.time_unit_label.layout.display = "none"
            self.time_type.layout.display = "none"
        
        # Update labels when time unit changes
        def _on_time_unit_change(change):
            if change['name'] == 'value':
                unit = change['new']
                self.exp_duration.description = f"Max timepoints ({unit})"
                self.min_track_length.description = f"Minimal length ({unit})"
                self.max_track_length.description = f"Maximal length ({unit})"
        
        self.time_type.observe(_on_time_unit_change, names='value')
        
        # ---- Run button ----
        self.btn_run = widgets.Button(
            description=f"Filter {self.cell_type} tracks & summarize",
            button_style="success",
            layout=widgets.Layout(width="fit-content")
        )
        self.btn_run.on_click(self._on_run_clicked)
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html],
            layout=widgets.Layout(align_items="center", gap="10px")
        )
        
        self.out_run = widgets.Output()
        
        # ---- UI ----
        layout_widgets = [
            self.section_title,
            self.en_exp_duration,
            self.row_exp,
            self.en_min_length,
            self.row_min,
            self.en_max_length,
            self.row_max,
            self.filter_t0_dead,
        ]
        
        # Add size filter row for organoids
        if self.category == "organoid" and self.row_size_t1 is not None:
            layout_widgets.append(self.row_size_t1)
        
        layout_widgets.extend([
            self.time_unit_label,
            self.time_type,
            self.run_row,
            self.out_run,
        ])
        
        self.ui = widgets.VBox(layout_widgets)
        
        # Choose the correct filter function based on category
        if self.category == "organoid":
            self._filter_tracks = filter_organoid_tracks
        else:
            self._filter_tracks = filter_cell_tracks
    
    def _lock(self, locked):
        """Lock/unlock all controls"""
        lock_widgets = [self.en_exp_duration, self.exp_duration, self.en_min_length,
                  self.min_track_length, self.en_max_length, self.max_track_length,
                  self.filter_t0_dead, self.time_type, self.btn_run]
        
        # Add organoid-specific widgets if they exist
        if self.category == "organoid" and self.min_size_t1 is not None:
            lock_widgets.append(self.min_size_t1)
        
        for w in lock_widgets:
            w.disabled = locked
    
    def _on_run_clicked(self, *_):
        """Run filtering when button clicked"""
        self._lock(True)
        self.out_run.clear_output()
        self.spinner_html.layout.display = None
        
        with self.out_run:
            try:
                print(f"▶️ Filtering {self.cell_type} tracks...")
                
                # Gather parameters
                output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
                metadata = self.metadata_loader.metadata
                
                exp_duration = int(self.exp_duration.value) if self.en_exp_duration.value else None
                min_track_length = int(self.min_track_length.value) if self.en_min_length.value else None
                max_track_length = int(self.max_track_length.value) if self.en_max_length.value else None
                
                # For organoids: size filter at t=1 instead of dead filter
                if self.category == "organoid":
                    filter_t0_dead = False
                    min_size_t1 = int(self.min_size_t1.value) if self.filter_t0_dead.value else None
                else:
                    # Immune/other: dead filter (only if dead channel exists)
                    filter_t0_dead = bool(self.filter_t0_dead.value) if self.has_dead else False
                    min_size_t1 = None
                
                time_type = str(self.time_type.value)
                plot_results = True  # Always generate plots
                
                # Persist to config
                cfg = self.metadata_loader.behav3d_parameters["track_filtering"][self.cell_type]
                cfg["exp_duration_enabled"] = self.en_exp_duration.value
                cfg["exp_duration"] = int(self.exp_duration.value)
                cfg["min_length_enabled"] = self.en_min_length.value
                cfg["min_track_length"] = int(self.min_track_length.value)
                cfg["max_length_enabled"] = self.en_max_length.value
                cfg["max_track_length"] = int(self.max_track_length.value)
                cfg["time_type"] = time_type
                cfg["plot_results"] = plot_results
                
                # Save category-specific configs
                if self.category == "organoid":
                    cfg["filter_min_size_t1"] = self.filter_t0_dead.value
                    cfg["min_size_t1"] = int(self.min_size_t1.value) if self.min_size_t1 else 500
                else:
                    cfg["filter_t0_dead"] = filter_t0_dead
                
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
                
                print(f"  output_dir = {output_dir}")
                if exp_duration is not None:
                    print(f"  exp_duration = {exp_duration}")
                if min_track_length is not None:
                    print(f"  min_track_length = {min_track_length}")
                if max_track_length is not None:
                    print(f"  max_track_length = {max_track_length}")
                if self.category == "organoid":
                    if min_size_t1 is not None:
                        print(f"  min_size_t1 = {min_size_t1}")
                elif self.has_dead:
                    print(f"  filter_t0_dead = {filter_t0_dead}")
                
                # Determine input path: use advanced features if available
                df_input_path = None
                if self._use_advanced_features:
                    df_input_path = str(self._advanced_features_path)
                    print(f"  Using ADVANCED features with active killing data")
                
                # Call adapted filter function (organoid or tcell)
                if self.category == "organoid":
                    self._filter_tracks(
                        metadata=metadata,
                        output_dir=output_dir,
                        exp_duration=exp_duration,
                        min_track_length=min_track_length,
                        max_track_length=max_track_length,
                        min_size=min_size_t1,
                        time_type=time_type,
                        cell_type=self.cell_type,
                        df_input_path=df_input_path
                    )
                    
                    # Summarize tracks after filtering (needed for behavioral analysis)
                    print("✅ Filtering done. Summarizing tracks...")
                    self.df_tracks_summ = summarize_track_features(
                        output_dir=str(output_dir),
                        cell_type=self.cell_type
                    )
                else:
                    self._filter_tracks(
                        metadata=metadata,
                        output_dir=output_dir,
                        exp_duration=exp_duration,
                        min_track_length=min_track_length,
                        max_track_length=max_track_length,
                        filter_t0_dead=filter_t0_dead,
                        cell_type=self.cell_type,
                        time_type=time_type,
                        plot_results=plot_results,
                        df_input_path=df_input_path
                    )
                
                # Summarize tracks after filtering (needed for behavioral analysis)
                print("✅ Filtering done. Summarizing tracks...")
                self.df_tracks_summ = summarize_track_features(
                    output_dir=str(output_dir),
                    cell_type=self.cell_type
                )
                
                print(f"✅ {self.cell_type} filtering complete!")
                
            except Exception:
                print(f"❌ Error while filtering {self.cell_type} tracks:")
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)


class MotileCellAnalysisPanel:
    """
    Generic behavioral analysis panel (DTW/UMAP/clustering) for ANY cell type.
    
    """
    
    def __init__(self, metadata_loader, cell_type):
        """
        Parameters
        ----------
        metadata_loader : object
            Metadata loader instance
        cell_type : str
            Name of the cell type (e.g., "tcell", "organoid", "macrophage")
        """
        from behav3d.utils import expand_column_patterns
        
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        # Use unified category detection
        self.category = detect_cell_type_category(self.cell_type, metadata_loader.metadata)
        
        # Load config
        params = dict(self.metadata_loader.behav3d_parameters or {})
        params.setdefault("analysis", {})
        if self.cell_type not in params["analysis"]:
            cat_cfg = _DEFAULT_CONFIG["analysis"].get(self.category, {})
            params["analysis"][self.cell_type] = deepcopy(cat_cfg)
        self._params = params
        self._panel_cfg = self._params["analysis"][self.cell_type]
        
        # Feature groups - expand wildcards based on ACTUAL DATA
        groups = deepcopy(behav3d_calculated_features)
        
        # Load a sample of the data to get actual column names
        feature_outdir = Path(self.output_dir, "analysis", self.cell_type, "track_features")
        active_killing_dir = Path(self.output_dir, "analysis", self.cell_type, "active_killing")
        
        df_tracks_path_filt = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")
        df_tracks_path_unfilt = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features.csv")
        # Check for advanced track features (with active killing data)
        df_advanced_path = Path(active_killing_dir, f"BEHAV3D_{self.cell_type}_advanced_track_features.csv")
        
        actual_columns = []
        self._use_advanced_features = False
        self._advanced_features_path = None
        
        # First, check if advanced track features exist (has active killing data)
        if df_advanced_path.exists():
            try:
                import pandas as pd
                df_sample = pd.read_csv(df_advanced_path, nrows=0)
                actual_columns = list(df_sample.columns)
                self._use_advanced_features = True
                self._advanced_features_path = df_advanced_path
                print(f"[MotileCellAnalysisPanel] Using ADVANCED features ({len(actual_columns)} columns) from {df_advanced_path.name} for {self.cell_type}")
            except Exception as e:
                print(f"[MotileCellAnalysisPanel] Failed to load advanced features: {e}")
        
        # Fallback to regular track features if no advanced features
        if not actual_columns:
            # Try filtered first, then unfiltered as fallback
            for df_tracks_path in [df_tracks_path_filt, df_tracks_path_unfilt]:
                if df_tracks_path.exists():
                    try:
                        # Read just the first row to get column names
                        import pandas as pd
                        df_sample = pd.read_csv(df_tracks_path, nrows=0)
                        actual_columns = list(df_sample.columns)
                        print(f"[MotileCellAnalysisPanel] Loaded {len(actual_columns)} columns from {df_tracks_path.name} for {self.cell_type}")
                        break  # Stop after first successful load
                    except Exception as e:
                        print(f"[MotileCellAnalysisPanel] Failed to load {df_tracks_path.name}: {e}")
                        pass
        
        if not actual_columns:
            print(f"[MotileCellAnalysisPanel] WARNING: No CSV found for {self.cell_type} - using template features only")
        
        # Check if death features actually exist in the data
        has_death_features = False
        if actual_columns:
            death_cols = {'mean_dead_dye', 'percentage_dead_mask', 'nr_dead_mask_pixels', 'increase_dead_mask', 'dead'}
            has_death_features = bool(death_cols.intersection(set(actual_columns)))
        
        # Remove death features if not present in actual data
        if not has_death_features:
            print(f"[MotileCellAnalysisPanel] No death features found in {self.cell_type} data - hiding death group")
            if "death" in groups:
                del groups["death"]
            if "intensity" in groups:
                groups["intensity"] = [f for f in groups["intensity"] if f != "mean_dead_dye"]
        
        # Check if active killing features exist (only present if Active Killing Analysis was run)
        has_active_killing = False
        if actual_columns:
            killing_cols = {'is_active_killing', 'killing_efficiency', 'n_killing_events_total'}
            has_active_killing = bool(killing_cols.intersection(set(actual_columns)))
        
        if not has_active_killing:
            if "active_killing" in groups:
                del groups["active_killing"]
        else:
            print(f"[MotileCellAnalysisPanel] Active killing features found for {self.cell_type}")
        
        # Expand wildcards in feature groups to show actual columns
        if actual_columns:
            expanded_groups = {}
            for group_name, patterns in groups.items():
                expanded_features = []
                for pattern in patterns:
                    if '*' in pattern or '?' in pattern or '[' in pattern:
                        # This is a wildcard - expand it
                        matches = expand_column_patterns(pattern, actual_columns)
                        # Filter out non-numeric columns
                        matches = [m for m in matches if m != 'orientation_vector' and not m.startswith('touching_')]
                        expanded_features.extend(matches)
                    else:
                        # Exact match - include if exists
                        if pattern in actual_columns:
                            expanded_features.append(pattern)
                # Only include groups that have at least one feature
                if expanded_features:
                    expanded_groups[group_name] = sorted(set(expanded_features))
            groups = expanded_groups
            print(f"[MotileCellAnalysisPanel] Expanded to {sum(len(v) for v in groups.values())} total features across {len(groups)} groups")
        
        # ---- Titles ----
        self.section_title = widgets.HTML(
            f'<div style="font-size:22px;font-weight:700;">{self.cell_type} behavioral analysis</div>'
        )
        self.sel_title = widgets.HTML(
            '<div style="font-size:20px;font-weight:700;">Select features for DTW:</div>'
        )
        self.norm_title = widgets.HTML(
            '<div style="font-size:20px;font-weight:700;">Normalize (z-score):</div>'
        )
        self.umap_title = widgets.HTML(
            '<div style="font-size:20px;font-weight:700;">UMAP settings:</div>'
        )
        
        # ---- Seed ----
        self.seed_widget = widgets.IntText(
            description="Seed",
            value=int(self._panel_cfg.get("seed", self._params.get("seed", 42))),
            style={"description_width": "80px"},
        )
        
        # ---- Grid layout helper ----
        def _grid_for(children, indent_px="24px", columns=3):
            return widgets.GridBox(
                children,
                layout=widgets.Layout(
                    grid_template_columns=" ".join(["max-content"] * columns),
                    grid_gap="6px 18px",
                    margin=f"0 0 0 {indent_px}",
                )
            )
        
        # ---- SELECTION (feature groups) ----
        self._group_rows = {}
        sel_group_boxes = []
        
        preset = list(self._panel_cfg.get("dtw_features_input", []))
        preset_set = set(preset)
        groups_enabled = dict(self._panel_cfg.get("dtw_feature_groups_enabled", {}))
        
        # Get all available features from expanded groups
        all_expanded_features = set()
        for feats in groups.values():
            all_expanded_features.update(feats)
        
        # Check if saved config is stale:
        # 1. Has wildcards (e.g., mean_intensity_*)
        # 2. Saved features don't exist in expanded set
        # 3. Saved config is too incomplete (< 20% of available features)
        has_wildcards_in_config = any('*' in f or '?' in f or '[' in f for f in preset)
        features_dont_exist = preset_set and not preset_set.intersection(all_expanded_features)
        too_incomplete = len(preset_set) > 0 and len(preset_set) < len(all_expanded_features) * 0.2
        
        config_is_stale = has_wildcards_in_config or features_dont_exist or too_incomplete
        
        # If stale config detected, reset to show nothing selected (user will select what they want)
        if config_is_stale:
            if too_incomplete:
                print(f"[MotileCellAnalysisPanel] Incomplete config for {self.cell_type} ({len(preset_set)}/{len(all_expanded_features)} features) - clearing selection")
            else:
                print(f"[MotileCellAnalysisPanel] Stale config detected for {self.cell_type} - clearing selection")
            preset_set.clear()
            # Update saved config to clear the stale data
            self._panel_cfg["dtw_features_input"] = []
            self._panel_cfg["z_normalize"] = {}
            # Don't auto-select features - let user choose
        elif not preset_set:
            # No saved config at all - start with enabled groups
            print(f"[MotileCellAnalysisPanel] No saved config for {self.cell_type} - using category defaults")
            for group_name, feats in groups.items():
                if groups_enabled.get(group_name, False):
                    preset_set.update(feats)
            if preset_set:
                print(f"[MotileCellAnalysisPanel] Initialized {len(preset_set)} features from enabled groups")
                # Update config with initialized features
                self._panel_cfg["dtw_features_input"] = list(preset_set)
            else:
                print(f"[MotileCellAnalysisPanel] No features auto-selected (all groups disabled in category defaults)")
        
        for group_name, feats in groups.items():
            child_cbs = [
                widgets.Checkbox(value=(f in preset_set), description=f, indent=True)
                for f in feats
            ]
            # Group checkbox: use saved state if available, otherwise check if any children selected
            g_init = groups_enabled.get(group_name, any(cb.value for cb in child_cbs)) if groups_enabled else any(cb.value for cb in child_cbs)
            
            gcb = widgets.Checkbox(value=g_init, indent=False)
            glabel = widgets.HTML(f"<b>{group_name}</b>")
            header = widgets.HBox([gcb, glabel])
            
            grid = _grid_for(child_cbs, columns=3)
            container = widgets.VBox([header, grid])
            
            self._group_rows[group_name] = {
                "group_cb": gcb, "child_cbs": child_cbs, "grid": grid, "container": container
            }
            sel_group_boxes.append(container)
        
        self.checkboxes_box = widgets.VBox(sel_group_boxes)
        
        # Wire selection events
        def make_group_handler(grp_name):
            def _on_group(change):
                if change["name"] != "value":
                    return
                val = bool(change["new"])
                row = self._group_rows[grp_name]
                for cb in row["child_cbs"]:
                    cb.unobserve(row["child_handler"], names="value")
                for cb in row["child_cbs"]:
                    cb.value = val
                for cb in row["child_cbs"]:
                    cb.observe(row["child_handler"], names="value")
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
        
        # Selection action buttons
        self.btn_select_all = widgets.Button(description="Select all")
        self.btn_clear_all = widgets.Button(description="Clear")
        self.btn_select_all.on_click(lambda *_: self._set_all(True))
        self.btn_clear_all.on_click(lambda *_: self._set_all(False))
        
        # ---- NORMALIZE (mirrors selection) ----
        self._norm_group_rows = {}
        norm_group_boxes = []
        
        for group_name, feats in groups.items():
            gcb = widgets.Checkbox(value=False, indent=False)
            glabel = widgets.HTML(f"<b>{group_name}</b>")
            header = widgets.HBox([gcb, glabel])
            
            child_cbs = []
            for f in feats:
                cb = widgets.Checkbox(
                    value=bool(self._panel_cfg.get("z_normalize", {}).get(f, False)),
                    description=f,
                    indent=True
                )
                cb.layout.visibility = 'hidden'
                child_cbs.append(cb)
            
            grid = _grid_for(child_cbs, columns=3)
            container = widgets.VBox([header, grid])
            
            self._norm_group_rows[group_name] = {
                "group_cb": gcb, "child_cbs": child_cbs, "grid": grid, "container": container
            }
            norm_group_boxes.append(container)
        
        # Normalize group handlers
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
        
        # Normalize action buttons
        self.norm_select_all_btn = widgets.Button(description="Select all")
        self.norm_clear_all_btn = widgets.Button(description="Clear")
        self.norm_select_all_btn.on_click(self._on_norm_select_all)
        self.norm_clear_all_btn.on_click(self._on_norm_clear_all)
        
        self.normalize_section = widgets.VBox([
            self.norm_title,
            widgets.HBox([self.norm_select_all_btn, self.norm_clear_all_btn],
                         layout=widgets.Layout(gap="8px")),
            widgets.VBox(norm_group_boxes)
        ])
        
        # ---- UMAP settings ----
        self.umap_distance_widget = widgets.FloatText(
            description="UMAP min_dist",
            value=float(self._panel_cfg.get("umap_min_dist", 0.1)),
            style={"description_width": "140px"}
        )
        self.umap_neighbors_widget = widgets.IntText(
            description="UMAP n_neighbors",
            value=int(self._panel_cfg.get("umap_n_neighbors", 15)),
            style={"description_width": "140px"}
        )
        self.num_clusters_widget = widgets.IntText(
            description="# clusters",
            value=int(self._panel_cfg.get("nr_of_clusters", 5)),
            style={"description_width": "140px"}
        )
        
        self.umap_box = widgets.VBox([
            self.umap_title,
            widgets.HBox([self.umap_distance_widget, self.umap_neighbors_widget, self.num_clusters_widget],
                         layout=widgets.Layout(flex_flow="row wrap", gap="12px"))
        ])
        
        # ---- Cluster percentage grouping selector ----
        self.groupby_title = widgets.HTML(
            '<div style="font-size:20px;font-weight:700;">Cluster % grid grouping:</div>'
        )
        self.groupby_description = widgets.HTML(
            '<div style="font-size:12px;color:#666;margin-bottom:8px;">'
            f'Select columns to group by in cluster percentage plots. '
            f'Rows will always be {self.cell_type}_line_condition values. '
            f'If none selected, only per-sample plots are generated.</div>'
        )
        
        # Build list of eligible columns for grouping
        # Eligible: exp_nr and all *_line_condition EXCEPT this cell type's own (EXCLUDE 'well')
        this_line_col = f"{self.cell_type}_line_condition"
        metadata = self.metadata_loader.metadata
        eligible_cols = []
        if metadata is not None:
            for col in metadata.columns:
                if col == 'exp_nr':  # 'well' is now excluded
                    eligible_cols.append(col)
                elif col.endswith('_line_condition') and col != this_line_col:
                    eligible_cols.append(col)
        
        # Get saved selection from config
        saved_groupby = list(self._panel_cfg.get("cluster_percentage_group_by", []))
        # Filter to only include columns that still exist
        saved_groupby = [c for c in saved_groupby if c in eligible_cols]
        
        self.groupby_selector = widgets.SelectMultiple(
            options=eligible_cols,
            value=saved_groupby,
            description="Group by:",
            style={"description_width": "80px"},
            layout=widgets.Layout(width="400px", height="100px")
        )
        
        self.groupby_box = widgets.VBox([
            self.groupby_title,
            self.groupby_description,
            self.groupby_selector
        ])
        
        # ---- Run button ----
        self.btn_run = widgets.Button(
            description=f"Run {self.cell_type} analysis",
            button_style="success",
            layout=widgets.Layout(width="300px")
        )
        self.btn_run.on_click(self._on_run_clicked)
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
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
            self.normalize_section,
            widgets.HTML("<hr>"),
            self.umap_box,
            widgets.HTML("<hr>"),
            self.groupby_box,
            widgets.HTML("<hr>"),
            self.run_row,
            self.out_run,
        ])
        
        # Internal state
        self.df_tracks_clustered = None
        
        # Store function reference
        self._expand_column_patterns = expand_column_patterns
        
        # Initial normalize visibility sync
        self._sync_normalize_visibility_batched(list(self._group_rows.keys()))
    
    def _selected_features(self):
        """Get list of selected feature patterns"""
        feats = []
        for row in self._group_rows.values():
            for cb in row["child_cbs"]:
                if cb.value:
                    feats.append(cb.description)
        return feats
    
    def _selected_normalize_patterns(self):
        """Get list of normalize patterns (only visible ones that are checked)"""
        pats = []
        for row in self._norm_group_rows.values():
            for cb in row["child_cbs"]:
                if cb.layout.visibility != 'hidden' and cb.value:
                    pats.append(cb.description)
        return pats
    
    def _visible_norm_children(self, grp_name):
        """Get visible normalize checkboxes for a group"""
        row = self._norm_group_rows[grp_name]
        return [cb for cb in row["child_cbs"] if cb.layout.visibility != 'hidden']
    
    def _sync_normalize_visibility_batched(self, groups_to_update):
        """Mirror selection to normalize section"""
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
            
            # Update group checkbox
            vis_children = self._visible_norm_children(grp_name)
            gh = nrow.get("group_handler")
            nrow["group_cb"].unobserve(gh, names="value")
            nrow["group_cb"].value = bool(vis_children) and all(cb.value for cb in vis_children)
            nrow["group_cb"].observe(gh, names="value")
            
            # Show/hide container
            nrow["container"].layout.display = None if any_visible else "none"
        
        # Show/hide whole section
        any_group_visible = any(
            self._norm_group_rows[g]["container"].layout.display != "none"
            for g in self._norm_group_rows
        )
        self.normalize_section.layout.display = None if any_group_visible else "none"
    
    def _set_all(self, value):
        """Select all or clear all features"""
        for grp_name, row in self._group_rows.items():
            row["group_cb"].unobserve(row["group_handler"], names="value")
            row["group_cb"].value = value
            row["group_cb"].observe(row["group_handler"], names="value")
            
            for cb in row["child_cbs"]:
                cb.unobserve(row["child_handler"], names="value")
            for cb in row["child_cbs"]:
                cb.value = value
            for cb in row["child_cbs"]:
                cb.observe(row["child_handler"], names="value")
        
        self._sync_normalize_visibility_batched(list(self._group_rows.keys()))
    
    def _on_norm_select_all(self, *_):
        """Select all visible normalize checkboxes"""
        for grp_name, row in self._norm_group_rows.items():
            for cb in row["child_cbs"]:
                if cb.layout.visibility != 'hidden':
                    cb.value = True
            
            vis_children = self._visible_norm_children(grp_name)
            gh = row.get("group_handler")
            row["group_cb"].unobserve(gh, names="value")
            row["group_cb"].value = bool(vis_children) and all(cb.value for cb in vis_children)
            row["group_cb"].observe(gh, names="value")
    
    def _on_norm_clear_all(self, *_):
        """Clear all visible normalize checkboxes"""
        for grp_name, row in self._norm_group_rows.items():
            for cb in row["child_cbs"]:
                if cb.layout.visibility != 'hidden':
                    cb.value = False
            
            gh = row.get("group_handler")
            row["group_cb"].unobserve(gh, names="value")
            row["group_cb"].value = False
            row["group_cb"].observe(gh, names="value")
    
    def _expand_patterns(self, patterns):
        """Expand wildcard patterns using actual CSV columns"""
        # Use advanced features path if available (has active killing columns)
        if self._use_advanced_features and self._advanced_features_path and self._advanced_features_path.exists():
            df = pd.read_csv(self._advanced_features_path, nrows=1)
            return self._expand_column_patterns(patterns, df.columns.tolist())
        
        # Load a sample CSV to get column names
        metadata = self.metadata_loader.metadata
        if metadata is None or len(metadata) == 0:
            return patterns
        
        # Get first available CSV
        csv_col = f"{self.cell_type}_combined_track_features_csv_path"
        if csv_col not in metadata.columns:
            return patterns
        
        csv_path = metadata[csv_col].dropna().iloc[0] if len(metadata[csv_col].dropna()) > 0 else None
        if csv_path is None or not Path(csv_path).exists():
            return patterns
        
        # Read columns
        df = pd.read_csv(csv_path, nrows=1)
        return self._expand_column_patterns(patterns, df.columns.tolist())
    
    def _lock(self, locked):
        """Lock/unlock all controls"""
        for row in self._group_rows.values():
            row["group_cb"].disabled = locked
            for cb in row["child_cbs"]:
                cb.disabled = locked
        
        for row in self._norm_group_rows.values():
            row["group_cb"].disabled = locked
            for cb in row["child_cbs"]:
                cb.disabled = locked
        
        lock_widgets = [self.seed_widget, self.btn_select_all, self.btn_clear_all,
                  self.umap_distance_widget, self.umap_neighbors_widget,
                  self.num_clusters_widget, self.norm_select_all_btn,
                  self.norm_clear_all_btn, self.groupby_selector, self.btn_run]
        
        for w in lock_widgets:
            w.disabled = locked
    
    def _on_run_clicked(self, *_):
        """Run analysis when button clicked"""
        self._lock(True)
        self.out_run.clear_output()
        self.spinner_html.layout.display = None
        
        with self.out_run:
            try:
                print(f"▶️ Running {self.cell_type} behavioral analysis...")
                
                seed = int(self.seed_widget.value)
                umap_min_dist = float(self.umap_distance_widget.value)
                umap_n_neighbors = int(self.umap_neighbors_widget.value)
                nr_of_clusters = int(self.num_clusters_widget.value)
                
                # Get cluster percentage grouping columns
                cluster_percentage_group_by = list(self.groupby_selector.value)
                
                # Get selected features (patterns)
                dtw_patterns = self._selected_features()
                if not dtw_patterns:
                    print("⚠️ Please select at least one DTW feature.")
                    return
                
                # Get normalize patterns (from checkboxes - includes contact features if checked)
                norm_patterns = self._selected_normalize_patterns()
                
                # Expand patterns using CSV columns
                columns_to_use = self._expand_patterns(dtw_patterns)
                columns_to_normalize = self._expand_patterns(norm_patterns)
                
                # Persist config
                self._panel_cfg["seed"] = seed
                self._panel_cfg["umap_min_dist"] = umap_min_dist
                self._panel_cfg["umap_n_neighbors"] = umap_n_neighbors
                self._panel_cfg["nr_of_clusters"] = nr_of_clusters
                self._panel_cfg["cluster_percentage_group_by"] = cluster_percentage_group_by
                self._panel_cfg["dtw_features_input"] = dtw_patterns
                self._panel_cfg["dtw_feature_groups_enabled"] = {
                    g: row["group_cb"].value for g, row in self._group_rows.items()
                }
                self._panel_cfg["z_normalize"] = {
                    cb.description: cb.value
                    for row in self._norm_group_rows.values()
                    for cb in row["child_cbs"]
                    if cb.layout.visibility != 'hidden'
                }
                self._panel_cfg["columns_to_normalize_input"] = norm_patterns
                
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(self._params, f, sort_keys=False)
                
                print(f"  output_dir = {self.output_dir}")
                print(f"  seed = {seed}")
                print(f"  UMAP: n_neighbors={umap_n_neighbors}, min_dist={umap_min_dist}")
                print(f"  clusters = {nr_of_clusters}")
                print(f"  cluster_percentage_group_by = {cluster_percentage_group_by}")
                print(f"  columns_to_use [{len(columns_to_use)}]: {columns_to_use}")
                print(f"  columns_to_normalize [{len(columns_to_normalize)}]: {columns_to_normalize}")
                
                Path(self.output_dir).mkdir(parents=True, exist_ok=True)
                import random
                random.seed(seed)
                
                # Determine which CSV to use for DTW analysis
                feature_outdir = Path(self.output_dir, "analysis", self.cell_type, "track_features")
                active_killing_dir = Path(self.output_dir, "analysis", self.cell_type, "active_killing")
                
                filtered_csv_path = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")
                summarized_csv_path = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features_summarized.csv")
                
                # ENFORCE: Summarized CSV is REQUIRED for DTW analysis (clustering needs it)
                if not summarized_csv_path.exists():
                    print(f"❌ ERROR: Summarized track features not found!")
                    print(f"   Expected: {summarized_csv_path}")
                    print(f"")
                    print(f"   ⚠️ You MUST run 'Filter {self.cell_type} tracks & summarize' before running behavioral analysis.")
                    print(f"   The summarization step is required for clustering.")
                    return
                
                # ENFORCE: Filtered CSV is REQUIRED (unfiltered tracks cause issues with DTW)
                if not filtered_csv_path.exists():
                    print(f"❌ ERROR: Filtered track features not found!")
                    print(f"   Expected: {filtered_csv_path}")
                    print(f"")
                    print(f"   ⚠️ You MUST run 'Filter {self.cell_type} tracks & summarize' before running behavioral analysis.")
                    print(f"   Filtering ensures tracks have equal lengths for Dynamic Time Warping.")
                    return
                
                df_tracks_path = filtered_csv_path
                
                # Check if filtered CSV has active killing columns
                if self._use_advanced_features and self._advanced_features_path:
                    try:
                        df_check = pd.read_csv(filtered_csv_path, nrows=0)
                        if 'is_active_killing' in df_check.columns or 'killing_efficiency' in df_check.columns:
                            print(f"  Using filtered track features WITH active killing data")
                        else:
                            print(f"  ⚠️ Note: Active killing features exist but filtering was run before Active Killing Analysis.")
                            print(f"     To include killing features in DTW, re-run 'Filter {self.cell_type} tracks & summarize'.")
                    except Exception:
                        pass
                else:
                    print(f"  Using filtered track features")
                
                # Call adapted tcell_analysis with cell_type parameter
                self.df_tracks_clustered = run_tcell_analysis(
                    cell_type=self.cell_type,
                    output_dir=self.output_dir,
                    df_tracks_path=df_tracks_path,
                    columns_to_use=columns_to_use,
                    columns_to_normalize=columns_to_normalize,
                    umap_minimal_distance=umap_min_dist,
                    umap_n_neighbors=umap_n_neighbors,
                    nr_of_clusters=nr_of_clusters,
                    cluster_percentage_group_by=cluster_percentage_group_by,
                    plot_results=True,
                    seed=seed
                )
                
                # Show preview
                try:
                    display(self.df_tracks_clustered.head(10))
                except Exception:
                    pass
                
                print(f"✅ {self.cell_type} analysis complete!")
                
            except Exception:
                print(f"❌ Error while running {self.cell_type} analysis:")
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)


class DeathDynamicsPanel:
    """
    Death dynamics analysis panel (organoid-specific).
    
    """
    
    def __init__(self, metadata_loader, cell_type):
        """
        Parameters
        ----------
        metadata_loader : object
            Metadata loader instance
        cell_type : str
            Name of the organoid cell type (e.g., "organoid")
        """
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        # Check if death features exist in the data
        feature_outdir = Path(self.output_dir, "analysis", self.cell_type, "track_features")
        df_tracks_path = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features.csv")
        
        self.has_death_features = False
        if df_tracks_path.exists():
            try:
                import pandas as pd
                df_sample = pd.read_csv(df_tracks_path, nrows=0)
                death_cols = {'mean_dead_dye', 'percentage_dead_mask', 'nr_dead_mask_pixels'}
                self.has_death_features = bool(death_cols.intersection(set(df_sample.columns)))
            except Exception:
                pass
        
        print(f"[DeathDynamicsPanel] Death features {'FOUND' if self.has_death_features else 'NOT FOUND'} for {self.cell_type}")
        
        # Load config
        params = dict(self.metadata_loader.behav3d_parameters or {})
        params.setdefault("death_dynamics", {})
        if self.cell_type not in params["death_dynamics"]:
            cat_cfg = _DEFAULT_CONFIG.get("death_dynamics", {}).get("organoid", {})
            params["death_dynamics"][self.cell_type] = deepcopy(cat_cfg)
        self._params = params
        self._panel_cfg = self._params["death_dynamics"][self.cell_type]
        
        # ---- Title ----
        self.section_title = widgets.HTML(
            f'<div style="font-size:22px;font-weight:700;">{self.cell_type} death dynamics</div>'
        )
        
        # If no death features, show warning message
        if not self.has_death_features:
            self.warning_html = widgets.HTML(
                value='<div style="color:#b00;font-size:14px;padding:10px;background:#fee;border:1px solid #fcc;border-radius:4px;">'
                      '<b>⚠️ Death Dynamics Not Available</b><br>'
                      f'No death features found for {self.cell_type}.<br>'
                      'Death channel was not present during feature extraction.'
                      '</div>'
            )
            self.ui = widgets.VBox([self.section_title, self.warning_html])
            return  # Don't create the rest of the UI
        
        # ---- Dead percentage threshold ----
        self.dead_perc_threshold = widgets.FloatText(
            description="Dead % threshold",
            value=float(self._panel_cfg.get("dead_perc_threshold", 0.02)),
            style={'description_width': '160px'},
            layout=widgets.Layout(width="220px")
        )
        self.dead_threshold_help = widgets.HTML(
            '<div style="font-size:12px;color:#666;margin-left:10px;">'
            '(percentage_dead_mask threshold to classify as dead)'
            '</div>'
        )
        
        # ---- Run button ----
        self.btn_run = widgets.Button(
            description=f"Run {self.cell_type} death dynamics",
            button_style="warning",
            layout=widgets.Layout(width="300px")
        )
        self.btn_run.on_click(self._on_run_clicked)
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html],
            layout=widgets.Layout(align_items="center", gap="10px")
        )
        
        self.out_run = widgets.Output()
        
        # ---- UI ----
        self.ui = widgets.VBox([
            self.section_title,
            widgets.HBox([self.dead_perc_threshold, self.dead_threshold_help], layout=widgets.Layout(align_items="center")),
            widgets.HTML("<hr>"),
            self.run_row,
            self.out_run,
        ])
        
        # Internal state
        self.death_data = None
    
    def _lock(self, locked):
        """Lock/unlock all controls"""
        for w in [self.dead_perc_threshold, self.btn_run]:
            w.disabled = locked
    
    def _on_run_clicked(self, *_):
        """Run death dynamics when button clicked"""
        self._lock(True)
        self.out_run.clear_output()
        self.spinner_html.layout.display = None
        
        with self.out_run:
            try:
                print(f"▶️ Running {self.cell_type} death dynamics...")
                
                dead_perc_threshold = float(self.dead_perc_threshold.value)
                plot_results = True  # Always generate plots
                
                # Persist config
                self._panel_cfg["dead_perc_threshold"] = dead_perc_threshold
                self._panel_cfg["plot_results"] = plot_results
                
                with self.metadata_loader.behav3d_parameters_path.open("w", encoding="utf-8") as f:
                    yaml.safe_dump(self._params, f, sort_keys=False)
                
                print(f"  output_dir = {self.output_dir}")
                print(f"  dead_perc_threshold = {dead_perc_threshold}")
                print(f"  plot_results = {plot_results}")
                
                Path(self.output_dir).mkdir(parents=True, exist_ok=True)
                
                # Call organoid death dynamics function
                run_organoid_analysis(
                    dead_perc_threshold=dead_perc_threshold,
                    output_dir=self.output_dir,
                    df_tracks_path=None,
                    org_type=self.cell_type,
                    metadata=self.metadata_loader.metadata  # Pass metadata so it can find dead_mask_path
                )
                
                print(f"✅ {self.cell_type} death dynamics complete!")
                
                # Show preview
                try:
                    display(self.death_data.head(10))
                except Exception:
                    pass
                
            except Exception:
                print(f"❌ Error while running {self.cell_type} death dynamics:")
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)

class InteractionAnalysisPanel:
    """
    Interaction analysis panel for organoids.
    
    Analyzes interactions from the organoid's point of view:
    - Cumulative interactions over time with selected cell types
    - Comparison between organoids that survive vs die
    - Per-sample breakdowns
    
    Outputs PDF plots to analysis/{organoid_type}/interaction_analysis/
    Each selected cell type gets its own separate PDF with statistics.
    
    Note: Uses death classification from Step 3 (Death Dynamics) if available.
    """
    
    def __init__(self, metadata_loader, cell_type):
        """
        Parameters
        ----------
        metadata_loader : object
            Metadata loader instance
        cell_type : str
            Name of the organoid cell type (e.g., "organoid")
        """
        from behav3d.utils import (
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata
        )
        
        self.metadata_loader = metadata_loader
        self.cell_type = str(cell_type).strip()
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())
        
        # Detect available cell types for interaction analysis
        md = self.metadata_loader.metadata
        self.immune_types = detect_immune_cell_types_from_metadata(md)
        self.other_types = detect_other_cell_types_from_metadata(md)
        self.available_cell_types = self.immune_types + self.other_types
        
        # Path to filtered data
        feature_outdir = Path(self.output_dir, "analysis", self.cell_type, "track_features")
        self.df_tracks_path = Path(feature_outdir, f"BEHAV3D_{self.cell_type}_combined_track_features_filtered.csv")
        
        # Load death threshold from config (saved by Step 3: Death Dynamics)
        self._load_death_threshold_from_config()
        
        # ---- Title ----
        self.section_title = widgets.HTML(
            f'<div style="font-size:22px;font-weight:700;">{self.cell_type} Interaction Analysis</div>'
        )
        
        # ---- Status message (will be updated dynamically) ----
        self.status_html = widgets.HTML(value="")
        
        # ---- Refresh button to check for new files ----
        self.btn_refresh = widgets.Button(
            description="🔄 Refresh",
            button_style="info",
            layout=widgets.Layout(width="100px")
        )
        self.btn_refresh.on_click(self._on_refresh_clicked)
        
        # ---- Cell type selector ----
        self.cell_type_selector_label = widgets.HTML(
            '<div style="font-weight:600;margin-bottom:5px;">Select cell types to analyze interactions with:</div>'
        )
        
        # Container for dynamically updated checkboxes
        self.cell_type_box = widgets.VBox([])
        self.cell_type_checkboxes = {}
        
        # ---- Death threshold info ----
        self.death_info_html = widgets.HTML(value="")
        
        # ---- Run button ----
        self.btn_run = widgets.Button(
            description=f"Run Interaction Analysis",
            button_style="warning",
            layout=widgets.Layout(width="250px")
        )
        self.btn_run.on_click(self._on_run_clicked)
        
        self.spinner_html = widgets.HTML(value=spinning_loader)
        self.spinner_html.layout.display = "none"
        
        self.run_row = widgets.HBox(
            [self.btn_run, self.spinner_html],
            layout=widgets.Layout(align_items="center", gap="10px")
        )
        
        self.out_run = widgets.Output()
        
        # ---- Output info ----
        self.output_info = widgets.HTML(
            '<div style="font-size:12px;color:#666;margin-top:10px;">'
            '<b>Generates for each selected cell type:</b><br>'
            '• PDF with all plots (cumulative overall, per-sample, alive vs dead overall, alive vs dead per-sample)<br>'
            '• CSV with summary statistics per sample'
            '</div>'
        )
        
        # ---- UI ----
        self.ui = widgets.VBox([
            self.section_title,
            widgets.HBox([self.status_html, self.btn_refresh], layout=widgets.Layout(align_items="center", gap="10px")),
            self.cell_type_selector_label,
            self.cell_type_box,
            self.death_info_html,
            self.output_info,
            widgets.HTML("<hr>"),
            self.run_row,
            self.out_run,
        ])
        
        # Initial refresh to populate available cell types
        self._refresh_data_status()
    
    def _load_death_threshold_from_config(self):
        """Load death threshold from Step 3's config if available"""
        self.dead_threshold = 0.02  # Default
        self.death_dynamics_ran = False
        
        try:
            params = dict(self.metadata_loader.behav3d_parameters or {})
            death_cfg = params.get("death_dynamics", {}).get(self.cell_type, {})
            if "dead_perc_threshold" in death_cfg:
                self.dead_threshold = float(death_cfg["dead_perc_threshold"])
                self.death_dynamics_ran = True
        except Exception:
            pass
    
    def _refresh_data_status(self):
        """
        Check for available filtered data and update the UI accordingly.
        This can be called at any time to refresh after Step 1 completes.
        """
        import pandas as pd
        
        # Reload death threshold from config
        self._load_death_threshold_from_config()
        
        self.has_data = self.df_tracks_path.exists()
        self.has_contact_columns = False
        self.contact_cell_types = []
        self.has_dead_column = False
        
        if self.has_data:
            try:
                df_sample = pd.read_csv(self.df_tracks_path, nrows=0)
                # Find which cell types have contact columns
                for ct in self.available_cell_types:
                    if f"{ct}_contact" in df_sample.columns:
                        self.contact_cell_types.append(ct)
                self.has_contact_columns = len(self.contact_cell_types) > 0
                
                # Check for death data
                self.has_dead_column = "dead" in df_sample.columns or "percentage_dead_mask" in df_sample.columns
            except Exception:
                pass
        
        # Update status message
        if not self.has_data:
            self.status_html.value = (
                '<div style="color:#b00;font-size:14px;padding:10px;background:#fee;border:1px solid #fcc;border-radius:4px;">'
                '<b>⚠️ Waiting for filtered data</b><br>'
                f'No filtered track features found for {self.cell_type}.<br>'
                'Run Step 1 (Track Filtering) first, then click Refresh.'
                '</div>'
            )
            self._set_controls_enabled(False)
        elif not self.has_contact_columns:
            self.status_html.value = (
                '<div style="color:#b00;font-size:14px;padding:10px;background:#fee;border:1px solid #fcc;border-radius:4px;">'
                '<b>⚠️ No Contact Data Found</b><br>'
                f'No contact columns found for {self.cell_type}.<br>'
                'Ensure contact features were calculated during feature extraction.'
                '</div>'
            )
            self._set_controls_enabled(False)
        else:
            self.status_html.value = (
                '<div style="color:#080;font-size:14px;padding:10px;background:#efe;border:1px solid #cfc;border-radius:4px;">'
                f'<b>✅ Ready</b> - Filtered data found for {self.cell_type}<br>'
                f'Available contact types: {", ".join(self.contact_cell_types)}'
                '</div>'
            )
            self._set_controls_enabled(True)
            self._update_cell_type_checkboxes()
        
        # Update death info message
        self._update_death_info()
    
    def _update_death_info(self):
        """Update the death threshold info display"""
        if self.has_dead_column:
            if self.death_dynamics_ran:
                self.death_info_html.value = (
                    '<div style="font-size:13px;padding:8px;background:#e3f2fd;border:1px solid #90caf9;border-radius:4px;margin-top:10px;">'
                    f'<b>📊 Death classification:</b> Using threshold from Step 3 (Dead % threshold = {self.dead_threshold})'
                    '</div>'
                )
            else:
                self.death_info_html.value = (
                    '<div style="font-size:13px;padding:8px;background:#fff3e0;border:1px solid #ffcc80;border-radius:4px;margin-top:10px;">'
                    f'<b>⚠️ Step 3 not run:</b> Will calculate death using default threshold ({self.dead_threshold})<br>'
                    '<i>Tip: Run Step 3 (Death Dynamics) first for consistent death classification.</i>'
                    '</div>'
                )
        else:
            self.death_info_html.value = (
                '<div style="font-size:13px;padding:8px;background:#fce4ec;border:1px solid #f48fb1;border-radius:4px;margin-top:10px;">'
                '<b>ℹ️ No death data:</b> Alive vs Dead plots will be skipped.'
                '</div>'
            )
    
    def _update_cell_type_checkboxes(self):
        """Update the cell type checkboxes based on available contact columns"""
        self.cell_type_checkboxes = {}
        checkbox_widgets = []
        
        for ct in self.contact_cell_types:
            category = "immune" if ct in self.immune_types else "other"
            cb = widgets.Checkbox(
                value=True,
                description=f"{ct} ({category})",
                indent=False,
                layout=widgets.Layout(width="200px")
            )
            self.cell_type_checkboxes[ct] = cb
            checkbox_widgets.append(cb)
        
        self.cell_type_box.children = checkbox_widgets
    
    def _set_controls_enabled(self, enabled: bool):
        """Enable or disable all controls"""
        self.btn_run.disabled = not enabled
    
    def _on_refresh_clicked(self, *_):
        """Handle refresh button click"""
        with self.out_run:
            print("🔄 Refreshing data status...")
        self._refresh_data_status()
        with self.out_run:
            if self.has_data and self.has_contact_columns:
                print("   ✅ Data found! Ready to run.")
            else:
                print("   ⚠️ Data not yet available.")
    
    def _lock(self, locked):
        """Lock/unlock all controls during processing"""
        self.btn_run.disabled = locked
        self.btn_refresh.disabled = locked
        for cb in self.cell_type_checkboxes.values():
            cb.disabled = locked
    
    def _on_run_clicked(self, *_):
        """Run interaction analysis when button clicked"""
        
        # First refresh to ensure we have latest data
        self._refresh_data_status()
        
        if not self.has_data or not self.has_contact_columns:
            with self.out_run:
                print("❌ Cannot run - no valid data found. Please run Step 1 (Track Filtering) first.")
            return
        
        self._lock(True)
        self.out_run.clear_output()
        self.spinner_html.layout.display = None
        
        with self.out_run:
            try:
                # Get selected cell types
                selected_cell_types = [
                    ct for ct, cb in self.cell_type_checkboxes.items() if cb.value
                ]
                
                if not selected_cell_types:
                    print("❌ Please select at least one cell type to analyze.")
                    return
                
                # Use the threshold from Step 3 (Death Dynamics)
                dead_threshold = self.dead_threshold
                
                print(f"▶️ Running {self.cell_type} Interaction Analysis...")
                print(f"   Analyzing interactions with: {', '.join(selected_cell_types)}")
                if self.has_dead_column:
                    if self.death_dynamics_ran:
                        print(f"   Using death threshold from Step 3: {dead_threshold}")
                    else:
                        print(f"   Using default death threshold: {dead_threshold}")
                else:
                    print(f"   No death data - skipping alive vs dead comparisons")
                print()
                
                # Create output directory
                results_dir = Path(self.output_dir, "analysis", self.cell_type, "interaction_analysis")
                results_dir.mkdir(parents=True, exist_ok=True)
                
                # Run the analysis - generates SEPARATE PDF per cell type
                # All plots are generated (no selection needed)
                self._run_interaction_analysis(
                    selected_cell_types=selected_cell_types,
                    dead_threshold=dead_threshold,
                    results_dir=results_dir
                )
                
            except Exception:
                print("❌ Error during Interaction Analysis:")
                traceback.print_exc()
            finally:
                self.spinner_html.layout.display = "none"
                self._lock(False)
    
    def _run_interaction_analysis(
        self,
        selected_cell_types: list,
        dead_threshold: float,
        results_dir: Path
    ):
        """
        Run the interaction analysis using the analysis module.
        Thin wrapper that calls behav3d.analysis.interaction_analysis.run_interaction_analysis()
        """
        from behav3d.analysis.interaction_analysis import run_interaction_analysis
        
        run_interaction_analysis(
            output_dir=self.output_dir,
            cell_type=self.cell_type,
            interacting_cell_types=selected_cell_types,
            dead_threshold=dead_threshold,
            df_tracks_path=str(self.df_tracks_path),
            show_plots=True,
        )


class BackprojectionPanel:
    def __init__(self, metadata_loader, cell_type=None):
        self.metadata_loader = metadata_loader
        self.output_dir = str(Path(self.metadata_loader.output_dir).expanduser())

        # ---- Detect ALL cell types from metadata ----
        from behav3d.utils import (
            detect_immune_cell_types_from_metadata,
            detect_organoid_types_from_metadata,
            detect_other_cell_types_from_metadata
        )
        
        md = self.metadata_loader.metadata
        if md is None:
            raise RuntimeError("metadata_loader.metadata must be loaded before creating BackprojectionPanel.")
        
        # Detect all cell types
        immune_types = detect_immune_cell_types_from_metadata(md)
        organoid_types = detect_organoid_types_from_metadata(md)
        other_types = detect_other_cell_types_from_metadata(md)
        
        self.all_cell_types = immune_types + organoid_types + other_types
        if not self.all_cell_types:
            raise RuntimeError("No cell types detected in metadata. Ensure metadata has columns with prefixes: im_, or_, ot_")
        
        # Set default cell type
        if cell_type is None:
            cell_type = self.all_cell_types[0]
        
        # ---- load groups ----
        try:
            from copy import deepcopy
            from behav3d.utils import has_dead_channel
            self._base_groups = deepcopy(behav3d_calculated_features)  # dict[group] -> list[str/pattern]
            
            # Check if dead channel exists and remove death features if not
            has_dead = has_dead_channel(md)
            if not has_dead:
                if "death" in self._base_groups:
                    del self._base_groups["death"]
                if "intensity" in self._base_groups:
                    self._base_groups["intensity"] = [f for f in self._base_groups["intensity"] if f != "mean_dead_dye"]
            
            # Dynamically add contact features for ALL detected cell types
            for ct in self.all_cell_types:
                if f"{ct}_contact" not in self._base_groups["contact"]:
                    self._base_groups["contact"].extend([
                        f"{ct}_contact",
                        f"{ct}_contact_pixels",
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

        # ---- cell type dropdown (dynamic for ALL detected cell types) ----
        # Create human-readable labels for cell types
        self._celltype_map = {ct.replace('_', ' ').title(): ct for ct in self.all_cell_types}
        
        # Get default from config or use first cell type
        inv_ct_map = {v: k for k, v in self._celltype_map.items()}
        ct_default = inv_ct_map.get(str(self._cfg.get("cell_type", cell_type)), list(self._celltype_map.keys())[0])
        
        self.celltype_dd = widgets.Dropdown(
            options=list(self._celltype_map.keys()),
            value=ct_default,
            description="Cell type",
            layout=widgets.Layout(width="360px")
        )
        self.celltype_dd.style.description_width = "80px"

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
        self.celltype_dd.observe(lambda ch: (self._on_celltype_changed(), self._update_status_for_mode()), names='value')

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
            self.celltype_dd,   # under sample (now dropdown instead of toggle buttons)
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
        cell_type = self._celltype_map[self.celltype_dd.value]
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

        # Load actual columns from CSV to expand patterns
        cell_type = self._celltype_map[self.celltype_dd.value]
        feature_outdir = Path(self.output_dir, "analysis", cell_type, "track_features")
        csv_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
        avail_cols = []
        if csv_path.exists():
            try:
                avail_cols = pd.read_csv(csv_path, nrows=1).columns.tolist()
            except Exception:
                avail_cols = []

        rows = {}
        boxes = []
        for group_name, feats in self._base_groups.items():
            # Expand patterns against available columns
            expanded_feats = []
            for f in feats:
                if any(ch in f for ch in "*?["):
                    matches = [c for c in avail_cols if fnmatch.fnmatchcase(str(c), f)]
                    expanded_feats.extend(matches if matches else [])
                else:
                    if f in avail_cols or not avail_cols:
                        expanded_feats.append(f)
            
            # De-duplicate while preserving order
            seen = set()
            expanded_feats = [x for x in expanded_feats if not (x in seen or seen.add(x))]
            
            if not expanded_feats:
                continue  # Skip empty groups
                
            child_cbs = [widgets.Checkbox(value=(f in preset_set), description=f, indent=True) for f in expanded_feats]
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

            rows[group_name] = {
                "group_cb": gcb, 
                "child_cbs": child_cbs, 
                "grid": grid, 
                "container": container,
                "group_handler": group_handler,
                "child_handler": child_handler
            }
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

        # Prepare target patterns per group: ALWAYS add mean_ prefix
        # The summarized CSV has mean_ prefix on ALL columns (e.g., mean_intensity_ch1 becomes mean_mean_intensity_ch1)
        def to_mean_pattern(name: str) -> str:
            return "mean_" + str(name)

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

            rows[gname] = {
                "group_cb": gcb, 
                "child_cbs": child_cbs, 
                "grid": grid, 
                "container": container,
                "group_handler": gh,
                "child_handler": ch
            }
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
        
        # Temporarily disable all handlers to prevent recursion
        for row in rows.values():
            # Unobserve group checkbox if handler exists
            if "group_handler" in row:
                try:
                    row["group_cb"].unobserve(row["group_handler"], names="value")
                except:
                    pass
            # Unobserve child checkboxes if handler exists
            if "child_handler" in row:
                for cb in row["child_cbs"]:
                    try:
                        cb.unobserve(row["child_handler"], names="value")
                    except:
                        pass
        
        # Set all values
        for row in rows.values():
            row["group_cb"].value = val
            for cb in row["child_cbs"]:
                cb.value = val
        
        # Re-observe handlers
        for row in rows.values():
            if "group_handler" in row:
                try:
                    row["group_cb"].observe(row["group_handler"], names="value")
                except:
                    pass
            if "child_handler" in row:
                for cb in row["child_cbs"]:
                    try:
                        cb.observe(row["child_handler"], names="value")
                    except:
                        pass

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
        if len(available_columns) == 0:
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
        for w in [self.sample_dd, self.celltype_dd, self.mode_tb, self.save_cb,
                  self.btn_select_all, self.btn_clear_all, self.btn_run, self.refresh_btn]:
            w.disabled = state
        for rows in (self._time_group_rows, self._mean_group_rows):
            for row in rows.values():
                row["group_cb"].disabled = state
                for cb in row["child_cbs"]:
                    cb.disabled = state

    def _persist(self, *, resolved=None):
        mode = self._mode_map[self.mode_tb.value]
        self._cfg["cell_type"] = self._celltype_map[self.celltype_dd.value]
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
                cell_type = self._celltype_map[self.celltype_dd.value]
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
                traceback.print_exc()
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