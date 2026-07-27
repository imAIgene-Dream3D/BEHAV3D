import os
import re
import threading
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
    detect_merged_cell_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    is_merged_celltype,
    multicolor_base_name,
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

class Debouncer:
    """Delays calling `func` until `wait` seconds pass with no new call.

    Each call cancels any pending invocation and reschedules, so rapid-fire
    calls (typing digit-by-digit, holding a spinner button) collapse into a
    single call once the value has settled.
    """
    def __init__(self, func, wait=1.0):
        self._func = func
        self._wait = wait
        self._timer = None
        self._lock = threading.Lock()

    def __call__(self, *args, **kwargs):
        with self._lock:
            if self._timer is not None:
                self._timer.cancel()
            self._timer = threading.Timer(self._wait, self._func, args=args, kwargs=kwargs)
            self._timer.daemon = True
            self._timer.start()


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
        "extra_channel_slots": 0,
        "labels_mode": "same_for_all",  # "same_for_all" or "per_sample"
        "channel_labels": {},  # {0: "organoid1", 1: "tcell", ...}
        "per_sample_channel_labels": {},  # {sample_name: {0: "organoid1", ...}, ...}
    },
    # Cellpose-SAM (cellpose v4, zero-shot). Channel selection is its own,
    # per-cell-type config here (channel_selection below) - separate from
    # classic Cellpose's per-channel "cellpose" block above, since Cellpose-SAM
    # supports selecting several raw channels per cell type.
    "cellpose_sam": {
        "model_name": "cpsam",
        "python_path": "",       # optional override for the sidecar interpreter
        "device": "auto",        # "auto" | "cpu" | "cuda:<n>"
        "force_cpu": False,
        "all_cell_types": False,
        "use_all_timepoints": True,
        "tp_start": 0,
        "tp_end": 0,
        "preview_sample": "",
        "preview_timepoint": 0,
        # Mirrors the official Cellpose GUI's segmentation settings, except
        # min_size/max_size_fraction - see DEFAULT_SAM_PARAMS in
        # behav3d.preprocessing.segmentation.cellpose_sam_prediction for why
        # those two deviate from cellpose's own defaults.
        "diameter": 0.0,             # 0 => auto
        "flow_threshold": 0.4,       # ignored by cellpose when do_3D is True
        "cellprob_threshold": 0.0,
        "niter": 0,
        "do_3D": True,
        "stitch_threshold": 0.0,
        "flow3D_smooth": 0.0,
        "batch_size": 8,
        "min_size": 1,
        "max_size_fraction": 1.0,
        "norm_percentile_low": 1.0,
        "norm_percentile_high": 99.0,
        "norm3D": True,
        "sharpen_radius": 0.0,
        "smooth_radius": 0.0,
        "tile_norm_blocksize": 0.0,
        "tile_norm_smooth3D": 0.0,
        "invert": False,
        # BEHAV3D post-process: drop labels occupying a single Z/Y/X slice.
        # Cellpose leaves flat fragments on the first/last Z slice, most often in
        # 2D+stitch mode.
        "drop_2d_segments": True,
        # Interactive size filter, native voxels, per cell type:
        # {"tcell": {"size_min": 10, "size_max": 0}}  (0 => bound disabled)
        "size_filter": {},
        # Which raw channels feed Cellpose-SAM for each cell type, uniform
        # across samples: {"tcell": [0], "organoid1": [1], "organoid2": [2, 3]}
        "channel_selection": {},
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
                "use_visual_features": False,
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
                "use_visual_features": False,
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
                "use_visual_features": False,
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
            "feature_scaling_preset": None,
            "min_track_length_enabled": False,
            "min_track_length": 1,
            "max_track_length_enabled": False,
            "max_track_length": 999999,
            "feature_dtw_napari_backprojection_workers": 4,
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
            "feature_scaling_preset": None,
            "min_track_length_enabled": False,
            "min_track_length": 1,
            "max_track_length_enabled": False,
            "max_track_length": 999999,
            "feature_dtw_napari_backprojection_workers": 4,
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
            "feature_scaling_preset": None,
            "min_track_length_enabled": False,
            "min_track_length": 1,
            "max_track_length_enabled": False,
            "max_track_length": 999999,
            "feature_dtw_napari_backprojection_workers": 4,
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
    },
    "behavioral_state_classification": {
        "defaults": {
            "hmm_n_states_mode": "fixed",
            "hmm_n_states": 4,
            "hmm_k_min": 2,
            "hmm_k_max": 8,
            "hmm_covariance_type": "full",
            "hmm_n_iter": 500,
            "hmm_tol": 0.001,
            "hmm_sticky": False,
            "hmm_stickiness_kappa": 8.0,
            "hmm_transmat_alpha": 1.0,
            "hmm_min_covar": 0.001,
            "hmm_feature_smoothing_window": 5,
            "hmm_smoothing_min_periods": 1,
            "hmm_start_offset": 1,
            "hmm_start_offset_fill_mode": "backfill",
            "hmm_backprojection_workers": 4,
            "hmm_window_features": ["net_displacement"],
            "hmm_window_features_window": 5,
            "hmm_log_scale_features": ["speed", "net_displacement"],
            "feature_quantile_capping_low_percentile": "",
            "feature_quantile_capping_high_percentile": "",
            "apply_hmm_deployment_artifact_path": "",
            "selected_features": ["elongation", "extent", "solidity", "sphericity", "speed"],
            "binary_features_to_group": ["any_immune_cell_contact", "any_organoid_contact", "dead", "is_active_killing"],
            "random_state": 12345,
        },
    },
    "behavioral_track_classification": {
        "defaults": {
            "behavioral_trajectory_size": 100,
            "n_clusters": 6,
            "n_per_cluster": 10,
            "random_state": 12345,
            "trajectory_trim_mode": "last",
            "linkage": "average",
            "missing_policy": "keep",
            "parallel": True,
            "save_distance_matrix": False,
            "window": "",
            "max_dist": "",
            "max_length_diff": "",
            "penalty": "",
            "psi": "",
            "make_overview_statebars": True,
            "make_backprojection_pdf": False,
            "make_backprojection_mp4": False,
            "classifier_n_estimators": 100,
            "classifier_min_samples_leaf": 1,
            "classifier_min_samples_split": 2,
            "classifier_max_features": "sqrt",
            "classifier_max_depth": "",
            "classifier_test_size": 0.2,
            "classifier_n_jobs": -1,
            "classifier_artifact_path": "",
        },
    },
    "viewer_display": {
        "channels": {},
    },
}

behav3d_calculated_features = {
    "morphology": [
        "nr_pixels", "volume", "bbox_volume", "extent", "solidity",
        "equivalent_diameter", "major_axis_length", "minor_axis_length",
        "axis1_length", "axis2_length", "axis3_length",
        "elongation", "surface_area", "sphericity", "convex_volume", "surface_to_volume_ratio",
    ],
    "movement": [
        "displacement", "cumulative_displacement", "displacement_from_origin",
        "mean_square_displacement", "speed", "mean_speed", "summed_displacement",
        "net_displacement", "straightness", "directional_persistence",
        "median_turning_angle", "fraction_reversed_movement",
    ],
    "intensity": [
        "mean_intensity_*", "mean_dead_dye", "q75_mean_intensity_*", "q90_mean_intensity_*",
    ],
    "death": [
        "percentage_dead_mask", "nr_dead_mask_pixels", "increase_dead_mask", "dead",
    ],
    "contact": [
        "*_contact", "*_contact_on_distance", "*_contact_pixels", "active_*_contact",
    ],
    "invasiveness": [
        "*_invasiveness", "*_invasiveness_perc", "any_org_invasiveness", "any_org_invasiveness_perc",
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
# ROW REORDER BUTTONS
# ===============================
def build_row_move_buttons(on_move_up, on_move_down):
    """Small stacked up/down buttons for reordering a row within a list of rows."""
    up = widgets.Button(
        icon="chevron-up", tooltip="Move up",
        layout=widgets.Layout(width="26px", height="17px", padding="0px"),
    )
    down = widgets.Button(
        icon="chevron-down", tooltip="Move down",
        layout=widgets.Layout(width="26px", height="17px", padding="0px"),
    )
    up.on_click(lambda _btn: on_move_up())
    down.on_click(lambda _btn: on_move_down())
    return widgets.VBox([up, down], layout=widgets.Layout(gap="0px"))


# ===============================
# PLOT ACTION BOX
# ===============================
def build_plot_box(title, description, run_row, settings=None, extra=None, output=None):
    """Build a collapsible plot box in the same visual style as the workflow-step accordions.

    Parameters
    ----------
    title : str
        Plot name, shown as the accordion header.
    description : str
        Plain-text "what does this plot?" explanation, shown once the box is expanded.
    run_row : widgets.Widget
        The always-visible run control (typically ``HBox([button, spinner])``).
    settings : list[widgets.Widget] or None
        Parameter widgets for this plot. Hidden behind a "Settings" toggle until clicked,
        rather than always shown inline.
    extra : list[widgets.Widget] or None
        Additional always-visible widgets (e.g. status/warning text) shown between the
        description and the run row.
    output : widgets.Output or None
        Output area appended at the end of the box.

    Returns
    -------
    widgets.Accordion
        A single-section, initially-collapsed accordion titled ``title``.
    """
    children = [widgets.HTML(f"<span style='color:#555;'>{description}</span>")]
    children.extend(extra or [])

    if settings:
        settings_toggle = widgets.ToggleButton(
            value=False,
            description="Settings",
            icon="cog",
            layout=widgets.Layout(width="110px"),
        )
        settings_panel = widgets.VBox(
            list(settings),
            layout=widgets.Layout(display="none", margin="6px 0 0 0"),
        )

        def _on_settings_toggle(change, panel=settings_panel):
            panel.layout.display = "" if change["new"] else "none"

        settings_toggle.observe(_on_settings_toggle, names="value")
        children.append(
            widgets.HBox([run_row, settings_toggle], layout=widgets.Layout(align_items="center", gap="8px"))
        )
        children.append(settings_panel)
    else:
        children.append(run_row)

    if output is not None:
        children.append(output)

    box = widgets.VBox(children)
    acc = widgets.Accordion(children=[box], selected_index=None)
    acc.set_title(0, title)
    return acc


# ===============================
# HELPER FUNCTIONS
# ===============================
def build_feature_distribution_figure(csv_path, feats, title=None, max_plots=36):
    """Build a matplotlib grid of per-feature histograms from a CSV.

    Shared by the napari plugin and the notebook log-scaling panels so the
    user can visually judge which features are skewed (and would benefit
    from log scaling) before applying it. Returns ``(figure, truncated)``
    where ``truncated`` is True when ``feats`` was capped at ``max_plots``.
    Returns ``(None, False)`` when there is nothing to plot.
    """
    from matplotlib.figure import Figure

    feats = [f for f in (feats or [])]
    if csv_path is None or not Path(csv_path).exists() or not feats:
        return None, False

    truncated = len(feats) > max_plots
    feats = feats[:max_plots]
    try:
        df = pd.read_csv(csv_path, usecols=feats, low_memory=False)
    except Exception:
        return None, False

    n = len(feats)
    ncols = min(4, n) or 1
    nrows = (n + ncols - 1) // ncols
    fig = Figure(figsize=(3.2 * ncols, 2.4 * nrows), tight_layout=True)
    for i, feat in enumerate(feats):
        ax = fig.add_subplot(nrows, ncols, i + 1)
        series = pd.to_numeric(df[feat], errors="coerce").dropna()
        if len(series):
            ax.hist(series.values, bins=40, color="#4C72B0", edgecolor="none")
        else:
            ax.text(0.5, 0.5, "no data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=8, color="#999")
        ax.set_title(str(feat), fontsize=8)
        ax.tick_params(labelsize=7)
    if title:
        fig.suptitle(str(title), fontsize=11)
    return fig, truncated


def build_track_length_distribution_figure(
    csv_path, group_col="TrackID", sample_col="sample_name", title=None,
):
    """Build a per-sample grid of track-length histograms from a track CSV.

    Shared by the Filtering tab's "See Track Length Distributions" preview
    and the cell-type Grouping dialog (so a freshly merged population can be
    sanity-checked the same way). Returns a matplotlib ``Figure``, or
    ``None`` if the CSV is missing/unreadable or lacks the required columns.
    """
    from matplotlib.figure import Figure

    if csv_path is None or not Path(csv_path).exists():
        return None
    try:
        df = pd.read_csv(csv_path, usecols=lambda c: c in {group_col, sample_col}, low_memory=False)
    except Exception:
        return None
    if group_col not in df.columns or sample_col not in df.columns:
        return None

    samples = sorted(df[sample_col].dropna().unique())
    if not samples:
        return None

    ncols = min(3, len(samples))
    nrows = (len(samples) + ncols - 1) // ncols
    fig = Figure(figsize=(5 * ncols, 4 * nrows), tight_layout=True)
    for i, sample in enumerate(samples):
        ax = fig.add_subplot(nrows, ncols, i + 1)
        lengths = df[df[sample_col] == sample].groupby(group_col).size()
        ax.hist(lengths.values, bins=50, color="skyblue", edgecolor="black")
        ax.set_title(str(sample))
        ax.set_xlabel("Track Length (timepoints)")
        ax.set_ylabel("Frequency")
    if title:
        fig.suptitle(str(title), fontsize=12)
    return fig


class UnitGroupManager:
    """ipywidgets counterpart of ``behav3d.napari._units.UnitGroupManager``.

    Manages the physical(µm)/pixel display of a group of ipywidgets
    number inputs (``BoundedFloatText``/``BoundedIntText``/``IntText``/
    ``FloatText``). Persistence stays native: :meth:`get_native` always
    returns the canonical pixel/voxel (or µm, for parameters whose native
    unit is already physical, e.g. tracking distances) value regardless of
    the currently displayed unit, so config round-trips unchanged.

    Conversions use :func:`behav3d.core.utils.resolution_from_metadata` /
    ``native_to_physical`` / ``physical_to_native``.
    """

    def __init__(self, metadata, default_physical=True):
        from behav3d.core.utils import resolution_from_metadata

        self.xy, self.z, self.valid = resolution_from_metadata(metadata)
        self.physical = bool(default_physical and self.valid)
        self._entries = []  # list of (widget, kind, native_unit)

        self.toggle = widgets.ToggleButtons(
            options=[("px", False), ("µm", True)],
            value=self.physical,
            style={"button_width": "50px"},
        )
        if not self.valid:
            self.toggle.disabled = True
        self.toggle.observe(self._on_toggle, names="value")

    # -- widget construction ---------------------------------------------
    def header_row(self, label="Units:"):
        """Return a small ``HBox`` row: ``label [px|µm]``."""
        children = []
        if label:
            children.append(widgets.HTML(f"<b>{label}</b>"))
        children.append(self.toggle)
        if not self.valid:
            children.append(widgets.HTML(
                "<i style='color:#999;font-size:10px;'>(no resolution in metadata)</i>"
            ))
        return widgets.HBox(children)

    # -- registration -------------------------------------------------------
    def register(self, widget, kind, native_value, native_unit="pixel"):
        """Register ``widget`` (a numeric ipywidgets input) holding a value
        of ``kind`` (``"distance"`` or ``"volume"``).

        ``native_value`` is the value in the parameter's *native* unit
        (what gets persisted/passed to the backend). ``native_unit`` is
        ``"pixel"`` (segmentation defaults) or ``"physical"`` (tracking,
        whose linking already runs in µm).
        """
        widget._unit_kind = kind
        widget._unit_native_unit = native_unit
        widget._unit_native = float(native_value)
        # Store the exact bound handler so later unobserve/observe calls
        # only ever touch this manager's own listener, never any other
        # "value" observer the panel itself may have registered.
        widget._unit_handler = lambda change, w=widget: self._on_edit(w)
        self._set_display(widget)
        self._update_description_suffix(widget)
        widget.observe(widget._unit_handler, names="value")
        self._entries.append(widget)

    def get_native(self, widget):
        """Return the canonical native value for ``widget``."""
        return float(getattr(widget, "_unit_native", widget.value))

    def set_native(self, widget, native_value):
        """Set ``widget``'s canonical native value and refresh its display."""
        widget._unit_native = float(native_value)
        self._set_display(widget)

    # -- internal -------------------------------------------------------
    def _display_to_native(self, value, kind, native_unit):
        from behav3d.core.utils import native_to_physical, physical_to_native

        if not self.valid:
            return float(value)
        if native_unit == "physical":
            if self.physical:
                return float(value)
            return native_to_physical(value, kind, self.xy, self.z)
        if self.physical:
            return physical_to_native(value, kind, self.xy, self.z)
        return float(value)

    def _native_to_display(self, native, kind, native_unit):
        from behav3d.core.utils import native_to_physical, physical_to_native

        if not self.valid:
            return float(native)
        if native_unit == "physical":
            if self.physical:
                return float(native)
            return physical_to_native(native, kind, self.xy, self.z)
        if self.physical:
            return native_to_physical(native, kind, self.xy, self.z)
        return float(native)

    def _on_edit(self, widget):
        widget._unit_native = self._display_to_native(
            widget.value, widget._unit_kind, widget._unit_native_unit
        )

    def _set_display(self, widget):
        disp = self._native_to_display(
            widget._unit_native, widget._unit_kind, widget._unit_native_unit
        )
        handler = getattr(widget, "_unit_handler", None)
        if handler is not None:
            widget.unobserve(handler, names="value")
        if isinstance(widget, (widgets.IntText, widgets.BoundedIntText)):
            widget.value = int(round(disp))
        else:
            widget.value = float(disp)
        if handler is not None:
            widget.observe(handler, names="value")

    def _update_description_suffix(self, widget):
        unit_label = "µm" if self.physical else ("px" if widget._unit_kind == "distance" else "vox")
        base = getattr(widget, "_unit_base_description", None)
        if base is None:
            base = widget.description
            widget._unit_base_description = base
        widget.description = f"{base} ({unit_label})" if base else f"({unit_label})"

    def _on_toggle(self, change):
        self.physical = bool(change["new"]) and self.valid
        for w in self._entries:
            self._set_display(w)
            self._update_description_suffix(w)


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
    merged_types = detect_merged_cell_types_from_metadata(metadata)

    if cell_type in merged_types:
        base_name = multicolor_base_name(cell_type)
        for candidate in organoid_types:
            if multicolor_base_name(candidate) == base_name:
                return 'organoid'
        for candidate in immune_types:
            if multicolor_base_name(candidate) == base_name:
                return 'immune'
        for candidate in other_types:
            if multicolor_base_name(candidate) == base_name:
                return 'other'
        return 'immune'
    
    if cell_type in organoid_types:
        return 'organoid'
    elif cell_type in immune_types:
        return 'immune'
    elif cell_type in other_types:
        return 'other'
    base_name = multicolor_base_name(cell_type)
    if base_name != cell_type:
        for candidate in organoid_types:
            if multicolor_base_name(candidate) == base_name:
                return 'organoid'
        for candidate in immune_types:
            if multicolor_base_name(candidate) == base_name:
                return 'immune'
        for candidate in other_types:
            if multicolor_base_name(candidate) == base_name:
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
        allow_zarr=False,
    ):
        super().__init__()
        self._mode = mode
        self._start_dir = start_dir
        self.allow_zarr = allow_zarr
        self.filter_pattern = filter_pattern
        
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
        
        if self.allow_zarr:
            self.toggle_zarr = widgets.ToggleButton(
                value=(self._mode == 'dir'),
                description='Zarr/Dir',
                tooltip='Toggle between File and Zarr/Dir selection mode',
                layout=widgets.Layout(width='auto', flex='0 0 auto')
            )
            self.toggle_zarr.observe(self._on_toggle_zarr, names='value')
            self.children = [self.text, self.btn, self.toggle_zarr, self.chooser_out]
        else:
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
                    filter_pattern=getattr(self, 'filter_pattern', None) if self._mode == 'file' else None
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

    def _on_toggle_zarr(self, change):
        self._mode = 'dir' if change['new'] else 'file'
        if self.chooser is not None:
            self.chooser_out.clear_output()
            self.chooser = FileChooser(
                self._start_dir,
                select_default=True,
                show_only_dirs=(self._mode == 'dir'),
                filter_pattern=getattr(self, 'filter_pattern', None) if self._mode == 'file' else None
            )
            self.chooser.register_callback(self._on_select)
            with self.chooser_out:
                display(self.chooser)

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
            description="File:", placeholder="Path to image (TIFF, TIF)",
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

        # External segmentation/tracking imports must be TIFF only.
        # (They are converted into per-sample .zarr and recorded in metadata.)
        if self._is_label_image and p.suffix.lower() not in (".tif", ".tiff"):
            self.info.value = (
                "<span style='color:red'>"
                f"Unsupported file type for external label import: {p.suffix}. "
                "Please provide a .tif or .tiff file."
                "</span>"
            )
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

        if self._is_label_image and src.suffix.lower() not in (".tif", ".tiff"):
            with self.out:
                print(
                    f"Unsupported file type for external label import: {src.suffix}. "
                    "Please provide a .tif or .tiff file."
                )
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
