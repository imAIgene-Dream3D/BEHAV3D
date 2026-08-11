from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import numpy as np
import pandas as pd
import scanpy as sc
import napari

from behav3d.analysis.behavior.state.visualization.backprojection import (
    _add_mapping_dock_widget,
    _align_labels_to_raw_shape_for_view,
    _build_state_mapping_text,
    _ensure_behavioral_state_backprojection_for_sample,
    _extract_state_code_color_map,
    _extract_state_label_map,
    _resolve_behavioral_state_image_path,
    _resolve_behavioral_states_h5ad_path,
    _resolve_raw_image_path,
    _resolve_tracked_image_path,
)
from behav3d.analysis.behavior.track.visualization.backprojection import (
    add_track_cluster_trajectory_layers,
)
from behav3d.analysis.behavior.track.visualization.plots.exemplar_coordinate_utils import (
    merge_coordinate_columns_into_obs,
    resolve_exemplar_positions_csv_path,
)
from behav3d.io.images import load_image

# NOTE: this script is intentionally self-contained rather than importing from
# test/show_behavioral_state_backprojection_split_channels.py. That script's bottom
# "# %%" interactive cell runs unconditionally at module import time (no
# `if __name__ == "__main__":` guard, by design, so it can be run cell-by-cell from
# an IDE) — importing functions from it would also execute that cell immediately.
# The raw/tracked/state resolution and split-channel layer setup below therefore
# mirror that script's `_resolve_inputs`, `build_split_channel_backprojection_payload`,
# and `launch_split_channel_backprojection_viewer` rather than reusing them directly.


# %%
# Interactive configuration
# Fill these in when running the script line by line from the IDE.
# Point this to the BEHAV3D results folder that contains `images/`, `analysis/`,
# and typically `metadata.csv`.
INTERACTIVE_BEHAV3D_FOLDER = r"F:\BHVD_BEHAV3D\BEHAV3D_python\runs\NatureBriefComm\LowDensity\behav3d"
# Optional explicit override if your BEHAV3D output directory differs.
INTERACTIVE_SAMPLE_NAME = "BHVD_SB1_Exp009_Img001"
INTERACTIVE_OUTPUT_DIR = None
INTERACTIVE_CELL_TYPE = "tcell"
INTERACTIVE_STATE_COL = "full_behavioral_cluster"
INTERACTIVE_RAW_IMAGE_PATH = None
INTERACTIVE_TRACKED_IMG_PATH = None
INTERACTIVE_STATE_IMG_PATH = None
INTERACTIVE_AUTO_CREATE_IF_MISSING = True
INTERACTIVE_REFRESH_IF_STALE = True
INTERACTIVE_VERBOSE = True
# Trajectory-specific options.
INTERACTIVE_SHOW_TRAJECTORIES = True
# Drop same-state runs shorter than this many timepoints (1 = keep every run,
# including single-frame ones).
INTERACTIVE_MIN_RUN_LENGTH = 1
# Passed straight to viewer.add_tracks; None lets each trajectory layer default
# to that state's own full time span so it never fades out.
INTERACTIVE_TAIL_LENGTH = None


def resolve_behav3d_output_dir(
    behav3d_folder: str | Path | None = None,
    output_dir: str | Path | None = None,
) -> Path:
    """Resolve the BEHAV3D output folder used by the backprojection helpers."""
    candidate = output_dir if output_dir not in {None, ""} else behav3d_folder
    if candidate in {None, ""}:
        raise ValueError(
            "Set INTERACTIVE_BEHAV3D_FOLDER to your BEHAV3D results folder "
            "or provide INTERACTIVE_OUTPUT_DIR explicitly."
        )

    resolved = Path(candidate).expanduser()
    if not resolved.exists():
        raise FileNotFoundError(f"BEHAV3D folder does not exist: '{resolved}'")

    expected_images = resolved / "images"
    expected_analysis = resolved / "analysis"
    if not expected_images.exists() or not expected_analysis.exists():
        raise FileNotFoundError(
            "Expected a BEHAV3D output folder containing 'images/' and 'analysis/', "
            f"but got '{resolved}'."
        )
    return resolved


def _resolve_inputs(
    sample_name: str,
    output_dir: Path,
    cell_type: str,
    state_col: str,
    raw_image_path: str | None,
    tracked_img_path: str | None,
    state_img_path: str | None,
    auto_create_if_missing: bool,
    refresh_if_stale: bool,
    verbose: bool,
) -> tuple[Path, Path, Path]:
    raw_path = Path(raw_image_path) if raw_image_path is not None else _resolve_raw_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        verbose=verbose,
    )
    if raw_path is None or not Path(raw_path).exists():
        raise FileNotFoundError(
            "Could not find raw image for sample "
            f"'{sample_name}'. Expected '{Path(output_dir, 'images', sample_name, f'{sample_name}.zarr')}' "
            "or '.zarr.zip'."
        )

    tracked_path = Path(tracked_img_path) if tracked_img_path is not None else _resolve_tracked_image_path(
        output_dir=output_dir,
        sample_name=sample_name,
        cell_type=cell_type,
        verbose=verbose,
    )
    if tracked_path is None or not Path(tracked_path).exists():
        raise FileNotFoundError(
            "Could not find tracked image for sample "
            f"'{sample_name}' and cell_type '{cell_type}'."
        )

    if state_img_path is not None:
        state_path = Path(state_img_path)
        if not state_path.exists():
            if auto_create_if_missing:
                state_path = _ensure_behavioral_state_backprojection_for_sample(
                    sample_name=sample_name,
                    output_dir=output_dir,
                    cell_type=cell_type,
                    state_col=state_col,
                    track_col="TrackID",
                    time_col="position_t",
                    sample_col="sample_name",
                    background_value=0,
                    enforce_time_coverage=True,
                    refresh_if_stale=refresh_if_stale,
                    verbose=verbose,
                )
            else:
                raise FileNotFoundError(
                    "Could not find behavioral-state image for sample "
                    f"'{sample_name}' and cell_type '{cell_type}'."
                )
    elif auto_create_if_missing:
        state_path = _ensure_behavioral_state_backprojection_for_sample(
            sample_name=sample_name,
            output_dir=output_dir,
            cell_type=cell_type,
            state_col=state_col,
            track_col="TrackID",
            time_col="position_t",
            sample_col="sample_name",
            background_value=0,
            enforce_time_coverage=True,
            refresh_if_stale=refresh_if_stale,
            verbose=verbose,
        )
    else:
        state_path = _resolve_behavioral_state_image_path(
            output_dir=output_dir,
            sample_name=sample_name,
            cell_type=cell_type,
            verbose=verbose,
        )
        if state_path is None or not Path(state_path).exists():
            raise FileNotFoundError(
                "Could not find behavioral-state image for sample "
                f"'{sample_name}' and cell_type '{cell_type}'."
            )

    return Path(raw_path), Path(tracked_path), Path(state_path)


# %%
# Pixel-coordinate + h5ad loading helpers
#

_PIXEL_POSITION_TRIPLETS = (
    ("pixel_position_z", "pixel_position_y", "pixel_position_x"),
    ("pixel_position_y", "pixel_position_x"),
)  # local copy — mirrors the constant in
# behav3d.analysis.behavior.track.visualization.backprojection; not imported from
# there to keep this script self-contained.


def _ensure_state_trajectory_pixel_coordinates(adata_full, output_dir: Path, cell_type: str):
    """Populate pixel coordinates for state trajectories, allowing either 3D or 2D
    pixel positions. Mirrors _ensure_trajectory_pixel_coordinates in
    test/show_track_cluster_backprojection_split_channels.py."""
    obs_cols = set(adata_full.obs.columns)
    pixel_triplet_3d = set(_PIXEL_POSITION_TRIPLETS[0])
    pixel_pair_2d = set(_PIXEL_POSITION_TRIPLETS[1])
    if pixel_triplet_3d.issubset(obs_cols) or pixel_pair_2d.issubset(obs_cols):
        return {
            "enriched": False,
            "csv_path": None,
            "added_columns": [],
            "coordinate_mode": "3d" if pixel_triplet_3d.issubset(obs_cols) else "2d",
        }

    positions_csv_path = resolve_exemplar_positions_csv_path(output_dir=output_dir, cell_type=cell_type)
    merge_info = merge_coordinate_columns_into_obs(
        adata=adata_full,
        positions_csv_path=positions_csv_path,
        sample_col="sample_name",
        track_col="TrackID",
        time_col="position_t",
    )
    obs_cols = set(adata_full.obs.columns)
    if pixel_triplet_3d.issubset(obs_cols):
        coordinate_mode = "3d"
    elif pixel_pair_2d.issubset(obs_cols):
        coordinate_mode = "2d"
    else:
        raise ValueError(
            "State trajectory overlays require pixel coordinates in adata_full.obs. "
            "Expected either ['pixel_position_z', 'pixel_position_y', 'pixel_position_x'] "
            "or ['pixel_position_y', 'pixel_position_x']. "
            f"Enrichment source: '{positions_csv_path}'."
        )
    return {
        "enriched": True,
        "csv_path": str(positions_csv_path),
        "added_columns": list(merge_info.get("added_columns", [])),
        "coordinate_mode": coordinate_mode,
    }


def _load_sample_state_obs(output_dir: Path, cell_type: str, sample_name: str, state_col: str, verbose: bool = True):
    """Load the full per-frame behavioral-state h5ad and filter it to one sample."""
    h5ad_path = _resolve_behavioral_states_h5ad_path(output_dir=output_dir, cell_type=cell_type)
    if not h5ad_path.exists():
        raise FileNotFoundError(
            "Cannot build state trajectories because the full behavioral-state h5ad is missing: "
            f"'{h5ad_path}'."
        )
    adata_full = sc.read_h5ad(h5ad_path)
    if not hasattr(adata_full, "obs"):
        raise ValueError(f"Loaded h5ad has no .obs: '{h5ad_path}'.")
    if state_col not in adata_full.obs.columns:
        raise ValueError(
            f"state_col '{state_col}' not found in '{h5ad_path}'.obs columns: "
            f"{list(adata_full.obs.columns)}."
        )

    sample_obs = adata_full.obs[
        adata_full.obs["sample_name"].astype("string") == str(sample_name)
    ].copy()
    if len(sample_obs) == 0:
        raise ValueError(
            f"Cannot build state trajectories because sample '{sample_name}' is absent in '{h5ad_path}'."
        )
    if bool(verbose):
        print(f"Loaded {len(sample_obs)} obs rows for sample '{sample_name}' from '{h5ad_path.name}'.")
    return adata_full, sample_obs


# %%
# Contiguous same-state trajectory segmentation
#

def _assign_state_run_segment_ids(obs: pd.DataFrame, track_col: str, time_col: str, state_col: str) -> pd.DataFrame:
    """Sort by (track, time) and assign a globally-unique integer id to each
    maximal run of rows that share the same track, the same state_col value,
    and strictly consecutive time_col values. A state change, a track change,
    or a time gap all start a new segment."""
    obs = obs.sort_values([track_col, time_col], kind="mergesort").reset_index(drop=True)
    track_vals = obs[track_col].to_numpy()
    time_vals = obs[time_col].to_numpy()
    state_vals = obs[state_col].astype("string").to_numpy()

    continues = np.zeros(len(obs), dtype=bool)
    if len(obs) > 1:
        continues[1:] = (
            (track_vals[1:] == track_vals[:-1])
            & (state_vals[1:] == state_vals[:-1])
            & ((time_vals[1:] - time_vals[:-1]) == 1)
        )
    obs = obs.copy()
    obs["__state_run_segment_id"] = np.cumsum(~continues) - 1
    return obs


def prepare_state_trajectory_data(
    sample_obs: pd.DataFrame,
    state_col: str,
    track_col: str = "TrackID",
    time_col: str = "position_t",
    min_run_length: int = 1,
) -> dict:
    """
    Build one array per behavioral-state label of contiguous same-state
    trajectory segments, ready for one viewer.add_tracks(...) call per label
    (see add_track_cluster_trajectory_layers in
    behav3d.analysis.behavior.track.visualization.backprojection, which is
    generic and reused as-is for this).

    Segments are runs of consecutive position_t frames within a single
    TrackID that share the same state_col value; a state change, a track
    change, or a time gap all start a new segment. Requires pixel-space
    position columns (pixel_position_z/y/x, or pixel_position_y/x for 2D) in
    sample_obs, for the same reason prepare_track_cluster_trajectory_data
    does (raw/labels layers carry no scale= transform).

    Parameters
    ----------
    sample_obs : pd.DataFrame
        adata_full.obs already filtered to one sample.
    min_run_length : int
        Drop segments with fewer than this many timepoints (default 1 = keep
        every run, including single-frame ones).

    Returns
    -------
    dict[str, np.ndarray]
        state label -> array of shape [N, 4] or [N, 5] with columns
        [segment_id, position_t, (position_z), position_y, position_x].
    """
    pos_triplet = None
    for candidate in _PIXEL_POSITION_TRIPLETS:
        if all(c in sample_obs.columns for c in candidate):
            pos_triplet = candidate
            break
    if pos_triplet is None:
        raise ValueError(
            "sample_obs is missing pixel-space position columns needed to build state "
            "trajectory segments (expected 'pixel_position_z'/'pixel_position_y'/"
            "'pixel_position_x', or 'pixel_position_y'/'pixel_position_x' for 2D data)."
        )

    needed_cols = [str(track_col), str(time_col), str(state_col)] + list(pos_triplet)
    missing = [c for c in needed_cols if c not in sample_obs.columns]
    if len(missing) > 0:
        raise ValueError(f"sample_obs missing required columns: {missing}")

    obs = sample_obs[needed_cols].copy()
    obs[track_col] = pd.to_numeric(obs[track_col], errors="coerce")
    obs[time_col] = pd.to_numeric(obs[time_col], errors="coerce")
    obs[state_col] = obs[state_col].astype("string").str.strip()
    obs = obs.dropna(subset=[track_col, time_col] + list(pos_triplet))
    obs = obs[obs[state_col].notna() & (obs[state_col].str.len() > 0)]
    if len(obs) == 0:
        return {}
    obs[track_col] = obs[track_col].astype(np.int64)
    obs[time_col] = obs[time_col].astype(np.int64)

    obs = _assign_state_run_segment_ids(obs, track_col=track_col, time_col=time_col, state_col=state_col)

    if int(min_run_length) > 1:
        run_sizes = obs.groupby("__state_run_segment_id")[time_col].transform("size")
        obs = obs[run_sizes >= int(min_run_length)]
        if len(obs) == 0:
            return {}

    trajectory_data = {}
    cols = ["__state_run_segment_id", time_col] + list(pos_triplet)
    for label, group in obs.groupby(state_col, observed=True, sort=False):
        trajectory_data[str(label)] = group[cols].to_numpy(dtype=np.float64, copy=True)
    return trajectory_data


# %%
# Payload / viewer construction
#

def build_state_trajectory_backprojection_payload(
    sample_name: str,
    output_dir: str | Path,
    cell_type: str,
    state_col: str = "full_behavioral_cluster",
    raw_image_path: str | None = None,
    tracked_img_path: str | None = None,
    state_img_path: str | None = None,
    auto_create_if_missing: bool = True,
    refresh_if_stale: bool = True,
    show_trajectories: bool = True,
    min_run_length: int = 1,
    tail_length: int | None = None,
    verbose: bool = True,
):
    sample_name = str(sample_name).strip()
    output_dir = Path(output_dir)
    cell_type = str(cell_type).strip()

    if len(sample_name) == 0:
        raise ValueError("sample_name is required")
    if len(cell_type) == 0:
        raise ValueError("cell_type is required")

    raw_path, tracked_path, state_path = _resolve_inputs(
        sample_name=sample_name,
        output_dir=output_dir,
        cell_type=cell_type,
        state_col=state_col,
        raw_image_path=raw_image_path,
        tracked_img_path=tracked_img_path,
        state_img_path=state_img_path,
        auto_create_if_missing=bool(auto_create_if_missing),
        refresh_if_stale=bool(refresh_if_stale),
        verbose=bool(verbose),
    )

    raw_img = load_image(raw_path)
    raw_shape = tuple(int(v) for v in raw_img.shape)
    if len(raw_shape) != 5:
        raise ValueError(
            "Raw image must have shape (T, C, Z, Y, X), "
            f"but got shape {raw_shape} from '{raw_path}'."
        )
    if int(raw_shape[1]) != 3:
        raise ValueError(
            "Raw image must contain exactly 3 channels on axis 1, "
            f"but got shape {raw_shape} from '{raw_path}'."
        )

    raw_channels = [raw_img[:, ch, ...] for ch in range(3)]

    tracked_img = load_image(tracked_path)
    state_img = load_image(state_path)
    raw_reference = raw_channels[0]

    tracked_img_view = _align_labels_to_raw_shape_for_view(
        tracked_img,
        raw_reference,
        layer_name="TrackID",
        verbose=verbose,
    )
    state_img_view = _align_labels_to_raw_shape_for_view(
        state_img,
        raw_reference,
        layer_name="behavioral_state_class",
        verbose=verbose,
    )

    label_map = _extract_state_label_map(state_path)
    mapping_text = _build_state_mapping_text(label_map)

    trajectory_data = {}
    code_colors = {}
    coordinate_info = None
    if bool(show_trajectories):
        code_colors = _extract_state_code_color_map(state_path)
        adata_full, _initial_sample_obs = _load_sample_state_obs(
            output_dir=output_dir,
            cell_type=cell_type,
            sample_name=sample_name,
            state_col=state_col,
            verbose=verbose,
        )
        coordinate_info = _ensure_state_trajectory_pixel_coordinates(
            adata_full=adata_full,
            output_dir=output_dir,
            cell_type=cell_type,
        )
        # merge_coordinate_columns_into_obs enriches adata_full.obs in place, so the
        # pre-enrichment sample_obs filtered inside _load_sample_state_obs may be missing
        # the pixel columns; re-filter from the (possibly enriched) adata_full.obs instead
        # of reusing _initial_sample_obs.
        sample_obs = adata_full.obs[
            adata_full.obs["sample_name"].astype("string") == str(sample_name)
        ].copy()
        trajectory_data = prepare_state_trajectory_data(
            sample_obs,
            state_col=state_col,
            track_col="TrackID",
            time_col="position_t",
            min_run_length=min_run_length,
        )

    return {
        "sample_name": sample_name,
        "output_dir": output_dir,
        "cell_type": cell_type,
        "state_col": state_col,
        "raw_path": raw_path,
        "tracked_path": tracked_path,
        "state_path": state_path,
        "raw_shape": raw_shape,
        "raw_channels": raw_channels,
        "tracked_img_view": tracked_img_view,
        "state_img_view": state_img_view,
        "label_map": label_map,
        "mapping_text": mapping_text,
        "trajectory_data": trajectory_data,
        "code_colors": code_colors,
        "show_trajectories": bool(show_trajectories),
        "min_run_length": int(min_run_length),
        "tail_length": tail_length,
        "coordinate_info": coordinate_info,
    }


def launch_state_trajectory_backprojection_viewer(payload, run: bool = True):
    viewer = napari.Viewer()
    for ch, channel_img in enumerate(payload["raw_channels"]):
        viewer.add_image(channel_img, name=f"raw_ch{ch}")
    viewer.add_labels(payload["tracked_img_view"], name="TrackID", visible=False)
    viewer.add_labels(
        payload["state_img_view"],
        name="behavioral_state_class",
        visible=True,
    )

    added_dock = _add_mapping_dock_widget(
        viewer=viewer,
        mapping_text=payload["mapping_text"],
        title="State Class Mapping",
    )
    if (not added_dock) and payload.get("mapping_text"):
        print(payload["mapping_text"])

    trajectory_layers = []
    if payload.get("show_trajectories") and len(payload.get("trajectory_data", {})) > 0:
        trajectory_layers = add_track_cluster_trajectory_layers(
            viewer,
            trajectory_data=payload["trajectory_data"],
            code_colors=payload["code_colors"],
            label_map=payload["label_map"],
            output_col=payload.get("state_col", "full_behavioral_cluster"),
            tail_length=payload.get("tail_length"),
            visible=True,
        )
        for trajectory_layer in trajectory_layers:
            trajectory_layer.blending = "opaque"
        print(
            f"Added {len(trajectory_layers)} state-trajectory layer(s) "
            f"(min_run_length={payload.get('min_run_length', 1)})."
        )
    elif payload.get("show_trajectories"):
        print("No state trajectory segments found for this sample; no trajectory layers added.")

    if payload.get("sample_name"):
        print(
            "Opened state-trajectory backprojection viewer for sample "
            f"'{payload['sample_name']}' with raw='{payload['raw_path'].name}' "
            f"shape={payload['raw_shape']}, "
            f"tracked='{payload['tracked_path'].name}' "
            f"shape={tuple(int(v) for v in payload['tracked_img_view'].shape)}, "
            f"states='{payload['state_path'].name}' "
            f"shape={tuple(int(v) for v in payload['state_img_view'].shape)}."
        )

    if bool(run):
        napari.run()
    return viewer


def show_behavioral_state_backprojection_trajectories(
    sample_name: str,
    output_dir: str | Path,
    cell_type: str,
    state_col: str = "full_behavioral_cluster",
    raw_image_path: str | None = None,
    tracked_img_path: str | None = None,
    state_img_path: str | None = None,
    auto_create_if_missing: bool = True,
    refresh_if_stale: bool = True,
    show_trajectories: bool = True,
    min_run_length: int = 1,
    tail_length: int | None = None,
    verbose: bool = True,
):
    payload = build_state_trajectory_backprojection_payload(
        sample_name=sample_name,
        output_dir=output_dir,
        cell_type=cell_type,
        state_col=state_col,
        raw_image_path=raw_image_path,
        tracked_img_path=tracked_img_path,
        state_img_path=state_img_path,
        auto_create_if_missing=auto_create_if_missing,
        refresh_if_stale=refresh_if_stale,
        show_trajectories=show_trajectories,
        min_run_length=min_run_length,
        tail_length=tail_length,
        verbose=verbose,
    )
    return launch_state_trajectory_backprojection_viewer(payload, run=True)


# %%
# Interactive usage example:
# 1. Fill the INTERACTIVE_* values above.
# 2. Run this cell to resolve the BEHAV3D folder and prepare the layers + trajectory data.
# 3. Inspect `payload` in the variable explorer if you want.
#
interactive_output_dir = resolve_behav3d_output_dir(
    behav3d_folder=INTERACTIVE_BEHAV3D_FOLDER,
    output_dir=INTERACTIVE_OUTPUT_DIR,
)

payload = build_state_trajectory_backprojection_payload(
    sample_name=INTERACTIVE_SAMPLE_NAME,
    output_dir=interactive_output_dir,
    cell_type=INTERACTIVE_CELL_TYPE,
    state_col=INTERACTIVE_STATE_COL,
    raw_image_path=INTERACTIVE_RAW_IMAGE_PATH,
    tracked_img_path=INTERACTIVE_TRACKED_IMG_PATH,
    state_img_path=INTERACTIVE_STATE_IMG_PATH,
    auto_create_if_missing=INTERACTIVE_AUTO_CREATE_IF_MISSING,
    refresh_if_stale=INTERACTIVE_REFRESH_IF_STALE,
    show_trajectories=INTERACTIVE_SHOW_TRAJECTORIES,
    min_run_length=INTERACTIVE_MIN_RUN_LENGTH,
    tail_length=INTERACTIVE_TAIL_LENGTH,
    verbose=INTERACTIVE_VERBOSE,
)


# %%
# Launch napari from a prepared payload:
#
viewer = launch_state_trajectory_backprojection_viewer(payload, run=True)
