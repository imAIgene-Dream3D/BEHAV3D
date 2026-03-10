import json
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path
from skimage.measure import regionprops_table
from tqdm import tqdm

import btrack
import btrack.io

from behav3d.io.images import load_image, get_filepath_stem
from behav3d.preprocessing.tracking import convert_segments_to_tracks

# ---------------------------------------------------------------------------
# Resolve config presets to bundled JSON paths
# ---------------------------------------------------------------------------
_MODELS_DIR = Path(__file__).parent / "models"

_PRESET_MAP = {
    "cell": _MODELS_DIR / "cell_config.json",
    "particle": _MODELS_DIR / "particle_config.json",
}


def _resolve_config(config_preset):
    """Return an absolute Path to a btrack JSON config file.

    Parameters
    ----------
    config_preset : str or Path
        Either a preset name ("cell", "particle") or a path to a custom
        JSON configuration file.
    """
    if isinstance(config_preset, Path):
        return config_preset
    key = str(config_preset).strip().lower()
    if key in _PRESET_MAP:
        return _PRESET_MAP[key]
    # Treat as a file path
    p = Path(config_preset)
    if p.exists():
        return p
    raise FileNotFoundError(
        f"btrack config not found: '{config_preset}'. "
        f"Use a preset name ({list(_PRESET_MAP.keys())}) or an existing JSON path."
    )


def _override_hypothesis_model(config_path, hypotheses=None, dist_thresh=None,
                                time_thresh=None):
    """Load a btrack legacy JSON config, optionally override hypothesis fields,
    and write a patched temp file so ``tracker.configure()`` can read it.

    The wrapper ``{"TrackerConfig": {...}}`` is preserved so btrack uses
    its legacy loader (PascalCase keys).

    Returns the path to the temporary file (caller is responsible for cleanup).
    """
    with open(config_path, "r") as f:
        cfg = json.load(f)

    # Work on the inner dict while keeping the wrapper intact.
    inner = cfg.get("TrackerConfig", cfg)
    hmodel = inner.get("HypothesisModel", {})
    if hypotheses is not None:
        hmodel["hypotheses"] = list(hypotheses)
    if dist_thresh is not None:
        hmodel["dist_thresh"] = int(dist_thresh)
    if time_thresh is not None:
        hmodel["time_thresh"] = int(time_thresh)
    inner["HypothesisModel"] = hmodel

    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False, encoding="utf-8"
    )
    json.dump(cfg, tmp)
    tmp.close()
    return Path(tmp.name)


# ---------------------------------------------------------------------------
# Single-image tracking
# ---------------------------------------------------------------------------
def btrack_image(
    segments=None,
    segments_path=None,
    tracked_img_outpath=None,
    tracked_csv_outpath=None,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    config_preset="cell",
    max_search_radius=100,
    update_method="EXACT",
    step_size=100,
    volume_padding=10,
    use_optimize=False,
    hypotheses=None,
    dist_thresh=None,
    time_thresh=None,
    return_trackimg=True,
):
    """Run btrack (Bayesian tracking) on a single segmentation array.

    Parameters
    ----------
    segments : ndarray, optional
        Segmentation label array (T, Z, Y, X) or (T, Y, X).
    segments_path : str or Path, optional
        Path to zarr / tiff segmentation (used if *segments* is None).
    tracked_img_outpath : Path, optional
        Where to save the tracked zarr image.
    tracked_csv_outpath : Path, optional
        Where to save the track CSV.
    element_size_x/y/z : float
        Physical pixel sizes for coordinate scaling.
    config_preset : str or Path
        ``"cell"``, ``"particle"``, or a path to a custom JSON.
    max_search_radius : int
        Maximum isotropic search radius (pixels).
    update_method : str
        ``"EXACT"`` or ``"APPROXIMATE"``.
    step_size : int
        Frames per tracking iteration batch.
    volume_padding : int
        Pixels added around data bounding-box for volume estimation.
    use_optimize : bool
        Whether to run the global hypothesis optimizer (Step 2).
    hypotheses : list of str, optional
        Override the hypothesis list in the config (e.g.
        ``["P_FP", "P_init", "P_term", "P_link"]``). ``None`` keeps
        the config-file default.
    dist_thresh : int, optional
        Override distance threshold for hypothesis generation.
    time_thresh : int, optional
        Override time threshold for hypothesis generation.
    return_trackimg : bool
        Whether to produce a tracked zarr image.
    """
    assert segments is not None or segments_path is not None, \
        "Either segments or segments_path must be provided"

    if segments_path is not None:
        segments_path = Path(segments_path)
    if segments is None:
        segments = load_image(segments_path)

    basename = get_filepath_stem(segments_path)
    if tracked_img_outpath is None:
        tracked_img_outpath = Path(segments_path.parent,
                                   f"{basename}_tracked.zarr")
    if tracked_csv_outpath is None:
        tracked_csv_outpath = Path(segments_path.parent,
                                   f"{basename}_tracks.csv")

    # ------------------------------------------------------------------
    # 1. Extract centroids per timepoint (same pattern as laptracking.py)
    # ------------------------------------------------------------------
    df_centroids = []
    for t, t_seg in enumerate(tqdm(segments, desc="Extracting centroids")):
        t_seg = np.asarray(t_seg)
        if t_seg.max() == 0:
            continue
        props = pd.DataFrame(
            regionprops_table(label_image=t_seg,
                              properties=["label", "centroid", "area"])
        )
        props["position_t"] = t
        df_centroids.append(props)

    if not df_centroids:
        print("WARNING: No objects found in segmentation — skipping btrack.")
        return

    df_centroids = pd.concat(df_centroids, ignore_index=True)
    df_centroids["position_z"] = df_centroids["centroid-0"] * element_size_z
    df_centroids["position_y"] = df_centroids["centroid-1"] * element_size_y
    df_centroids["position_x"] = df_centroids["centroid-2"] * element_size_x

    # Keep pixel positions for output
    df_centroids.rename(columns={
        "centroid-0": "pixel_position_z",
        "centroid-1": "pixel_position_y",
        "centroid-2": "pixel_position_x",
    }, inplace=True)

    # ------------------------------------------------------------------
    # 2. Build btrack objects from centroids
    # ------------------------------------------------------------------
    obj_cols = ["position_t", "position_x", "position_y", "position_z", "label"]
    objects_arr = df_centroids[obj_cols].to_numpy()
    objects = btrack.io.objects_from_array(
        objects_arr,
        default_keys=["t", "x", "y", "z", "label"],
    )

    # ------------------------------------------------------------------
    # 3. Configure & run tracker
    # ------------------------------------------------------------------
    config_path = _resolve_config(config_preset)

    # Load config, optionally overriding hypothesis params
    # Resolve the config file path, writing a patched temp file when needed.
    need_override = use_optimize and (
        hypotheses is not None or dist_thresh is not None or time_thresh is not None
    )
    if need_override:
        active_config_path = _override_hypothesis_model(
            config_path,
            hypotheses=hypotheses,
            dist_thresh=dist_thresh,
            time_thresh=time_thresh,
        )
        _delete_active_config = True
    else:
        active_config_path = config_path
        _delete_active_config = False

    try:
        with btrack.BayesianTracker() as tracker:
            tracker.configure(str(active_config_path))

            tracker.append(objects)

            # Auto-compute volume from data bounds + padding.
            # Do NOT clamp to 0 — negative bounds are valid and ensure all
            # objects have padding on every side, which the hypothesis model needs.
            x_coords = [o.x for o in objects]
            y_coords = [o.y for o in objects]
            z_coords = [o.z for o in objects]
            tracker.volume = (
                (min(x_coords) - volume_padding, max(x_coords) + volume_padding),
                (min(y_coords) - volume_padding, max(y_coords) + volume_padding),
                (min(z_coords) - volume_padding, max(z_coords) + volume_padding),
            )

            # Set update method
            if update_method.upper() == "APPROXIMATE":
                from btrack.constants import BayesianUpdates
                tracker.update_method = BayesianUpdates.APPROXIMATE
                tracker.max_search_radius = max_search_radius
            else:
                tracker.max_search_radius = max_search_radius

            # Step 1 — Kalman filter linking
            tracker.track(step_size=step_size)

            # Step 2 — Global hypothesis optimizer (optional)
            if use_optimize:
                tracker.optimize()

            tracks = tracker.tracks
    finally:
        if _delete_active_config:
            active_config_path.unlink(missing_ok=True)

    print(f"btrack produced {len(tracks)} tracks")

    # ------------------------------------------------------------------
    # 4. Convert tracks → BEHAV3D standard CSV
    # ------------------------------------------------------------------
    # Use track.refs (original object indices into objects_arr) for lookup.
    # This is reliable because btrack.refs always points back to the exact
    # input object, regardless of Kalman smoothing on the returned positions.
    # Dummy observations inserted by btrack to fill gaps have refs == -1
    # and must be skipped — they have no corresponding segmented object and
    # their position is (0, 0, 0), which caused the "radiating from corner"
    # artefact seen in the track visualisation.
    idx_to_label = df_centroids["label"].to_numpy()          # shape (N,)
    idx_to_px_x  = df_centroids["pixel_position_x"].to_numpy()
    idx_to_px_y  = df_centroids["pixel_position_y"].to_numpy()
    idx_to_px_z  = df_centroids["pixel_position_z"].to_numpy()
    idx_to_pos_x = df_centroids["position_x"].to_numpy()
    idx_to_pos_y = df_centroids["position_y"].to_numpy()
    idx_to_pos_z = df_centroids["position_z"].to_numpy()
    idx_to_pos_t = df_centroids["position_t"].to_numpy()

    rows = []
    for track in tracks:
        tid = int(track.ID) + 1  # 1-based; 0 reserved for background
        refs = np.asarray(track.refs, dtype=int)
        t_arr = np.asarray(track.t, dtype=int)

        for i, ref in enumerate(refs):
            if ref < 0:
                # Dummy object inserted by btrack to bridge a gap — skip.
                continue
            rows.append({
                "TrackID": tid,
                "SegmentID": int(idx_to_label[ref]),
                "position_t": int(idx_to_pos_t[ref]),
                "position_x": float(idx_to_pos_x[ref]),
                "position_y": float(idx_to_pos_y[ref]),
                "position_z": float(idx_to_pos_z[ref]),
                "pixel_position_x": float(idx_to_px_x[ref]),
                "pixel_position_y": float(idx_to_px_y[ref]),
                "pixel_position_z": float(idx_to_px_z[ref]),
            })

    df_tracks = pd.DataFrame(rows)
    df_tracks = df_tracks[[
        "TrackID", "SegmentID", "position_t",
        "position_x", "position_y", "position_z",
        "pixel_position_x", "pixel_position_y", "pixel_position_z",
    ]]
    df_tracks.to_csv(tracked_csv_outpath, sep=",", index=False)

    if return_trackimg:
        print("Overwriting the original segment IDs with the tracked IDs")
        convert_segments_to_tracks(
            df_tracks=df_tracks,
            segments=segments,
            outpath=tracked_img_outpath,
        )


# ---------------------------------------------------------------------------
# Batch wrapper (mirrors run_laptracking)
# ---------------------------------------------------------------------------
def run_btracking(
    metadata,
    output_dir,
    cell_type,
    config_preset="cell",
    max_search_radius=100,
    update_method="EXACT",
    step_size=100,
    volume_padding=10,
    use_optimize=False,
    hypotheses=None,
    dist_thresh=None,
    time_thresh=None,
    return_trackimg=True,
    overwrite=False,
    log_callback=None,
    **kwargs,
):
    """Run btrack (Bayesian tracking) on all samples for a given cell type.

    Parameters
    ----------
    metadata : pd.DataFrame
        DataFrame containing sample information.
    output_dir : str or Path
        Root output directory.
    cell_type : str
        Name of the cell type to track.
    config_preset : str or Path
        ``"cell"``, ``"particle"``, or path to custom JSON.
    max_search_radius : int
        Max isotropic search radius (pixels).
    update_method : str
        ``"EXACT"`` or ``"APPROXIMATE"``.
    step_size : int
        Frames per tracking batch.
    volume_padding : int
        Pixels padding around data bounding-box.
    use_optimize : bool
        Run global hypothesis optimizer (Step 2).
    hypotheses : list of str, optional
        Override hypothesis list (None keeps config default).
    dist_thresh : int, optional
        Override distance threshold.
    time_thresh : int, optional
        Override time threshold.
    return_trackimg : bool
        Whether to save tracked zarr image.
    overwrite : bool
        Whether to overwrite existing results.
    log_callback : callable, optional
        Function to log messages (for GUI integration).
    """
    _log = log_callback or print

    for idx, sample in metadata.iterrows():
        sample_name = sample["sample_name"]
        _log(f"Tracking sample: {sample_name}")

        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name,
                                  cell_type)

        # Find the correct prefixed column (or_, im_, ot_)
        segments_col = None
        for prefix in ["or", "im", "ot"]:
            col_name = f"{prefix}_{cell_type}_segments_image_path"
            if col_name in sample.index and pd.notna(sample[col_name]):
                segments_col = col_name
                break

        if segments_col is None:
            segments_col = f"{cell_type}_segments_image_path"

        segments_path = sample[segments_col]
        tracked_img_outpath = Path(
            tracked_img_outdir,
            f"{sample_name}_{cell_type}_tracked.zarr",
        )
        tracked_csv_outpath = Path(
            tracked_csv_outdir,
            f"{sample_name}_{cell_type}_tracks.csv",
        )

        tracked_img_outdir.mkdir(parents=True, exist_ok=True)
        tracked_csv_outdir.mkdir(parents=True, exist_ok=True)

        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]

        if (not tracked_csv_outpath.exists()
                or not tracked_img_outpath.exists()
                or overwrite):
            btrack_image(
                segments_path=segments_path,
                tracked_img_outpath=tracked_img_outpath,
                tracked_csv_outpath=tracked_csv_outpath,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                config_preset=config_preset,
                max_search_radius=max_search_radius,
                update_method=update_method,
                step_size=step_size,
                volume_padding=volume_padding,
                use_optimize=use_optimize,
                hypotheses=hypotheses,
                dist_thresh=dist_thresh,
                time_thresh=time_thresh,
                return_trackimg=return_trackimg,
            )
        else:
            _log("Tracking already exists... Provide overwrite=True to "
                 "overwrite... Loading existing tracking data")

        # Update metadata with prefixed column names
        if (segments_col is not None
                and segments_col.startswith(("or_", "im_", "ot_"))):
            prefix = segments_col.split("_")[0]
            img_col = f"{prefix}_{cell_type}_tracks_image_path"
            csv_col = f"{prefix}_{cell_type}_tracks_csv_path"
        else:
            img_col = f"{cell_type}_tracks_image_path"
            csv_col = f"{cell_type}_tracks_csv_path"

        for col in [img_col, csv_col]:
            if col not in metadata.columns or metadata[col].dtype != object:
                metadata[col] = metadata.get(
                    col, pd.Series(dtype=object)
                ).astype(object)

        metadata.at[idx, img_col] = str(tracked_img_outpath)
        metadata.at[idx, csv_col] = str(tracked_csv_outpath)

    return metadata
