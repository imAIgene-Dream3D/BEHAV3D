import json
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path
from skimage.measure import regionprops_table
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool

import btrack
import btrack.io

from behav3d.io.images import _ensure_zarr, load_image, get_filepath_stem
from behav3d.preprocessing.tracking import convert_segments_to_tracks

# ---------------------------------------------------------------------------
# Resolve config presets to bundled JSON paths
# ---------------------------------------------------------------------------
_MODELS_DIR = Path(__file__).parent / "models"

_PRESET_MAP = {
    "cell": _MODELS_DIR / "cell_config.json",
    "particle": _MODELS_DIR / "particle_config.json",
}

_SUPPORTED_BTRACK_VISUAL_FEATURES = (
    "area",
    "major_axis_length",
    "minor_axis_length",
    "intensity_mean",
)

_BTRACK_VISUAL_FEATURES = (
    "area",
    "major_axis_length",
    "minor_axis_length",
    "intensity_mean",
)


def _run_parallel_with_fallback(fn, args_list, n_workers):
    try:
        with ProcessPoolExecutor(max_workers=n_workers) as ex:
            return list(tqdm(ex.map(fn, args_list), total=len(args_list)))
    except BrokenProcessPool:
        raise RuntimeError(
            f"Worker processes crashed while using {n_workers} workers "
            f"(likely out of memory). Please lower the number of workers "
            f"and try again."
        )


def _get_btrack_visual_feature_spec(use_visual_features=False):
    if not use_visual_features:
        return {
            "selected_features": [],
            "regionprops_properties": [],
            "feature_columns": [],
            "needs_raw_image": False,
        }

    selected_features = list(dict.fromkeys(_BTRACK_VISUAL_FEATURES))
    invalid_features = [
        feature for feature in selected_features
        if feature not in _SUPPORTED_BTRACK_VISUAL_FEATURES
    ]
    if invalid_features:
        supported = ", ".join(_SUPPORTED_BTRACK_VISUAL_FEATURES)
        invalid = ", ".join(invalid_features)
        raise ValueError(
            f"Unsupported btrack visual feature(s): {invalid}. "
            f"Supported features are: {supported}."
        )

    regionprops_properties = []
    feature_columns = []
    needs_raw_image = False
    for feature in selected_features:
        if feature == "intensity_mean":
            needs_raw_image = True
        else:
            feature_columns.append(feature)
        regionprops_properties.append(feature)

    return {
        "selected_features": selected_features,
        "regionprops_properties": regionprops_properties,
        "feature_columns": feature_columns,
        "needs_raw_image": needs_raw_image,
    }


def _resolve_btrack_visual_feature_columns(df, visual_feature_spec):
    missing_cols = [
        col for col in visual_feature_spec["feature_columns"]
        if col not in df.columns
    ]
    if missing_cols:
        available = ", ".join(df.columns)
        missing = ", ".join(missing_cols)
        raise ValueError(
            f"Requested visual feature columns are missing: {missing}. "
            f"Available columns are: {available}."
        )

    intensity_cols = []
    if "intensity_mean" in visual_feature_spec["selected_features"]:
        intensity_cols = sorted(
            [col for col in df.columns if col.startswith("intensity_mean-")],
            key=lambda col: int(col.split("-")[-1]),
        )
        if not intensity_cols:
            available = ", ".join(df.columns)
            raise ValueError(
                "Requested visual feature 'intensity_mean' but no "
                f"'intensity_mean-*' columns were extracted. Available columns are: {available}."
            )

    feature_columns = []
    for feature in visual_feature_spec["selected_features"]:
        if feature == "intensity_mean":
            feature_columns.extend(intensity_cols)
        else:
            feature_columns.append(feature)
    return feature_columns


def _sanitize_feature_columns(df, feature_cols):
    if not feature_cols:
        return df
    df = df.copy()
    for col in feature_cols:
        if col in df.columns:
            # Pandas can hand back a read-only NumPy view here, so force a writable copy
            # before replacing non-finite values.
            values = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float, copy=True)
            values[~np.isfinite(values)] = 0.0
            df[col] = values
    return df


def _extract_btrack_props_single_timepoint(
    t,
    t_seg,
    t_raw=None,
    visual_feature_spec=None,
):
    t_seg = np.asarray(t_seg)
    if t_seg.max() == 0:
        return None

    if visual_feature_spec is None:
        visual_feature_spec = _get_btrack_visual_feature_spec(False)

    props_kwargs = {
        "label_image": t_seg,
        "properties": ["label", "centroid"],
    }

    if visual_feature_spec["regionprops_properties"]:
        props_kwargs["properties"] += visual_feature_spec["regionprops_properties"]

    if visual_feature_spec["needs_raw_image"]:
        if t_raw is None:
            raise ValueError(
                "Visual-feature btrack requires raw_image_path so mean intensities "
                "can be extracted from the raw image."
            )
        t_raw = np.asarray(t_raw)
        if t_raw.ndim != 4:
            raise ValueError(
                f"Raw image timepoint should have shape (C, Z, Y, X), but "
                f"got {t_raw.shape} at t={t}."
            )
        if tuple(t_raw.shape[1:]) != tuple(t_seg.shape):
            raise ValueError(
                "Raw image and segmentation spatial shapes must match for "
                "motion+visual btrack. "
                f"Got raw {t_raw.shape[1:]} and segmentation {t_seg.shape} "
                f"at t={t}."
            )
        props_kwargs["intensity_image"] = np.moveaxis(t_raw, 0, -1)

    props = pd.DataFrame(regionprops_table(**props_kwargs))
    props["position_t"] = t
    return props


def _extract_btrack_props_single_timepoint_from_paths(args):
    t, segments_path, raw_image_path, visual_feature_spec = args
    segments = load_image(segments_path)
    t_seg = np.asarray(segments[t])
    t_raw = None
    if visual_feature_spec["needs_raw_image"]:
        raw_image = load_image(raw_image_path)
        t_raw = np.asarray(raw_image[t])
    return _extract_btrack_props_single_timepoint(
        t=t,
        t_seg=t_seg,
        t_raw=t_raw,
        visual_feature_spec=visual_feature_spec,
    )


def _extract_btrack_objects_dataframe(
    segments=None,
    segments_path=None,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    raw_image=None,
    raw_image_path=None,
    use_visual_features=False,
    n_workers=1,
):
    n_workers = max(1, int(n_workers or 1))
    df_objects = []
    visual_feature_spec = _get_btrack_visual_feature_spec(use_visual_features)

    if visual_feature_spec["needs_raw_image"] and raw_image is None and raw_image_path is None:
        raise ValueError(
            "Visual-feature btrack requires raw_image_path so mean intensities "
            "can be extracted from the raw image."
        )

    can_parallel = (
        n_workers > 1
        and segments_path is not None
        and (not visual_feature_spec["needs_raw_image"] or raw_image_path is not None)
    )

    if can_parallel:
        segments_ref = load_image(segments_path)
        timepoints = int(segments_ref.shape[0])
        args_list = [
            (t, segments_path, raw_image_path, visual_feature_spec)
            for t in range(timepoints)
        ]
        results = _run_parallel_with_fallback(
            _extract_btrack_props_single_timepoint_from_paths,
            args_list,
            n_workers,
        )
        df_objects = [props for props in results if props is not None and not props.empty]
    else:
        if segments is None:
            segments = load_image(segments_path)
        if visual_feature_spec["needs_raw_image"] and raw_image is None:
            raw_image = load_image(raw_image_path)

        if visual_feature_spec["needs_raw_image"]:
            if raw_image.ndim != 5:
                raise ValueError(
                    f"Raw image should have 5 dimensions (T, C, Z, Y, X), but "
                    f"{raw_image.ndim} found."
                )
            if len(raw_image) != len(segments):
                raise ValueError(
                    "Raw image and segmentation must have the same number of timepoints "
                    f"for motion+visual btrack: got {len(raw_image)} and {len(segments)}."
                )

        for t, t_seg in enumerate(tqdm(segments, desc="Extracting objects")):
            t_raw = None
            if visual_feature_spec["needs_raw_image"]:
                t_raw = np.asarray(raw_image[t])
            props = _extract_btrack_props_single_timepoint(
                t=t,
                t_seg=t_seg,
                t_raw=t_raw,
                visual_feature_spec=visual_feature_spec,
            )
            if props is not None and not props.empty:
                df_objects.append(props)

    if not df_objects:
        return None, []

    df_objects = pd.concat(df_objects, ignore_index=True)
    df_objects["position_z"] = df_objects["centroid-0"] * element_size_z
    df_objects["position_y"] = df_objects["centroid-1"] * element_size_y
    df_objects["position_x"] = df_objects["centroid-2"] * element_size_x
    df_objects.rename(columns={
        "centroid-0": "pixel_position_z",
        "centroid-1": "pixel_position_y",
        "centroid-2": "pixel_position_x",
    }, inplace=True)
    visual_feature_cols = _resolve_btrack_visual_feature_columns(
        df_objects,
        visual_feature_spec,
    )
    df_objects = _sanitize_feature_columns(df_objects, visual_feature_cols)

    return df_objects, visual_feature_cols


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
    raw_image=None,
    raw_image_path=None,
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
    n_workers=1,
    use_visual_features=None,
    return_trackimg=True,
):
    """Run btrack (Bayesian tracking) on a single segmentation array.

    Parameters
    ----------
    segments : ndarray, optional
        Segmentation label array (T, Z, Y, X) or (T, Y, X).
    segments_path : str or Path, optional
        Path to zarr / tiff segmentation (used if *segments* is None).
    raw_image : ndarray, optional
        Raw image array in BEHAV3D order (T, C, Z, Y, X).
    raw_image_path : str or Path, optional
        Path to raw image used for visual features when enabled.
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
    n_workers : int, optional
        Number of parallel workers for regionprops extraction. Path-backed
        inputs use multiprocessing when ``n_workers > 1``.
    use_visual_features : bool, optional
        Whether to use raw-image-derived visual features during linking.
        ``None`` is treated the same as ``False``.
    return_trackimg : bool
        Whether to produce a tracked zarr image.
    """
    assert segments is not None or segments_path is not None, \
        "Either segments or segments_path must be provided"

    if segments_path is not None:
        segments_path = Path(segments_path)
    if raw_image_path is not None:
        raw_image_path = Path(raw_image_path)
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
    # Resolve config and whether this run should use visual updates
    # ------------------------------------------------------------------
    config_path = _resolve_config(config_preset)
    use_visual_updates = bool(use_visual_features)
    visual_feature_spec = _get_btrack_visual_feature_spec(use_visual_updates)
    if visual_feature_spec["needs_raw_image"] and raw_image is None:
        if raw_image_path is None:
            raise ValueError(
                "Visual-feature btrack requires raw_image_path in metadata."
            )
        raw_image = load_image(raw_image_path)

    # ------------------------------------------------------------------
    # 1. Extract detections and optional visual features
    # ------------------------------------------------------------------
    df_centroids, visual_feature_cols = _extract_btrack_objects_dataframe(
        segments=segments,
        segments_path=segments_path,
        element_size_x=element_size_x,
        element_size_y=element_size_y,
        element_size_z=element_size_z,
        raw_image=raw_image,
        raw_image_path=raw_image_path,
        use_visual_features=use_visual_updates,
        n_workers=n_workers,
    )

    if df_centroids is None:
        print("WARNING: No objects found in segmentation — skipping btrack.")
        return

    # ------------------------------------------------------------------
    # 2. Build btrack objects from detections
    # ------------------------------------------------------------------
    obj_cols = ["position_t", "position_x", "position_y", "position_z", "label"]
    objects_arr = df_centroids[obj_cols].to_numpy()
    objects = btrack.io.objects_from_array(
        objects_arr,
        default_keys=["t", "x", "y", "z", "label"],
    )
    if visual_feature_cols:
        feature_records = df_centroids[visual_feature_cols].to_dict(orient="records")
        for obj, props in zip(objects, feature_records):
            obj.properties = props

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
            if visual_feature_cols:
                tracker.features = visual_feature_cols

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
            tracking_updates = ["motion", "visual"] if visual_feature_cols else ["motion"]
            tracker.track(
                step_size=step_size,
                tracking_updates=tracking_updates,
            )

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
            n_workers=n_workers,
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
    n_workers=1,
    use_visual_features=None,
    return_trackimg=True,
    overwrite=False,
    log_callback=None,
    progress_cb=None,
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
    n_workers : int
        Number of parallel workers for regionprops extraction.
    use_visual_features : bool, optional
        Whether to use raw-image-derived visual features during linking.
        ``None`` is treated the same as ``False``.
    return_trackimg : bool
        Whether to save tracked zarr image.
    overwrite : bool
        Whether to overwrite existing results.
    log_callback : callable, optional
        Function to log messages (for GUI integration).
    progress_cb : callable, optional
        Called as ``progress_cb(current, total, label)`` from inside the
        per-sample loop so a GUI can drive a progress bar.  ``None`` is a
        no-op (default).  Notebook callers omit this kwarg.
    """
    _log = log_callback or print

    total_samples = len(metadata)
    for i, (idx, sample) in enumerate(metadata.iterrows()):
        sample_name = sample["sample_name"]
        if progress_cb is not None:
            try:
                progress_cb(i, total_samples, f"{cell_type} / {sample_name}")
            except Exception:
                pass
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

        segments_path = Path(sample[segments_col])
        _ensure_zarr(segments_path, label=f"Segments for '{sample_name}'")

        raw_image_path = sample.get("raw_image_path", None)
        if raw_image_path:
            _ensure_zarr(Path(raw_image_path), label=f"Raw image for '{sample_name}'")
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
                raw_image_path=raw_image_path,
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
                n_workers=n_workers,
                use_visual_features=use_visual_features,
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

    if progress_cb is not None:
        try:
            progress_cb(total_samples, total_samples, f"{cell_type} done")
        except Exception:
            pass
    return metadata
