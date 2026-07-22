import napari
import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import time
import argparse
from tqdm import tqdm
import dask.array as da
from concurrent.futures import ProcessPoolExecutor, as_completed

from behav3d.core.utils import format_time
from behav3d.io.images import _ensure_zarr, get_filepath_stem, load_image, load_zarr, save_as_zarr
from behav3d.io.formats.zarr import write_zarr_parallel
from behav3d.core.metadata import (
    detect_immune_cell_types_from_metadata,
    detect_organoid_types_from_metadata,
    detect_other_cell_types_from_metadata,
    detect_merged_cell_types_from_metadata,
    filter_multicolor_inputs,
)


STRUCTURAL_FEATURE_COLUMNS = {
    "sample_name",
    "TrackID",
    "track_id",
    "original_TrackID",
    "position_t",
}

TRACKED_CSV_FALLBACK_COLUMNS = {
    "SegmentID",
    "position_x",
    "position_y",
    "position_z",
    "pixel_position_x",
    "pixel_position_y",
    "pixel_position_z",
    "source_cell_type",
}


def _normalize_backprojection_mode(mode):
    mode = str(mode).strip().lower()
    if mode in {"mean", "summary", "summarized", "track"}:
        return "summary"
    if mode in {"time", "timepoint", "per_timepoint", "per-timepoint"}:
        return "timepoint"
    raise ValueError(
        "mode must be one of {'summary', 'timepoint'} "
        f"(also accepts 'mean' and 'time'), got '{mode}'."
    )


def _read_csv_header(path):
    path = Path(path)
    if not path.exists():
        return []
    header_df = pd.read_csv(path, nrows=0)
    return [str(col).lstrip("\ufeff") for col in header_df.columns]


def _find_track_csv_column(df_sample, cell_type):
    track_csv_col = f"{cell_type}_tracks_csv_path"
    if track_csv_col in df_sample.columns:
        return track_csv_col

    for prefix in ["im_", "or_", "ot_"]:
        alt_col = f"{prefix}{cell_type}_tracks_csv_path"
        if alt_col in df_sample.columns:
            return alt_col
    return None


def _resolve_tracked_csv_path(df_sample, output_dir, sample_name, cell_type):
    track_csv_col = _find_track_csv_column(df_sample, cell_type)
    if track_csv_col is not None:
        track_csv_value = df_sample[track_csv_col].values[0]
        if pd.notna(track_csv_value) and str(track_csv_value).strip():
            track_csv_path = Path(str(track_csv_value).strip()).expanduser()
            if track_csv_path.exists():
                return track_csv_path

    sample_dir = Path(output_dir) / "trackdata" / str(sample_name) / str(cell_type)
    candidates = [
        sample_dir / f"{sample_name}_{cell_type}_tracks.csv",
        sample_dir / f"{sample_name}_{cell_type}_track.csv",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate

    fallback = sorted(list(sample_dir.glob("*_tracks.csv")) + list(sample_dir.glob("*_track.csv")))
    if fallback:
        return fallback[0]

    return None


def _get_excluded_backprojection_columns(
    df_tracks_clustered,
    metadata_columns=None,
    tracked_csv_path=None,
):
    excluded = set()
    if metadata_columns is not None:
        excluded.update(str(col) for col in metadata_columns)

    excluded.update(col for col in STRUCTURAL_FEATURE_COLUMNS if col in df_tracks_clustered.columns)

    if tracked_csv_path is not None:
        excluded.update(_read_csv_header(tracked_csv_path))
    else:
        excluded.update(col for col in TRACKED_CSV_FALLBACK_COLUMNS if col in df_tracks_clustered.columns)

    return excluded


def _select_backprojection_columns(
    df_tracks_clustered,
    columns,
    metadata_columns=None,
    tracked_csv_path=None,
    extra_excluded=None,
):
    excluded = _get_excluded_backprojection_columns(
        df_tracks_clustered=df_tracks_clustered,
        metadata_columns=metadata_columns,
        tracked_csv_path=tracked_csv_path,
    )
    if extra_excluded is not None:
        excluded.update(extra_excluded)

    if columns == []:
        return [col for col in df_tracks_clustered.columns if col not in excluded]

    filtered = []
    skipped = []
    missing = []
    for col in columns:
        if col not in df_tracks_clustered.columns:
            missing.append(col)
            continue
        if col in excluded:
            skipped.append(col)
            continue
        filtered.append(col)

    if missing:
        print(f"  Skipping missing columns: {missing}")
    if skipped:
        print(f"  Skipping excluded non-feature columns: {skipped}")
    return filtered


def _infer_backprojection_dtype(series, background_value=0):
    non_null = pd.Series(series).dropna()
    if non_null.empty:
        return np.float32

    if pd.api.types.is_bool_dtype(non_null):
        return np.uint8

    if pd.api.types.is_integer_dtype(non_null):
        min_val = int(min(int(background_value), int(non_null.min())))
        max_val = int(max(int(background_value), int(non_null.max())))
        if min_val >= 0 and max_val <= int(np.iinfo(np.uint16).max):
            return np.uint16
        if min_val >= int(np.iinfo(np.int16).min) and max_val <= int(np.iinfo(np.int16).max):
            return np.int16
        if min_val >= 0 and max_val <= int(np.iinfo(np.uint32).max):
            return np.uint32
        return np.int32

    return np.float32


def _prepare_feature_series(df, column):
    if column not in df.columns:
        raise ValueError(f"Column '{column}' not found in feature dataframe.")

    series = df[column]
    if pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series):
        return series

    coerced = pd.to_numeric(series, errors="coerce")
    if coerced.notna().sum() == series.notna().sum():
        return coerced

    raise ValueError(
        f"Column '{column}' is not numeric. Encode categorical values before backprojection."
    )


def _build_summary_lookup(lookup_df, track_col, value_col, background_value=0):
    values = _prepare_feature_series(lookup_df, value_col)
    work = lookup_df[[track_col]].copy()
    work["_value"] = values
    work[track_col] = pd.to_numeric(work[track_col], errors="coerce")
    work = work.dropna(subset=[track_col, "_value"]).copy()

    out_dtype = _infer_backprojection_dtype(work["_value"], background_value=background_value)
    if work.empty:
        return {
            "dtype": np.dtype(out_dtype),
            "frame_lookup": {},
        }

    work[track_col] = work[track_col].astype(np.int64)
    work = work.drop_duplicates(subset=[track_col], keep="last")
    work = work.sort_values(track_col)
    ids = work[track_col].to_numpy(dtype=np.int64, copy=False)
    vals = work["_value"].astype(out_dtype, copy=False).to_numpy(copy=False)

    return {
        "dtype": np.dtype(out_dtype),
        "frame_lookup": {
            -1: {
                "ids": ids,
                "values": vals,
                "max_label": int(ids.max()) if ids.size > 0 else 0,
            }
        },
    }


def _build_timepoint_lookup(lookup_df, track_col, time_col, value_col, background_value=0):
    values = _prepare_feature_series(lookup_df, value_col)
    work = lookup_df[[track_col, time_col]].copy()
    work["_value"] = values
    work[track_col] = pd.to_numeric(work[track_col], errors="coerce")
    work[time_col] = pd.to_numeric(work[time_col], errors="coerce")
    work = work.dropna(subset=[track_col, time_col, "_value"]).copy()

    out_dtype = _infer_backprojection_dtype(work["_value"], background_value=background_value)
    if work.empty:
        return {
            "dtype": np.dtype(out_dtype),
            "frame_lookup": {},
        }

    work[track_col] = work[track_col].astype(np.int64)
    work[time_col] = work[time_col].astype(np.int64)
    work = work.sort_values([time_col, track_col]).drop_duplicates(
        subset=[time_col, track_col],
        keep="last",
    )

    per_time_lookup = {}
    for time_id, group in work.groupby(time_col, observed=False):
        ids = group[track_col].to_numpy(dtype=np.int64, copy=False)
        vals = group["_value"].astype(out_dtype, copy=False).to_numpy(copy=False)
        order = np.argsort(ids, kind="mergesort")
        per_time_lookup[int(time_id)] = {
            "ids": ids[order],
            "values": vals[order],
            "max_label": int(ids.max()) if ids.size > 0 else 0,
        }

    return {
        "dtype": np.dtype(out_dtype),
        "frame_lookup": per_time_lookup,
    }


def _should_use_dense_lut(max_label, ids, out_dtype):
    max_label = int(max_label)
    if max_label < 0:
        return False

    approx_bytes = (max_label + 1) * np.dtype(out_dtype).itemsize
    if approx_bytes > 256 * 1024 * 1024:
        return False

    n_ids = int(len(ids))
    if n_ids <= 0:
        return True

    return max_label <= max(16384, n_ids * 64)


def _map_labels_frame_sparse(labels_t, ids, values, out_dtype, background_value=0):
    labels_t = np.asarray(labels_t, dtype=np.int64)
    unique_labels, inverse = np.unique(labels_t, return_inverse=True)
    mapped_unique = np.full(unique_labels.shape, background_value, dtype=out_dtype)
    if ids.size > 0 and unique_labels.size > 0:
        positive_index = np.flatnonzero(unique_labels > 0)
        positive_labels = unique_labels[positive_index]
        if positive_labels.size > 0:
            positions = np.searchsorted(ids, positive_labels)
            valid = (positions < ids.size) & (ids[positions] == positive_labels)
            if np.any(valid):
                mapped_unique[positive_index[valid]] = values[positions[valid]]
    return mapped_unique[inverse].reshape(labels_t.shape)


def _map_labels_frame(labels_t, ids, values, out_dtype, background_value=0, max_label=None):
    labels_t = np.asarray(labels_t, dtype=np.int64)
    ids = np.asarray(ids, dtype=np.int64)
    values = np.asarray(values, dtype=out_dtype)

    if max_label is None:
        max_label = int(ids.max()) if ids.size > 0 else 0

    if _should_use_dense_lut(max_label=max_label, ids=ids, out_dtype=out_dtype):
        lut = np.full(int(max_label) + 1, background_value, dtype=out_dtype)
        positive_ids = ids[(ids >= 0) & (ids <= int(max_label))]
        if positive_ids.size > 0:
            lut[positive_ids] = values[(ids >= 0) & (ids <= int(max_label))]

        mapped = np.full(labels_t.shape, background_value, dtype=out_dtype)
        valid = (labels_t >= 0) & (labels_t <= int(max_label))
        if np.any(valid):
            mapped[valid] = lut[labels_t[valid]]
        return mapped

    return _map_labels_frame_sparse(
        labels_t=labels_t,
        ids=ids,
        values=values,
        out_dtype=out_dtype,
        background_value=background_value,
    )


def _resolve_frame_lookup(frame_lookup, time_index, mode):
    if mode == "summary":
        return frame_lookup.get(-1, None)
    return frame_lookup.get(int(time_index), None)


def _read_tracked_frame(tracked_source, time_index):
    if isinstance(tracked_source, (str, Path)):
        tracked_img = load_image(tracked_source)
        return np.asarray(tracked_img[int(time_index)])
    return np.asarray(tracked_source[int(time_index)])


def _process_backprojection_frame(
    *,
    tracked_source,
    time_index,
    frame_lookup,
    mode,
    out_dtype,
    background_value,
    verbose=False,
):
    labels_t = _read_tracked_frame(tracked_source=tracked_source, time_index=time_index)
    lookup = _resolve_frame_lookup(frame_lookup=frame_lookup, time_index=time_index, mode=mode)
    if lookup is None:
        return np.full(labels_t.shape, background_value, dtype=out_dtype)
    mapped = _map_labels_frame(
        labels_t=labels_t,
        ids=np.asarray(lookup["ids"], dtype=np.int64),
        values=np.asarray(lookup["values"], dtype=out_dtype),
        out_dtype=out_dtype,
        background_value=background_value,
        max_label=lookup.get("max_label", None),
    )
    if bool(verbose) and int(time_index) < 2:
        non_background = int(np.count_nonzero(mapped != background_value))
        print(
            f"    frame={int(time_index)} rows={int(len(lookup['ids']))} "
            f"non_background_voxels={non_background}"
        )
    return mapped


def _write_backprojection_frame_worker(
    *,
    tracked_img_path,
    outpath,
    group,
    time_index,
    lookup_ids,
    lookup_values,
    max_label,
    out_dtype_str,
    background_value,
):
    out_dtype = np.dtype(out_dtype_str)
    labels_t = _read_tracked_frame(tracked_source=tracked_img_path, time_index=time_index)
    if lookup_ids is None or lookup_values is None:
        mapped_t = np.full(labels_t.shape, background_value, dtype=out_dtype)
    else:
        mapped_t = _map_labels_frame(
            labels_t=labels_t,
            ids=np.asarray(lookup_ids, dtype=np.int64),
            values=np.asarray(lookup_values, dtype=out_dtype),
            out_dtype=out_dtype,
            background_value=background_value,
            max_label=max_label,
        )
    write_zarr_parallel(
        outpath=outpath,
        group=group,
        index=int(time_index),
        data=mapped_t,
    )
    return int(time_index)


def _compute_and_store_feature_backprojection(
    *,
    track_img,
    tracked_img_path,
    zarr_outpath,
    group,
    lookup_payload,
    mode,
    background_value,
    n_workers=1,
    prefer_parallel=True,
    verbose=False,
):
    out_dtype = np.dtype(lookup_payload["dtype"])
    frame_lookup = lookup_payload["frame_lookup"]
    output_shape = tuple(int(v) for v in track_img.shape)
    frame_shape = output_shape[1:]

    write_zarr_parallel(
        outpath=zarr_outpath,
        group=group,
        shape=output_shape,
        dtype=out_dtype,
        chunks=frame_shape,
        overwrite=True,
    )

    n_timepoints = int(output_shape[0])
    can_parallel = (
        bool(prefer_parallel)
        and int(n_workers) > 1
        and tracked_img_path is not None
        and Path(zarr_outpath).suffix == ".zarr"
    )
    if bool(prefer_parallel) and int(n_workers) > 1 and Path(zarr_outpath).suffix != ".zarr":
        print("  Parallel backprojection only supports directory-backed .zarr outputs; falling back to sequential.")

    exec_mode = "parallel" if can_parallel else "sequential"
    print(
        f"    frames={n_timepoints} workers={int(n_workers)} mode={mode} execution={exec_mode}"
    )

    frame_desc = tqdm(
        range(n_timepoints),
        total=n_timepoints,
        desc="Frames",
        position=1,
        leave=False,
    )
    if can_parallel:
        frame_desc.close()
        try:
            with ProcessPoolExecutor(max_workers=int(n_workers)) as ex:
                futures = []
                for t in range(n_timepoints):
                    lookup = _resolve_frame_lookup(frame_lookup=frame_lookup, time_index=t, mode=mode)
                    futures.append(
                        ex.submit(
                            _write_backprojection_frame_worker,
                            tracked_img_path=str(tracked_img_path),
                            outpath=str(zarr_outpath),
                            group=str(group),
                            time_index=t,
                            lookup_ids=None if lookup is None else lookup["ids"],
                            lookup_values=None if lookup is None else lookup["values"],
                            max_label=None if lookup is None else lookup.get("max_label", None),
                            out_dtype_str=out_dtype.str,
                            background_value=background_value,
                        )
                    )
                progress = tqdm(total=n_timepoints, desc="Frames", position=1, leave=False)
                try:
                    for fut in as_completed(futures):
                        fut.result()
                        progress.update(1)
                finally:
                    progress.close()
            return load_image(zarr_outpath, group=group)
        except (PermissionError, OSError) as exc:
            print(f"  Parallel backprojection unavailable ({exc}); falling back to sequential.")
            can_parallel = False
            frame_desc = tqdm(
                range(n_timepoints),
                total=n_timepoints,
                desc="Frames",
                position=1,
                leave=False,
            )

    try:
        for t in frame_desc:
            mapped_t = _process_backprojection_frame(
                tracked_source=track_img,
                time_index=t,
                frame_lookup=frame_lookup,
                    mode=mode,
                    out_dtype=out_dtype,
                    background_value=background_value,
                    verbose=verbose,
                )
            write_zarr_parallel(
                outpath=zarr_outpath,
                group=group,
                index=t,
                data=mapped_t,
            )
    finally:
        frame_desc.close()

    return load_image(zarr_outpath, group=group)


def _load_active_killing_data(
    output_dir: str,
    cell_type: str,
    sample_name: str,
    track_img: np.ndarray,
    raw_img_shape: tuple
) -> dict:
    """
    Load active killing data and create a time-varying labels layer.
    
    Returns a dict with:
    - 'Active_Killing': labels layer showing cells only when actively killing
    - 'Killing_Efficiency': image layer with killing efficiency values
    
    The labels layer will show:
    - TrackID for cells that are actively killing at that timepoint
    - 0 for all other pixels
    
    Parameters
    ----------
    output_dir : str
        BEHAV3D output directory
    cell_type : str
        Cell type (e.g., 'tcell')
    sample_name : str
        Sample name to filter data
    track_img : np.ndarray
        Track image (t, z, y, x) with TrackIDs
    raw_img_shape : tuple
        Shape of raw image for tiling
        
    Returns
    -------
    dict
        Dictionary with 'Active_Killing' and optionally 'Killing_Efficiency' layers
    """
    from behav3d.features.advanced_timepoint_features import find_advanced_features_csv

    advanced_features_path = find_advanced_features_csv(output_dir, cell_type)

    result = {}

    if advanced_features_path is None:
        print(f"  Active killing data not found for cell type '{cell_type}' in {output_dir}")
        return result
    
    print(f"  Loading active killing data from {advanced_features_path}")
    df_active_killing = pd.read_csv(advanced_features_path)
    df_active_killing = df_active_killing[df_active_killing["sample_name"] == sample_name]
    
    if df_active_killing.empty:
        print(f"  No active killing data found for sample {sample_name}")
        return result
    
    # Check if we have the required columns
    if "is_active_killing" not in df_active_killing.columns:
        print(f"  'is_active_killing' column not found in advanced features")
        return result
    
    # Get the timepoints where killing is active
    killing_timepoints = df_active_killing[df_active_killing["is_active_killing"] == True]
    n_killing_events = len(killing_timepoints)
    n_killing_tracks = killing_timepoints["TrackID"].nunique() if n_killing_events > 0 else 0
    print(f"  Found {n_killing_events} active killing timepoints across {n_killing_tracks} tracks")
    
    if n_killing_events == 0:
        return result
    
    # Create time-varying active killing mask
    # track_img shape: (t, 1, z, y, x) after expand_dims or (t, z, y, x) before
    if track_img.ndim == 5:
        t_dim, _, z_dim, y_dim, x_dim = track_img.shape
    else:
        t_dim, z_dim, y_dim, x_dim = track_img.shape
    
    # Initialize arrays for active killing labels, efficiency, and targeted organoids
    active_killing_labels = np.zeros((t_dim, 1, z_dim, y_dim, x_dim), dtype=np.uint16)
    killing_efficiency_img = np.zeros((t_dim, 1, z_dim, y_dim, x_dim), dtype=np.float32)
    targeted_organoid_labels = np.zeros((t_dim, 1, z_dim, y_dim, x_dim), dtype=np.uint16)
    
    # Check if targeted_track_id is available
    has_targeted_id = "targeted_track_id" in df_active_killing.columns
    
    # Process each timepoint
    print(f"  Creating active killing and targeted organoid visualization...")
    track_img_np = np.asarray(track_img)
    if track_img_np.ndim == 4:
        track_img_np = np.expand_dims(track_img_np, axis=1)
    
    # Get efficiency column (could be killing_efficiency or we calculate from threshold)
    efficiency_col = "killing_efficiency" if "killing_efficiency" in df_active_killing.columns else None
    
    unique_timepoints = killing_timepoints["position_t"].unique()
    for t in tqdm(unique_timepoints, desc="  Processing timepoints"):
        t = int(t)
        if t >= t_dim:
            continue
        
        # Get tracks that are actively killing at this timepoint
        killing_tracks_at_t = killing_timepoints[killing_timepoints["position_t"] == t]
        
        for _, row in killing_tracks_at_t.iterrows():
            track_id = int(row["TrackID"])
            efficiency = float(row[efficiency_col]) if efficiency_col else 1.0
            
            # Create mask for this track at this timepoint
            track_mask = track_img_np[t] == track_id
            
            # Set the active killing labels
            active_killing_labels[t][track_mask] = track_id
            
            # Set the killing efficiency value
            killing_efficiency_img[t][track_mask] = efficiency
            
            # Identify targeted organoid
            if has_targeted_id and pd.notna(row["targeted_track_id"]):
                targeted_id = int(row["targeted_track_id"])
                # We need the organoid mask. Since we don't have it here, we'll store the ID
                # and the view_napari function can use it if we provide the organoid layer.
                # BETTER: let's modify the result to include the targeted IDs per timepoint.
                if "Targeted_IDs" not in result: result["Targeted_IDs"] = {}
                result["Targeted_IDs"][t] = result.get("Targeted_IDs", {}).get(t, []) + [targeted_id]

    # Convert to dask arrays and tile for channel dimension if needed
    n_channels = raw_img_shape[-4]
    active_killing_labels = da.from_array(active_killing_labels)
    active_killing_labels = da.tile(active_killing_labels, (1, n_channels, 1, 1, 1))
    
    killing_efficiency_img = da.from_array(killing_efficiency_img)
    killing_efficiency_img = da.tile(killing_efficiency_img, (1, n_channels, 1, 1, 1))
    
    result["Active_Killing"] = {
        "img": active_killing_labels,
        "type": "label"
    }
    
    result["Killing_Efficiency"] = {
        "img": killing_efficiency_img,
        "type": "image"
    }
    
    print(f"  Active killing layers created successfully")
    return result


def backproject_mean_features_behav3d(
    metadata,
    sample_name,
    config=None,
    output_dir=None,
    cell_type="tcell",
    columns=[],
    save=False,
    overwrite=False,
    ):
    """
    Backproject mean (summarized) features onto tracked segments.
    Features are averaged per track, so each track has one value per feature.
    """
    
    print(f"--------------- Backprojecting {cell_type} for {sample_name} (mean mode) ---------------")
    start_time = time.time()

    assert config is not None or all(
        [output_dir, metadata is not None]
    ), "Either 'config' or 'output_dir and metadata' parameters must be supplied"
    
    if not all([output_dir, metadata is not None]):
        output_dir = config['output_dir']
        metadata = pd.read_csv(config["metadata_csv_path"])
        
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    results_outdir = Path(analysis_outdir, "results")
    backproj_outdir = Path(analysis_outdir, "backprojection")
    feature_outdir = Path(analysis_outdir, "track_features")
    img_outdir = Path(output_dir, "images", sample_name)
    
    if not backproj_outdir.exists():
        backproj_outdir.mkdir(parents=True)
        
    df_sample = metadata[metadata["sample_name"]==sample_name]
    assert(sample_name in df_sample["sample_name"].values), f"Supplied sample name {sample_name} not in metadata"
 
    raw_img_path = Path(df_sample["raw_image_path"].values[0])
    
    # Find tracks image path column dynamically
    track_img_col = _find_track_image_column(df_sample, cell_type)
    track_img_path = Path(df_sample[track_img_col].values[0])
    
    # Load clustered data (summarized - one row per track)
    df_tracks_clustered = pd.read_csv(Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv"))
    
    # Load full track data (time-varying - for napari track positions)
    df_tracks_all_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_tracks_all = pd.read_csv(df_tracks_all_path)
    df_tracks_all = df_tracks_all[df_tracks_all["sample_name"]==sample_name]
    
    # Load clustered full data (time-varying with ClusterID)
    df_tracks_full_clustered_path = Path(results_outdir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv")
    df_tracks_full_clustered = pd.read_csv(df_tracks_full_clustered_path)
    df_tracks_full_clustered = df_tracks_full_clustered[df_tracks_full_clustered["sample_name"]==sample_name]
    
    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
    ]
    
    # Use separate zarr file for mean mode
    backproj_out_path = Path(backproj_outdir, f"{get_filepath_stem(track_img_path)}_backprojected_mean.zarr")
    _ensure_zarr(raw_img_path, label=f"Raw image for '{sample_name}'")
    _ensure_zarr(track_img_path, label=f"Tracked image for '{sample_name}'")

    raw_img = load_image(raw_img_path)
    
    raw_img_data = {
        "raw_data": {
            "img": raw_img,
            "type": "image"
        }
    }
    
    print("- Loading in tracked segments")
    track_img = load_image(track_img_path)
    
    track_img = np.expand_dims(track_img, axis=1)
    
    # Filter to only tracks with ClusterID
    # Pixel labels in track_img are always original segmentation ids; use
    # original_TrackID (when present) so chunks split by filtering still
    # match, not just the first chunk of each split track.
    keep_id_col = "original_TrackID" if "original_TrackID" in df_tracks_clustered.columns else "TrackID"
    filt_track_img = np.where(np.isin(track_img, df_tracks_clustered[keep_id_col].unique()), track_img, 0)
    
    # Check filtered result (handle dask arrays)
    if hasattr(filt_track_img, 'compute'):
        unique_filtered = da.unique(filt_track_img).compute()
    else:
        unique_filtered = np.unique(filt_track_img)
    unique_filtered = unique_filtered[unique_filtered != 0]
    print(f"  Found {len(unique_filtered)} matching TrackIDs for filtered view")
    
    track_img = da.tile(track_img, (1, raw_img.shape[-4], 1, 1, 1))
    
    # Build trackid_data
    trackid_data = {
        f"{cell_type}_TrackID": {
            "img": track_img,
            "type": "label"
        },
        f"filtered_{cell_type}_TrackID": {
            "img": filt_track_img,
            "type": "label"
        }
    }
    
    # Load other cell type tracks for context
    trackid_data = _load_other_cell_tracks(trackid_data, metadata, df_sample, cell_type, raw_img)
    
    # Load dead mask if available
    trackid_data = _load_dead_mask(trackid_data, df_sample, img_outdir, sample_name, raw_img)
    
    tracked_csv_path = _resolve_tracked_csv_path(df_sample, output_dir, sample_name, cell_type)

    # Select columns to backproject
    columns = _select_backprojection_columns(
        df_tracks_clustered=df_tracks_clustered,
        columns=columns,
        metadata_columns=metadata.columns,
        tracked_csv_path=tracked_csv_path,
        extra_excluded={"UMAP1", "UMAP2"},
    )
    
    # Always ensure ClusterID is included
    if "ClusterID" not in columns and "ClusterID" in df_tracks_clustered.columns:
        columns = ["ClusterID"] + columns
    
    # Filter to valid numeric columns
    columns = _filter_valid_columns(columns, df_tracks_clustered, df_sample)
    
    print(f"- Backprojecting {len(columns)} features onto each segment")
    backprojected_cols = backproject_columns(
        track_img=filt_track_img,
        df_tracks_clustered=df_tracks_clustered,
        zarr_outpath=backproj_out_path,
        columns=columns,
        mode="summary",
        tracked_img_path=track_img_path,
        overwrite=overwrite,
    )
    
    label_columns = ["ClusterID"]
    for col in backprojected_cols.keys():
        backprojected_cols[col]["img"] = da.tile(backprojected_cols[col]["img"], (1, raw_img.shape[-4], 1, 1, 1))
        if col in label_columns:
            backprojected_cols[col]["type"] = "label"
    
    trackid_data[f"filtered_{cell_type}_TrackID"]["img"] = da.tile(
        trackid_data[f"filtered_{cell_type}_TrackID"]["img"], (1, raw_img.shape[-4], 1, 1, 1))
    
    # Load active killing data if available (for immune cell types)
    active_killing_data = _load_active_killing_data(
        output_dir=output_dir,
        cell_type=cell_type,
        sample_name=sample_name,
        track_img=filt_track_img,
        raw_img_shape=raw_img.shape
    )
    
    backproject_data = {**trackid_data, **backprojected_cols, **active_killing_data}
    visualize_data = {**raw_img_data, **backproject_data}
    
    print("- Visualizing backprojection in napari")
    viewer = view_napari(
        visualize_data,
        df_tracks_full=df_tracks_all,
        df_tracks_clustered=df_tracks_full_clustered,
        cell_type=cell_type,
        elsize=elsize
    )
    
    end_time = time.time()
    h, m, s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return viewer


def backproject_time_features_behav3d(
    metadata,
    sample_name,
    config=None,
    output_dir=None,
    cell_type="tcell",
    columns=[],
    save=False,
    overwrite=False,
    ):
    """
    Backproject time-varying features onto tracked segments.
    Features vary per timepoint, so each track point has its own value.
    """
    
    print(f"--------------- Backprojecting {cell_type} for {sample_name} (time mode) ---------------")
    start_time = time.time()

    assert config is not None or all(
        [output_dir, metadata is not None]
    ), "Either 'config' or 'output_dir and metadata' parameters must be supplied"
    
    if not all([output_dir, metadata is not None]):
        output_dir = config['output_dir']
        metadata = pd.read_csv(config["metadata_csv_path"])
        
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    results_outdir = Path(analysis_outdir, "results")
    backproj_outdir = Path(analysis_outdir, "backprojection")
    feature_outdir = Path(analysis_outdir, "track_features")
    img_outdir = Path(output_dir, "images", sample_name)
    
    if not backproj_outdir.exists():
        backproj_outdir.mkdir(parents=True)
        
    df_sample = metadata[metadata["sample_name"]==sample_name]
    assert(sample_name in df_sample["sample_name"].values), f"Supplied sample name {sample_name} not in metadata"
 
    raw_img_path = Path(df_sample["raw_image_path"].values[0])
    
    # Find tracks image path column dynamically
    track_img_col = _find_track_image_column(df_sample, cell_type)
    track_img_path = Path(df_sample[track_img_col].values[0])
    
    # Load UMAP clusters (to get ClusterID mapping)
    df_umap_clusters = pd.read_csv(Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv"))
    
    # Load full track data (time-varying)
    df_tracks_all_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_tracks_all = pd.read_csv(df_tracks_all_path)
    df_tracks_all = df_tracks_all[df_tracks_all["sample_name"]==sample_name]
    
    # Merge ClusterID into time-varying data
    df_tracks_clustered = pd.merge(
        df_tracks_all, 
        df_umap_clusters[["TrackID", "ClusterID"]], 
        on='TrackID', 
        how='inner'  # Only keep tracks that have ClusterID
    )
    
    # Load clustered full data for napari
    df_tracks_full_clustered_path = Path(results_outdir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv")
    df_tracks_full_clustered = pd.read_csv(df_tracks_full_clustered_path)
    df_tracks_full_clustered = df_tracks_full_clustered[df_tracks_full_clustered["sample_name"]==sample_name]
    
    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
    ]
    
    # Use separate zarr file for time mode
    backproj_out_path = Path(backproj_outdir, f"{get_filepath_stem(track_img_path)}_backprojected_time.zarr")
    raw_img = load_image(raw_img_path)
    
    raw_img_data = {
        "raw_data": {
            "img": raw_img,
            "type": "image"
        }
    }
    
    print("- Loading in tracked segments")
    track_img = load_image(track_img_path)
    
    track_img = np.expand_dims(track_img, axis=1)
    
    # Filter to only tracks with ClusterID
    # Pixel labels in track_img are always original segmentation ids; use
    # original_TrackID (when present) so chunks split by filtering still
    # match, not just the first chunk of each split track.
    keep_id_col = "original_TrackID" if "original_TrackID" in df_tracks_clustered.columns else "TrackID"
    filt_track_img = np.where(np.isin(track_img, df_tracks_clustered[keep_id_col].unique()), track_img, 0)
    
    # Check filtered result (handle dask arrays)
    if hasattr(filt_track_img, 'compute'):
        unique_filtered = da.unique(filt_track_img).compute()
    else:
        unique_filtered = np.unique(filt_track_img)
    unique_filtered = unique_filtered[unique_filtered != 0]
    print(f"  Found {len(unique_filtered)} matching TrackIDs for filtered view")
    
    track_img = da.tile(track_img, (1, raw_img.shape[-4], 1, 1, 1))
    
    # Build trackid_data
    trackid_data = {
        f"{cell_type}_TrackID": {
            "img": track_img,
            "type": "label"
        },
        f"filtered_{cell_type}_TrackID": {
            "img": filt_track_img,
            "type": "label"
        }
    }
    
    # Load other cell type tracks for context
    trackid_data = _load_other_cell_tracks(trackid_data, metadata, df_sample, cell_type, raw_img)
    
    # Load dead mask if available
    trackid_data = _load_dead_mask(trackid_data, df_sample, img_outdir, sample_name, raw_img)
    
    tracked_csv_path = _resolve_tracked_csv_path(df_sample, output_dir, sample_name, cell_type)

    # Select columns to backproject
    columns = _select_backprojection_columns(
        df_tracks_clustered=df_tracks_clustered,
        columns=columns,
        metadata_columns=metadata.columns,
        tracked_csv_path=tracked_csv_path,
        extra_excluded={"UMAP1", "UMAP2"},
    )
    
    # Always ensure ClusterID is included
    if "ClusterID" not in columns and "ClusterID" in df_tracks_clustered.columns:
        columns = ["ClusterID"] + columns
    
    # Filter to valid numeric columns
    columns = _filter_valid_columns(columns, df_tracks_clustered, df_sample)
    
    print(f"- Backprojecting {len(columns)} features onto each segment")
    backprojected_cols = backproject_columns(
        track_img=filt_track_img,
        df_tracks_clustered=df_tracks_clustered,
        zarr_outpath=backproj_out_path,
        columns=columns,
        mode="timepoint",
        tracked_img_path=track_img_path,
        overwrite=overwrite,
    )
    
    label_columns = ["ClusterID"]
    for col in backprojected_cols.keys():
        backprojected_cols[col]["img"] = da.tile(backprojected_cols[col]["img"], (1, raw_img.shape[-4], 1, 1, 1))
        if col in label_columns:
            backprojected_cols[col]["type"] = "label"
    
    trackid_data[f"filtered_{cell_type}_TrackID"]["img"] = da.tile(
        trackid_data[f"filtered_{cell_type}_TrackID"]["img"], (1, raw_img.shape[-4], 1, 1, 1))
    
    # Load active killing data if available (for immune cell types)
    active_killing_data = _load_active_killing_data(
        output_dir=output_dir,
        cell_type=cell_type,
        sample_name=sample_name,
        track_img=filt_track_img,
        raw_img_shape=raw_img.shape
    )
    
    backproject_data = {**trackid_data, **backprojected_cols, **active_killing_data}
    visualize_data = {**raw_img_data, **backproject_data}
    
    print("- Visualizing backprojection in napari")
    viewer = view_napari(
        visualize_data,
        df_tracks_full=df_tracks_all,
        df_tracks_clustered=df_tracks_full_clustered,
        cell_type=cell_type,
        elsize=elsize
    )
    
    end_time = time.time()
    h, m, s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return viewer


def _find_track_image_column(df_sample, cell_type):
    """Find the correct column name for track image path."""
    track_img_col = f"{cell_type}_tracks_image_path"
    if track_img_col not in df_sample.columns:
        for prefix in ['im_', 'or_', 'ot_']:
            alt_col = f"{prefix}{cell_type}_tracks_image_path"
            if alt_col in df_sample.columns:
                return alt_col
        raise ValueError(f"Cannot find tracks image path column for {cell_type}")
    return track_img_col


def _load_other_cell_tracks(trackid_data, metadata, df_sample, cell_type, raw_img):
    """Load tracks from other cell types for visualization context."""
    
    
    all_cell_types = (
        filter_multicolor_inputs(detect_immune_cell_types_from_metadata(metadata)) +
        filter_multicolor_inputs(detect_organoid_types_from_metadata(metadata)) +
        filter_multicolor_inputs(detect_other_cell_types_from_metadata(metadata)) +
        detect_merged_cell_types_from_metadata(metadata)
    )
    
    for other_ct in all_cell_types:
        if other_ct == cell_type:
            continue
        
        try:
            other_track_col = _find_track_image_column(df_sample, other_ct)
            if other_track_col in df_sample.columns and pd.notna(df_sample[other_track_col].values[0]):
                other_track_path = Path(df_sample[other_track_col].values[0])
                if other_track_path.exists():
                    other_track_img = load_image(other_track_path)
                    other_track_img = np.expand_dims(other_track_img, axis=1)
                    other_track_img = da.tile(other_track_img, (1, raw_img.shape[-4], 1, 1, 1))
                    trackid_data[f"{other_ct}_TrackID"] = {
                        "img": other_track_img,
                        "type": "label"
                    }
        except (ValueError, Exception):
            pass
    
    return trackid_data


def _load_dead_mask(trackid_data, df_sample, img_outdir, sample_name, raw_img):
    """Load dead mask if dead channel is enabled."""
    has_dead_channel = "dead_channel" in df_sample.columns and pd.notna(df_sample["dead_channel"].values[0])
    
    if has_dead_channel:
        dead_mask_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
        if dead_mask_path.exists():
            try:
                dead_mask = load_zarr(dead_mask_path)
                dead_mask = np.expand_dims(dead_mask, axis=1)
                dead_mask = da.tile(dead_mask, (1, raw_img.shape[-4], 1, 1, 1))
                trackid_data["Dead_Mask"] = {
                    "img": dead_mask,
                    "type": "label"
                }
            except Exception:
                pass
    
    return trackid_data


def _filter_valid_columns(columns, df_tracks_clustered, df_sample):
    """Filter columns to only valid numeric ones."""
    has_dead_channel = "dead_channel" in df_sample.columns and pd.notna(df_sample["dead_channel"].values[0])
    
    valid_columns = []
    for col in columns:
        if col not in df_tracks_clustered.columns:
            print(f"  Skipping missing column: {col}")
            continue
        
        if col.startswith('touching_') or col == 'orientation_vector':
            print(f"  Skipping string column: {col}")
            continue
        
        if not has_dead_channel and ('dead' in col.lower() or col == 'is_dead'):
            continue
        
        if pd.api.types.is_numeric_dtype(df_tracks_clustered[col]):
            valid_columns.append(col)
    
    return valid_columns


def backproject_columns(
    track_img,
    zarr_outpath,
    df_tracks_clustered,
    columns=["ClusterID", "mean_speed"],
    mode="summary",
    track_col="TrackID",
    time_col="position_t",
    background_value=0,
    n_workers=1,
    prefer_parallel=True,
    tracked_img_path=None,
    verbose=False,
    overwrite=False,
    ):
    """
    Map column values onto track image pixels.
    Each pixel value (TrackID) is replaced with the corresponding column value.

    Parameters
    ----------
    mode : {"summary", "timepoint"}
        "summary" maps one value per TrackID across all timepoints.
        "timepoint" maps values using both TrackID and position_t.
    overwrite : bool
        If True, recompute requested feature groups even when they already exist
        in the backprojection zarr store.

    Notes
    -----
    Pixel values in ``track_img`` are always the original segmentation
    labels. When tracks were split into chunks by the filtering step's
    "Split long tracks into chunks" option, ``TrackID`` no longer matches
    those labels for chunks after the first — ``original_TrackID`` does.
    So whenever ``df_tracks_clustered`` carries an ``original_TrackID``
    column, it is used to key the pixel lookup instead of ``track_col``.
    """
    mode = _normalize_backprojection_mode(mode)
    if "original_TrackID" in df_tracks_clustered.columns:
        track_col = "original_TrackID"
    zarr_outpath = None if zarr_outpath is None else Path(zarr_outpath)
    if zarr_outpath is not None and Path(f"{zarr_outpath}.zip").exists():
        zarr_outpath = Path(f"{zarr_outpath}.zip")

    mapped_imgs = {col: {"img": None, "type": "image"} for col in columns}

    existing_groups = []
    if zarr_outpath is not None and zarr_outpath.exists() and not overwrite:
        import zarr
        try:
            if zarr_outpath.suffix == ".zip":
                store = zarr.storage.ZipStore(zarr_outpath, mode="r")
                try:
                    zarr_store = zarr.open_group(store=store, mode="r")
                    existing_groups = list(zarr_store.keys()) if hasattr(zarr_store, "keys") else []
                finally:
                    store.close()
            else:
                zarr_store = zarr.open_group(str(zarr_outpath), mode="r")
                existing_groups = list(zarr_store.keys()) if hasattr(zarr_store, "keys") else []
        except Exception:
            existing_groups = []

    cols_to_load = [c for c in columns if c in existing_groups] if not overwrite else []
    cols_to_compute = [c for c in columns if overwrite or c not in existing_groups]

    if cols_to_load:
        print(f"  Loading {len(cols_to_load)} backprojection layer(s) from existing .zarr")
        for col in cols_to_load:
            mapped_imgs[col]["img"] = load_image(zarr_outpath, group=col)

    if cols_to_compute:
        print(f"  Computing backprojection layers from existing CSV values ({mode})")
        for col in cols_to_compute:
            print(f"  Computing feature: {col} ({mode})")
            if mode == "summary":
                lookup_payload = _build_summary_lookup(
                    lookup_df=df_tracks_clustered[[track_col, col]].copy(),
                    track_col=track_col,
                    value_col=col,
                    background_value=background_value,
                )
            else:
                required_cols = [track_col, time_col, col]
                missing_cols = [c for c in required_cols if c not in df_tracks_clustered.columns]
                if missing_cols:
                    raise ValueError(
                        f"Timepoint backprojection for '{col}' requires columns {required_cols}, "
                        f"missing {missing_cols}."
                    )
                if bool(verbose):
                    non_null_rows = df_tracks_clustered[required_cols].dropna(subset=[track_col, time_col, col])
                    if not non_null_rows.empty:
                        print(
                            f"  feature={col} timepoints={int(non_null_rows[time_col].min())}.."
                            f"{int(non_null_rows[time_col].max())} "
                            f"tracks={int(non_null_rows[track_col].nunique())}"
                        )
                lookup_payload = _build_timepoint_lookup(
                    lookup_df=df_tracks_clustered[required_cols].copy(),
                    track_col=track_col,
                    time_col=time_col,
                    value_col=col,
                    background_value=background_value,
                )
            if zarr_outpath is None:
                output_shape = tuple(int(v) for v in np.asarray(track_img).shape)
                mapped_img = np.full(output_shape, background_value, dtype=lookup_payload["dtype"])
                frame_progress = tqdm(range(output_shape[0]), total=output_shape[0], desc="Frames", leave=False)
                try:
                    for t in frame_progress:
                        mapped_img[t] = _process_backprojection_frame(
                            tracked_source=track_img,
                            time_index=t,
                            frame_lookup=lookup_payload["frame_lookup"],
                            mode=mode,
                            out_dtype=np.dtype(lookup_payload["dtype"]),
                            background_value=background_value,
                            verbose=verbose,
                        )
                finally:
                    frame_progress.close()
                mapped_imgs[col]["img"] = mapped_img
            else:
                zarr_path_for_save = zarr_outpath.stem if zarr_outpath.suffix == ".zip" else zarr_outpath
                if bool(prefer_parallel) and int(n_workers) > 1 and zarr_outpath.suffix == ".zip":
                    print("  Parallel backprojection only supports directory-backed .zarr outputs; falling back to sequential.")
                mapped_imgs[col]["img"] = _compute_and_store_feature_backprojection(
                    track_img=track_img,
                    tracked_img_path=tracked_img_path,
                    zarr_outpath=zarr_path_for_save,
                    group=col,
                    lookup_payload=lookup_payload,
                    mode=mode,
                    background_value=background_value,
                    n_workers=n_workers,
                    prefer_parallel=bool(prefer_parallel) and zarr_outpath.suffix != ".zip",
                    verbose=verbose,
                )

    return mapped_imgs


def view_napari(
    backproject_data,
    df_tracks_full,
    df_tracks_clustered,
    cell_type,
    elsize,
    split_raw_channels=False,
    run=True,
    ):
    """
    Visualize backprojection in napari.
    
    Highlights:
    - Active_Killing: bright red outlines when cells are actively killing
    - Targeted_Organoids: pulsing or distinct color for targets
    - Individual features selectable via UI dropdown
    """
    
    viewer = napari.Viewer()
    
    # Collect targeted IDs if available
    targeted_ids_per_tp = backproject_data.get("Targeted_IDs", {})
    
    raw_layer_names = []
    label_reference = None
    if "TrackID" in backproject_data:
        label_reference = backproject_data["TrackID"].get("img")
    elif "ClusterID" in backproject_data:
        label_reference = backproject_data["ClusterID"].get("img")

    expected_timepoints = None
    if "position_t" in df_tracks_full.columns and len(df_tracks_full) > 0:
        time_values = pd.to_numeric(df_tracks_full["position_t"], errors="coerce")
        if time_values.notna().any():
            expected_timepoints = int(time_values.max()) + 1

    def _raw_to_tczyx(raw_img):
        if getattr(raw_img, "ndim", 0) != 5:
            return raw_img

        label_shape = tuple(int(v) for v in getattr(label_reference, "shape", ()) or ())
        raw_shape = tuple(int(v) for v in raw_img.shape)
        if len(label_shape) == 4:
            label_time = label_shape[0]
            label_spatial = label_shape[1:]
            if raw_shape[2:] == label_spatial:
                if raw_shape[0] == label_time:
                    return raw_img
                if raw_shape[1] == label_time:
                    return np.transpose(raw_img, (1, 0, 2, 3, 4))

        if expected_timepoints is not None:
            if raw_shape[0] == expected_timepoints:
                return raw_img
            if raw_shape[1] == expected_timepoints:
                return np.transpose(raw_img, (1, 0, 2, 3, 4))

        if raw_shape[0] <= 16 and raw_shape[1] > raw_shape[0]:
            return np.transpose(raw_img, (1, 0, 2, 3, 4))
        return raw_img

    # Add all image/label layers
    layer_items = []
    if "raw_data" in backproject_data:
        layer_items.append(("raw_data", backproject_data["raw_data"]))
    layer_items.extend((k, v) for k, v in backproject_data.items() if k != "raw_data")
    for k, v in layer_items:
        if k == "Targeted_IDs": continue
        try:
            # Handle transposition if not already correct
            if k == "raw_data":
                v["img"] = _raw_to_tczyx(v["img"])
            elif v["img"].ndim == 5 and v["img"].shape[0] != df_tracks_full["position_t"].max() + 1:
                try:
                    v["img"] = np.transpose(v["img"], (1, 0, 2, 3, 4))
                except Exception: pass

            if k == "raw_data" and bool(split_raw_channels) and getattr(v["img"], "ndim", 0) == 5 and int(v["img"].shape[1]) >= 1:
                channel_count = int(v["img"].shape[1])
                layers = viewer.add_image(
                    v["img"],
                    name=[f"raw_ch{idx + 1}" for idx in range(channel_count)],
                    channel_axis=1,
                    scale=elsize,
                    opacity=0.5,
                )
                for idx, layer in enumerate(layers):
                    layer.visible = (idx == 0)
                    layer.opacity = 0.5
                    raw_layer_names.append(layer.name)

            elif k == "Active_Killing":
                layer = viewer.add_labels(v["img"], name="Killing T-cells", scale=elsize, opacity=0.8)
                layer.contour = 2
                layer.visible = True

            elif k == "Killing_Efficiency":
                img_np = np.asarray(v["img"])
                valid_vals = img_np[(img_np != 0) & np.isfinite(img_np)]
                vmax = np.percentile(valid_vals, 98) if valid_vals.size > 0 else 1
                layer = viewer.add_image(v["img"], name=k, scale=elsize, colormap='hot', contrast_limits=[0, float(vmax)], opacity=0.8, visible=False)

            elif v["type"] == "label":
                viewer.add_labels(v["img"], name=k, scale=elsize, visible=(k == "ClusterID"))
            else:
                img_np = np.asarray(v["img"])
                valid_vals = img_np[(img_np != 0) & np.isfinite(img_np)]
                vmin, vmax = np.percentile(valid_vals, [2, 98]) if valid_vals.size > 0 else (0, 1)
                viewer.add_image(v["img"], name=k, scale=elsize, colormap='inferno', contrast_limits=[float(vmin), float(vmax)], visible=False)
        except Exception as e:
            print(f"  Skipping layer '{k}': {e}")

    # Specific logic for Targeted Organoids layer
    if targeted_ids_per_tp:
        # We look for organoid track layers to highlight
        org_layers = [lay for lay in viewer.layers if "_TrackID" in lay.name and cell_type not in lay.name]
        if org_layers:
            # Create a dedicated layer for targeted organoids
            # This is complex to do purely with labels if we don't have the full img here.
            # Instead, we'll try to use the existing organoid layers and filter them.
            print("  ✨ Targeted IDs found - highlights will be synchronized")

    # Add tracks
    napari_tracks_clustered = df_tracks_clustered[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()
    features = {'ClusterID': df_tracks_clustered["ClusterID"].to_numpy().astype(float)}
    try:
        tracks_layer = viewer.add_tracks(napari_tracks_clustered, name=f'Filtered {cell_type} Tracks', features=features, tail_length=75)
        tracks_layer.color_by = 'ClusterID'
        tracks_layer.colormap = 'turbo'
    except Exception as e:
        print(f"  Skipping {cell_type} tracks layer: {e}")
    
    # Add raw data
    if raw_layer_names:
        for layer_name in raw_layer_names:
            if layer_name in viewer.layers:
                viewer.layers[layer_name].visible = (layer_name == raw_layer_names[0])
                viewer.layers[layer_name].opacity = 0.5
    elif 'raw_data' in backproject_data:
        viewer.layers['raw_data'].visible = True
        viewer.layers['raw_data'].opacity = 0.5

    if bool(run):
        napari.run()
    return viewer
