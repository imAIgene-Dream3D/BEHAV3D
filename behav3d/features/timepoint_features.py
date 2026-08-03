### TODO Perhaps set a "static" speed based on quantiles of mean_speed to define static and actively interacting cells
### However this excludes a lot of cells to actively interact, now we take the fastest of the contacting T cells as actively interacting

### TODO WHy is tcell_cntact TRUE everywhere?
"""
This script calculates the features of tracks for BEHAV3D analysis.

Flly flexible - works with ANY cell types defined in metadata using prefixes:
- or_ = Organoid types 
- im_ = Immune types 
- ot_ = Other types 

-------------------------------------
--------------- INPUT ---------------
-------------------------------------

.csv containing the following columns:
- TrackID       (The ID of the track a segments belongs to)
- SegmentID     (The unique ID of the segment of a specific timepoint)
- position_t    (The timepoint of the segment)
- position_z
- position_y
- position_x

-------------------------------------
--------------- OUTPUT --------------
-------------------------------------

# Features of a track at each timepoint per sample in the metadata csv (.csv)
- See "FEATURES TRACKS" below

# Combined summarized features for each track for all samples in metadata csv (.csv)
- Merged DataFrame with all features across samples

-------------------------------------
---------- FEATURES TRACKS ---------- 
-------------------------------------

DYNAMIC CONTACT FEATURES (generated for ALL cell types in metadata):
- {cell_type}_contact               (bool) - Pixel-based contact (1.73 diagonal)
- {cell_type}_contact_on_distance   (bool) - Real distance-based contact
- touching_{cell_type}s             (str)  - Comma-separated TrackIDs
- any_organoid_contact              (bool) - Any organoid-type pixel-based contact
- any_organoid_contact_on_distance  (bool) - Any organoid-type real distance-based contact
- any_immune_cell_contact           (bool) - Any immune-type pixel-based contact
- any_immune_cell_contact_on_distance (bool) - Any immune-type real distance-based contact
- active_{cell_type}_contact        (bool) - Active interaction (works for any cell type)

INVASIVENESS FEATURES (immune vs organoids):
- {org_type}_invasiveness           (bool) - True if >= 50% surface in contact
- {org_type}_invasiveness_perc      (float) - Percentage (0-100) of surface in contact
- any_org_invasiveness              (bool) - True if invasive against any organoid type
- any_org_invasiveness_perc         (float) - Max invasiveness percentage across organoid types

MORPHOLOGY FEATURES:
- nr_pixels, volume, bbox_volume, elongation, extent, equivalent_diameter
- major_axis_length, minor_axis_length, surface_area, sphericity
- convex_volume, orientation_vector

MOVEMENT FEATURES:
- displacement
- cumulative_displacement
- displacement_from_origin
- mean_square_displacement
- speed
- mean_speed
- interpolated
- time

INTENSITY FEATURES:
- mean_intensity_ch{N}  (for each channel in raw image)
- mean_dead_dye         (if dead_channel specified in metadata)

DEATH FEATURES:
- percentage_dead_mask
- nr_dead_mask_pixels
- increase_dead_mask
- dead                  (bool flag set when threshold crossed)

-------------------------------------
------ CONTACT FEATURE DETAILS ------
-------------------------------------

### {cell_type}_contact
- True/False
Per segment, uses a pixel-based threshold (1.73 = diagonal distance) without
physical spacing. More lenient threshold for pixel-touching.

### {cell_type}_contact_on_distance  
- True/False
Creates a ZYX cutout with extended range around the segment border and
calculates a distance transform using physical spacing (µm). Any other cell
within "contact_threshold" µm counts as contacting.

### touching_{cell_type}s
- String (comma-separated TrackIDs)
Contains TrackIDs of all contacting cells of this type. Empty string if none.
Self-contact excluded when calculating from the same cell type.

### any_organoid_contact / any_organoid_contact_on_distance
- True/False
Aggregate contact flags across all detected organoid types. These columns are
separate from per-type columns so a literal cell type named "organoid" does not
collide with the aggregate output.

### any_immune_cell_contact / any_immune_cell_contact_on_distance
- True/False
Aggregate contact flags across all detected immune cell types.

### active_{cell_type}_contact
- True/False
Identifies actively interacting cells vs passively contacted cells. When multiple
cells of the same type contact each other, ranks by mean_speed over rolling window
(default 10 timepoints). Only the cell with highest speed is marked as "active" interaction.
Works for any cell type (tcell, macro, nk, organoid, etc.).

### mean_dead_dye
- Float
The mean intensity of the dead dye inside of each segment per timepoint. Calculated based on the 
channel of the dead dye. Supplied as an index of the channel image supplied as "dead_channel" 

### displacement
- Float
The displacement of a segment in a track at each timepoint.

### cumulative_displacement
- Float
The cumulative displacement of a segment in a track at each timepoint. This is the sum of all 
displacement values of that timepoint and all timepoints before it

### displacement_from_origin
- Float
The displacement of a segment at a specific timepoint compared to the position of the first timepoint
of that specific track

### mean_square_displacement
- Float
The mean_square_displacement (MSD) is a measure of the deviation of the position of a particle with 
respect to a reference position over time. It can tell you for example if a cell is moving due to random motion
or some outside force is directing its movement.

### speed
- Float
The speed is the same as "displacement", but now normalized to um/h

### mean_speed
- Float
The mean_speed is a rolling window over previous timepoints (defined by "rolling_meanspeed_window", default 10)
That average these timepoints to get a better indication of the cells actual speed. Used in 
calculating the "active_{cell_type}_contact"

### time
- Float
The time is "position_t" but then normalized to the defined supplied "time_interval" to get real-time
time of each position_t. Converted in the code to hours

### interpolated
- True/False
As some timepoint may be missing in a track, missing values are interpolated from the existing data
This indicates if a timepoint is actually found in the data or interpolated by this script

"""

import numpy as np
import pandas as pd

import warnings

from scipy.ndimage import distance_transform_edt, find_objects
from scipy.spatial.distance import cdist

from skimage.measure import regionprops_table
import argparse
import yaml
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

import seaborn as sns

import ast
import math
import time
import traceback
from behav3d.analysis import smooth_value_over_time
from behav3d.core.utils import get_current_time, format_time, convert_time, convert_distance
from behav3d.core.metadata import is_multicolor_celltype, resolve_immune_track_paths_for_sample
from behav3d.io.images import load_image, load_image_timepoint, convert_raw_file_to_zarr
from tqdm import tqdm
from datetime import datetime

import multiprocessing
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool


def _run_parallel_with_fallback(fn, args_list, n_workers, use_processes=True):
    """
    Run fn over args_list in parallel with n_workers.

    By default uses ThreadPoolExecutor (shared memory, no fork overhead).
    Set use_processes=True for ProcessPoolExecutor with spawn context
    (full CPU parallelism but higher memory per worker).
    """
    if use_processes:
        try:
            ctx = multiprocessing.get_context("spawn")
            with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
                return list(tqdm(ex.map(fn, args_list), total=len(args_list)))
        except BrokenProcessPool:
            raise RuntimeError(
                f"Worker processes crashed while using {n_workers} workers "
                f"(likely out of memory). Please lower the number of workers "
                f"and try again."
            )
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        return list(tqdm(ex.map(fn, args_list), total=len(args_list)))


def _normalize_merge_keys(df, keys=("TrackID", "position_t")):
    """Normalize merge keys to nullable integer dtype for robust joins."""
    for key in keys:
        if key in df.columns:
            df[key] = pd.to_numeric(df[key], errors="coerce").astype("Int64")
    return df


def _segment_touches_border(mask=None, *, slices=None, stack_shape=None):
    """Return True when any foreground voxel touches the 3D image border.

    Can be called with either:
    - mask: full boolean mask (legacy)
    - slices + stack_shape: bounding-box slices from find_objects (fast path)
    """
    if slices is not None and stack_shape is not None:
        for sl, dim_size in zip(slices, stack_shape):
            if sl.start == 0 or sl.stop >= dim_size:
                return True
        return False
    if mask is None or mask.ndim != 3:
        return False
    return bool(
        np.any(mask[0, :, :])
        or np.any(mask[-1, :, :])
        or np.any(mask[:, 0, :])
        or np.any(mask[:, -1, :])
        or np.any(mask[:, :, 0])
        or np.any(mask[:, :, -1])
    )


def _voxel_surface_area(mask, voxel_spacing):
    """Estimate surface area from exposed voxel faces in physical units."""
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 3 or not np.any(mask):
        return np.nan

    vz, vy, vx = map(float, voxel_spacing)
    padded = np.pad(mask, 1, mode="constant", constant_values=False)

    z_faces = np.count_nonzero(padded[1:, :, :] != padded[:-1, :, :])
    y_faces = np.count_nonzero(padded[:, 1:, :] != padded[:, :-1, :])
    x_faces = np.count_nonzero(padded[:, :, 1:] != padded[:, :, :-1])

    return float((z_faces * vy * vx) + (y_faces * vz * vx) + (x_faces * vz * vy))


def _principal_axis_metrics(coords, voxel_spacing):
    """Compute spacing-aware principal-axis lengths and orientation from voxel centers."""
    coords = np.asarray(coords, dtype=float)
    if coords.size == 0:
        nan_vec = [np.nan, np.nan, np.nan]
        return {
            "axis_lengths": np.array([np.nan, np.nan, np.nan], dtype=float),
            "orientation_vector": nan_vec,
            "oblateness": np.nan,
            "prolateness": np.nan,
        }

    spacing = np.asarray(voxel_spacing, dtype=float)
    pts = coords * spacing

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        center = pts.mean(axis=0, keepdims=True)
        X = pts - center
        _, _, Vt = np.linalg.svd(X, full_matrices=False)
        V = Vt.T
        proj = X @ V
        center_lengths = np.ptp(proj, axis=0)
        support_widths = np.abs(V).T @ spacing
        lengths = center_lengths + support_widths

        if (
            V.shape != (3, 3)
            or lengths.shape != (3,)
            or (not np.all(np.isfinite(V)))
            or (not np.all(np.isfinite(lengths)))
        ):
            raise FloatingPointError("invalid principal-axis geometry")

        order = np.argsort(lengths)[::-1]
        axis_lengths = lengths[order]
        V_sorted = V[:, order]

        a, b, c = axis_lengths
        major_vec = V_sorted[:, 0]
        major_norm = float(np.linalg.norm(major_vec))
        orientation_vector = [np.nan, np.nan, np.nan]
        if np.isfinite(major_norm) and major_norm > 0:
            orientation_vector = (major_vec / major_norm).tolist()

        oblateness = np.nan
        prolateness = np.nan
        if np.isfinite(a) and a > 0:
            if np.isfinite(c):
                oblateness = float(1.0 - (c / a))
            if np.isfinite(b):
                prolateness = float(1.0 - (b / a))

    return {
        "axis_lengths": axis_lengths.astype(float),
        "orientation_vector": orientation_vector,
        "oblateness": oblateness,
        "prolateness": prolateness,
    }


def calculate_smoothed_death_columns(df_tracks, window_size):
    death_columns = [
        "percentage_dead_mask",
        "nr_dead_mask_pixels",
        "increase_dead_mask",
    ]
    
    for column in death_columns:
        smoothed_column = f"smoothed_{column}"
        if column not in df_tracks.columns:
            df_tracks[smoothed_column] = np.nan
            continue

        if window_size is None:
            df_tracks[smoothed_column] = df_tracks[column]
        else:
            df_tracks[smoothed_column] = smooth_value_over_time(
                df_tracks,
                column=column,
                window_size=window_size,
                min_periods=1,
                groupby=["TrackID", "sample_name"],
            )

    return df_tracks


def rerun_death_classification(
    output_dir,
    cell_type,
    new_threshold,
    threshold_column="percentage_dead_mask",
):
    """Re-apply :func:`calculate_death` to an existing combined feature CSV
    with a new ``dead_mask_percentage_threshold``.

    This is a fast operation: it does **not** re-run any image-based
    feature extraction. It assumes ``percentage_dead_mask`` (or whichever
    ``threshold_column`` is supplied) is already present in the combined
    CSV produced by a prior full :func:`run_feature_extraction` call.

    The same file is rewritten in place. The per-sample track CSVs under
    ``trackdata/<sample>/<cell_type>/`` are left untouched; downstream
    analysis reads the combined CSV in
    ``analysis/<cell_type>/track_features/`` only.

    Parameters
    ----------
    output_dir : str or pathlib.Path
        Base output directory used by BEHAV3D.
    cell_type : str
        Cell type name (must match the directory used by feature
        extraction).
    new_threshold : float
        New value to threshold ``percentage_dead_mask`` (or
        ``threshold_column``) against.
    threshold_column : str, optional
        Column in the combined CSV used by :func:`calculate_death`.
        Defaults to ``"percentage_dead_mask"``.

    Returns
    -------
    pathlib.Path
        Path to the combined CSV that was rewritten.

    Raises
    ------
    FileNotFoundError
        If the combined CSV does not exist.
    KeyError
        If the required ``threshold_column`` is missing from the CSV.
    """
    output_dir = Path(output_dir)
    combined_csv = (
        output_dir
        / "analysis"
        / cell_type
        / "track_features"
        / f"BEHAV3D_{cell_type}_combined_track_features.csv"
    )
    if not combined_csv.exists():
        raise FileNotFoundError(
            f"Combined feature CSV not found for cell type '{cell_type}': "
            f"{combined_csv}. Run a full feature extraction first."
        )

    print(
        f"--------------- Re-running death classification for {cell_type} "
        f"with threshold={new_threshold} ---------------"
    )
    df = pd.read_csv(combined_csv)
    if threshold_column not in df.columns:
        raise KeyError(
            f"Column '{threshold_column}' not present in {combined_csv}. "
            "Cannot recompute the 'dead' classification without it. "
            "Re-run full feature extraction (with 'death' selected) first."
        )

    if "dead" in df.columns:
        df = df.drop(columns=["dead"])

    df = calculate_death(
        df,
        threshold=float(new_threshold),
        threshold_column=threshold_column,
    )
    df.to_csv(combined_csv, index=False)
    print(f"Re-wrote {combined_csv} with updated 'dead' column.")
    return combined_csv


def run_feature_extraction(
    metadata, 
    config=None,
    output_dir=None,
    cell_type="tcell",
    features_choice=["movement", "intensity", "morphology", "contact", "death"],
    imaris=False,
    dead_mask_percentage_threshold=None,
    contact_threshold=None,
    rolling_meanspeed_window=10,
    overwrite=False,
    n_workers=1,
    progress_cb=None,
    ):
    assert config is not None or output_dir is not None, "Either 'config' or 'output_dir' must be supplied"
    
    if output_dir is None:
        output_dir = config['output_dir']

    df_all_tracks=pd.DataFrame()

    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")

    _total_samples = len(metadata)
    for _i, (_, sample_metadata) in enumerate(metadata.iterrows()):
        if progress_cb is not None:
            try:
                progress_cb(_i, _total_samples, f"{cell_type} / {sample_metadata['sample_name']}")
            except Exception:
                pass

        print(f"--------------- Processing {cell_type}: {sample_metadata['sample_name']} ---------------")
        if is_multicolor_celltype(cell_type):
            raise ValueError(
                f"Feature extraction does not run on per-channel multicolor type '{cell_type}'. "
                "Please run it on the merged/grouped output (for example '*_merged' or '*_grouped')."
            )

        start_time = time.time()

        sample_name = sample_metadata['sample_name']

        try:
            distance_unit=sample_metadata['distance_unit']
            # dead_dye_threshold=sample_metadata['dead_dye_threshold']
            
            # Sometimes excel saves the encoding for µm differently, the following lines converts
            # other variants of µm to ones comparable in this code
            # The two written "μm" have different formatting
            if distance_unit=='_m':
                distance_unit ='μm'
            if distance_unit=='µm':
                distance_unit="μm"
            
            element_size_x=sample_metadata['pixel_distance_xy']
            element_size_y=sample_metadata['pixel_distance_xy'] 
            element_size_z=sample_metadata['pixel_distance_z']
            
            #Convert elekment size to um and hours, default settings for behav3d
            element_size_x = convert_distance(element_size_x, distance_unit)
            element_size_y = convert_distance(element_size_y, distance_unit)
            element_size_z = convert_distance(element_size_z, distance_unit)
            
            time_interval = sample_metadata['time_interval']
            time_unit = sample_metadata['time_unit']
            
            # For Imaris: Collect all cell-type-specific contact thresholds from metadata
            contact_thresholds = {}
            if imaris:
                for col in sample_metadata.index:
                    if col.endswith('_contact_threshold'):
                        # Extract cell type from column name (e.g., 'tcell_contact_threshold' -> 'tcell')
                        cell_type_name = col.replace('_contact_threshold', '')
                        if pd.notna(sample_metadata[col]):
                            contact_thresholds[cell_type_name] = sample_metadata[col]

            dead_channel=sample_metadata.get('dead_channel', None)
            
            print("###### Running track feature calculation")
            img_outdir = Path(output_dir, "images", sample_name)
            if not img_outdir.exists():
                img_outdir.mkdir(parents=True)

            track_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
            track_intermediate_outdir= Path(track_outdir, "intermediate_results")
            
            raw_image_path = sample_metadata['raw_image_path']
            
            # Dynamically find the current cell type's segments path
            current_cell_segments_path = None
            for prefix in ['or', 'im', 'ot']:
                col_name = f"{prefix}_{cell_type}_tracks_image_path"
                if col_name in sample_metadata.index and pd.notna(sample_metadata[col_name]):
                    current_cell_segments_path = sample_metadata[col_name]
                    break
            
            if current_cell_segments_path is None:
                print(f"⚠️  {cell_type} channel is empty for sample {sample_name} (no segments/tracks path in metadata) – skipping.")
                continue
            
            if not Path(current_cell_segments_path).exists():
                print(f"⚠️  {cell_type} channel is empty for sample {sample_name} (file not found: {current_cell_segments_path}) – skipping.")
                continue
            
            # Dynamically collect ALL organoid types' paths (for contact calculation)
            # Only include paths that actually exist on disk.
            organoid_segments_paths = {}
            for col in sample_metadata.index:
                if col.startswith('or_') and col.endswith('_tracks_image_path'):
                    parts = col.split('_')
                    if len(parts) >= 4:
                        organoid_type = '_'.join(parts[1:-3])
                        if is_multicolor_celltype(organoid_type):
                            continue
                        if pd.notna(sample_metadata[col]) and Path(sample_metadata[col]).exists():
                            organoid_segments_paths[organoid_type] = sample_metadata[col]
            
            # Dynamically collect ALL immune cell types' paths (for contact
            # calculation and dead-mask cleaning). Shares its detection logic
            # with the death-threshold preview via
            # resolve_immune_track_paths_for_sample so the two can never
            # disagree about which types count as immune.
            immune_segments_paths = {
                immune_type: path_val
                for immune_type, path_val in resolve_immune_track_paths_for_sample(sample_metadata).items()
                if Path(path_val).exists()
            }
            
            # Dynamically collect ALL other cell types' paths (for contact calculation)
            other_segments_paths = {}
            for col in sample_metadata.index:
                if col.startswith('ot_') and col.endswith('_tracks_image_path'):
                    parts = col.split('_')
                    if len(parts) >= 4:
                        other_type = '_'.join(parts[1:-3])
                        if is_multicolor_celltype(other_type):
                            continue
                        if pd.notna(sample_metadata[col]) and Path(sample_metadata[col]).exists():
                            other_segments_paths[other_type] = sample_metadata[col]

            # Old: Construct dead_mask_path from output directory (kept for reference)
            # dead_mask_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
            
            # Dead mask path is metadata-driven only (no hidden fallbacks).
            dead_mask_path = None
            if 'dead_mask_path' in sample_metadata.index and pd.notna(sample_metadata['dead_mask_path']):
                dead_mask_path = Path(sample_metadata['dead_mask_path'])

            needs_dead_mask = ("death" in features_choice) and (dead_channel is not None and pd.notna(dead_channel))
            if needs_dead_mask:
                if dead_mask_path is None:
                    print(
                        f"⚠️  {cell_type} in sample {sample_name} requires dead_mask_path in metadata "
                        "(death feature requested with dead_channel), but it is missing – skipping."
                    )
                    continue
                if not dead_mask_path.exists():
                    print(
                        f"⚠️  {cell_type} in sample {sample_name} dead_mask_path file not found: "
                        f"{dead_mask_path} – skipping."
                    )
                    continue
            
            print(f"{get_current_time()} - Converting all input files to .zarr for memory efficiency...")
            axis_order = None
            if "dimension_order" in sample_metadata.index and pd.notna(sample_metadata.get("dimension_order")):
                val = sample_metadata["dimension_order"]
                if isinstance(val, str) and 2 <= len(val) <= 5:
                    axis_order = val.upper()

            seg_path = Path(current_cell_segments_path)
            if seg_path.suffix != ".zarr" and not str(seg_path).endswith(".zarr.zip"):
                print(
                    f"⚠️  {cell_type} in sample {sample_name} uses a tracked image that is not .zarr: {seg_path}. "
                    "Please run Tracking or Import external tracking first."
                )
                continue
            current_cell_segments_path = seg_path

            raw_path = Path(raw_image_path)
            raw_zarr = raw_path if raw_path.suffix == ".zarr" else Path(img_outdir, f"{sample_name}.zarr")
            convert_raw_file_to_zarr(
                path=raw_path,
                outpath=raw_zarr,
                axis_order=axis_order,
                overwrite=overwrite,
            )
            raw_image_path = raw_zarr

            if not track_outdir.exists():
                track_outdir.mkdir(parents=True)
            if not track_intermediate_outdir.exists():
                track_intermediate_outdir.mkdir()
            if not analysis_outdir.exists():
                analysis_outdir.mkdir(parents=True)
            if not feature_outdir.exists():
                feature_outdir.mkdir(parents=True)
            
            print(f"{get_current_time()} - Loading in tracks csv...")
            # Find the correct prefixed column (or_, im_, ot_)
            tracks_csv_col = None
            for prefix in ['or', 'im', 'ot']:
                col_name = f"{prefix}_{cell_type}_tracks_csv_path"
                if col_name in sample_metadata.index and pd.notna(sample_metadata[col_name]):
                    tracks_csv_col = col_name
                    break
            
            if tracks_csv_col is None:
                # Fallback to old non-prefixed format for backward compatibility
                tracks_csv_col = f"{cell_type}_tracks_csv_path"
            
            df_tracks_path = sample_metadata.get(tracks_csv_col)
            if df_tracks_path is None or (isinstance(df_tracks_path, float) and pd.isna(df_tracks_path)):
                print(f"⚠️  {cell_type} channel is empty for sample {sample_name} (no tracks CSV path) – skipping.")
                continue
            
            df_tracks_path = Path(df_tracks_path)
            if not df_tracks_path.exists():
                print(f"⚠️  Tracks CSV not found for {cell_type} in sample {sample_name}: {df_tracks_path} – skipping.")
                continue
            
            df_tracks=pd.read_csv(df_tracks_path, sep=",")
            # Adding a sample name for later combination of multiple track experiments
            df_tracks['sample_name']=sample_name

            if df_tracks.empty or df_tracks['TrackID'].nunique() == 0:
                print(f"⚠️  No tracks found for {cell_type} in sample {sample_name} – skipping.")
                continue
            
            print(f"{get_current_time()} - Generalizing the units of position and time to um and hours")
            df_tracks = generalize_units_of_track_features(
                df_tracks=df_tracks,
                distance_unit=distance_unit,
                time_interval=time_interval,
                time_unit=time_unit
                )
            time_interval=convert_time(time_interval, time_unit)
            time_unit = "h"
            
            required_features = {"morphology", "intensity", "contact", "death"}

            if any(feature in features_choice for feature in required_features):
                ### Calculate image based features (per timepoint)
                print(f"{get_current_time()} - Calculating single-timepoint image-based features...")
                df_intensity_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_intensity.csv")
                df_top_quantile_intensity_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_top_quantile_intensity.csv")
                df_contacts_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_contact.csv")
                df_dead_mask_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_dead_mask.csv")
                df_morphology_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_morphology.csv")
                
                if imaris:
                    df_tracks=calculate_imaris_track_features(
                        df_tracks=df_tracks,
                        cell_type=cell_type,
                        contact_threshold=contact_threshold,
                        contact_thresholds=contact_thresholds,
                        distance_unit=distance_unit,
                    )
                else:
                    df_tracks=calculate_image_based_track_features(
                        df_tracks=df_tracks,
                        cell_type=cell_type,
                        features_choice=features_choice,
                        df_dead_mask_outpath=df_dead_mask_outpath,
                        df_morphology_outpath=df_morphology_outpath,
                        df_intensity_outpath=df_intensity_outpath,
                        df_top_quantile_intensity_outpath=df_top_quantile_intensity_outpath,
                        df_contacts_outpath=df_contacts_outpath,
                        dead_channel=dead_channel,
                        contact_threshold=contact_threshold,
                        element_size_x=element_size_x,
                        element_size_y=element_size_y,
                        element_size_z=element_size_z,
                        raw_image_path=raw_image_path,
                        dead_mask_path=dead_mask_path,
                        current_cell_segments_path=current_cell_segments_path,
                        organoid_segments_paths=organoid_segments_paths,
                        immune_segments_paths=immune_segments_paths,
                        other_segments_paths=other_segments_paths,
                        overwrite=overwrite,
                        n_workers=n_workers
                    )
            else:
                print(f"{get_current_time()} - Skipping image-based features as none requested in features_choice (morphology, intensity, contact, death)")

            if "intensity" in features_choice:
                intensity_cols = [c for c in df_tracks.columns if c.startswith("mean_intensity_ch")]
                if intensity_cols:
                    print(f"{get_current_time()} - Calculating fold change over baseline for intensity features...")
                    df_tracks = calculate_fold_change_over_baseline(df_tracks, intensity_cols)

                top_quantile_cols = [
                    c for c in df_tracks.columns
                    if c.startswith("q") and "_mean_intensity_ch" in c
                    and c.split("_mean_intensity_ch", 1)[0][1:].isdigit()
                ]
                if top_quantile_cols:
                    print(f"{get_current_time()} - Calculating fold change over track minimum for top-quantile (hot-pixel) intensity features...")
                    df_tracks = calculate_fold_change_over_track_min(df_tracks, top_quantile_cols)

                    print(f"{get_current_time()} - Calculating previous-frame fold change / flicker score for top-quantile intensity features...")
                    df_tracks = calculate_fold_change_over_previous_frame(df_tracks, top_quantile_cols)

            print(f"{get_current_time()} - Interpolating missing timepoints based on time interval")
            # As sometimes 1 or several timepoints are missing in a track, interpolate these missing rows
            # Values are interpolated linearly, forward filled or left blank based on the column
            # More explanation within the function
            df_tracks= interpolate_missing_positions(df_tracks)

            if "movement" in features_choice:
                print(f"{get_current_time()} - Calculating movement features...")
                df_tracks=calculate_movement_features(
                    df_tracks,
                    time_interval = time_interval,
                    rolling_meanspeed_window=rolling_meanspeed_window
                    )
                df_tracks = df_tracks.sort_values(['TrackID', 'position_t'])
            else:
                print(f"{get_current_time()} - Skipping movement features as not requested in features_choice")
            
            if "death" in features_choice:
                print("Smoothing death features with a rolling window of 10 minutes...")
                def get_death_smoothing_window_size(time_interval, time_unit, smoothing_minutes=10):
                    interval_minutes = convert_time(time_interval, time_unit, convert_to="m")
                    if interval_minutes <= 0 or interval_minutes > smoothing_minutes:
                        return None
                    return max(1, int(round(smoothing_minutes / interval_minutes)))
                death_smoothing_window_size = get_death_smoothing_window_size(time_interval, time_unit)
                df_tracks = calculate_smoothed_death_columns(
                    df_tracks,
                    window_size=death_smoothing_window_size,
                )
                
            if "death" in features_choice:
                # Calculate death features if threshold specified and dead_channel exists
                if dead_mask_percentage_threshold is not None and dead_channel is not None and pd.notna(dead_channel):
                    print(f"{get_current_time()} - Calculating cell death based on dead_mask_percentage_threshold {dead_mask_percentage_threshold}")
                    df_tracks = calculate_death(df_tracks, threshold=dead_mask_percentage_threshold, threshold_column="percentage_dead_mask")
                
            # if "contact" in features_choice:    
            #     # Calculate active contact for ALL cell types that have touching columns
            #     for col in df_tracks.columns:
            #         if col.startswith('touching_') and col.endswith('s'):
            #             # Extract target cell type from column name
            #             target_type = col[len('touching_'):-1]
            #             contact_col = f'{target_type}_contact'
            #             if contact_col in df_tracks.columns:
            #                 print(f"{get_current_time()} - Calculating active contact features for {cell_type} vs {target_type}...")
            #                 df_tracks = calculate_active_contact_features(
            #                     df_tracks,
            #                     cell_type=target_type
            #                 )
            #             else:
            #                 print(f"{get_current_time()} - Skipping active contact for {target_type}: no {contact_col} column found")
           
                
            tracks_out_path = Path(track_outdir, f"{sample_name}_{cell_type}_track_features.csv")
            print(f"{get_current_time()} - Writing output to {tracks_out_path}")
            df_tracks.to_csv(tracks_out_path, sep=",", index=False)
            
            df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
            df_all_tracks = pd.concat([df_all_tracks, df_tracks])
            
            end_time = time.time()
            h,m,s = format_time(start_time, end_time)
            print(f"###### DONE - elapsed time: {h}:{m:02}:{s:02}\n")

        except Exception:
            traceback.print_exc()
            print(f"⚠️  Error processing {cell_type} for sample {sample_name} – skipping to next sample.\n")
            continue
        
    feature_outdir.mkdir(parents=True, exist_ok=True)
    all_tracks_out_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_all_tracks.to_csv(all_tracks_out_path, index=False)

    stale_steps = []
    filtered_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_filtered.csv")
    if filtered_path.exists():
        stale_steps.append("Filtering")
    active_killing_dir = Path(analysis_outdir, "active_killing")
    if active_killing_dir.exists() and any(active_killing_dir.rglob(f"BEHAV3D_{cell_type}_advanced_track_features.csv")):
        stale_steps.append("Active Killing")
    if stale_steps:
        print(
            f"{get_current_time()} - WARNING: Feature Extraction for '{cell_type}' was just rerun. "
            f"Existing {' and '.join(stale_steps)} results for this cell type were computed from the "
            f"previous track features and are now stale. Re-run {' and then '.join(stale_steps)} to refresh them."
        )

    if progress_cb is not None:
        try:
            progress_cb(_total_samples, _total_samples, f"{cell_type} done")
        except Exception:
            pass
    return(df_all_tracks)       

def calculate_image_based_track_features(
    df_tracks,
    cell_type,
    dead_channel,
    
    features_choice,
    
    contact_threshold,
    element_size_x,
    element_size_y,
    element_size_z,
    #Paths
   
    dead_mask_path,
    current_cell_segments_path,
    organoid_segments_paths={},
    immune_segments_paths={},
    other_segments_paths={},
    raw_image_path = "",
    
    df_dead_mask_outpath="",
    df_morphology_outpath="",
    df_intensity_outpath="",
    df_top_quantile_intensity_outpath="",
    df_contacts_outpath="",
    # Overwrite/redo df_intensity_outpath and df_contacts_outpath
    dead_mask_percentage_threshold=None,
    n_workers=1,
    overwrite=False,
    ):
    """
    This code calculates the various features for each timepoint in a track for each 
    separate experiment.
    
    Fully flexible - works with any cell type from metadata.
    
    Output:
    - A .csv file containing all timepoints of all TrackIDs and their time-related 
      features
    """
    segments_path = current_cell_segments_path
    has_dead_mask = dead_mask_path is not None
    merge_keys = ["TrackID", "position_t"]
    df_tracks = _normalize_merge_keys(df_tracks, keys=merge_keys)
    
    if "morphology" in features_choice:
        print(f"{get_current_time()} - Calculating morphology features...")
        morph_dtypes = {
            "TrackID": "Int64",
            "position_t": "Int64",
            "nr_pixels": "Int64",
            "volume": float,
            "bbox_volume": float,
            "elongation": float,
            "extent": float,
            "equivalent_diameter": float,
            "major_axis_length": float,
            "minor_axis_length": float,
            "axis1_length": float,
            "axis2_length": float,
            "axis3_length": float,
            "oblateness": float,
            "prolateness": float,
            "surface_area": float,
            "sphericity": float,
            "convex_volume": float,
            "surface_to_volume_ratio": float,
            "orientation_vector": object,
            "border_touching_segment": "boolean",
        }
        
        if df_morphology_outpath.exists() and not overwrite:
            print("Morphology calculation .csv already exists. Loading in morphology information...")
            df_morphology = pd.read_csv(df_morphology_outpath)
            if "orientation_vector" in df_morphology.columns:
                df_morphology["orientation_vector"] = df_morphology["orientation_vector"].apply(
                    lambda v: ast.literal_eval(v) if isinstance(v, str) and v.startswith("[") else v
                )
        else:
            df_morphology=calculate_morphology_features(
                segments_path=segments_path,
                n_workers=n_workers,
                voxel_spacing=(element_size_z, element_size_y, element_size_x)
            )
            if df_morphology_outpath != "":
                df_morphology.to_csv(df_morphology_outpath, sep=",", index=False)  

        for col, dtype in morph_dtypes.items():
            if col in df_morphology.columns:
                try:
                    df_morphology[col] = df_morphology[col].astype(dtype)
                except Exception:
                    pass

        df_morphology = _normalize_merge_keys(df_morphology, keys=merge_keys)
        df_tracks = pd.merge(df_tracks, df_morphology, how="left", on=merge_keys)
    
    else:
        print(f"{get_current_time()} - Skipping full morphology features, only computing fast nr_pixels & volume")
        df_basic = calculate_basic_morphology(
            segments_path=segments_path,
            voxel_spacing=(element_size_z, element_size_y, element_size_x),
            n_workers=n_workers,
        )

        basic_dtypes = {
            "TrackID": "Int64",
            "position_t": "Int64",
            "nr_pixels": "Int64",
            "volume": float,
        }
        df_basic = df_basic.astype(basic_dtypes, copy=False)
        df_basic = _normalize_merge_keys(df_basic, keys=merge_keys)
        df_tracks = pd.merge(df_tracks, df_basic, how="left", on=merge_keys)

    if "intensity" in features_choice:
        print(f"{get_current_time()} - Calculating channel and especially death dye intensities...")

        if df_intensity_outpath.exists() and not overwrite:
            print("Intensity calculation .csv already exists. Loading in intensity information...")
            df_intensity = pd.read_csv(df_intensity_outpath)
        else:
            df_intensity=calculate_segment_intensity(
                segments=segments_path,
                intensity_image=raw_image_path,
                n_workers=n_workers,
            )
            if dead_channel is not None and pd.notna(dead_channel):
                df_intensity = df_intensity.rename(columns={f"mean_intensity_ch{dead_channel}":"mean_dead_dye"})
            
            if df_intensity_outpath != "":
                df_intensity.to_csv(df_intensity_outpath, sep=",", index=False)
        df_intensity = _normalize_merge_keys(df_intensity, keys=merge_keys)
        df_tracks = pd.merge(df_tracks, df_intensity, how="left", on=merge_keys)

        print(f"{get_current_time()} - Calculating top-quartile/decile (hot-pixel) intensities...")
        df_top_quantile_intensity = None
        if df_top_quantile_intensity_outpath.exists() and not overwrite:
            df_cached = pd.read_csv(df_top_quantile_intensity_outpath)
            if any(c.startswith("q90_mean_intensity_ch") for c in df_cached.columns):
                print("Top-quartile/decile intensity calculation .csv already exists. Loading in top-quantile intensity information...")
                df_top_quantile_intensity = df_cached
            else:
                print("Cached top-quantile intensity .csv predates the top-decile (q90) feature - recalculating...")
        if df_top_quantile_intensity is None:
            df_top_quantile_intensity = calculate_segment_top_quantile_intensity(
                segments=segments_path,
                intensity_image=raw_image_path,
                n_workers=n_workers,
            )
            if df_top_quantile_intensity_outpath != "":
                df_top_quantile_intensity.to_csv(df_top_quantile_intensity_outpath, sep=",", index=False)
        df_top_quantile_intensity = _normalize_merge_keys(df_top_quantile_intensity, keys=merge_keys)
        df_tracks = pd.merge(df_tracks, df_top_quantile_intensity, how="left", on=merge_keys)
    else:
        print(f"{get_current_time()} - Skipping intensity features as not requested in features_choice")
    
    if "death" in features_choice:
        if dead_channel is not None and pd.notna(dead_channel) and has_dead_mask:
            print(f"{get_current_time()} - Calculating number of dead mask pixels...")
            if df_dead_mask_outpath.exists() and not overwrite:
                print("Dead mask calculation .csv already exists. Loading in dead mask calculation information...")
                df_dead_mask = pd.read_csv(df_dead_mask_outpath)
            else:
                # Clean the dead mask by zeroing voxels that fall inside any
                # immune segment, but only when the current cell type is NOT
                # itself an immune type. This prevents immune cells trespassing
                # through organoids (or other passive structures) from inflating
                # the target's `percentage_dead_mask`.
                exclude_immune_from_dead = (
                    cell_type not in immune_segments_paths
                    and bool(immune_segments_paths)
                )
                if exclude_immune_from_dead:
                    print(
                        f"{get_current_time()} - Cleaning dead mask: zeroing voxels "
                        f"inside immune segments {list(immune_segments_paths.keys())} "
                        f"before computing {cell_type} death features."
                    )

                df_dead_mask=calculate_dead_mask(
                    segments=segments_path,
                    dead_mask=dead_mask_path,
                    n_workers=n_workers,
                    immune_segments_paths=immune_segments_paths if exclude_immune_from_dead else None,
                )
                if df_dead_mask_outpath != "":
                    df_dead_mask.to_csv(df_dead_mask_outpath, sep=",", index=False)
            df_dead_mask = _normalize_merge_keys(df_dead_mask, keys=merge_keys)
            df_tracks = pd.merge(df_tracks, df_dead_mask, how="left", on=merge_keys)
        else:
            if not has_dead_mask and dead_channel is not None:
                print(f"{get_current_time()} - Skipping death mask calculations: dead_mask_path not found in metadata or file doesn't exist")
            else:
                print(f"{get_current_time()} - Skipping death mask calculations: no dead_channel specified in metadata")
    else:
        print(f"{get_current_time()} - Skipping dead mask calculations as not requested in features_choice")
        
    if "contact" in features_choice:
        print(f"{get_current_time()} - Calculating contacts between {cell_type} and all other cell types...")
        print(f"Using a contact threshold of {contact_threshold} um")
        # Calculate contacts dynamically for all cell types found in metadata
        # The resulting DataFrame will have columns for ALL cell types (e.g., organoid1_contact, macro_contact, etc.)
        
        # Determine if invasiveness should be calculated (immune cells + invasiveness feature requested)
        should_calculate_invasiveness = "invasiveness" in features_choice and cell_type in immune_segments_paths.keys()
        
        cache_valid = df_contacts_outpath.exists() and not overwrite
        if cache_valid and should_calculate_invasiveness:
            cached_cols = pd.read_csv(df_contacts_outpath, nrows=0).columns
            if not any(c.endswith("_invasiveness_perc") for c in cached_cols):
                print(
                    "Contact .csv already exists but predates invasiveness "
                    "calculation. Recomputing contact information..."
                )
                cache_valid = False

        if cache_valid:
            print("Contact .csv already exists. Loading in contact information...")
            df_contacts = pd.read_csv(df_contacts_outpath)

        else:
            df_contacts=calculate_contact_features(
                current_cell_segments_path=current_cell_segments_path,
                organoid_segments_paths=organoid_segments_paths,
                immune_segments_paths=immune_segments_paths,
                other_segments_paths=other_segments_paths,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                contact_threshold=contact_threshold,
                calculate_from=cell_type,
                n_workers=n_workers,
                calculate_invasiveness=should_calculate_invasiveness
            )
            if df_contacts_outpath != "":
                df_contacts.to_csv(df_contacts_outpath, sep=",", index=False)
        df_contacts = _normalize_merge_keys(df_contacts, keys=merge_keys)
        df_tracks = pd.merge(df_tracks, df_contacts, how="left", on=merge_keys)
        if "invasiveness" in features_choice:
            perc_cols = [
                c for c in df_tracks.columns
                if c.endswith("_invasiveness_perc") and c != "any_org_invasiveness_perc"
            ]
            if perc_cols:
                bases = [c[: -len("_invasiveness_perc")] for c in perc_cols]
                for base, perc_col in zip(bases, perc_cols):
                    bool_col = f"{base}_invasiveness"
                    if bool_col not in df_tracks.columns:
                        df_tracks[bool_col] = np.nan
                    mask = df_tracks[perc_col].notna()
                    df_tracks.loc[mask, bool_col] = (
                        df_tracks.loc[mask, bool_col]
                        .fillna(df_tracks.loc[mask, perc_col] >= 50.0)
                        .astype(bool)
                    )
                if "any_org_invasiveness_perc" not in df_tracks.columns:
                    df_tracks["any_org_invasiveness_perc"] = df_tracks[perc_cols].max(axis=1, skipna=True)
                bool_cols = [
                    f"{base}_invasiveness"
                    for base in bases
                    if f"{base}_invasiveness" in df_tracks.columns
                ]
                if bool_cols:
                    if "any_org_invasiveness" not in df_tracks.columns:
                        df_tracks["any_org_invasiveness"] = np.nan
                    any_mask = df_tracks["any_org_invasiveness"].isna()
                    df_tracks.loc[any_mask, "any_org_invasiveness"] = (
                        df_tracks.loc[any_mask, bool_cols].fillna(False).any(axis=1)
                    )
    else:
        print(f"{get_current_time()} - Skipping contact calculations as not requested in features_choice")
        
    df_tracks = df_tracks.sort_values(
        by=["sample_name", "TrackID", "position_t"], 
        ascending=[True, True, True]
    )    
    new_order = ["sample_name"] + [c for c in df_tracks.columns.tolist() if c != "sample_name"]
    df_tracks = df_tracks[new_order]
    
    return(df_tracks)

def generalize_units_of_track_features(
    df_tracks,
    distance_unit,
    time_interval,
    time_unit
    ):
    # Converting the time and distance values to a default unit to allow comparison 
    # with differently provided units (defaults to µm and hours)
    print(f"{get_current_time()} - Converting distance and time unit to default um and hours...")   
    df_tracks["position_z"]=df_tracks["position_z"].apply(convert_distance, args=(distance_unit,))
    df_tracks["position_y"]=df_tracks["position_y"].apply(convert_distance, args=(distance_unit,))
    df_tracks["position_x"]=df_tracks["position_x"].apply(convert_distance, args=(distance_unit,))
    
    # Calculate relative time, where each track begins at timepoint 1
    # Use transform so the TrackID column is preserved (groupby+apply+reset_index(drop=True)
    # silently drops the TrackID column in newer pandas versions)
    min_time_per_track = df_tracks.groupby('TrackID')['position_t'].transform('min')
    df_tracks['relative_time'] = df_tracks['position_t'] - min_time_per_track + 1

    df_tracks["time"]=df_tracks["position_t"]*time_interval
    df_tracks["time"]=df_tracks["time"].apply(convert_time, args=(time_unit))
    
    # Make default distance unit um over μm, as encoding of μm between python 
    # and e.g Excel can lead to formatting errors
    df_tracks["distance_unit"] = "um"
    df_tracks["time_unit"] = "h"
    return(df_tracks)
    
def calculate_active_contact_features(df_tracks, cell_type):
    """
    Calculate active contact features for any cell type.
    
    Determines which cells are "actively" interacting vs "passively" contacted.
    When multiple cells of the same type are in contact, the one with the highest
    mean_speed is considered to be actively engaging.
    
    Args:
        df_tracks: DataFrame with track features including mean_speed and touching_{cell_type}s columns
        cell_type: The cell type to calculate active contact
    Returns:
        DataFrame with added active_{cell_type}_contact column
    """
    touching_col = f'touching_{cell_type}s'
    contact_col = f'{cell_type}_contact'
    active_contact_col = f'active_{cell_type}_contact'
    
    # Check if mean_speed column exists (requires movement features to be calculated first)
    if 'mean_speed' not in df_tracks.columns:
        print(f"{get_current_time()} - Warning: 'mean_speed' column not found. Skipping active contact calculation.")
        print(f"    Make sure 'movement' is included in features_choice to enable active contact detection.")
        # Set all active contacts to False since we can't determine activity without speed
        df_tracks[active_contact_col] = False
        return df_tracks
    
    print(f"{get_current_time()} - Determining active contact of {cell_type}")

    # --- Step 1: Explode touching cells column ---
    df_explode = (
        df_tracks[['TrackID', 'position_t', touching_col]]
        .dropna()
        .assign(**{touching_col: lambda d: d[touching_col].astype(str).str.split(',')})
        .explode(touching_col)
    )
    df_explode = df_explode[df_explode[touching_col].str.strip() != '']
    df_explode[touching_col] = pd.to_numeric(df_explode[touching_col], errors='coerce').astype('Int64')

    # --- Step 2: Get mean_speed of each touching cell ---
    speed_map = df_tracks[['TrackID', 'position_t', 'mean_speed']].rename(
        columns={'TrackID': touching_col, 'mean_speed': 'touching_speed'}
    )
    df_explode = df_explode.merge(speed_map, on=[touching_col, 'position_t'], how='left')

    # --- Step 3: Aggregate max touching speed for each cell ---
    max_touching_speed = df_explode.groupby(['TrackID', 'position_t'])['touching_speed'].max()

    # --- Step 4: Compare own speed with max of touching cells ---
    df_tracks = df_tracks.set_index(['TrackID', 'position_t'])
    df_tracks['max_touching_speed'] = max_touching_speed
    df_tracks['max_touching_speed'] = df_tracks['max_touching_speed'].fillna(-1)

    # Active contact = in contact AND moving faster than all touching cells
    df_tracks[active_contact_col] = (
        (df_tracks[contact_col]) & 
        (df_tracks['mean_speed'] >= df_tracks['max_touching_speed'])
    )

    df_tracks.reset_index(inplace=True)
    df_tracks.drop(columns=['max_touching_speed'], inplace=True)

    return df_tracks




def calculate_imaris_track_features(
    df_tracks,
    cell_type,
    contact_threshold,
    distance_unit,
    ):
    """
    Calculate contact features from Imaris-extracted statistics.
    Flexible - works with any cell type from metadata.
    
    Imaris provides pre-calculated distance columns in the CSV.
    This function thresholds those distances to determine contacts.
    
    Expected Imaris distance columns in df_tracks:
    - For organoid types: 'organoid_distance', 'organoid1_distance', 'organoid2_distance', etc.
    - For immune types: 'tcell_distance', 'macro_distance', etc.
    - For complementary contacts: 'complementary_{cell_type}_distance'
    
    Output:
    - DataFrame with contact columns for all detected cell types
    """
    
    print("Performing feature calculation from Imaris processing..")
    print(f"Using a contact threshold of {contact_threshold}{distance_unit}")
    
    # Detect all distance columns in the DataFrame
    distance_cols = [col for col in df_tracks.columns if col.endswith('_distance')]
    
    # Calculate contacts for each distance column
    for dist_col in distance_cols:
        # Extract cell type name from column (e.g., 'organoid_distance' -> 'organoid')
        other_cell_type = dist_col.replace('_distance', '').replace('complementary_', '')
        
        # Skip if this is the same cell type (will be handled by same-type contact below)
        if other_cell_type == cell_type:
            continue
        
        # Create contact column
        contact_col = f"{other_cell_type}_contact"
        df_tracks[contact_col] = df_tracks[dist_col] <= contact_threshold
        
        print(f"{get_current_time()} - Calculated contact with {other_cell_type} (From Imaris)")
    
    # Calculate same-type contacts (e.g., tcell-to-tcell, macro-to-macro)
    # This is based on centroid distances between cells of the same type
    print(f"{get_current_time()} - Calculating same-type contacts ({cell_type}-to-{cell_type})... (From Imaris)")
    
    grouped = df_tracks.groupby('position_t')
    new_dfs = []
    
    for group_name, group_df in grouped:
        positions = group_df[['position_x', 'position_y', 'position_z']].values
        distances = cdist(positions, positions)
        np.fill_diagonal(distances, np.inf)
        distances_mask = distances <= contact_threshold
        
        same_type_contacts_list = []
        for i, row in enumerate(distances_mask):
            contacts = np.where(row)[0].tolist()
            contact_ids = [group_df.reset_index().loc[idx]["TrackID"] for idx in contacts]
            same_type_contacts_list.append(contact_ids)
        
        group_df[f'touching_{cell_type}s'] = same_type_contacts_list
        group_df[f'{cell_type}_contact'] = group_df[f'touching_{cell_type}s'].apply(lambda x: len(x) > 0)
        
        # Check for complementary distance column (if Imaris calculated it)
        complementary_col = f'complementary_{cell_type}_distance'
        if complementary_col in df_tracks.columns:
            group_df[f'{cell_type}_contact'] = group_df[complementary_col] <= contact_threshold
            group_df[f'touching_{cell_type}s'].append("unknown")
        
        group_df[f'touching_{cell_type}s'] = group_df[f'touching_{cell_type}s'].apply(
            lambda x: ",".join(map(str, x)) if isinstance(x, list) and len(x) > 0 else None
        )
        
        new_dfs.append(group_df)
    
    df_tracks = pd.concat(new_dfs, ignore_index=True)
    return(df_tracks)

def calculate_death(
    df_tracks,
    threshold,
    threshold_column="mean_dead_dye"
    ):
    # print(f"- Calculating cell death based on defined dead_dye_threshold {dead_dye_threshold}")
    df_tracks["dead"] = False

    # ``threshold`` is always on the same scale as ``threshold_column``.
    # For ``percentage_dead_mask`` that means a fraction (0.0-1.0), matching
    # how that column is computed (mean of a 0/1 mask). Callers (napari's
    # ``_threshold_as_fraction`` / the notebook's fraction-scale widget)
    # already convert the percent-scale UI value to a fraction before it
    # reaches this function - do not convert again here.

    # For any cell crossing the dead_dye_threshold, set the cell to dead. Any timepoint after this timepoint are
    # Also set to dead, even if the mean dead dye intensity goes under the threshold again
    for track_id in df_tracks["TrackID"].unique():
        track_df = df_tracks[df_tracks["TrackID"] == track_id]
        track_df_reset = track_df.reset_index(drop=True)
        threshold_indices = track_df_reset.reset_index(drop=True)[track_df_reset[threshold_column] >= threshold].index
        
        if not threshold_indices.empty:
            first_threshold_index = threshold_indices.min()
            df_tracks.loc[track_df.index[first_threshold_index:], "dead"] = True
    return df_tracks
           
def interpolate_missing_positions(
    df_tracks,
    cols_to_copy=None,  # Will auto-detect contact columns + metadata
    cols_to_interpolate=None,
    col_to_none = [
        "SegmentID",
    ]
    ):
    """
    As not every track has a segment at every timepoint, interpolate the missing values of
    the missing timepoints.
    Returns df_tracks unchanged when there are no rows to process.
    
    Fully dynamic - automatically detects ALL contact columns from dataframe.
    
    It interpolates various columns in different ways:
    -   Interpolates the numerical columns of [cols_to_interpolate] such as speed using linear
        interpolation
    -   Copies the columns of [cols_to_copy] using a forward fill from the last non-interpolated
        row of each TrackID. Binary columns are always copied from the previous timepoint,
        which keeps dynamically added contact/death/border flags stable across gaps.
    -   Puts None in any column not specified, such as SegmentID, as no actual segment exists
    """

    # Nothing to interpolate – return as-is to avoid empty-concat errors downstream.
    if df_tracks.empty:
        return df_tracks

    base_copy_cols = set(cols_to_copy or [])
    base_copy_cols.update(
        {
            "sample_name",
            "TrackID",
            "distance_unit",
            "time_unit",
            "orientation_vector",
            "principal_axes",
            "well",
            "exp_nr",
        }
    )
    base_copy_cols.update(
        {
            col
            for col in df_tracks.columns
            if col.endswith("_line_condition")
        }
    )

    def _is_binary_series(series):
        if pd.api.types.is_bool_dtype(series):
            return True

        numeric = pd.to_numeric(series.dropna(), errors="coerce")
        numeric = numeric[np.isfinite(numeric)]
        if numeric.empty:
            return False

        unique_values = set(pd.unique(numeric.astype(float)))
        return unique_values.issubset({0.0, 1.0})

    def _should_copy_series(series):
        if series.name in base_copy_cols:
            return True
        if not pd.api.types.is_numeric_dtype(series):
            return True
        non_na = series.dropna()
        if non_na.nunique(dropna=True) <= 1:
            return True
        return _is_binary_series(series)

    def _copy_previous_values(series, missing_mask):
        values = series.copy()
        previous_value = pd.NA
        has_previous = False

        def _has_observed_value(value):
            if value is pd.NA or value is None:
                return False
            if isinstance(value, (list, tuple, dict, set, np.ndarray)):
                return True
            try:
                return not bool(pd.isna(value))
            except TypeError:
                return True

        for idx, value in values.items():
            if _has_observed_value(value):
                previous_value = value
                has_previous = True
            elif bool(missing_mask.loc[idx]) and has_previous:
                values.at[idx] = previous_value

        return values

    grouped_df = df_tracks.groupby('TrackID')

    def interpolate_group(group):
        group = group.sort_values("position_t").copy()
        group["interpolated"] = False

        min_time = int(group["position_t"].min())
        max_time = int(group["position_t"].max())
        all_times = pd.Index(np.arange(min_time, max_time + 1), name="position_t")

        if group["position_t"].duplicated().any():
            group = group.drop_duplicates(subset="position_t", keep="first")

        df_interpolated = group.set_index("position_t").reindex(all_times).reset_index()
        existing_rows = df_interpolated["TrackID"].notna()
        missing_rows = ~existing_rows
        df_interpolated["interpolated"] = ~existing_rows

        explicit_interpolate_cols = set(cols_to_interpolate or [])
        none_cols = {col for col in col_to_none if col in df_interpolated.columns}
        copy_cols = {
            col
            for col in group.columns
            if col not in none_cols and col != "position_t" and _should_copy_series(group[col])
        }
        copy_cols.update(
            {
                col
                for col in df_interpolated.columns
                if col in base_copy_cols and col not in none_cols and col != "position_t"
            }
        )

        if explicit_interpolate_cols:
            interpolate_cols = [
                col
                for col in explicit_interpolate_cols
                if col in df_interpolated.columns
                and col not in copy_cols
                and col not in none_cols
                and col != "position_t"
                and pd.api.types.is_numeric_dtype(group[col])
            ]
        else:
            interpolate_cols = [
                col
                for col in group.columns
                if col not in copy_cols
                and col not in none_cols
                and col != "position_t"
                and pd.api.types.is_numeric_dtype(group[col])
            ]

        for col in interpolate_cols:
            series = pd.to_numeric(group.set_index("position_t")[col], errors="coerce")
            series = series.reindex(all_times)
            df_interpolated[col] = series.interpolate(method="index", limit_direction="both").to_numpy()

        for col in copy_cols:
            filled = _copy_previous_values(df_interpolated[col], missing_rows)
            df_interpolated.loc[missing_rows, col] = filled.loc[missing_rows]

        for col in none_cols:
            df_interpolated.loc[missing_rows, col] = pd.NA

        df_interpolated["interpolated"] = df_interpolated["interpolated"].astype("boolean")
        assert len(all_times) == len(df_interpolated), (
            f"Length of expected nr of timepoints ({len(all_times)}) is not the same as "
            f"resulting timepoints ({df_interpolated})"
        )
        return df_interpolated

    df_interpolated = pd.concat([interpolate_group(group) for _, group in grouped_df])
    df_interpolated = df_interpolated.reset_index(drop=True)
    return(df_interpolated)
        
def calculate_movement_features(
    df_tracks, 
    time_interval,
    rolling_meanspeed_window=5
    ):
    """
    Calculates various movement features for each timepoint of a track
    """

    def calculate_displacement(track_coordinates):
        """calculate the displacement per timepoint compared to previous timepoint"""
        track_relative_pos = np.diff(track_coordinates,axis=0,prepend=track_coordinates[[0]])
        displacement=np.apply_along_axis(np.linalg.norm, 1, track_relative_pos)
        return (displacement)
    def calculate_displacement_from_origin(track_coordinates):
        """calculate the displacement to the first timepoint"""
        displacement_from_origin=np.apply_along_axis(np.linalg.norm, 1, track_coordinates)
        return (displacement_from_origin)
    def compute_MSD(track_coordinates):
        nr_rows=len(track_coordinates)
        msd_values=np.zeros(nr_rows)
        for i in range(nr_rows):
            squared_displacements = np.sum((track_coordinates[:i+1] - track_coordinates[i])**2, axis=1)
            msd_values[i] = np.mean(squared_displacements)
        return msd_values
    
    def calculate_directional_persistence(track_coordinates, eps=1e-12):
        """
        Calculate directional persistence per timepoint compared to previous timepoint
        """
        steps = np.diff(track_coordinates, axis=0, prepend=track_coordinates[[0]])
        step_norms = np.apply_along_axis(np.linalg.norm, 1, steps)

        N = len(track_coordinates)
        persistence = np.zeros(N, dtype=float)

        # persistence needs previous step, so start at t=2
        for t in range(2, N):
            a = steps[t-1]
            b = steps[t]
            denom = step_norms[t-1] * step_norms[t]
            if denom > eps:
                persistence[t] = np.dot(a, b) / denom
            else:
                persistence[t] = 0.0
        persistence = np.clip(persistence, -1.0, 1.0)
        return persistence
    
    ## split by unique trackID2 and process
    df_tracks_processed = []
    for track in df_tracks['TrackID'].unique():
        df_track = df_tracks[df_tracks['TrackID'] == track ].sort_values(by="position_t").reset_index(drop=True)
        df_track_pos = df_track[['position_x', 'position_y', 'position_z']] ## select the data of interest
        
        # convert to array
        track_array = df_track_pos.to_numpy()
        track_array_rel = track_array - track_array[0]
        displacement = calculate_displacement(track_array_rel)
        displacement_from_origin = calculate_displacement_from_origin(track_array_rel)
        cumulative_displacement = np.cumsum(displacement, axis = 0)
        mean_square_displacement=compute_MSD(track_array_rel)
        directional_persistence = calculate_directional_persistence(track_array_rel)

        # combine
        df_computed = pd.DataFrame({
            'displacement': displacement, 
            'cumulative_displacement': cumulative_displacement, 
            'displacement_from_origin': displacement_from_origin, 
            'mean_square_displacement':mean_square_displacement,
            'directional_persistence':directional_persistence
            })
        
        df_result= pd.concat([df_track_pos,df_computed], axis=1)
        df_result = pd.concat([df_result, df_track[["position_t", "SegmentID", "TrackID"]]], axis=1)

        df_result['speed'] = df_result["displacement"]/time_interval
        # Calculate the mean speed (default um/h) over the last {rolling_meanspeed_window} timepoints
        # df_result['mean_speed'] = df_result.groupby('TrackID')['speed'].apply(lambda x: x.iloc[1:].rolling(window=rolling_meanspeed_window, min_periods=1).mean()).reset_index(0, drop=True)
        # df_result['mean_speed'] = df_result['mean_speed'].fillna(0)
   
        df_tracks_processed.append(df_result)
    df_tracks_processed = pd.concat(df_tracks_processed)
    df_tracks_processed=pd.merge(df_tracks, df_tracks_processed, how="left")
    return(df_tracks_processed)

def calculate_cell_distance(
    current_cell_segments, 
    target_cell_segments, 
    element_size_x, 
    element_size_y, 
    element_size_z,
    current_cell_type="cell",
    target_cell_type="target"
    ):
    """
    Calculate distance from current cell type to target cell type.
    Fully flexible - works with any cell type combination.
    
    Args:
        current_cell_segments: Segments of the cell type being analyzed
        target_cell_segments: Segments of the target cell type
        element_size_x, element_size_y, element_size_z: Physical spacing in µm
        current_cell_type: Name of current cell type (for column naming)
        target_cell_type: Name of target cell type (for column naming)
    
    Returns:
        DataFrame with distance measurements for each cell at each timepoint
    """
    df_dist = []
    for t, current_stack in tqdm(enumerate(current_cell_segments), total=len(current_cell_segments)):
        target_stack = target_cell_segments[t,:,:,:]
        mask_target = np.ma.masked_where(target_stack==0, target_stack)
        dist_target = distance_transform_edt(mask_target.mask)
        real_dist_target = distance_transform_edt(
            mask_target.mask,
            sampling=[element_size_z, element_size_y, element_size_x]
        )
        properties_pix = pd.DataFrame(
            regionprops_table(
                label_image=current_stack, 
                intensity_image=dist_target, 
                properties=['label', 'intensity_min']
            )
        )
        properties_pix = properties_pix.rename(columns={"intensity_min": f"pix_distance_{target_cell_type}s"})
        properties_real = pd.DataFrame(
            regionprops_table(
                label_image=current_stack, 
                intensity_image=real_dist_target, 
                properties=['label', 'intensity_min']
            )
        )
        properties_real = properties_real.rename(columns={"intensity_min": f"real_distance_{target_cell_type}s"})
        properties = pd.merge(properties_pix, properties_real, how="left")
        properties["position_t"] = t
        df_dist.append(properties)
    df_dist = pd.concat(df_dist)
    df_dist = df_dist.rename(columns={"label": "TrackID"})
    df_dist[f"pix_{target_cell_type}_contact"] = df_dist[f"pix_distance_{target_cell_type}s"] <= 1.73
    return(df_dist)


def _calculate_contact_single_timepoint(args):
    """
    Calculate contacts between current cell type and ALL other cell types.
    Fully flexible - works with any combination of cell types.
    Optionally also calculates invasiveness (for immune cells against organoids).
    """
    (
        t,
        current_cell_segments_path,
        organoid_segments_paths,
        immune_segments_paths,
        other_segments_paths,
        element_size_x,
        element_size_y,
        element_size_z,
        contact_threshold,
        calculate_from,
        calculate_invasiveness
    ) = args

    # Validate element sizes before any division
    for name, value in [("element_size_x (pixel_distance_xy)", element_size_x),
                        ("element_size_y (pixel_distance_xy)", element_size_y),
                        ("element_size_z (pixel_distance_z)",  element_size_z)]:
        if value == 0 or value is None:
            raise ValueError(
                f"Element size '{name}' is {value!r}. "
                "Check that 'pixel_distance_xy' and 'pixel_distance_z' are set correctly "
                "in the metadata for all samples."
            )

    # Load current cell type's segments for this timepoint
    current_segments = load_image_timepoint(current_cell_segments_path, t)

    # Load all organoid types' segments
    organoid_segments_dict = {}
    for org_type, org_path in organoid_segments_paths.items():
        organoid_segments_dict[org_type] = load_image_timepoint(org_path, t)

    # Load all immune cell types' segments
    immune_segments_dict = {}
    for immune_type, immune_path in immune_segments_paths.items():
        immune_segments_dict[immune_type] = load_image_timepoint(immune_path, t)

    # Load all other cell types' segments
    other_segments_dict = {}
    for other_type, other_path in other_segments_paths.items():
        other_segments_dict[other_type] = load_image_timepoint(other_path, t)
    
    df_contacts = []
    segment_ids = np.unique(current_segments)

    # If there are no foreground segments at this timepoint, return an empty
    # DataFrame so that downstream concatenation can proceed without error.
    foreground_ids = [sid for sid in segment_ids if sid != 0]
    if not foreground_ids:
        return pd.DataFrame(columns=["TrackID", "position_t"])
    
    for segment_id in foreground_ids:
        
        stack_max_z, stack_max_y, stack_max_x = current_segments.shape
        seg_locs = np.argwhere(current_segments == segment_id)
        min_z, min_y, min_x = seg_locs.min(axis=0)
        max_z, max_y, max_x = seg_locs.max(axis=0)
        
        z_ext = 2 * math.ceil(contact_threshold / element_size_z)
        y_ext = 2 * math.ceil(contact_threshold / element_size_y)
        x_ext = 2 * math.ceil(contact_threshold / element_size_x)
        
        slicer = (
            slice(max(0, min_z - z_ext), min(stack_max_z, max_z + z_ext + 1)),
            slice(max(0, min_y - y_ext), min(stack_max_y, max_y + y_ext + 1)),
            slice(max(0, min_x - x_ext), min(stack_max_x, max_x + x_ext + 1))
        )
        
        seg_cutout = current_segments[slicer]
        
        real_distances = distance_transform_edt(
            seg_cutout != segment_id,
            sampling=[element_size_z, element_size_y, element_size_x]
        )
        pix_distances = distance_transform_edt(seg_cutout != segment_id)
        
        contact_data = {
            'TrackID': segment_id,
            'position_t': t,
        }
        any_org_contact = []
        any_org_contact_on_distance = []
        any_immune_contact = []
        any_immune_contact_on_distance = []
        
        # Calculate contacts with ALL organoid types
        for org_type, org_segments in organoid_segments_dict.items():
            org_cutout = org_segments[slicer]
            
            org_contacts = [
                str(x) for x in np.unique(org_cutout[real_distances <= contact_threshold]) if x != 0
            ]
            # If calculating from this organoid type, exclude self-contact
            if calculate_from == org_type:
                org_contacts = [x for x in org_contacts if x != str(segment_id)]
            
            real_org_contact = len(org_contacts) > 0
            
            pix_org_contacts = [
                str(x) for x in np.unique(org_cutout[pix_distances <= 1.73]) if x != 0
            ]
            if calculate_from == org_type:
                pix_org_contacts = [x for x in pix_org_contacts if x != str(segment_id)]
            pix_org_contact = len(pix_org_contacts) > 0
            
            # Add columns for this organoid type
            contact_data[f'{org_type}_contact'] = pix_org_contact
            contact_data[f'{org_type}_contact_on_distance'] = real_org_contact
            contact_data[f'touching_{org_type}s'] = ",".join(org_contacts) if real_org_contact else ""
            any_org_contact.append(pix_org_contact)
            any_org_contact_on_distance.append(real_org_contact)
        
        # Calculate contacts with ALL immune cell types
        for immune_type, immune_segments in immune_segments_dict.items():
            immune_cutout = immune_segments[slicer]
            
            immune_contacts = [
                str(x) for x in np.unique(immune_cutout[real_distances <= contact_threshold]) if x != 0
            ]
            # If calculating from this immune type, exclude self-contact
            if calculate_from == immune_type:
                immune_contacts = [x for x in immune_contacts if x != str(segment_id)]
            
            real_immune_contact = len(immune_contacts) > 0
            
            pix_immune_contacts = [
                str(x) for x in np.unique(immune_cutout[pix_distances <= 1.73]) if x != 0
            ]
            if calculate_from == immune_type:
                pix_immune_contacts = [x for x in pix_immune_contacts if x != str(segment_id)]
            pix_immune_contact = len(pix_immune_contacts) > 0
            
            # Add columns for this immune type
            contact_data[f'{immune_type}_contact'] = pix_immune_contact
            contact_data[f'{immune_type}_contact_on_distance'] = real_immune_contact
            contact_data[f'touching_{immune_type}s'] = ",".join(immune_contacts) if real_immune_contact else ""
            any_immune_contact.append(pix_immune_contact)
            any_immune_contact_on_distance.append(real_immune_contact)
        
        # Calculate contacts with ALL other cell types
        for other_type, other_segments in other_segments_dict.items():
            other_cutout = other_segments[slicer]
            
            other_contacts = [
                str(x) for x in np.unique(other_cutout[real_distances <= contact_threshold]) if x != 0
            ]
            # If calculating from this other type, exclude self-contact
            if calculate_from == other_type:
                other_contacts = [x for x in other_contacts if x != str(segment_id)]
            
            real_other_contact = len(other_contacts) > 0
            
            pix_other_contacts = [
                str(x) for x in np.unique(other_cutout[pix_distances <= 1.73]) if x != 0
            ]
            if calculate_from == other_type:
                pix_other_contacts = [x for x in pix_other_contacts if x != str(segment_id)]
            pix_other_contact = len(pix_other_contacts) > 0
            
            # Add columns for this other type
            contact_data[f'{other_type}_contact'] = pix_other_contact
            contact_data[f'{other_type}_contact_on_distance'] = real_other_contact
            contact_data[f'touching_{other_type}s'] = ",".join(other_contacts) if real_other_contact else ""

        contact_data['any_organoid_contact'] = any(any_org_contact)
        contact_data['any_organoid_contact_on_distance'] = any(any_org_contact_on_distance)
        contact_data['any_immune_cell_contact'] = any(any_immune_contact)
        contact_data['any_immune_cell_contact_on_distance'] = any(any_immune_contact_on_distance)
        
        # Calculate invasiveness for immune cells (only if requested and organoids exist)
        if calculate_invasiveness and organoid_segments_dict:
            # Define "surface" as pixels within 2 µm of cell boundary
            surface_threshold = 2.0  # µm
            surface_mask = real_distances <= surface_threshold
            total_surface_pixels = np.sum(surface_mask)
            
            if total_surface_pixels > 0:
                any_invasiveness_list = []
                invasiveness_perc_list = []
                
                # Calculate invasiveness for each organoid type
                for org_type, org_segments in organoid_segments_dict.items():
                    org_cutout = org_segments[slicer]
                    
                    # Count surface pixels in contact with this organoid type
                    org_contact_mask = (real_distances <= contact_threshold) & (org_cutout != 0)
                    contacted_surface_pixels = np.sum(org_contact_mask)
                    
                    # Calculate percentage
                    invasiveness_perc = (contacted_surface_pixels / total_surface_pixels) * 100.0
                    
                    # Boolean: invasive if >= 50% of surface is in contact
                    invasiveness_bool = invasiveness_perc >= 50.0
                    
                    contact_data[f'{org_type}_invasiveness'] = invasiveness_bool
                    contact_data[f'{org_type}_invasiveness_perc'] = invasiveness_perc
                    any_invasiveness_list.append(invasiveness_bool)
                    invasiveness_perc_list.append(invasiveness_perc)
                
                # Add aggregate: True if invasive against ANY organoid type
                contact_data['any_org_invasiveness'] = any(any_invasiveness_list)
                contact_data['any_org_invasiveness_perc'] = max(invasiveness_perc_list) if invasiveness_perc_list else 0.0
            else:
                # No surface detected, set invasiveness to False/0
                for org_type in organoid_segments_dict.keys():
                    contact_data[f'{org_type}_invasiveness'] = False
                    contact_data[f'{org_type}_invasiveness_perc'] = 0.0
                contact_data['any_org_invasiveness'] = False
                contact_data['any_org_invasiveness_perc'] = 0.0

        df_contacts.append(pd.DataFrame([contact_data]))

    # It is still possible (though unlikely) that no contact data is produced;
    # in that case, return an empty-but-valid DataFrame.
    if not df_contacts:
        return pd.DataFrame(columns=["TrackID", "position_t"])

    return pd.concat(df_contacts, ignore_index=True)

def calculate_contact_features(
    current_cell_segments_path,
    organoid_segments_paths={},
    immune_segments_paths={},
    other_segments_paths={},
    contact_threshold=None,
    element_size_x=None,
    element_size_y=None,
    element_size_z=None,
    calculate_from=None,
    n_workers=1,
    calculate_invasiveness=False
    ):
    """
    Flexible contact calculation - works with any cell types from metadata.
    
    Args:
        current_cell_segments_path: Path to the current cell type being analyzed
        organoid_segments_paths: Dict of {organoid_type: path} for all organoid types
        immune_segments_paths: Dict of {immune_type: path} for all immune types
        other_segments_paths: Dict of {other_type: path} for all other types
        calculate_from: The cell type we're calculating features for (required for self-contact exclusion)
        calculate_invasiveness: If True, calculate invasiveness features for immune cells against organoids
        
    Returns:
        DataFrame of contact annotations for the current cell type with ALL other cell types.
        If calculate_invasiveness=True, also includes *_invasiveness and *_invasiveness_perc columns.
    """
    current_segments = load_image(current_cell_segments_path)
    timepoints = current_segments.shape[0]

    args_list = [
        (
            t,
            current_cell_segments_path,
            organoid_segments_paths,
            immune_segments_paths,
            other_segments_paths,
            element_size_x,
            element_size_y,
            element_size_z,
            contact_threshold,
            calculate_from,
            calculate_invasiveness
        )
        for t in range(timepoints)
    ]

    if n_workers > 1:
        results = _run_parallel_with_fallback(_calculate_contact_single_timepoint, args_list, n_workers)
    else:
        results = [_calculate_contact_single_timepoint(args) for args in tqdm(args_list)]

    return pd.concat(results, ignore_index=True)

def _calculate_segment_intensity_single_timepoint(args):
    t, segments_path, intensity_image_path, calculation = args
    tcell_stack = load_image_timepoint(segments_path, t)
    intensity_stack = load_image_timepoint(intensity_image_path, t).transpose(1, 2, 3, 0)

    properties = pd.DataFrame(
        regionprops_table(
            label_image=tcell_stack,
            intensity_image=intensity_stack,
            properties=["label", f"intensity_{calculation}"],
        )
    )
    if properties.empty:
        return pd.DataFrame(columns=["label", "position_t"])

    properties["position_t"] = t
    return properties


def calculate_segment_intensity(segments, intensity_image, calculation="mean", n_workers=1):
    """
    Calculates the intensity of a specific marker features for each segment.
    The calculation can be the minimum, maximum, mean or median
    """
    assert calculation in ["min", "max", "mean", "median"]
    if isinstance(segments, (str, Path)) and isinstance(intensity_image, (str, Path)):
        timepoints = int(load_image(segments).shape[0])
        args_list = [(t, segments, intensity_image, calculation) for t in range(timepoints)]
        if n_workers > 1:
            results = _run_parallel_with_fallback(_calculate_segment_intensity_single_timepoint, args_list, n_workers)
        else:
            results = [_calculate_segment_intensity_single_timepoint(args) for args in tqdm(args_list)]
    else:
        intensity_image = np.transpose(intensity_image, axes=[0, 2, 3, 4, 1])
        results = []
        for t, (tcell_stack, intensity_stack) in tqdm(enumerate(zip(segments, intensity_image)), total=len(segments)):
            tcell_stack = np.asarray(tcell_stack)
            intensity_stack = np.asarray(intensity_stack)
            properties = pd.DataFrame(
                regionprops_table(
                    label_image=tcell_stack,
                    intensity_image=intensity_stack,
                    properties=["label", f"intensity_{calculation}"],
                )
            )
            if properties.empty:
                continue
            properties["position_t"] = t
            results.append(properties)

    if not results:
        return pd.DataFrame(columns=["TrackID", "position_t"])

    df_intensity = pd.concat(results, ignore_index=True)
    if df_intensity.empty:
        return pd.DataFrame(columns=["TrackID", "position_t"])

    # Define a dictionary mapping old column names to new column names using regex
    column_mapping = {}
    for i in range(df_intensity.shape[1]):
        old_col_name = f'intensity_mean-{i}'
        new_col_name = f'mean_intensity_ch{i}'
        column_mapping[old_col_name] = new_col_name
    column_mapping["label"]="TrackID"
    df_intensity=df_intensity.rename(columns=column_mapping)
    return(df_intensity)


def _calculate_segment_top_quantile_intensity_single_timepoint(args):
    t, segments_path, intensity_image_path, quantiles = args
    label_stack = np.asarray(load_image_timepoint(segments_path, t))
    intensity_stack = np.asarray(load_image_timepoint(intensity_image_path, t)).transpose(1, 2, 3, 0)

    labels = np.unique(label_stack)
    labels = labels[labels != 0]
    n_channels = intensity_stack.shape[-1]

    rows = []
    for label in labels:
        mask = label_stack == label
        row = {"label": int(label), "position_t": t}
        for c in range(n_channels):
            values = intensity_stack[..., c][mask]
            for q in quantiles:
                col_name = f"q{int(round(q * 100))}_mean_intensity_ch{c}"
                if values.size == 0:
                    row[col_name] = np.nan
                    continue
                cutoff = np.percentile(values, q * 100)
                top_values = values[values >= cutoff]
                row[col_name] = float(top_values.mean()) if top_values.size else np.nan
        rows.append(row)

    if not rows:
        return pd.DataFrame(columns=["label", "position_t"])
    return pd.DataFrame(rows)


def calculate_segment_top_quantile_intensity(segments, intensity_image, quantiles=(0.75, 0.9), n_workers=1):
    """Per-segment, per-timepoint mean intensity of the pixels at/above each
    given quantile of that segment's own pixel intensity distribution (e.g.
    quantile=0.75 -> mean of the brightest 25% of pixels, quantile=0.9 ->
    mean of the brightest 10% of pixels), for each channel.

    Unlike ``calculate_segment_intensity``'s whole-segment mean, this tracks
    a segment's brightest sub-region rather than being diluted by the rest
    of the cell body - the "hot-pixel mean" feature. Each quantile is
    computed from the same per-segment pixel values in a single pass and
    produces its own column, named ``q{int(quantile*100)}_mean_intensity_ch{N}``
    (e.g. ``q75_mean_intensity_ch0``, ``q90_mean_intensity_ch0``).
    """
    if isinstance(quantiles, (int, float)):
        quantiles = (quantiles,)
    quantiles = tuple(quantiles)
    assert len(quantiles) > 0
    for q in quantiles:
        assert 0.0 <= q < 1.0
    timepoints = int(load_image(segments).shape[0])
    args_list = [(t, segments, intensity_image, quantiles) for t in range(timepoints)]
    if n_workers > 1:
        results = _run_parallel_with_fallback(_calculate_segment_top_quantile_intensity_single_timepoint, args_list, n_workers)
    else:
        results = [_calculate_segment_top_quantile_intensity_single_timepoint(args) for args in tqdm(args_list)]

    if not results:
        return pd.DataFrame(columns=["TrackID", "position_t"])
    df_top_quantile = pd.concat(results, ignore_index=True)
    if df_top_quantile.empty:
        return pd.DataFrame(columns=["TrackID", "position_t"])
    df_top_quantile = df_top_quantile.rename(columns={"label": "TrackID"})
    return df_top_quantile


def calculate_relative_increase(df, column, nr_timepoints_back, groupby="TrackID"):
    df = df.sort_values(by=[groupby, "position_t"]).copy()
    def relative_increase(group):
        values = group[column].values
        increases = []
        for i in range(len(group)):
            if i < nr_timepoints_back:
                prev = values[0]
            else:
                prev = values[i - nr_timepoints_back]
            increases.append(values[i] - prev)
        return pd.Series(increases, index=group.index)

    return df.groupby(groupby).apply(relative_increase, include_groups=False).reset_index(drop=True).values

# def calculate_relative_increase(df, column, nr_timepoints_back, groupby):
#     """
#     Calculate the percentage increase of a column compared to the value
#     """
#     df = df.sort_values(by=[groupby, "position_t"]).copy()

#     def compute_increase(group):
#         values = group[column].values
#         increases = []
#         for i in range(len(values)):
#             if i >= nr_timepoints_back:
#                 ref = values[i - nr_timepoints_back]
#             else:
#                 ref = values[0]
#             current = values[i]
#             if ref == 0:
#                 increase = float('inf') if current > 0 else 0.0
#             else:
#                 increase = (current - ref) / ref
#             increases.append(increase)
#         return pd.Series(increases, index=group.index)

#     # Use `transform`-like behavior by resetting the index before recombining
#     result = df.groupby(groupby).apply(compute_increase)
#     if isinstance(result.index, pd.MultiIndex):
#         result.index = result.index.droplevel(0)

#     return result

def calculate_fold_change_over_baseline(df, intensity_cols, baseline_percentile=10, groupby=("sample_name", "TrackID")):
    """Per-track fold change of each intensity column over that track's own baseline.

    The baseline is a robust low percentile of the track's own intensity
    history (not its minimum, so a single noisy low dip doesn't set an
    artificially low denominator, and not its mean/std, which would collapse
    for a constant, non-flickering track and amplify ordinary noise into a
    false signal). Dividing by each track's own baseline level makes cells
    with different absolute intensities comparable, while a constant track
    stays at fold-change ~1 throughout.
    """
    df = df.sort_values(list(groupby) + ["position_t"]).copy()
    grouped = df.groupby(list(groupby), sort=False, observed=True)
    for col in intensity_cols:
        if col not in df.columns:
            continue
        baseline = grouped[col].transform(lambda s: np.nanpercentile(s, baseline_percentile))
        baseline = baseline.where(baseline > 0, np.nan)
        df[f"{col}_fold_change"] = df[col] / baseline
    return df


def calculate_fold_change_over_track_min(df, intensity_cols, groupby=("sample_name", "TrackID"), baseline_n_lowest=5):
    """Per-track fold change of each intensity column over that track's own
    low-end baseline: the median of the ``baseline_n_lowest`` lowest values
    ever observed for that track (not the literal single minimum, and not a
    percentile). Intended for hot-pixel-mean columns (e.g.
    ``q75_mean_intensity_ch{N}`` / ``q90_mean_intensity_ch{N}``, the mean of
    a segment's brightest quartile/decile of pixels) where the caller wants a
    strict low baseline without letting one anomalously dim frame (noise,
    segmentation glitch) single-handedly set the denominator for the whole
    track, as a literal min would.
    """
    df = df.sort_values(list(groupby) + ["position_t"]).copy()
    grouped = df.groupby(list(groupby), sort=False, observed=True)
    n_lowest = int(baseline_n_lowest)
    for col in intensity_cols:
        if col not in df.columns:
            continue
        baseline = grouped[col].transform(lambda s: s.nsmallest(n_lowest).median())
        baseline = baseline.where(baseline > 0, np.nan)
        df[f"{col}_fold_change"] = df[col] / baseline
    return df


def calculate_fold_change_over_previous_frame(
    df, intensity_cols, groupby=("sample_name", "TrackID"),
    baseline_n_lowest=5, floor=1.0, ceiling=3.0,
):
    """Per-frame fold change of each intensity column over its own previous
    timepoint within the same track (a frame-to-frame ratio), meant to catch
    abrupt brief jumps ("flickers") that a whole-track baseline dilutes by
    comparing every frame against the track's overall dim history instead of
    against what it looked like one frame ago.

    Guarded against an anomalously low previous frame: if the previous
    value is below the track's own low-end baseline (the median of its
    ``baseline_n_lowest`` lowest values, see ``calculate_fold_change_over_track_min``),
    that baseline is used as the denominator instead, so a single near-zero
    noisy dip can't blow up the ratio. The first timepoint of each track
    (no previous frame) also uses this baseline as its denominator.

    Produces two columns per input column:
      - ``{col}_fold_change_prevframe``: the raw, unclipped ratio.
      - ``{col}_flicker_score``: the same ratio clipped to ``[floor, ceiling]``
        (decreases collapsed to ``floor``, jumps at/above ``ceiling`` treated
        as equally flicker-worthy) - a bounded score for thresholding.
    """
    df = df.sort_values(list(groupby) + ["position_t"]).copy()
    grouped = df.groupby(list(groupby), sort=False, observed=True)
    n_lowest = int(baseline_n_lowest)
    for col in intensity_cols:
        if col not in df.columns:
            continue
        track_floor = grouped[col].transform(lambda s: s.nsmallest(n_lowest).median())
        prev = grouped[col].shift(1)
        baseline = prev.where(prev >= track_floor, track_floor)
        baseline = baseline.where(baseline > 0, np.nan)
        ratio = df[col] / baseline
        df[f"{col}_fold_change_prevframe"] = ratio
        df[f"{col}_flicker_score"] = ratio.clip(lower=floor, upper=ceiling)
    return df


def _zero_dead_mask_under_segments(dead_mask_t, immune_arrays_t):
    """Return a copy of ``dead_mask_t`` with voxels set to 0 wherever any
    immune segment at the same timepoint is non-zero.

    Parameters
    ----------
    dead_mask_t : array-like
        Binary (0/1) dead mask volume for a single timepoint, shape (Z, Y, X).
    immune_arrays_t : dict[str, array-like]
        Mapping of immune cell type name -> immune segment volume for the same
        timepoint. Shape must match ``dead_mask_t``; mismatched arrays are
        skipped with a warning.

    Returns
    -------
    np.ndarray
        Cleaned dead mask (same shape and dtype as input).
    """
    cleaned = np.asarray(dead_mask_t).copy()
    for name, im_t in immune_arrays_t.items():
        im_t = np.asarray(im_t)
        if im_t.shape != cleaned.shape:
            warnings.warn(
                f"[death-mask cleaning] Shape mismatch for immune '{name}' "
                f"({im_t.shape} vs dead mask {cleaned.shape}) - skipping."
            )
            continue
        cleaned[im_t > 0] = 0
    return cleaned
            
def _calculate_dead_mask_single_timepoint(args):
    t, segments_path, dead_mask_path, immune_segments_paths = args
    tcell_stack = load_image_timepoint(segments_path, t)
    dead_mask_stack = load_image_timepoint(dead_mask_path, t)

    if immune_segments_paths:
        immune_arrays_t = {name: load_image_timepoint(path, t) for name, path in immune_segments_paths.items()}
        dead_mask_stack = _zero_dead_mask_under_segments(dead_mask_stack, immune_arrays_t)

    properties = pd.DataFrame(
        regionprops_table(
            label_image=tcell_stack,
            intensity_image=dead_mask_stack,
            properties=["label", "num_pixels", "intensity_mean"],
        )
    )
    if properties.empty:
        return pd.DataFrame(
            columns=["TrackID", "position_t", "percentage_dead_mask", "nr_dead_mask_pixels"]
        )

    properties["position_t"] = t
    properties.rename(columns={"intensity_mean": "percentage_dead_mask"}, inplace=True)
    properties["nr_dead_mask_pixels"] = properties["num_pixels"] * properties["percentage_dead_mask"]
    properties = properties.rename(columns={"label": "TrackID"})
    return properties[["TrackID", "position_t", "percentage_dead_mask", "nr_dead_mask_pixels"]]


def calculate_dead_mask(segments, dead_mask, n_workers=1, immune_segments_paths=None):
    """
    Calculates the intensity of a specific marker features for each segment.
    The calculation can be the minimum, maximum, mean or median

    Parameters
    ----------
    segments : array-like
        Label image (T, Z, Y, X) for the cell type whose death features are
        being computed.
    dead_mask : array-like
        Binary dead mask (T, Z, Y, X). Each voxel is 1 where the death
        classifier marked the pixel as positive.
    immune_segments_paths : dict[str, str | Path] | None
        Optional mapping of immune cell type name -> path to immune segment
        array (T, Z, Y, X). When provided, the per-timepoint dead mask is
        cleaned via :func:`_zero_dead_mask_under_segments` before regionprops,
        so that immune cells trespassing through ``segments`` do not contribute
        to the ``percentage_dead_mask`` of the target segment.
    """
    if isinstance(segments, (str, Path)) and isinstance(dead_mask, (str, Path)):
        timepoints = int(load_image(segments).shape[0])
        args_list = [(t, segments, dead_mask, immune_segments_paths) for t in range(timepoints)]
        if n_workers > 1:
            results = _run_parallel_with_fallback(_calculate_dead_mask_single_timepoint, args_list, n_workers)
        else:
            results = [_calculate_dead_mask_single_timepoint(args) for args in tqdm(args_list)]
    else:
        results = []
        for t, (tcell_stack, dead_mask_stack) in tqdm(enumerate(zip(segments, dead_mask)), total=len(segments)):
            tcell_stack = np.asarray(tcell_stack)
            dead_mask_stack = np.asarray(dead_mask_stack)
            properties = pd.DataFrame(
                regionprops_table(
                    label_image=tcell_stack,
                    intensity_image=dead_mask_stack,
                    properties=["label", "num_pixels", "intensity_mean"],
                )
            )
            if properties.empty:
                continue
            properties["position_t"] = t
            properties.rename(columns={"intensity_mean": "percentage_dead_mask"}, inplace=True)
            properties["nr_dead_mask_pixels"] = properties["num_pixels"] * properties["percentage_dead_mask"]
            properties = properties.rename(columns={"label": "TrackID"})
            results.append(properties[["TrackID", "position_t", "percentage_dead_mask", "nr_dead_mask_pixels"]])

    # If no segments/dead-mask pixels were found across all timepoints,
    # return an empty but well-formed DataFrame so downstream code can continue.
    if not results:
        return pd.DataFrame(
            columns=[
                "TrackID",
                "position_t",
                "percentage_dead_mask",
                "nr_dead_mask_pixels",
                "increase_dead_mask",
            ]
        )

    df_intensity = pd.concat(results, ignore_index=True)
    if df_intensity.empty:
        # Should not happen because of the guard above, but keep this
        # as an extra safety net.
        df_intensity["increase_dead_mask"] = pd.Series(dtype=float)
        return df_intensity

    df_intensity = df_intensity.sort_values(by=["TrackID", "position_t"]).reset_index(drop=True)

    rel_increase = calculate_relative_increase(
        df=df_intensity,
        column="nr_dead_mask_pixels",
        nr_timepoints_back=10,
        groupby="TrackID",
    )

    # Ensure we always assign a 1D array/Series, even if calculation returns empty.
    rel_increase = np.asarray(rel_increase).reshape(-1,) if np.asarray(rel_increase).ndim > 1 else rel_increase
    df_intensity["increase_dead_mask"] = rel_increase
    return df_intensity

def _basic_counts_single_timepoint(args):
    """
    Worker: compute nr_pixels and volume for one timepoint using np.bincount.
    """
    t, segments_path, voxel_volume = args
    seg_t = load_image_timepoint(segments_path, t)
    labels = seg_t.astype(np.int64).ravel()

    if labels.size == 0:
        return pd.DataFrame(columns=["TrackID", "position_t", "nr_pixels", "volume"])

    # Count voxels per label (0 is background)
    counts = np.bincount(labels)
    if counts.size <= 1:
        # no foreground labels
        return pd.DataFrame(columns=["TrackID", "position_t", "nr_pixels", "volume"])

    # Build per-object table (skip background index 0)
    track_ids = np.nonzero(counts[1:] > 0)[0] + 1
    nr_pixels = counts[track_ids]
    volumes = nr_pixels.astype(np.float64) * float(voxel_volume)

    return pd.DataFrame(
        {
            "TrackID": track_ids.astype(np.int64),
            "position_t": np.full(track_ids.shape[0], t, dtype=np.int64),
            "nr_pixels": nr_pixels.astype(np.int64),
            "volume": volumes.astype(np.float64),
        }
    )

def calculate_basic_morphology(segments_path, voxel_spacing=(1.0, 1.0, 1.0), n_workers=8):
    """
    Fast calculation of nr_pixels and volume per object per timepoint.

    Parameters:
        segments_path: path to segments .zarr
        voxel_spacing: (z, y, x) physical spacing
        n_workers: number of processes

    Returns:
        pd.DataFrame with columns: TrackID, position_t, nr_pixels, volume
    """
    segments = load_image(segments_path)
    timepoints = int(segments.shape[0])
    vz, vy, vx = map(float, voxel_spacing)
    voxel_volume = vz * vy * vx

    args_list = [(t, segments_path, voxel_volume) for t in range(timepoints)]
    if n_workers and n_workers > 1:
        results = _run_parallel_with_fallback(_basic_counts_single_timepoint, args_list, int(n_workers))
    else:
        results = [
            _basic_counts_single_timepoint(args)
            for args in tqdm(
                args_list, 
                # desc="Fast morphology (nr_pixels & volume)"
                )
        ]

    if not results:
        return pd.DataFrame(columns=["TrackID", "position_t", "nr_pixels", "volume"])

    return pd.concat(results, ignore_index=True)

def _calculate_morphology_single_timepoint(args):
    """
    Function to compute morphology features for a single timepoint.

    Features calculated:
    - volume
    - bbox_volume
    - extent (orientation-invariant, principal-axis-based)
    - solidity
    - equivalent_diameter
    - major_axis_length
    - minor_axis_length
    - elongation
    - surface_area
    - sphericity
    - convex_volume
    - orientation_vector
    - oblateness
    - axis_length_a
    - axis_length_b
    - axis_length_c
    
    """
    t, segments_path, voxel_spacing = args
    stack = load_image_timepoint(segments_path, t)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*convex hull image.*")
        warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide")

        properties = pd.DataFrame(
            regionprops_table(
                label_image=stack,
                properties=[
                    "label",
                    "num_pixels",
                    "area",
                    "bbox_area",
                    "solidity",
                    "equivalent_diameter",
                    "inertia_tensor",
                    "inertia_tensor_eigvals",
                ],
                spacing=voxel_spacing
            )
        )

    properties.rename(columns={
        "label": "TrackID",
        "area": "volume",
        "bbox_area": "bbox_volume",
        "num_pixels": "nr_pixels",
    }, inplace=True)

    properties["position_t"] = t

    surface_areas = []
    sphericities = []
    convex_volumes = []
    axis_length_a_list = []
    axis_length_b_list = []
    axis_length_c_list = []
    major_axis_lengths = []
    minor_axis_lengths = []
    elongations = []
    extents = []
    surface_to_volume_ratios = []
    orientation_vectors = []
    oblateness_list = []
    prolateness_list = []
    border_touching_segments = []

    label_slices = find_objects(stack)
    stack_shape = stack.shape

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide")
        warnings.filterwarnings("ignore", message=".*convex hull image.*")

        tensor_columns = [f"inertia_tensor-{i}-{j}" for i in range(3) for j in range(3)]
        properties_by_track = properties.set_index("TrackID", drop=False)

        for region_label in properties["TrackID"]:
            slices = label_slices[region_label - 1]
            cropped = stack[slices]
            mask = (cropped == region_label)
            coords = np.argwhere(mask)
            row = properties_by_track.loc[region_label]
            volume = float(row["volume"]) if np.isfinite(row["volume"]) else np.nan
            solidity = float(row["solidity"]) if np.isfinite(row["solidity"]) else np.nan

            border_touching_segments.append(bool(
                _segment_touches_border(slices=slices, stack_shape=stack_shape)
            ))

            surface_area = _voxel_surface_area(mask, voxel_spacing)
            surface_areas.append(surface_area)

            if np.isfinite(surface_area) and surface_area > 0 and np.isfinite(volume) and volume > 0:
                sphericities.append(float((np.pi ** (1.0 / 3.0)) * ((6.0 * volume) ** (2.0 / 3.0)) / surface_area))
                surface_to_volume_ratios.append(float(surface_area / volume))
            else:
                sphericities.append(np.nan)
                surface_to_volume_ratios.append(np.nan)

            if np.isfinite(solidity) and solidity > 0 and np.isfinite(volume) and volume > 0:
                convex_volumes.append(float(volume / solidity))
            else:
                convex_volumes.append(np.nan)

            try:
                metrics = _principal_axis_metrics(coords, voxel_spacing)
            except Exception:
                metrics = {
                    "axis_lengths": np.array([np.nan, np.nan, np.nan], dtype=float),
                    "orientation_vector": [np.nan, np.nan, np.nan],
                    "oblateness": np.nan,
                    "prolateness": np.nan,
                }

            axis_lengths = np.asarray(metrics["axis_lengths"], dtype=float)
            axis1 = float(axis_lengths[0]) if axis_lengths.shape == (3,) else np.nan
            axis2 = float(axis_lengths[1]) if axis_lengths.shape == (3,) else np.nan
            axis3 = float(axis_lengths[2]) if axis_lengths.shape == (3,) else np.nan

            if (not np.all(np.isfinite(axis_lengths))) or np.any(axis_lengths <= 0):
                try:
                    tensor = row[tensor_columns].to_numpy(dtype=float).reshape(3, 3)
                    if np.all(np.isfinite(tensor)):
                        eigvals, eigvecs = np.linalg.eigh(tensor)
                        order = np.argsort(eigvals)[::-1]
                        major_vec = eigvecs[:, order[0]]
                        major_norm = float(np.linalg.norm(major_vec))
                        if np.isfinite(major_norm) and major_norm > 0:
                            metrics["orientation_vector"] = (major_vec / major_norm).tolist()
                except Exception:
                    pass

            axis_length_a_list.append(axis1)
            axis_length_b_list.append(axis2)
            axis_length_c_list.append(axis3)
            major_axis_lengths.append(axis1)
            minor_axis_lengths.append(axis3)
            orientation_vectors.append(metrics["orientation_vector"])
            oblateness_list.append(metrics["oblateness"])
            prolateness_list.append(metrics["prolateness"])

            if np.isfinite(axis1) and np.isfinite(axis2) and axis2 > 0:
                elongations.append(float(axis1 / axis2))
            else:
                elongations.append(np.nan)

            obb_volume = axis1 * axis2 * axis3
            if np.isfinite(obb_volume) and obb_volume > 0 and np.isfinite(volume) and volume > 0:
                extents.append(float(volume / obb_volume))
            else:
                extents.append(np.nan)

    properties["surface_area"] = surface_areas
    properties["sphericity"] = sphericities
    properties["convex_volume"] = convex_volumes
    properties["surface_to_volume_ratio"] = surface_to_volume_ratios
    properties["axis1_length"] = axis_length_a_list
    properties["axis2_length"] = axis_length_b_list
    properties["axis3_length"] = axis_length_c_list
    properties["major_axis_length"] = major_axis_lengths
    properties["minor_axis_length"] = minor_axis_lengths
    properties["elongation"] = elongations
    properties["extent"] = extents
    properties["orientation_vector"] = orientation_vectors
    properties["oblateness"] = oblateness_list
    properties["prolateness"] = prolateness_list
    properties["border_touching_segment"] = border_touching_segments

    # Drop heavy intermediate columns
    columns_to_drop = [col for col in properties.columns
                       if col.startswith("inertia_tensor")
                       or col.startswith("inertia_tensor_eigvals")
                       or col.startswith("moments_central")]
    properties = properties.drop(columns=columns_to_drop)

    return properties

def calculate_morphology_features(segments_path, voxel_spacing=(1.0, 1.0, 1.0), n_workers=8):
    """
    Calculates morphological features (volume, shape, sphericity, etc.) for 3D segments.
    
    Parameters:
        segments_path: path to segments .zarr
        voxel_spacing: Tuple of (z, y, x) spacing in physical units (e.g., µm). Default is isotropic.
    
    Returns:
        pd.DataFrame with one row per object per time point.
    """
    segments=load_image(segments_path)
    timepoints = segments.shape[0]
    args_list = [(t, segments_path, voxel_spacing) for t in range(timepoints)]
    if n_workers > 1:
        results = _run_parallel_with_fallback(_calculate_morphology_single_timepoint, args_list, n_workers)
    else:
        results = [_calculate_morphology_single_timepoint(args) for args in tqdm(args_list)]

    return pd.concat(results, ignore_index=True)

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
#     parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
#     args = parser.parse_args()
#     with open(args.config, "r") as parameters:
#         config=yaml.load(parameters, Loader=yaml.SafeLoader)
#     metadata = pd.read_csv(config["metadata_csv_path"])
#     run_behav3d_feature_extraction(config, metadata)
