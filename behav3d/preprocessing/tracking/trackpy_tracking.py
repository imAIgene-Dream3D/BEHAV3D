from pathlib import Path
import sys
import numpy as np
import pandas as pd
import shutil
from skimage.measure import regionprops, regionprops_table
from scipy.ndimage import distance_transform_edt
from behav3d.io.images import load_image, save_as_zarr
from behav3d.preprocessing.tracking import convert_segments_to_tracks
from typing import Optional, Tuple, Sequence

try:
    import trackpy as tp
except ImportError:
    raise ImportError(
        "trackpy is required for tracking. "
        "Install via: pip install trackpy"
    )
tp.quiet()

# ---------------------------------------------------------------------------
# Raise TrackPy's subnet size limit from the default (30) to 300.
# Dense datasets (e.g. 100×100 µm FOV) can generate subnets with >30 points,
# causing SubnetOversizeException.  300 is a reasonable trade-off between
# allowing real-world data and avoiding factorial-time blowup.
# ---------------------------------------------------------------------------
_BEHAV3D_MAX_SUB_NET_SIZE = 300
try:
    import trackpy.linking as _tpl
    _ORIGINAL_MAX = getattr(_tpl.Linker, 'MAX_SUB_NET_SIZE', 30)
    _tpl.Linker.MAX_SUB_NET_SIZE = _BEHAV3D_MAX_SUB_NET_SIZE
except Exception:
    _ORIGINAL_MAX = 30

def mask_from_one_hot_encoding_3D(one_hot):
    mask = np.zeros((one_hot.shape[1:]), dtype='int')
    for i in range(one_hot.shape[0]):
        mask[one_hot[i]==True] = i
    return mask

def one_hot_encoding_3D(img):
    object_labels=np.unique(img)
    one_hot = np.zeros((len(object_labels), img.shape[0], img.shape[1], img.shape[2]), dtypef='bool')
    for i, unique_value in enumerate(object_labels):
        one_hot[i,...][img == unique_value] = True
    return one_hot

def calculate_ious(outputs, targets):
    outputs_2D = np.max(outputs, axis=1)  # Project along Z axis
    targets_2D = np.max(targets, axis=1)
    ious = np.zeros((outputs.shape[0], targets.shape[0]), dtype='float32')
    for i in range(outputs.shape[0]):
        for j in range(targets.shape[0]):
            if (np.logical_and(outputs_2D[i], targets_2D[j]).sum() > 0):
                ious[i, j] = iou_score(outputs[i], targets[j])
    return ious
        
# Segmentation metrics
def iou_score(outputs, targets):
    intersection = np.logical_and(outputs, targets).sum()
    union = outputs.sum() + targets.sum() - intersection
    if union != 0:
        iou = intersection / union
    else:
        iou = 0.0  # No overlap and no objects
    return iou

def track_by_overlap_3d(prev_boxes, curr_boxes, iou_threshold=0.3):
    iou_matrix = calculate_ious(prev_boxes, curr_boxes)
    iou_matched = iou_matrix > iou_threshold
    curr_inds = np.where(np.sum(iou_matched, axis=0) > 1)[0]
    for j in range(curr_inds.size):
        row_inds = np.where(iou_matched[:, curr_inds[j]] == True)[0]
        reg_max = np.argmax(iou_matrix[row_inds, curr_inds[j]])
        to_zero = np.setdiff1d(np.arange(row_inds.size), reg_max)
        iou_matrix[row_inds[to_zero], curr_inds[j]] = 0

    row_ind, col_ind = np.where(iou_matrix > iou_threshold)
    matches = np.concatenate((row_ind[:, None], col_ind[:, None]), axis=1)
    return matches



def run_tracking_by_cc(segments, backpropagate=True, iou_threshold=0.1):
    tracked_segments = np.zeros_like(segments)
    labels_prev = -1 * np.ones(segments[0].max() + 1)
    for t in range(segments.shape[0] - 1):
        print(f"Tracking frame {t}")
        _, indices = distance_transform_edt(segments[t] == 0, return_indices=True)
        voronoi = segments[t][tuple(indices)]
        prev = one_hot_encoding_3D(voronoi)
        curr = one_hot_encoding_3D(segments[t + 1])[1:]
        matches = track_by_overlap_3d(prev, curr, iou_threshold=iou_threshold)
        prev = one_hot_encoding_3D(segments[t])[1:]
        ts_prev = tracked_segments[t]
        ts_curr = tracked_segments[t + 1]
        if t == 0:
            for k in range(len(matches)):
                ts_prev[prev[matches[k][0]]] = k + 1
                labels_prev[matches[k][0]] = k + 1
        labels_curr = -1 * np.ones(curr.shape[0])
        for k in range(len(matches)):
            if labels_prev[matches[k][0]] >= 0:
                ts_curr[curr[matches[k][1]]] = labels_prev[matches[k][0]]
                labels_curr[matches[k][1]] = labels_prev[matches[k][0]]
            else:
                new_label = tracked_segments.max() + 1
                ts_prev[prev[matches[k][0]]] = new_label
                labels_prev[matches[k][0]] = new_label
                ts_curr[curr[matches[k][1]]] = new_label
                labels_curr[matches[k][1]] = new_label
        labels_prev = labels_curr

    if backpropagate:
        tracked_segments = run_tracking_by_cc(np.flip(tracked_segments, axis=0), backpropagate=False, iou_threshold=iou_threshold)
        return np.flip(tracked_segments, axis=0)
    else:
        return tracked_segments

def run_trackpy_tracking(
    segments,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    search_range=31,
    memory=2,
    adaptive_stop=10.0,
    adaptive_step=0.95, 
    return_trackimg=True,
    ):
    
    df_centroids = []
    for t, object_stack in enumerate(segments):
        properties=pd.DataFrame(regionprops_table(label_image=object_stack, properties=['label', f'centroid']))
        properties["position_t"]=t
        df_centroids.append(properties)
    df_centroids = pd.concat(df_centroids)
    df_centroids["position_z"]=df_centroids["centroid-0"]*element_size_z
    df_centroids["position_y"]=df_centroids["centroid-1"]*element_size_y
    df_centroids["position_x"]=df_centroids["centroid-2"]*element_size_x
    # Tracking

    df_tracks = tp.link(
        df_centroids, 
        search_range, 
        memory=memory, 
        adaptive_stop=adaptive_stop, 
        adaptive_step=adaptive_step,
        pos_columns=['position_z','position_y','position_x'],
        t_column='position_t'
        )
    df_tracks=df_tracks.rename(columns={'particle':'track_id'})
    df_tracks = df_tracks.reset_index()
    df_tracks["track_id"]+=1
    if return_trackimg:
        tracked_img = np.zeros_like(segments)

        for t, t_seg in enumerate(segments):
            t_df_tracks = df_tracks[df_tracks["position_t"]==t]
            for _, row in t_df_tracks.iterrows():
                # print(row["label"], row["track_id"], (tracked_img[t]==row["label"]).any())
                tracked_img[t][t_seg==row["label"]] = row["track_id"]

        return df_tracks, tracked_img
    else:
        return df_tracks
    


def run_trackpy_tracking_generic(
    metadata,
    output_dir,
    cell_type,
    overwrite=False,
    search_range=31,
    memory=2,
    adaptive_stop=10.0,
    adaptive_step=0.95,
    n_workers=1,
    return_trackimg=True,
    log_callback=None,
    progress_cb=None,
    **kwargs
):
    """Run trackpy tracking on any cell type.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        DataFrame containing sample information
    output_dir : str or Path
        Root output directory
    cell_type : str
        Name of the cell type to track 
    overwrite : bool
        Whether to overwrite existing tracking results
    search_range : float
        Maximum distance cells can move between frames (in physical units)
    memory : int
        Number of frames objects can disappear before being considered lost
    adaptive_stop : float
        Stop adaptive search when correlation coefficient falls below this value
    adaptive_step : float
        Reduce search range by this factor when adaptive search fails
    return_trackimg : bool
        Whether to return and save tracked image
    progress_cb : callable, optional
        Called as ``progress_cb(current, total, label)`` from inside the
        per-sample loop so a GUI can drive a progress bar.  ``None`` is a
        no-op (default).  Notebook callers omit this kwarg.
    **kwargs : dict
    """
    output_dir = Path(output_dir)

    total_samples = len(metadata)
    for i, (idx, sample) in enumerate(metadata.iterrows()):
        sample_name = sample["sample_name"]
        if progress_cb is not None:
            try:
                progress_cb(i, total_samples, f"{cell_type} / {sample_name}")
            except Exception:
                pass
        print(f"\nProcessing {sample_name}...")
        
        # Set up paths
        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        
        # Find the correct prefixed column (or_, im_, ot_)
        segments_col = None
        for prefix in ['or', 'im', 'ot']:
            col_name = f"{prefix}_{cell_type}_segments_image_path"
            if col_name in sample.index and pd.notna(sample[col_name]):
                segments_col = col_name
                break
        
        if segments_col is None:
            raise ValueError(f"No segments found for {cell_type} in sample {sample_name}. Expected prefixed column (or_/im_/ot_).")
        
        segments_path = sample[segments_col]
        tracked_img_outpath = Path(tracked_img_outdir, f"{sample_name}_{cell_type}_tracked.zarr")
        tracked_csv_outpath = Path(tracked_csv_outdir, f"{sample_name}_{cell_type}_tracks.csv")
    
        if not tracked_img_outdir.exists():
            tracked_img_outdir.mkdir(parents=True)
        if not tracked_csv_outdir.exists():
            tracked_csv_outdir.mkdir(parents=True)
            
        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]
        
        if (
            (
                not tracked_csv_outpath.exists() or 
                not tracked_img_outpath.exists()
            ) or overwrite
            ):
            
            print(f"Loading segments from {segments_path}")
            segments = load_image(segments_path)
            print(f"Running trackpy tracking on {len(segments)} frames")
            
            df_centroids = []
            for t, object_stack in enumerate(segments):
                object_stack = np.asarray(object_stack)
                properties = pd.DataFrame(regionprops_table(label_image=object_stack, properties=['label', f'centroid']))
                properties["position_t"] = t
                df_centroids.append(properties)
            df_centroids = pd.concat(df_centroids)
            df_centroids["position_z"] = df_centroids["centroid-0"]*element_size_z
            df_centroids["position_y"] = df_centroids["centroid-1"]*element_size_y
            df_centroids["position_x"] = df_centroids["centroid-2"]*element_size_x

            # print("Running trackpy tracking...")
            # Tracking
            from trackpy.linking.utils import SubnetOversizeException

            if _BEHAV3D_MAX_SUB_NET_SIZE != _ORIGINAL_MAX:
                _msg = (f"[TrackPy] MAX_SUB_NET_SIZE raised from {_ORIGINAL_MAX} "
                        f"to {_BEHAV3D_MAX_SUB_NET_SIZE}. Large subnets may be slow.")
                print(_msg, file=sys.stderr)
                if log_callback:
                    log_callback(f"⚠️ {_msg}")

            try:
                df_tracks = tp.link(
                    df_centroids,
                    search_range=search_range,
                    memory=memory,
                    adaptive_stop=adaptive_stop,
                    adaptive_step=adaptive_step,
                    pos_columns=['position_z', 'position_y', 'position_x'],
                    t_column='position_t',
                )
            except SubnetOversizeException as e:
                _warn = (f"[WARNING] Numba linking hit subnet limit: {e}\n"
                         f"Falling back to recursive (slower) linking strategy...")
                print(_warn, file=sys.stderr)
                if log_callback:
                    log_callback(f"⚠️ {_warn}")

                df_tracks = tp.link(
                    df_centroids,
                    search_range=search_range,
                    memory=memory,
                    adaptive_stop=adaptive_stop,
                    adaptive_step=adaptive_step,
                    pos_columns=['position_z', 'position_y', 'position_x'],
                    t_column='position_t',
                    link_strategy='recursive'
                )
                 
            # Rename columns to match BEHAV3D convention
            df_tracks = df_tracks.reset_index()
            df_tracks.rename(
                columns={
                    'particle': 'TrackID',
                    'label': 'SegmentID',
                    'centroid-0': "pixel_position_z",
                    'centroid-1': "pixel_position_y",
                    'centroid-2': "pixel_position_x",
                }, 
                inplace=True
            )
            df_tracks["TrackID"] += 1
            
            # Select only needed columns
            df_tracks = df_tracks[[
                "TrackID", "SegmentID", "position_t", "position_x", "position_y", "position_z",
                "pixel_position_x", "pixel_position_y", "pixel_position_z"
            ]]
            
            print(f"Saving tracks CSV to {tracked_csv_outpath}")
            df_tracks.to_csv(tracked_csv_outpath, sep=",", index=False)
            
            if return_trackimg:
                print(f"Converting segments to tracks and saving to {tracked_img_outpath}")
                tracked_img_outpath = convert_segments_to_tracks(
                    df_tracks=df_tracks,
                    segments=segments,
                    outpath=tracked_img_outpath,
                    n_workers=n_workers,
                )
        else:
            print("Tracking already exists... Provide overwrite=True to overwrite... Loading existing tracking data")
        
        # Update metadata with prefixed column names
        if segments_col is not None and segments_col.startswith(('or_', 'im_', 'ot_')):
            # Use the same prefix as the segments column
            prefix = segments_col.split('_')[0]
            img_col = f"{prefix}_{cell_type}_tracks_image_path"
            csv_col = f"{prefix}_{cell_type}_tracks_csv_path"
        else:
            # Fallback to old non-prefixed format
            img_col = f"{cell_type}_tracks_image_path"
            csv_col = f"{cell_type}_tracks_csv_path"
            
        # Ensure columns are object dtype to avoid LossySetitemError
        for col in [img_col, csv_col]:
            if col not in metadata.columns or metadata[col].dtype != object:
                metadata[col] = metadata.get(col, pd.Series(dtype=object)).astype(object)
                
        metadata.at[idx, img_col] = str(tracked_img_outpath)
        metadata.at[idx, csv_col] = str(tracked_csv_outpath)

    if progress_cb is not None:
        try:
            progress_cb(total_samples, total_samples, f"{cell_type} done")
        except Exception:
            pass
    return metadata


def run_trackpy_tracking_with_method(
    metadata,
    output_dir,
    cell_type,
    overwrite=False,
    search_range=31,
    memory=2,
    adaptive_stop=10.0,
    adaptive_step=0.95,
    tracking_method='trackpy_post',
    **kwargs
):
    """Run tracking on any cell type with specified tracking method.
        
    Parameters
    ----------
    metadata : pd.DataFrame
        DataFrame containing sample information
    output_dir : str or Path
        Root output directory
    cell_type : str
        Name of the cell type to track (e.g., 'organoid1', 'organoid2', 'tcell', 'macro', 'tum')
    overwrite : bool
        Whether to overwrite existing tracking results
    search_range : float
        Maximum distance cells can move between frames (in physical units)
    memory : int
        Number of frames objects can disappear before being considered lost
    adaptive_stop : float
        Stop adaptive search when correlation coefficient falls below this value
    adaptive_step : float
        Reduce search range by this factor when adaptive search fails
    tracking_method : str
        Either "trackpy_post" or "cc3D"
    **kwargs : dict
    """
    output_dir = Path(output_dir)
    
    for idx, sample in metadata.iterrows():
        sample_name = sample["sample_name"]
        print(f"\nProcessing {sample_name}...")
        
        # Set up paths
        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        
        # Find the correct prefixed column (or_, im_, ot_)
        segments_col = None
        for prefix in ['or', 'im', 'ot']:
            col_name = f"{prefix}_{cell_type}_segments_image_path"
            if col_name in sample.index and pd.notna(sample[col_name]):
                segments_col = col_name
                break
        
        if segments_col is None:
            raise ValueError(f"No segments found for {cell_type} in sample {sample_name}. Expected prefixed column (or_/im_/ot_).")
        
        segments_path = sample[segments_col]
        tracked_img_outpath = Path(tracked_img_outdir, f"{sample_name}_{cell_type}_tracked.zarr")
        tracked_csv_outpath = Path(tracked_csv_outdir, f"{sample_name}_{cell_type}_tracks.csv")
    
        if not tracked_img_outdir.exists():
            tracked_img_outdir.mkdir(parents=True)
        if not tracked_csv_outdir.exists():
            tracked_csv_outdir.mkdir(parents=True)
        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]
        segments = np.asarray(load_image(segments_path))

        if (
            (
                not tracked_csv_outpath.exists() or 
                not tracked_img_outpath.exists()
            ) or overwrite
            ):
            
            print(f"Running {tracking_method} tracking...")
            if tracking_method == 'trackpy_post':
                df_tracks, tracked_img = run_trackpy_tracking(
                    segments=segments,
                    element_size_x=element_size_x,
                    element_size_y=element_size_y,
                    element_size_z=element_size_z,
                    search_range=search_range,
                    memory=memory,
                    adaptive_stop=adaptive_stop,
                    adaptive_step=adaptive_step,
                    return_trackimg=True,
                )
                # Rename and reorder columns to your required format:
                # [TrackID, SegmentID, position_t, position_x, position_y, position_z, pixel_position_x, pixel_position_y, pixel_position_z]

                df_tracks = df_tracks.rename(columns={
                    'track_id': 'TrackID',
                    'label': 'SegmentID'
                })
                df_tracks['pixel_position_z'] = df_tracks['centroid-0']
                df_tracks['pixel_position_y'] = df_tracks['centroid-1']
                df_tracks['pixel_position_x'] = df_tracks['centroid-2']
                df_tracks['position_z'] = df_tracks['centroid-0']*element_size_z        
                df_tracks['position_y'] = df_tracks['centroid-1']*element_size_y
                df_tracks['position_x'] = df_tracks['centroid-2']*element_size_x
                df_tracks = df_tracks[[
                    'TrackID',
                    'SegmentID',
                    'position_t',
                    'position_x',
                    'position_y',
                    'position_z',
                    'pixel_position_x',
                    'pixel_position_y',
                    'pixel_position_z',
                ]]
                print(f"Found {len(df_tracks['TrackID'].unique())} tracks")
            elif tracking_method == 'cc3D':
                tracked_img = run_tracking_by_cc(
                    segments=segments,
                    backpropagate=True,
                    iou_threshold=0.01,
                )
                print(f"Found {len(np.unique(tracked_img)) - 1} tracks")  # -1 for background
                df_tracks = pd.DataFrame()  # Placeholder
                
            else:
                raise ValueError(f"Tracking method {tracking_method} not implemented.")

            # Verify tracked_img values fit in uint16
            max_val = tracked_img.max()
            if max_val >= 65535:
                print(f"[WARNING] Track IDs exceed uint16 range (max={max_val})")
                
            # Save outputs
            # Clean up any existing tracking files if we're overwriting
            if tracked_img_outpath.exists():
                shutil.rmtree(tracked_img_outpath)
                
            print(f"Saving tracked image to {tracked_img_outpath}")
            save_as_zarr(tracked_img.astype(np.uint16), tracked_img_outpath)
            
            if not df_tracks.empty:
                print(f"Saving tracks CSV to {tracked_csv_outpath}")
                
                df_tracks.to_csv(tracked_csv_outpath, index=False)
        else:
            print("Tracking already exists... Provide overwrite=True to overwrite... Loading existing tracking data")

        # Update metadata with prefixed column names
        if segments_col is not None and segments_col.startswith(('or_', 'im_', 'ot_')):
            # Use the same prefix as the segments column
            prefix = segments_col.split('_')[0]
            img_col = f"{prefix}_{cell_type}_tracks_image_path"
            csv_col = f"{prefix}_{cell_type}_tracks_csv_path"
        
        # Ensure columns are object dtype
        for col in [img_col, csv_col]:
            if col not in metadata.columns or metadata[col].dtype != object:
                metadata[col] = metadata.get(col, pd.Series(dtype=object)).astype(object)
                
        metadata.at[idx, img_col] = str(tracked_img_outpath)
        metadata.at[idx, csv_col] = str(tracked_csv_outpath)

    return metadata
