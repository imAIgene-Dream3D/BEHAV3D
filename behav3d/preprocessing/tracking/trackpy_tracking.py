from pathlib import Path
import numpy as np
import pandas as pd
import shutil
from skimage.measure import regionprops, regionprops_table
from scipy.ndimage import distance_transform_edt
from behav3d.utils.fileio import load_image, save_as_zarr

try:
    import trackpy as tp
except ImportError:
    raise ImportError(
        "trackpy is required for tracking. "
        "Install via: pip install trackpy"
    )



def mask_from_one_hot_encoding_3D(one_hot):
    mask = np.zeros((one_hot.shape[1:]), dtype='int')
    for i in range(one_hot.shape[0]):
        mask[one_hot[i]==True] = i
    return mask

def one_hot_encoding_3D(img):
    object_labels=np.unique(img)
    one_hot = np.zeros((len(object_labels), img.shape[0], img.shape[1], img.shape[2]), dtype='bool')
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
    df_tracks = tp.link(df_centroids, search_range, memory=memory, adaptive_stop=adaptive_stop, adaptive_step=adaptive_step,pos_columns=['position_z','position_y','position_x'],t_column='position_t')
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
    


def run_organoid_trackpy_tracking(
    metadata,
    output_dir,
    overwrite=False,
    cell_type="tcell",
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    search_range=31,
    memory=2,
    adaptive_stop=10.0,
    adaptive_step=0.95,
    tracking_method='trackpy_post',
):
    """Run tracking on segmented objects.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        DataFrame containing sample information
    output_dir : str or Path
        Root output directory
    cell_type : str
        Either "tcell" or "organoid"
    tracking_method : str
        Either "trackpy_post" or "cc3D"
    """
    output_dir = Path(output_dir)
    
    for idx, sample in metadata.iterrows():
        sample_name = sample["sample_name"]
        print(f"\nProcessing {sample_name}...")
        
        # Set up paths
        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        
        segments_path = sample[f"{cell_type}_segments_image_path"]
        tracked_img_outpath = Path(tracked_img_outdir, f"{sample_name}_{cell_type}_tracked.zarr")
        tracked_csv_outpath = Path(tracked_csv_outdir, f"{sample_name}_{cell_type}_tracks.csv")
    
        if not tracked_img_outdir.exists():
            tracked_img_outdir.mkdir(parents=True)
        if not tracked_csv_outdir.exists():
            tracked_csv_outdir.mkdir(parents=True)
        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]
        segments = load_image(segments_path)

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

        # Update metadata with BEHAV3D standard paths
        metadata.at[idx, f"{cell_type}_tracks_image_path"] = str(tracked_img_outpath)
        metadata.at[idx, f"{cell_type}_tracks_csv_path"] = str(tracked_csv_outpath)

    return metadata