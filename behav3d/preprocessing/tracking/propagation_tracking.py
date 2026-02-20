from behav3d.io.images import load_image, append_to_zarr, save_as_zarr
from behav3d.preprocessing import dilate_mask
from behav3d.preprocessing.segmentation import segment_size_filter
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv
import pandas as pd

from skimage.segmentation import watershed
import numpy as np
import pandas as pd
from pathlib import Path
import zarr
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import shutil
import dask.array as da

# def propagate_segmentation(args):
#     zarr_path, t, dilation_nr_pixels, segment_size_min = args
#     t_seg = np.asarray(load_image(zarr_path)[t])
#     seg_prev_tp = np.asarray(load_image(zarr_path)[max(0, t-1)])
#     mask = t_seg!=0
#     seg_prev_tp[mask==0]=0
#     # seeds = keep_largest_connected_components(seeds)
#     new_seg = watershed(mask, markers = seg_prev_tp, mask=mask)
#     mask_dilated = dilate_mask(mask, nr_pixels=dilation_nr_pixels)
#     new_seg = watershed(mask_dilated, markers = new_seg, mask=mask_dilated)
#     new_seg[mask==0]=0
#     new_seg = segment_size_filter(new_seg, size_min=segment_size_min)
#     return new_seg

def _run_single_timepoint_propagation(
    t_seg,
    seg_prev_tp,
    dilation_nr_pixels=2,
    segment_size_min=100,
    ):

    mask = t_seg != 0
    seg_prev_tp[mask==0]=0
    # seeds = keep_largest_connected_components(seeds)
    new_seg = watershed(mask, markers = seg_prev_tp, mask=mask)
    mask_dilated = dilate_mask(mask, nr_pixels=dilation_nr_pixels)
    new_seg = watershed(mask_dilated, markers = new_seg, mask=mask_dilated)
    new_seg[mask==0]=0
    new_seg = segment_size_filter(new_seg, size_min=segment_size_min)
    return(new_seg)
    # return new_seg

def propagate_tracks(
    segments_path,
    tracked_img_outpath,
    tracked_csv_outpath,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    dilation_nr_pixels=2,
    segment_size_min=100,
    **kwargs
    ):
    
    seg = load_image(segments_path)
    
    if tracked_img_outpath.exists():
        shutil.rmtree(tracked_img_outpath)
    
    seg_prev_tp = seg[0]
    for t, t_seg in tqdm(enumerate(seg), total=seg.shape[0]):
        t_seg = np.asarray(t_seg)
        
        if t==0:
            t_tracked_seg = np.asarray(t_seg)
        else:
            t_tracked_seg = _run_single_timepoint_propagation(
                t_seg,
                seg_prev_tp,
                dilation_nr_pixels=dilation_nr_pixels,
                segment_size_min=segment_size_min,
            )
        seg_prev_tp = t_tracked_seg.copy()
        t_tracked_seg = np.expand_dims(t_tracked_seg, axis=0)
        append_to_zarr(t_tracked_seg, tracked_img_outpath)
        
    df_tracks = convert_tracked_image_to_csv(
            img_path=tracked_img_outpath,
            outpath=tracked_csv_outpath,
            element_size_x=element_size_x,
            element_size_y=element_size_y,
            element_size_z=element_size_z
        )
    # Save the tracked image as zarr
def run_propagation_tracking(
    metadata,
    output_dir,
    cell_type,
    overwrite=False,
    dilation_nr_pixels=2,
    segment_size_min=100,
    **kwargs
    ):
    """Run propagation-based tracking on any cell type.
    
    This method propagates segmentation from one timepoint to the next using
    watershed with the previous timepoint's segmentation as markers.
    Typically used for large, slow-moving objects like organoids.
    
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
    dilation_nr_pixels : int
        Number of pixels to dilate the mask for watershed propagation
    segment_size_min : int
        Minimum segment size in voxels (smaller segments are filtered out)
    **kwargs : dict
    """
    for idx, sample in metadata.iterrows():
        sample_name=sample['sample_name']
        print(f"Tracking sample: {sample_name}")

        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]
        
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
            # Fallback to old non-prefixed format for backward compatibility
            segments_col = f"{cell_type}_segments_image_path"
        
        segments_path = sample.get(segments_col)
        
        # Check if segments_path is valid
        if pd.isna(segments_path) or segments_path is None:
            print(f"Warning: No segmentation data found for {sample_name}, {cell_type}. Skipping tracking.")
            continue
            
        segments_path = Path(segments_path)
        if not segments_path.exists():
            print(f"Warning: Segmentation file not found: {segments_path}. Skipping tracking.")
            continue
            
        tracked_img_outpath = Path(tracked_img_outdir, f"{sample_name}_{cell_type}_tracked.zarr")
        tracked_csv_outpath = Path(tracked_csv_outdir, f"{sample_name}_{cell_type}_tracks.csv")
    
        if not tracked_img_outdir.exists():
            tracked_img_outdir.mkdir(parents=True)
        if not tracked_csv_outdir.exists():
            tracked_csv_outdir.mkdir(parents=True)
        if (
            (
                not tracked_csv_outpath.exists() or 
                not tracked_img_outpath.exists()
            ) or overwrite
            ):
            propagate_tracks(
                segments_path=segments_path,
                tracked_img_outpath=tracked_img_outpath,
                tracked_csv_outpath=tracked_csv_outpath,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                dilation_nr_pixels=dilation_nr_pixels,
                segment_size_min=segment_size_min,
            )
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