from behav3d.utils.fileio import load_image, append_to_zarr, save_as_zarr
from behav3d.utils.preprocessing import dilate_mask
from behav3d.utils.segmentation import segment_size_filter
from behav3d.utils.tracking import convert_tracked_image_to_csv
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
    cell_type="organoid",
    overwrite=False,
    **kwargs
    ):
    for idx, sample in metadata.iterrows():
        sample_name=sample['sample_name']
        print(f"Tracking sample: {sample_name}")

        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]
        
        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        
        segments_path = sample[f"{cell_type}_segments_image_path"]
        
        # Check if segments_path is valid
        if pd.isna(segments_path) or segments_path is None:
            print(f"Warning: No segmentation data found for {sample_name}. Skipping tracking.")
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
            )
        else:
            print("Tracking already exists... Provide overwrite=True to overwrite... Loading existing tracking data")
        
        metadata.at[idx, f"{cell_type}_tracks_image_path"] = str(tracked_img_outpath)
        metadata.at[idx, f"{cell_type}_tracks_csv_path"] = str(tracked_csv_outpath)
        
    return metadata