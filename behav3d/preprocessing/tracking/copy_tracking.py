from behav3d.utils.fileio import load_image, append_to_zarr
from behav3d.utils.preprocessing import dilate_mask
from behav3d.utils.segmentation import segment_size_filter
from skimage.segmentation import watershed
import numpy as np
import pandas as pd
from pathlib import Path
import zarr


def propagate_segmentation(args):
    zarr_path, t, dilation_nr_pixels, segment_size_min = args
    t_seg = zarr.open(zarr_path, mode='r+')[t]
    seg_prev_tp = zarr.open(zarr_path, mode='r+')[max(0, t-1)]
    mask = t_seg!=0
    seg_prev_tp[mask==0]=0
    # seeds = keep_largest_connected_components(seeds)
    new_seg = watershed(mask, markers = seg_prev_tp, mask=mask)
    mask_dilated = dilate_mask(mask, nr_pixels=dilation_nr_pixels)
    new_seg = watershed(mask_dilated, markers = new_seg, mask=mask_dilated)
    new_seg[mask==0]=0
    new_seg = segment_size_filter(new_seg, size_min=segment_size_min)
    append_to_zarr(new_seg, tracked_img_outpath)
    
def copytrack_image(
    segments_path,
    tracked_img_outpath,
    tracked_csv_outpath,
    dilation_nr_pixels=2,
    segment_size_min=100,
    **kwargs
    ):
    
    seg = load_image(segments_path)
    seg_prev_tp = seg[0]
    
    for t, t_seg in enumerate(seg):
        t_seg = np.asarray(t_seg)
        seeds = seg_prev_tp.copy()
        mask = t_seg!=0
        seeds[mask==0]=0
        # seeds = keep_largest_connected_components(seeds)
        new_seg = watershed(mask, markers = seeds, mask=mask)
        mask_dilated = dilate_mask(mask, nr_pixels=dilation_nr_pixels)
        new_seg = watershed(mask_dilated, markers = new_seg, mask=mask_dilated)
        new_seg[mask==0]=0
        new_seg = segment_size_filter(new_seg, size_min=segment_size_min)
        append_to_zarr(new_seg, tracked_img_outpath)

def run_tcell_laptracking(
    metadata,
    output_dir,
    cell_type="tcell",
    overwrite=False,
    **kwargs
    ):
    for idx, sample in metadata.iterrows():
        sample_name=sample['sample_name']
        print(f"Tracking sample: {sample_name}")
        
        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        
        segments_path = sample[f"{cell_type}_segments_image_path"]
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
            copytrack_image(
                segments_path=segments_path,
                tracked_img_outpath=tracked_img_outpath,
                tracked_csv_outpath=tracked_csv_outpath,
            )
        else:
            print("Tracking already exists... Provide overwrite=True to overwrite... Loading existing tracking data")
        
        metadata.at[idx, f"{cell_type}_tracks_image_path"] = str(tracked_img_outpath)
        metadata.at[idx, f"{cell_type}_tracks_csv_path"] = str(tracked_csv_outpath)
        
    return metadata