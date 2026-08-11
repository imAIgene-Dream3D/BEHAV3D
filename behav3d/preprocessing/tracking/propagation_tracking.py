from behav3d.io.images import load_image, append_to_zarr
from behav3d.preprocessing import dilate_mask
from behav3d.preprocessing.segmentation import segment_size_filter
from behav3d.preprocessing.tracking import convert_tracked_image_to_csv

from skimage.segmentation import watershed
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import shutil

def _watershed_from_markers(mask, markers, dilation_nr_pixels=2):
    """Two-pass watershed propagation: seed the raw mask, then close small
    gaps by re-seeding a dilated mask with the first pass's result.
    """
    new_seg = watershed(mask, markers=markers, mask=mask)
    mask_dilated = dilate_mask(mask, nr_pixels=dilation_nr_pixels)
    new_seg = watershed(mask_dilated, markers=new_seg, mask=mask_dilated)
    new_seg[mask==0]=0
    return new_seg

def _run_single_timepoint_propagation(
    t_seg,
    seg_prev_tp,
    dilation_nr_pixels=2,
    segment_size_min=100,
    ):

    mask = t_seg != 0
    seg_prev_tp[mask==0]=0
    new_seg = _watershed_from_markers(mask, seg_prev_tp, dilation_nr_pixels)
    new_seg = segment_size_filter(new_seg, size_min=segment_size_min)
    return(new_seg)

def _resolve_segments_column(sample, cell_type):
    for prefix in ['or', 'im', 'ot']:
        col_name = f"{prefix}_{cell_type}_segments_image_path"
        if col_name in sample.index and pd.notna(sample[col_name]):
            return col_name, sample[col_name]
    legacy_col = f"{cell_type}_segments_image_path"
    return legacy_col, sample.get(legacy_col)

def _resolve_tracking_output_columns(cell_type, segments_col):
    if segments_col and segments_col.startswith(('or_', 'im_', 'ot_')):
        prefix = segments_col.split('_', 1)[0]
        return (
            f"{prefix}_{cell_type}_tracks_image_path",
            f"{prefix}_{cell_type}_tracks_csv_path",
        )
    return (
        f"{cell_type}_tracks_image_path",
        f"{cell_type}_tracks_csv_path",
    )

def _ensure_metadata_output_columns(metadata, *cols):
    for col in cols:
        if col not in metadata.columns:
            metadata[col] = pd.Series(index=metadata.index, dtype=object)
        elif metadata[col].dtype != object:
            metadata[col] = metadata[col].astype(object)

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
    all_organoids=False,
    progress_cb=None,
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
    progress_cb : callable, optional
        Called as ``progress_cb(current, total, label)`` from inside the
        per-sample loop so a GUI can drive a progress bar.  ``None`` is
        a no-op (default).  Notebook / queue callers omit this kwarg and
        the function behaves identically to its pre-existing version.
    **kwargs : dict
    """
    if all_organoids:
        from behav3d.preprocessing.tracking.propagation_tracking_all_organoids import (
            run_propagation_tracking_all_organoids,
        )

        return run_propagation_tracking_all_organoids(
            metadata=metadata,
            output_dir=output_dir,
            overwrite=overwrite,
            dilation_nr_pixels=dilation_nr_pixels,
            segment_size_min=segment_size_min,
            progress_cb=progress_cb,
            **kwargs,
        )

    total_samples = len(metadata)
    for i, (idx, sample) in enumerate(metadata.iterrows()):
        sample_name=sample['sample_name']
        if progress_cb is not None:
            try:
                progress_cb(i, total_samples, f"{cell_type} / {sample_name}")
            except Exception:
                pass
        print(f"Tracking sample: {sample_name}")

        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]
        
        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        
        segments_col, segments_path = _resolve_segments_column(sample, cell_type)
        
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
        
        img_col, csv_col = _resolve_tracking_output_columns(cell_type, segments_col)
        _ensure_metadata_output_columns(metadata, img_col, csv_col)
        metadata.at[idx, img_col] = str(tracked_img_outpath)
        metadata.at[idx, csv_col] = str(tracked_csv_outpath)

    if progress_cb is not None:
        try:
            progress_cb(total_samples, total_samples, f"{cell_type} done")
        except Exception:
            pass
    return metadata
