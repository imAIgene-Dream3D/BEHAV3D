import pandas as pd
import numpy as np
from pathlib import Path
from skimage.measure import regionprops_table
from laptrack import LapTrack
from tqdm import tqdm
from behav3d.io.images import load_image, append_to_zarr, get_filepath_stem
from behav3d.preprocessing.tracking import convert_segments_to_tracks

def laptrack_image(
    segments=None,
    segments_path=None,
    tracked_img_outpath=None,
    tracked_csv_outpath=None,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    track_cost_cutoff=25**2, 
    gap_closing_cost_cutoff=60**2,
    gap_closing_max_frame_count=3,
    merging_cost_cutoff=False,
    splitting_cost_cutoff=False,
    return_trackimg=True,
    
    ):
    
    assert segments is not None or segments_path is not None, "Either segments or segments_path must be provided"
    if segments_path is not None:
        segments_path = Path(segments_path)
    
    if segments is None:
        segments = load_image(segments_path)
    
    basename = get_filepath_stem(segments_path)
    if tracked_img_outpath is None:
        tracked_img_outpath = Path(segments_path.parent, f"{basename}_tracked.zarr")
        
    if tracked_csv_outpath is None:
        tracked_csv_outpath = Path(segments_path.parent, f"{basename}_tracks.csv")        
          
    df_centroids = []
    for t, tcell_stack in enumerate(segments):
        tcell_stack = np.asarray(tcell_stack)
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, properties=['label', f'centroid']))
        properties["position_t"]=t
        df_centroids.append(properties)
    df_centroids = pd.concat(df_centroids)
    df_centroids["position_z"]=df_centroids["centroid-0"]*element_size_z
    df_centroids["position_y"]=df_centroids["centroid-1"]*element_size_y
    df_centroids["position_x"]=df_centroids["centroid-2"]*element_size_x

    laptracker = LapTrack(
        track_cost_cutoff=track_cost_cutoff, 
        gap_closing_cost_cutoff=gap_closing_cost_cutoff,
        gap_closing_max_frame_count=gap_closing_max_frame_count,
        merging_cost_cutoff=merging_cost_cutoff,
        splitting_cost_cutoff=splitting_cost_cutoff,
        )
    
    df_tracks, df_splits, df_merges = laptracker.predict_dataframe(
        df_centroids.copy(),
        coordinate_cols=["position_z", "position_y", "position_x"],
        frame_col="position_t",
        only_coordinate_cols=False,
    )
    df_tracks = df_tracks.reset_index()
    
    df_tracks.rename(
            columns={
                    'track_id': 'TrackID',
                    'label': 'SegmentID',
                    'centroid-0': "pixel_position_z",
                    'centroid-1': "pixel_position_y",
                    'centroid-2': "pixel_position_x",
                }, 
            inplace=True
        )
    df_tracks["TrackID"]+=1
    # select only the columns we need
    df_tracks = df_tracks[["TrackID", "SegmentID", "position_t", "position_x", "position_y", "position_z", "pixel_position_x", "pixel_position_y", "pixel_position_z"]]
    
    df_tracks.to_csv(tracked_csv_outpath, sep=",", index=False)
    
    if return_trackimg:
        print("Overwriting the original segments IDs with the tracked IDs")
        tracked_img_outpath = convert_segments_to_tracks(
            df_tracks=df_tracks,
            segments=segments,
            outpath=tracked_img_outpath
        )
   
def run_laptracking(
    metadata,
    output_dir,
    cell_type,
    track_cost_cutoff=45**2, 
    gap_closing_cost_cutoff=60**2,
    gap_closing_max_frame_count=3,
    merging_cost_cutoff=False,
    splitting_cost_cutoff=False,
    return_trackimg=True,
    overwrite=False,
    **kwargs
    ):
    """Run LAP tracking on any cell type.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        DataFrame containing sample information
    output_dir : str or Path
        Root output directory
    cell_type : str
        Name of the cell type to track
    track_cost_cutoff : float
        Maximum cost for linking detections between consecutive frames
    gap_closing_cost_cutoff : float
        Maximum cost for closing gaps (linking across missing frames)
    gap_closing_max_frame_count : int
        Maximum number of frames a track can be missing
    merging_cost_cutoff : float or False
        Maximum cost for merging tracks (False to disable)
    splitting_cost_cutoff : float or False
        Maximum cost for splitting tracks (False to disable)
    return_trackimg : bool
        Whether to save tracked image
    overwrite : bool
        Whether to overwrite existing tracking results
    **kwargs : dict
     """
    for idx, sample in metadata.iterrows():
        sample_name=sample['sample_name']
        print(f"Tracking sample: {sample_name}")
        
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
            laptrack_image(
                segments_path=segments_path,
                tracked_img_outpath=tracked_img_outpath,
                tracked_csv_outpath=tracked_csv_outpath,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                track_cost_cutoff=track_cost_cutoff, 
                gap_closing_cost_cutoff=gap_closing_cost_cutoff,
                gap_closing_max_frame_count=gap_closing_max_frame_count,
                merging_cost_cutoff=merging_cost_cutoff,
                splitting_cost_cutoff=splitting_cost_cutoff,
                return_trackimg=return_trackimg
            )
        else:
            print("Tracking already exists... Provide overwrite=True to overwrite... Loading existing tracking data")
        
        # Update metadata with prefixed column names
        if segments_col is not None and segments_col.startswith(('or_', 'im_', 'ot_')):
            # Use the same prefix as the segments column
            prefix = segments_col.split('_')[0]
            metadata.at[idx, f"{prefix}_{cell_type}_tracks_image_path"] = str(tracked_img_outpath)
            metadata.at[idx, f"{prefix}_{cell_type}_tracks_csv_path"] = str(tracked_csv_outpath)
        else:
            # Fallback to old non-prefixed format
            metadata.at[idx, f"{cell_type}_tracks_image_path"] = str(tracked_img_outpath)
            metadata.at[idx, f"{cell_type}_tracks_csv_path"] = str(tracked_csv_outpath)
        
    return metadata