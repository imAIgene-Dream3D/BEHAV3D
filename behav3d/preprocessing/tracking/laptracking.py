import shutil
import pandas as pd
import numpy as np
from pathlib import Path
from skimage.measure import regionprops_table
from laptrack import LapTrack
from tqdm import tqdm
from behav3d.utils.fileio import load_image, append_to_zarr, get_filepath_stem

def laptrack_image(
    segments=None,
    segments_path=None,
    trackimg_outdir=None,
    trackcsv_outdir=None,
    basename=None,
    element_size_x=1,
    element_size_y=1,
    element_size_z=1,
    track_cost_cutoff=45**2, 
    gap_closing_cost_cutoff=60**2,
    gap_closing_max_frame_count=3,
    merging_cost_cutoff=False,
    splitting_cost_cutoff=False,
    return_trackimg=True
    ):
    
    assert segments is not None or segments_path is not None, "Either segments or segments_path must be provided"
    if segments is None:
        segments = load_image(segments_path)
    
    if trackimg_outdir is None:
        trackimg_outdir = Path(segments_path).parent
    if trackcsv_outdir is None:
        trackcsv_outdir = Path(segments_path).parent
    
    if basename is None and segments_path is not None:
        basename = get_filepath_stem(segments_path)
          
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
    
    df_tracks_outpath = Path(trackcsv_outdir, f"{basename}_tracks.csv")
    df_tracks.to_csv(df_tracks_outpath, sep=",", index=False)
    
    if return_trackimg:
        print("Overwriting the original segments IDs with the tracked IDs")
        # tracked_img = np.zeros_like(segments)
        tracked_img_outpath = Path(trackimg_outdir, f"{basename}_tracked.zarr")
        if tracked_img_outpath.exists():
            # Remove the file if it already exists
            tracked_img_outpath.unlink()
        for t, t_seg in tqdm(enumerate(segments), total=len(segments)):
            t_seg = np.asarray(t_seg)
            tracked_img = np.zeros_like(t_seg)
            t_df_tracks = df_tracks[df_tracks["position_t"]==t]
            for _, row in t_df_tracks.iterrows():
                # print(row["SegmentID"], row["TrackID"], (tracked_img[t]==row["SegmentID"]).any())
                tracked_img[t_seg==row["SegmentID"]] = row["TrackID"]

            tracked_img = np.expand_dims(tracked_img, axis=0)
            append_to_zarr(
                img=tracked_img, 
                outpath=tracked_img_outpath
                )
        shutil.make_archive(tracked_img_outpath, "zip", tracked_img_outpath)
        shutil.rmtree(tracked_img_outpath)
        tracked_img_outpath = Path(f"{tracked_img_outpath}.zip")
                
        return df_tracks_outpath, tracked_img_outpath
    else:
        return df_tracks_outpath
        
        # return df_tracks, tracked_img
    # else:
        # return df_tracks
   
def run_tcell_laptracking(
    metadata,
    output_dir,
    track_cost_cutoff=45**2, 
    gap_closing_cost_cutoff=60**2,
    gap_closing_max_frame_count=3,
    merging_cost_cutoff=False,
    splitting_cost_cutoff=False,
    return_trackimg=True,
    overwrite=False,
    **kwargs
    ):
     for idx, sample in metadata.iterrows():
        sample_name=sample['sample_name']
        print(f"Tracking sample: {sample_name}")
        
        tracked_img_outpath = Path(output_dir, "images", sample_name)
        tracked_csv_outpath = Path(output_dir, "trackdata", sample_name)
        if not tracked_img_outpath.exists():
            tracked_img_outpath.mkdir(parents=True)
        if not tracked_csv_outpath.exists():
            tracked_csv_outpath.mkdir(parents=True)
        
        element_size_x = sample["pixel_distance_xy"]
        element_size_y = sample["pixel_distance_xy"]
        element_size_z = sample["pixel_distance_z"]
        
        basename = f"{sample_name}_tcell"
        df_tracks_outpath, tracked_img_outpath = laptrack_image(
            segments_path=sample["tcell_segments_path"],
            basename=basename,
            trackimg_outdir=tracked_img_outpath,
            trackcsv_outdir=tracked_csv_outpath,
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
        
        metadata.at[idx, "tcell_segments_path"] = str(tracked_img_outpath)
        metadata.at[idx, "tcell_tracks_csv"] = str(tracked_img_outpath)
        
        return metadata