# https://imagej.net/plugins/trackmate/scripting/trackmate-detectors-trackers-keys
# TODO https://github.com/yfukai/laptrack
import sys
import os
import pandas as pd
import numpy as np
import psutil
from pathlib import Path
import yaml
from tifffile import imread, imwrite
from skimage.measure import regionprops_table
import argparse
import time
import math
import time
from laptrack import LapTrack

def run_lap_tracking(
    segments,
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
    
    
    df_centroids = []
    for t, tcell_stack in enumerate(segments):
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
                }, 
            inplace=True
        )
    
    # select only the columns we need
    df_tracks = df_tracks[["position_t", "position_x", "position_y", "position_z", "TrackID", "SegmentID"]]
    
    if return_trackimg:
        tracked_img = np.zeros_like(segments)
        
        for t, t_seg in enumerate(segments):
            print(t)
            t_df_tracks = df_tracks[df_tracks["position_t"]==t]
            for _, row in t_df_tracks.iterrows():
                # print(row["label"], row["track_id"], (tracked_img[t]==row["label"]).any())
                tracked_img[t][t_seg==row["SegmentID"]] = row["TrackID"]
        
        return df_tracks, tracked_img
    else:
        return df_tracks

def track_file(path):
    start_time = time.time() 
    
    path = Path(path)
    filename = Path(path).stem
    
    print(f"--------------- Tracking {filename} ---------------")
    segments = imread(path)
    df_tracks, tracked_img = run_lap_tracking(
        segments=segments,
        
    )
    
    tcell_tracked_csv_out_path= Path(path.parent, f"{filename}_tracked.csv")
    tcell_tracked_img_out_path= Path(path.parent, f"{filename}_tracked.tiff")
    
    # df_tracks.to_csv(tcell_tracked_csv_out_path)
    # imwrite(tcell_tracked_img_out_path, tracked_img)
    
    df_tracks = pd.read_csv(tcell_tracked_csv_out_path)
    df_tracks.rename(
            columns={
                    'position_t': 'Time',
                    'position_x': 'Position X',
                    'position_y': 'Position Y',
                    'position_z': 'Position Z',
                    'track_id': 'TrackID',
                }, 
            inplace=True
        )
    df_tracks["Length X"]=10
    df_tracks["Length Y"]=10
    df_tracks["Length Z"]=10
    df_tracks["parentID"]=0
    df_tracks = df_tracks[["Time", "Position X", "Position Y", "Position Z", "TrackID", "Length X", "Length Y", "Length Z", "parentID"]]
    
    df_path = Path("/Volumes/T7_Sam/SMI_SmartMicroscopy/TrackingTimeIntervalTest/lap_tracking", tcell_tracked_csv_out_path.name)
    df_tracks.to_csv(df_path, index=False)
    
    print(f"Tracking took {time.time() - start_time:.2f} s")

def main():
    path = "/Volumes/T7_Sam/SMI_SmartMicroscopy/BEHAV3Dsetup/segmentation/unet_predictions/FUNC_EV1_Exp080_Img001_13T-ctrl_upscaled_T0.tiff"
    track_file(path)
    segments_path = Path("/Volumes/T7_Sam/SMI_SmartMicroscopy/TrackingTimeIntervalTest/segmentation/downsampling")
    segmented_files = [x.resolve() for x in segments_path.glob(f'*.tif')]

    for interval_img_path in segmented_files:
        start_time = time.time() 
        
        filename = Path(interval_img_path).stem
        print(f"--------------- Tracking {filename} ---------------")
        # segments = imread(interval_img_path)
        # df_tracks, tracked_img = run_lap_tracking(
        #     segments=segments,
            
        # )
        
        tcell_tracked_csv_out_path= Path(interval_img_path.parent, f"{filename}_tracked.csv")
        tcell_tracked_img_out_path= Path(interval_img_path.parent, f"{filename}_tracked.tiff")
        
        # df_tracks.to_csv(tcell_tracked_csv_out_path)
        # imwrite(tcell_tracked_img_out_path, tracked_img)
        
        df_tracks = pd.read_csv(tcell_tracked_csv_out_path)
        df_tracks.rename(
                columns={
                        'position_t': 'Time',
                        'position_x': 'Position X',
                        'position_y': 'Position Y',
                        'position_z': 'Position Z',
                        'track_id': 'TrackID',
                    }, 
                inplace=True
            )
        df_tracks["Length X"]=10
        df_tracks["Length Y"]=10
        df_tracks["Length Z"]=10
        df_tracks["parentID"]=0
        df_tracks = df_tracks[["Time", "Position X", "Position Y", "Position Z", "TrackID", "Length X", "Length Y", "Length Z", "parentID"]]
        
        df_path = Path("/Volumes/T7_Sam/SMI_SmartMicroscopy/TrackingTimeIntervalTest/lap_tracking", tcell_tracked_csv_out_path.name)
        df_tracks.to_csv(df_path, index=False)
        
        print(f"Tracking took {time.time() - start_time:.2f} s")

