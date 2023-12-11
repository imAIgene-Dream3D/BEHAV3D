### TODO Perhaps set a "static" speed based on quantiles of mean_speed to define static and actively interacting cells
### However this excludes a lot of cells to actively interact, now we take the fastest of the contacting T cells as actively interacting

"""
This script calculates the features of tracks for BEHAV3D analysis.

-------------------------------------
--------------- INPUT ---------------
-------------------------------------


.csv containing the following columns:
- TrackID       (The ID of the track a segments belongs to)
- SegmentID     (The unique ID of the segment of a specific timepoint)
- position_t    (The timepoint of the segment)
- position_z
- position_y
- position_z

-------------------------------------
--------------- OUTPUT --------------
-------------------------------------

# Features of a track at each timepoint per sample in the metadata csv (.csv)
- See "FEATURES TRACKS"

# Combined summarized features for each track for all samples in metadata csv (.csv)
- 

-------------------------------------
---------- FEATURES TRACKS ---------- 
-------------------------------------
- organoid_contact
- organoid_contact_pixels
- touching_organoids
- tcell_contact
- tcell_contact_pixels
- touching_tcells
- active_tcell_interaction
- mean_dead_dye
- displacement
- cumulative_displacement
- displacement_from_origin
- mean_square_displacement
- speed
- mean_speed
- interpolated
- time

### organoid_contact
- True/False
Per segment, creates a zyx cutout of the segment with a range of pixels around 
the segment border and calculates a distance transform from the T cell border. 
Any other segment inside a range specified by "contact_threshold" counts as a contacting organoid

### organoid_contact_pixels
- True/False
Same as "organoid_contact", but a contact is now specified as anotehr segment touching the segment
based on pixels without taking pixel_distances into account.

### touching_organoids
- String (List separated by ",")
These are the TrackIDs of touching organoids, separated by ",". 
NaN if none are touching

### tcell_contact
- True/False
Same as "organoid_contact", but now checks for touching T cells

### tcell_contact_pixels
- True/False
Same as "organoid_contact_pixels", but now checks for touching T cells

### touching_tcells
- String (List separated by ",")
Same as "touching_organoids", but now checks for touching T cells

### active_tcell_interaction
- True/False
For cell interaction we can consider the following:
When two cells interact it is often the one cell moves and interacts with another one that is static
In this case one might consider that only one motile cell is actively interacting and the other cells
are just passively interacting. To determine when a cell is actively interacting we measure for each 
cell what was its mean_speed over the last "rolling_meanspeed_window" timepoints. We then rank the 
T cells that are touching on mean_speed and only the one with the highest speed is labeled as an 
active tcellinteraction.

### mean_dead_dye
- Float
The mean intensity of the dead dye inside of each segment per timepoint. Calculated based on the 
channel of the dead dye. Supplied as an index of the channel image supplied as "dead_dye_channel" 

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
calculating the "active_tcell_interaction"

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
from scipy.ndimage import distance_transform_edt
from scipy.spatial.distance import cdist
from skimage.measure import regionprops_table
import argparse
import yaml
from tifffile import imread
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from plotnine import (
    ggplot, 
    aes, 
    geom_bar, 
    geom_violin, 
    geom_jitter, 
    facet_grid, 
    labs, 
    theme, 
    element_text, 
    theme_minimal,
    theme_bw,
    scale_x_continuous
)
import math
import time
from behav3d import format_time
 
def run_behav3d_feature_extraction(config, metadata, cell_type="tcells"):
   
    df_all_tracks=calculate_track_features(config, metadata, cell_type=cell_type)
    
    # Calculate features for each timepoint in a track
    df_all_tracks_filt = filter_tracks(df_all_tracks, config, metadata)
    
    # Calculate features for each timepoint in a track
    df_summarized_tracks=summarize_track_features(df_all_tracks_filt, config, metadata) 
         
def calculate_track_features(
    metadata, 
    config=None,
    output_dir=None,
    cell_type="tcells",
    imaris=False,
    ):
    """
    This code calculates the various features for each timepoint in a track for each 
    separateexperiment
    
    This codes works with either:
    - The generated segments and tracks from the BEHAV3D preprocessing modules
    - Imaris extracted statistics
    
    Output:
    - A .csv file containing all timepoints of all TrackIDs and their time-related 
      features
    """

    assert config is not None or output_dir is not None, "Either 'config' or 'output_dir' must be supplied"
    
    if output_dir is None:
        output_dir = config['output_dir']
        if "imaris" in config.keys():
            imaris = config["imaris"]
        else:
            imaris = False
            
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)
        
    rolling_meanspeed_window=10
    
    df_all_tracks=pd.DataFrame()
    
    for _, sample_metadata in metadata.iterrows():
        print(f"--------------- Processing {cell_type}: {sample_metadata['sample_name']} ---------------")
        start_time = time.time()
        sample_name = sample_metadata['sample_name']
        element_size_x=sample_metadata['pixel_distance_xy']
        element_size_y=sample_metadata['pixel_distance_xy'] 
        element_size_z=sample_metadata['pixel_distance_z']
        distance_unit=sample_metadata['distance_unit']
        dead_dye_threshold=sample_metadata['dead_dye_threshold']
        contact_threshold = sample_metadata["contact_threshold"]
        
        img_outdir = Path(output_dir, "images", sample_name)
        track_outdir = Path(output_dir, "trackdata", sample_name)
        analysis_outdir = Path(output_dir, "analysis", cell_type)
        
        if not img_outdir.exists():
            img_outdir.mkdir(parents=True)
        if not track_outdir.exists():
            track_outdir.mkdir(parents=True)
        if not analysis_outdir.exists():
            analysis_outdir.mkdir(parents=True)
        # Sometimes excel saves the encoding for µm differently, the following lines converts
        # other variants of µm to ones comparable in this code
        # The two written "μm" have different formatting
        if distance_unit=='_m':
            distance_unit ='μm'
        if distance_unit=='µm':
            distance_unit="μm"
        time_interval = sample_metadata['time_interval']
        time_unit = sample_metadata['time_unit']

        print("###### Running track feature calculation")
        
        if imaris:
            print("Performing feature calculation from Imaris processing..")
        
        # ### Set specific parameters based on Imaris or BEHAV3D processing
        # if imaris:
        #     tcell_contact_threshold = sample_metadata["tcell_contact_threshold"]
        #     organoid_contact_threshold = sample_metadata["organoid_contact_threshold"]
        # else:
        raw_image_path = sample_metadata['raw_image_path']
        organoid_segments_path = sample_metadata['organoid_segments_path']
        tcell_segments_path = sample_metadata['tcell_segments_path']
        
        red_lym_channel=sample_metadata['dead_dye_channel']
        contact_threshold = sample_metadata["contact_threshold"]

        print("- Loading in tracks csv...")
        ### Load in the specified track csv
        # if sample_metadata["tcell_tracks_csv"]=="" or sample_metadata["tcell_tracks_csv"]==None or math.isnan(sample_metadata["tcell_tracks_csv"]):
        #     df_tracks_path = Path(track_outdir, f"{sample_name}_{cell_type}_tracks.csv")
        # else:
        df_tracks_path = sample_metadata["tcell_tracks_csv"]
        df_tracks=pd.read_csv(df_tracks_path, sep=",")
        
        ### Calculate organoid distance, t cell distance and dead dye mean form Imaris or BEHAV3D processing
        if imaris and cell_type=="tcells":
            print("- Calculating contact with organoids and other T cells... (From Imaris)")
            # Threshold the distance to organoid based on the supplied "organoid_contact_threshold"
            # This distance is calculated before in Imaris and supplied as a separate channel and
            # Extracted as a statistic (...Intensity_Min_Ch<#>_img=<#>.csv)
            df_tracks["organoid_contact"]=df_tracks["organoid_distance"]<=contact_threshold
            df_tracks["tcell_contact"]=df_tracks["tcell_distance"]<=contact_threshold
        else:
            print("- Calculating contact with organoids and other T cells...")
            print(f"Using a contact threshold of {contact_threshold}{distance_unit}")
            
            # Load in the images containing the organoid segments and T cell segments
            # organoid_segments_path=Path(img_outdir, f"{sample_name}_organoids_tracked.tiff")
            organoid_segments=imread(organoid_segments_path)
            
            # tcell_segments_path=Path(img_outdir, f"{sample_name}_tcells_tracked.tiff")
            tcell_segments=imread(tcell_segments_path)
            
            # Calculate the contact od each T cell with an organoid or a T cell
            # Explanation on how in the function itself
            df_contacts=calculate_organoid_and_tcell_contact(
                tcell_segments=tcell_segments,
                organoid_segments=organoid_segments,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                contact_threshold=contact_threshold,
                calculate_from=cell_type
            )
            df_tracks = pd.merge(df_tracks, df_contacts, how="left")
            
            print("- Calculating channel and especially death dye intensities...")
            intensity_image = imread(raw_image_path)
            df_intensity=calculate_segment_intensity(
                tcell_segments=tcell_segments,
                intensity_image=intensity_image
            )
            if red_lym_channel is not None:
                df_intensity = df_intensity.rename(columns={f"mean_intensity_ch{red_lym_channel}":"mean_dead_dye"})
            df_tracks = pd.merge(df_tracks, df_intensity, how="left")
        
        # As sometimes 1 or several timepoints are missing in a track, interpolate these missing rows
        # Values are interpolated linearly, forward filled or left blank based on the column
        # More explanation within the function
        print("- Interpolating missing timepoints based on time interval")
        df_tracks= interpolate_missing_positions(
            df_tracks
        )
        
        print(f"- Calculating cell death based on defined dead_dye_threshold {dead_dye_threshold}")
        df_tracks["dead"] = False
        
        # For any cell crossing the dead_dye_threshold, set the cell to dead. Any timepoint after this timepoint are
        # Also set to dead, even if the mean dead dye intensity goes under the threshold again
        for track_id in df_tracks["TrackID"].unique():
            track_df = df_tracks[df_tracks["TrackID"] == track_id]
            track_df_reset = track_df.reset_index(drop=True)
            threshold_indices = track_df_reset.reset_index(drop=True)[track_df_reset["mean_dead_dye"] >= dead_dye_threshold].index
            
            if not threshold_indices.empty:
                first_threshold_index = threshold_indices.min()
                df_tracks.loc[track_df.index[first_threshold_index:], "dead"] = True
        
        print("- Converting distance and time unit to default um and hours...")
        
        # Converting the time and distance values to a default unit to allow comparison 
        # with differently provided units (defaults to µm and hours)
        def convert_distance(distance, distance_unit):
            distance_conversions={
                "nm":1000,
                "μm":1,
                "um":1,
                "mm":0.001
            }
            assert distance_unit in list(distance_conversions.keys()), f"distance unit needs to be one of: {list(distance_conversions.keys())}, is {distance_unit}"
            distance = distance/distance_conversions[distance_unit]
            return(distance)
        
        def convert_time(time_interval, time_unit):
            time_conversions={
                "s": 3600,
                "m": 60,
                "h": 1
            }
            assert time_unit in time_conversions.keys(), f"time unit needs to be one of: {time_conversions.keys()}, is {time_unit}"
            time_interval = time_interval/time_conversions[time_unit]
            return(time_interval)
        
        df_tracks["position_z"]=df_tracks["position_z"].apply(convert_distance, args=(distance_unit,))
        df_tracks["position_y"]=df_tracks["position_y"].apply(convert_distance, args=(distance_unit,))
        df_tracks["position_x"]=df_tracks["position_x"].apply(convert_distance, args=(distance_unit,))
        
        # Calculate relative time, where each track begins at timepoint 1
        def calculate_relative_time(group):
            min_position = group['position_t'].min()
            group['relative_time'] = group['position_t'].sub(min_position).add(1)
            return group

        df_tracks = df_tracks.groupby('TrackID').apply(calculate_relative_time).reset_index(drop=True)

        df_tracks["time"]=df_tracks["position_t"]*time_interval
        df_tracks["time"]=df_tracks["time"].apply(convert_time, args=(time_unit))
        time_interval=convert_time(time_interval, time_unit)
        contact_threshold = convert_distance(contact_threshold, distance_unit)

        element_size_x = convert_distance(element_size_x, distance_unit)
        element_size_y = convert_distance(element_size_y, distance_unit)
        element_size_z = convert_distance(element_size_z, distance_unit)
        
        # Make default distance unit um over μm, as encoding of μm between python 
        # and e.g Excel can lead to formatting errors
        
        df_tracks["distance_unit"] = "um"
        df_tracks["time_unit"] = "h"

        # Calculate various movement features such as speed and mean square displacement of the tracks  
        print("- Calculating movement features...")
        df_tracks=calculate_movement_features(
            df_tracks,
            time_interval = time_interval,
            rolling_meanspeed_window=rolling_meanspeed_window
            )
        df_tracks = df_tracks.sort_values(['TrackID', 'position_t'])
        
        if cell_type=="tcells" and not imaris:
            print("- Determining active contact of T cells")
            # Determining if a T cell is actively interacting with another T cell based on speed
            # More explanation at the top of this code
            df_tracks['list_touching_tcells'] = df_tracks['touching_tcells'].apply(lambda x: [] if pd.isna(x) else list(map(int, x.split(','))))
            active_interaction = []
            for _, row in df_tracks.iterrows():
                if row['tcell_contact']:
                    max_mean_speed = max(row['mean_speed'], df_tracks.loc[df_tracks['SegmentID'].isin(row['list_touching_tcells']), 'mean_speed'].max())
                    active_interaction.append(row['mean_speed'] == max_mean_speed)
                else:
                    active_interaction.append(False)
            df_tracks=df_tracks.drop('list_touching_tcells', axis=1)      
            df_tracks['active_tcell_contact'] = active_interaction
        
        print("- Perform z-normalization on selected feature columns")
        df_tracks=normalize_track_features(
            df_tracks, 
            columns=[
                "mean_square_displacement",
                "speed",
                "mean_dead_dye"
            ]
        )
        
        df_tracks=df_tracks.sort_values(by=["TrackID", "relative_time"])
        
        tracks_out_path = Path(track_outdir, f"{sample_name}_{cell_type}_track_features.csv")
        print(f"- Writing output to {tracks_out_path}")
        df_tracks.to_csv(tracks_out_path, sep=",", index=False)
        
        # Adding a sample name for later combination of multiple track experiments
        df_tracks['sample_name']=sample_name
        df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
        df_all_tracks = pd.concat([df_all_tracks, df_tracks])
        
        end_time = time.time()
        h,m,s = format_time(start_time, end_time)
        print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    
    all_tracks_out_path = Path(feature_outdir, f"BEHAV3D_combined_track_features.csv")
    df_all_tracks.to_csv(all_tracks_out_path, index=False) 
    return(df_all_tracks)

def filter_tracks(
    df_all_tracks,
    metadata,
    config=None,
    output_dir=None,
    tcell_exp_duration=None,
    tcell_min_track_length=None,
    tcell_max_track_length=None,
    cell_type="tcells"
    ):
    """
    This code filters tracks based on supplied parameters in the config.yml
    
    Filtering is based on:
    - Maximum experiment length (tcell_exp_duration)    
    - Minimum track length (tcell_min_track_length)
    - Tracks starting at timepoint 1 with a dead dye mean over the dead_dye_threshold (dead_dye_threshold)
    
    Additonally, all tracks are cut down to:
    - Maximum track length (tcell_max_track_length)
    
    Output:
    - A .csv file containing filtered tracks from all experiments
    """
    
    assert config is not None or all(
        [output_dir, tcell_exp_duration, tcell_min_track_length, tcell_max_track_length]
    ), "Either 'config' or 'output_dir, tcell_exp_duration, tcell_min_track_length and tcell_max_track_length' parameters must be supplied"
    
    start_time = time.time()
    
    print(f"--------------- Filtering tracks ---------------")
    
    if not all([output_dir, tcell_exp_duration, tcell_min_track_length, tcell_max_track_length]):
        output_dir = config['output_dir']
        tcell_exp_duration = config['tcell_exp_duration']
        tcell_min_track_length = config['tcell_min_track_length']
        tcell_max_track_length = config['tcell_max_track_length']
    
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    qc_outdir = Path(analysis_outdir, "quality_control")
    
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)
    if not qc_outdir.exists():
        qc_outdir.mkdir(parents=True)
     
    group_cols = ['TrackID', 'sample_name', 'organoid_line', 'tcell_line', 'exp_nr', 'well']
    df_all_tracks_filt = pd.merge(df_all_tracks, metadata, how="left", on="sample_name")

    # Function to count the number of unique tracks in the DataFrame
    def count_tracks(df_all_tracks, col_name="nr_tracks", df_track_counts=None):
        nr_tracks=df_all_tracks.groupby([
            'sample_name', 'organoid_line', 'tcell_line', 'exp_nr', 'well']
            ).agg(nr_tracks=pd.NamedAgg(column='TrackID', aggfunc='nunique')).reset_index()
        nr_tracks=nr_tracks.rename(columns={"nr_tracks":col_name})
        if df_track_counts is None:
            return(nr_tracks)
        else:
            return(pd.merge(df_track_counts, nr_tracks, how="left"))
    
    # Counting the nr of tracks before filtering
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_before_filtering")
    
    # Filtering the tracks based on the total experimental duration
    # Any timepoint after this will be filtered out 
    df_all_tracks_filt=df_all_tracks_filt.drop(df_all_tracks_filt[df_all_tracks_filt["position_t"]>tcell_exp_duration].index)
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_exp_duration", df_track_counts=df_track_counts)

    # Filtering out tracks under specific track length and cutting them down to specified max track length
    df_all_tracks_filt=df_all_tracks_filt.groupby(group_cols).filter(lambda group: len(group) >= tcell_min_track_length).reset_index(drop=True)
    df_all_tracks_filt=df_all_tracks_filt.groupby(group_cols).apply(lambda group: group.iloc[:tcell_max_track_length]).reset_index(drop=True)
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_min_track_length", df_track_counts=df_track_counts)

    # Plot the number of cells having contact with another T cell and Organoid for analysis
    # of the set contact_threshold
    plot_tcell_touching = plot_touching_nontouching_distribution(
        df_all_tracks_filt, 
        contact_column="tcell_contact"
        )
    plot_org_touching = plot_touching_nontouching_distribution(
        df_all_tracks_filt, 
        contact_column="organoid_contact"
        )
    
    # Plot the distribution of dead dye intensity of all timepoints and at timepoint 1
    # Can be used to aid in the choice of dead_dye_threshold
    dead_dye_distr=plot_dead_dye_distribution(df_all_tracks_filt)
    dead_dye_distr_t0=plot_dead_dye_distribution(df_all_tracks_filt[df_all_tracks_filt["relative_time"]==1])
    
    # Filter out all T cells that are dead based on the threshold at the first timepoint of a track
    dead_t0 = df_all_tracks_filt[
        (df_all_tracks_filt["relative_time"]==1) & 
        (df_all_tracks_filt["dead"])
        ][["TrackID","sample_name"]]
    df_all_tracks_filt=df_all_tracks_filt[~df_all_tracks_filt.set_index(['TrackID', 'sample_name']).index.isin(dead_t0.set_index(['TrackID', 'sample_name']).index)]   
    df_track_counts=count_tracks(df_all_tracks_filt, col_name="nr_tracks_dead_t1", df_track_counts=df_track_counts)

    # TODO create filter plot to see how much is filtered
    
    # Write the filtered tracks to a .csv
    filt_tracks_out_path = Path(feature_outdir, f"BEHAV3D_combined_track_features_filtered.csv")
    print(f"- Writing filtered tracks to {filt_tracks_out_path}")
    df_all_tracks_filt.to_csv(filt_tracks_out_path, sep=",", index=False)
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_all_tracks_filt)

def normalize_track_features(
    df_tracks,
    columns = [
        "mean_square_displacement",
        "speed",
        "mean_dead_dye"
    ]
    ):
    
    for column_name in columns:
        df_tracks.loc[:, f'z_{column_name}'] = df_tracks[column_name].transform(lambda x: (x - x.mean()) / x.std())
    return (df_tracks)

def summarize_track_features(
    df_tracks,
    config=None,
    output_dir=None,
    imaris=False,
    cell_type="tcells"
    ):
    """
    This code calculates summarized features (e.g. mean speed of the whole track) 
    for each TrackID for every experiment specified in the provided metadata.csv
    
    Output:
    - A .csv file containing all tracks from all experiments with their track-summarized features
    """
    
    assert config is not None or output_dir is not None, "Either 'config' or 'output_dir' must be supplied"
    
    start_time = time.time()
            
    print(f"--------------- Summarizing track features ---------------")
    
    if output_dir is None:
        output_dir = config['output_dir']
        if "imaris" in config.keys():
            imaris = config["imaris"]
        else:
            imaris = False
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    feature_outdir = Path(analysis_outdir, "track_features")
    qc_outdir = Path(analysis_outdir, "quality_control")
    
    if not analysis_outdir.exists():
        analysis_outdir.mkdir(parents=True)
    if not feature_outdir.exists():
        feature_outdir.mkdir(parents=True)
    if not qc_outdir.exists():
        qc_outdir.mkdir(parents=True)

    # Calculate mean values of track features over the whole track
    grouped_df_tracks=df_tracks.groupby(['sample_name','TrackID'])
    df_summarized_tracks = grouped_df_tracks.size().reset_index(name="track_length")
    df_summarized_tracks['mean_dead_dye'] = grouped_df_tracks['mean_dead_dye'].mean().reset_index()["mean_dead_dye"]
    df_summarized_tracks['mean_MSD'] =  grouped_df_tracks['mean_square_displacement'].mean().reset_index()['mean_square_displacement']
    df_summarized_tracks['mean_speed'] =  grouped_df_tracks['speed'].mean().reset_index()['speed']
    df_summarized_tracks['mean_organoid_contact'] =  grouped_df_tracks['organoid_contact'].mean().reset_index()['organoid_contact']
    df_summarized_tracks['mean_tcell_contact'] =  grouped_df_tracks['tcell_contact'].mean().reset_index()['tcell_contact']
    df_summarized_tracks['mean_displacement'] =  grouped_df_tracks['displacement'].mean().reset_index()['displacement']
    df_summarized_tracks['dies'] =  grouped_df_tracks['dead'].any().reset_index()["dead"]
    
    # For some values, take the maximum of the track such as "displacement_from_origin"
    df_summarized_tracks['displacement_from_origin'] =  grouped_df_tracks['displacement_from_origin'].last().reset_index()['displacement_from_origin']
    df_summarized_tracks['cumulative_displacement'] =  grouped_df_tracks['cumulative_displacement'].last().reset_index()['cumulative_displacement']
    
    if not imaris:
        # Calculate for the contact that occurs, what percentage has been active contact
        # As it only takes points of contact, this can mean the mean contact is 1% (0.01)
        # while the active contact can then still be 100% (1.0)
        def calculate_active_contact_when_contact(group):
            if group['tcell_contact'].any():
                return group[group['tcell_contact']]['active_tcell_contact'].mean()
            else:
                return 0
        df_summarized_tracks['active_tcell_contact'] = grouped_df_tracks.apply(calculate_active_contact_when_contact).reset_index(drop=True)

    df_trackinfo = df_tracks[['TrackID', 'sample_name','well', 'exp_nr', 'organoid_line', 'tcell_line']].drop_duplicates()
    df_summarized_tracks = pd.merge(df_trackinfo, df_summarized_tracks, how="left")
    # Write the summarized features to a .csv
    summ_tracks_out_path = Path(feature_outdir, f"BEHAV3D_combined_track_features_summarized.csv")
    print(f"- Writing summarized tracks to {summ_tracks_out_path}")
    df_summarized_tracks.to_csv(summ_tracks_out_path, sep=",", index=False)
    
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return(df_summarized_tracks)

def plot_dead_dye_distribution(
    df_tracks
    ):
    """
    Create a violin plot with an underlying scatterplot that provides an
    overview of the mean dead dye intensity per segment at the first timepoint
    in each experiment
    """
    df_time1 = df_tracks[df_tracks["relative_time"]==1]
    figure = (
        ggplot(df_time1, aes(x='time', y='mean_dead_dye')) +
        geom_jitter() +
        geom_violin(aes(fill='sample_name')) +
        facet_grid('~sample_name') +
        scale_x_continuous(breaks=[]) +
        theme_bw()
    )   
    return(figure)

def plot_touching_nontouching_distribution(
    df_tracks,
    contact_column='organoid_contact',
    ):
    """
    Create a barplot that provides an overview of how many cells make contact with
    other organoids/cells
    """
    figure = (
        ggplot(df_tracks, aes(x=contact_column, fill='tcell_line')) +
        geom_bar(position='dodge') +
        labs(x='Contact', y='Count', title='Touching vs. Non_touching organoids', fill='T cell line') +
        facet_grid('sample_name ~ organoid_line', scales='free') +
        theme_minimal() +
        theme(
            strip_text_x=element_text(angle=0, ha='center'),
            strip_text_y=element_text(angle=0, va='center')# Rotate row facet titles horizontally
        )   
    )
    # Display the combined plot
    # print(figure)
    return(figure)
    
def interpolate_missing_positions(
    df_tracks,
    cols_to_copy=[
        "organoid_contact",
        "organoid_contact_pixels",
        "touching_organoids",
        "tcell_contact",
        "tcell_contact_pixels",
        "touching_tcells"
        ],
    cols_to_interpolate=[
        "position_t", 
        "position_z", 
        "position_y", 
        "position_x",
        "mean_dead_dye"
        ]
    ):
    """
    As not every track has a segment at every timepoint, interpolate the missing values of
    the missing timepoints
    
    It interpolates various columns in different ways:
    -   Interpolates the numerical columns of [cols_to_interpolate] such as speed using linear
        interpolation
    -   Copies the columns of [cols_to_copy] using a forward fill from the last non-interpolated
        row of each TrackID
    -   Puts None in any column not specified, such as SegmentID, as no actual segment exists
    """
     # Interpolate missing timepoints so each calculation takes the same intervals
    grouped_df = df_tracks.groupby('TrackID')
    def interpolate_group(group, cols_to_interpolate, cols_to_copy):
        # group=group.set_index('time', drop=False)
        group["interpolated"]=False
        min_time = group['position_t'].min()
        max_time = group['position_t'].max()
        all_times = [time for time in list(np.arange(min_time, max_time, 1))]+[max_time]
        missing_times = pd.DataFrame({'position_t':[x for x in all_times if x not in group['position_t'].tolist()]})
        missing_times["TrackID"]=group['TrackID'].unique()[0]
        
        df_interpolated=group.copy()
        df_interpolated = pd.concat([group, missing_times], ignore_index=True).sort_values(by="position_t")

        for col in cols_to_interpolate:
            df_interpolated[col] = np.interp(df_interpolated['position_t'], group['position_t'], group[col])
        
        cols_to_copy = [col for col in cols_to_copy if col in df_interpolated.columns]
        # Apply forward-fill only to the newly added rows
        newly_added_rows = df_interpolated.loc[df_interpolated['interpolated'].isna()]
        for col in cols_to_copy:
            for idx in newly_added_rows.index:
                previous_idx = idx - 1
                while previous_idx >= 0 and df_interpolated.loc[previous_idx, 'interpolated']:
                    previous_idx -= 1
                if previous_idx >= 0:
                    df_interpolated.loc[idx, col] = df_interpolated.loc[previous_idx, col]
        
        df_interpolated["interpolated"] = df_interpolated["interpolated"].fillna(True) 
        assert(len(all_times)==len(df_interpolated)), f"Length of expected nr of timepoints ({len(all_times)}) is not the same as resulting timepoints ({df_interpolated})"
        return df_interpolated
    df_interpolated = pd.concat([interpolate_group(group, cols_to_interpolate, cols_to_copy) for _, group in grouped_df])
    df_interpolated = df_interpolated.reset_index(drop=True)
    return(df_interpolated)
        
def calculate_movement_features(
    df_tracks, 
    time_interval,
    rolling_meanspeed_window=10
    ):
    """
    Calculates various movement features for each timepoint of a track
    """
    ## Convert the coordinates to time series
    
    #TODO Angleness/directionality: How much does it move in a single direction
    # Calculate by calcualting standard deviation of angle changes ?
    
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

        # combine
        df_computed = pd.DataFrame({
            'displacement': displacement, 
            'cumulative_displacement': cumulative_displacement, 
            'displacement_from_origin': displacement_from_origin, 
            'mean_square_displacement':mean_square_displacement
            })
        df_result= pd.concat([df_track_pos,df_computed], axis=1)
        df_result = pd.concat([df_result, df_track[["position_t", "SegmentID", "TrackID"]]], axis=1)

        df_result['speed'] = df_result["displacement"]/time_interval
        # Calculate the mean speed (default um/h) over the last {rolling_meanspeed_window} timepoints
        df_result['mean_speed'] = df_result.groupby('TrackID')['speed'].apply(lambda x: x.iloc[1:].rolling(window=rolling_meanspeed_window, min_periods=1).mean()).reset_index(0, drop=True)
        df_result['mean_speed'] = df_result['mean_speed'].fillna(0)
   
        df_tracks_processed.append(df_result)
    df_tracks_processed = pd.concat(df_tracks_processed)
    df_tracks_processed=pd.merge(df_tracks, df_tracks_processed, how="left")
    return(df_tracks_processed)

def calculate_organoid_distance(
    tcell_segments, 
    organoid_segments, 
    element_size_x, 
    element_size_y, 
    element_size_z
    ):
    df_dist_organoid = []
    for t, tcell_stack in enumerate(tcell_segments):
        org_stack = organoid_segments[t,:,:,:]
        mask_org= np.ma.masked_where(org_stack==0, org_stack)
        dist_org=distance_transform_edt(mask_org.mask)
        real_dist_org=distance_transform_edt(
            mask_org.mask,
            sampling=[element_size_z, element_size_y, element_size_x]
            )
        properties_pix=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=dist_org, properties=['label', 'intensity_min']))
        properties_pix=properties_pix.rename(columns={"intensity_min":"pix_distance_organoids"})
        properties_real=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=real_dist_org, properties=['label', 'intensity_min']))
        properties_real=properties_real.rename(columns={"intensity_min":"real_distance_organoids"})
        properties=pd.merge(properties_pix,properties_real, how="left")
        properties["position_t"]=t
        df_dist_organoid.append(properties)
    df_dist_organoid = pd.concat(df_dist_organoid)
    df_dist_organoid=df_dist_organoid.rename(columns={"label":"TrackID"})
    df_dist_organoid["pix_organoid_contact"] =  df_dist_organoid["pix_distance_organoids"] <= 1.73
    return(df_dist_organoid)

def calculate_organoid_and_tcell_contact(
    tcell_segments,
    organoid_segments,
    element_size_x, 
    element_size_y, 
    element_size_z,
    contact_threshold,
    calculate_from="tcells"
    ):
    """
    Calculates contact with organoids by looping through
    the segments and cutting out a small area around it.
    It then calculates the euclidian distance of every pixel
    outside the segment and sees if any other segments are in this
    area.
    
    organoid_contact: 
    Based on the provided 'contact_threshold' value. Any other segment
    in this distance is seen as contacting
    
    organoid_contact_pixel:
    Based on direct pixel contact of one t cell and another, not influenced
    by the element_sizes
    """
    df_contacts= []    
    for t, tcell_stack in enumerate(tcell_segments):
        org_stack = organoid_segments[t,:,:,:]
        if calculate_from=="tcells":
            segments_stack=tcell_stack
        elif calculate_from=="organoids":
            segments_stack=org_stack
        else:
            print(f"calculate_from has to be either 'tcells' or 'organoids', is {calculate_from}")
            return
        segment_ids=np.unique(segments_stack)
        for segment_id in segment_ids:
            if segment_id==0:
                continue
            
            stack_max_z, stack_max_y, stack_max_x = segments_stack.shape
            seg_locs = np.argwhere(segments_stack==segment_id)
            min_z, min_y, min_x = seg_locs.min(axis=0)
            max_z, max_y, max_x = seg_locs.max(axis=0)
            
            z_ext = 2*math.ceil(contact_threshold / element_size_z)
            y_ext = 2*math.ceil(contact_threshold / element_size_y)
            x_ext = 2*math.ceil(contact_threshold / element_size_x)
            tcell_cutout = tcell_stack[
                max(0, min_z-z_ext):min(stack_max_z, max_z+z_ext+1),
                max(0, min_y-y_ext):min(stack_max_y, max_y+y_ext+1),
                max(0, min_x-x_ext):min(stack_max_x, max_x+x_ext+1),
                ]
            org_cutout = org_stack[
                max(0, min_z-z_ext):min(stack_max_z, max_z+z_ext+1),
                max(0, min_y-y_ext):min(stack_max_y, max_y+y_ext+1),
                max(0, min_x-x_ext):min(stack_max_x, max_x+x_ext+1),
                ]
            seg_cutout = segments_stack[
                max(0, min_z-z_ext):min(stack_max_z, max_z+z_ext+1),
                max(0, min_y-y_ext):min(stack_max_y, max_y+y_ext+1),
                max(0, min_x-x_ext):min(stack_max_x, max_x+x_ext+1),
                ]
            
            real_distances=distance_transform_edt(
                seg_cutout!=segment_id,
                sampling=[element_size_z, element_size_y, element_size_x]
                )
            pix_distances=distance_transform_edt(
                seg_cutout!=segment_id
                )
             
            organoid_contacts = [str(x) for x in np.unique(org_cutout[real_distances<=contact_threshold]) if x!=0]
            real_organoid_contact = len(organoid_contacts)>0
            pix_organoid_contacts = [str(x) for x in np.unique(org_cutout[pix_distances<= 1.73]) if x!=0]
            pix_organoid_contact = len(pix_organoid_contacts)>0
            
            tcell_contacts = [str(x) for x in np.unique(tcell_cutout[real_distances<=contact_threshold]) if x not in [0, segment_id]]
            real_tcell_contact = len(tcell_contacts)>0
            pix_tcell_contacts = [str(x) for x in np.unique(tcell_cutout[pix_distances<= 1.73]) if x not in [0, segment_id]]
            pix_tcell_contact = len(pix_tcell_contacts)>0
            
            if real_tcell_contact:
                touching_tcells = ",".join(tcell_contacts)
            else:
                touching_tcells=None
            if real_organoid_contact:
                touching_organoids = ",".join(organoid_contacts)
            else:
                touching_organoids = None
                
            df_contacts.append(pd.DataFrame([{
                'TrackID': segment_id, 
                'position_t': t,
                'organoid_contact': real_organoid_contact, 
                'organoid_contact_pixels': pix_organoid_contact,
                'touching_organoids':touching_organoids,
                'tcell_contact': real_tcell_contact, 
                'tcell_contact_pixels': pix_tcell_contact,
                'touching_tcells': touching_tcells
            }]))
    
    df_contacts=pd.concat(df_contacts)
    return(df_contacts)

def calculate_segment_intensity(tcell_segments, intensity_image, calculation="mean"):
    """
    Calculates the intensity of a specific marker features for each segment.
    The calculation can be the minimum, maximum, mean or median
    """
    assert calculation in ["min", "max", "mean", "median"]
    
    intensity_image=np.transpose(intensity_image, axes=[1,2,3,4,0])
    df_intensity = []
    for t, (tcell_stack, intensity_stack) in enumerate(zip(tcell_segments, intensity_image)):
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=intensity_stack, properties=['label', f'intensity_{calculation}']))
        properties["position_t"]=t
        df_intensity.append(properties)
    df_intensity = pd.concat(df_intensity)
    
    # Define a dictionary mapping old column names to new column names using regex
    column_mapping = {}
    for i in range(df_intensity.shape[1]):
        old_col_name = f'intensity_mean-{i}'
        new_col_name = f'mean_intensity_ch{i}'
        column_mapping[old_col_name] = new_col_name
    column_mapping["label"]="TrackID"
    df_intensity=df_intensity.rename(columns=column_mapping)
    return(df_intensity)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    metadata = pd.read_csv(config["metadata_csv_path"])
    run_behav3d_feature_extraction(config, metadata)