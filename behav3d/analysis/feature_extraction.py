### TODO Perhaps set a "static" speed based on quantiles of mean_speed to define static and actively interacting cells
### However this excludes a lot of cells to actively interact, now we take the fastest of the contacting T cells as actively interacting

### TODO WHy is tcell_cntact TRUE everywhere?
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
channel of the dead dye. Supplied as an index of the channel image supplied as "dead_channel" 

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

import warnings

from scipy.ndimage import distance_transform_edt
from scipy.spatial.distance import cdist
from scipy.spatial import ConvexHull

from skimage.measure import regionprops_table, marching_cubes, mesh_surface_area
import argparse
import yaml
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

import seaborn as sns

import ast
import math
import time
from behav3d.utils import get_current_time, format_time, convert_time, convert_distance
from behav3d.utils.fileio import load_image, convert_input_files_to_zarr
from tqdm import tqdm
from datetime import datetime

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

def run_feature_extraction(
    metadata, 
    config=None,
    output_dir=None,
    cell_type="tcell",
    imaris=False,
    rolling_meanspeed_window=10,
    overwrite=False,
    n_workers=1
    ):
    assert config is not None or output_dir is not None, "Either 'config' or 'output_dir' must be supplied"
    
    if output_dir is None:
        output_dir = config['output_dir']

    df_all_tracks=pd.DataFrame()
    
    for _, sample_metadata in metadata.iterrows():
        print(f"--------------- Processing {cell_type}: {sample_metadata['sample_name']} ---------------")
        start_time = time.time()

        sample_name = sample_metadata['sample_name']
        
        distance_unit=sample_metadata['distance_unit']
        dead_dye_threshold=sample_metadata['dead_dye_threshold']
        
        # Sometimes excel saves the encoding for µm differently, the following lines converts
        # other variants of µm to ones comparable in this code
        # The two written "μm" have different formatting
        if distance_unit=='_m':
            distance_unit ='μm'
        if distance_unit=='µm':
            distance_unit="μm"
        
        element_size_x=sample_metadata['pixel_distance_xy']
        element_size_y=sample_metadata['pixel_distance_xy'] 
        element_size_z=sample_metadata['pixel_distance_z']
        
        #Convert elekment size to um and hours, default settings for behav3d
        element_size_x = convert_distance(element_size_x, distance_unit)
        element_size_y = convert_distance(element_size_y, distance_unit)
        element_size_z = convert_distance(element_size_z, distance_unit)
        
        time_interval = sample_metadata['time_interval']
        time_unit = sample_metadata['time_unit']
        
        if imaris:
            tcell_contact_threshold = sample_metadata["tcell_contact_threshold"]
            organoid_contact_threshold = sample_metadata["organoid_contact_threshold"]
        else:
            contact_threshold = sample_metadata["contact_threshold"]
        dead_channel=sample_metadata['dead_channel']
        
        print("###### Running track feature calculation")
        img_outdir = Path(output_dir, "images", sample_name)
        if not img_outdir.exists():
            img_outdir.mkdir(parents=True)

        track_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        track_intermediate_outdir= Path(track_outdir, "intermediate_results")
        analysis_outdir = Path(output_dir, "analysis", cell_type)
        feature_outdir = Path(analysis_outdir, "track_features")
        
        raw_image_path = sample_metadata['raw_image_path']
        organoid_segments_path = sample_metadata['organoid_tracks_image_path']
        tcell_segments_path = sample_metadata['tcell_tracks_image_path']
        dead_mask_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
        
        print(f"{get_current_time()} - Converting all input files to .zarr for memory efficiency...")
        tcell_segments_path, organoid_segments_path, raw_image_path = convert_input_files_to_zarr(
            tcell_segments_path=tcell_segments_path,
            organoid_segments_path=organoid_segments_path,
            raw_image_path=raw_image_path,
            output_dir=img_outdir,
            overwrite=overwrite
        )

        if not track_outdir.exists():
            track_outdir.mkdir(parents=True)
        if not track_intermediate_outdir.exists():
            track_intermediate_outdir.mkdir()
        if not analysis_outdir.exists():
            analysis_outdir.mkdir(parents=True)
        if not feature_outdir.exists():
            feature_outdir.mkdir(parents=True)
        
        print(f"{get_current_time()} - Loading in tracks csv...")
        if cell_type=="tcell":
            df_tracks_path = sample_metadata["tcell_tracks_csv_path"]
        elif cell_type=="organoid":
            df_tracks_path = sample_metadata["organoid_tracks_csv_path"]
        else:
            raise ValueError(f"Unknown cell type: {cell_type}. Expected 'tcell' or 'organoid'.")
        
        df_tracks=pd.read_csv(df_tracks_path, sep=",")
        
        print(f"{get_current_time()} - Generalizing the units of position and time to um and hours")
        df_tracks = generalize_units_of_track_features(
            df_tracks=df_tracks,
            distance_unit=distance_unit,
            time_interval=time_interval,
            time_unit=time_unit
            )
        time_interval=convert_time(time_interval, time_unit)
        time_unit = "h"
        
        ### Calculate image based features (per timepoint)
        print(f"{get_current_time()} - Calculating single-timepoint image-based features...")
        df_intensity_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_intensity.csv")
        df_contacts_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_contact.csv")
        df_dead_mask_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_dead_mask.csv")
        df_morphology_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_morphology.csv")
        
        if imaris:
            df_tracks=calculate_imaris_track_features(
                df_tracks=df_tracks,
                organoid_contact_threshold=organoid_contact_threshold,
                tcell_contact_threshold=tcell_contact_threshold,
                distance_unit=distance_unit,
            )
        else:
            df_tracks=calculate_image_based_track_features(
                df_tracks=df_tracks,
                cell_type=cell_type,
                df_dead_mask_outpath=df_dead_mask_outpath,
                df_morphology_outpath=df_morphology_outpath,
                df_intensity_outpath=df_intensity_outpath,
                df_contacts_outpath=df_contacts_outpath,
                dead_channel=dead_channel,
                contact_threshold=contact_threshold,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                raw_image_path=raw_image_path,
                dead_mask_path=dead_mask_path,
                organoid_segments_path=organoid_segments_path,
                tcell_segments_path=tcell_segments_path,
                overwrite=overwrite,
                n_workers=n_workers
            )
        
        
        print(f"{get_current_time()} - Interpolating missing timepoints based on time interval")
        # As sometimes 1 or several timepoints are missing in a track, interpolate these missing rows
        # Values are interpolated linearly, forward filled or left blank based on the column
        # More explanation within the function
        df_tracks= interpolate_missing_positions(df_tracks)

        print(f"{get_current_time()} - Calculating movement features...")
        df_tracks=calculate_movement_features(
            df_tracks,
            time_interval = time_interval,
            rolling_meanspeed_window=rolling_meanspeed_window
            )
        df_tracks = df_tracks.sort_values(['TrackID', 'position_t'])
        
        if cell_type=="tcell":
            print(f"{get_current_time()} - Calculating T-cell specific features...")
            df_tracks = calculate_tcell_specific_track_features(
                df_tracks,
                dead_dye_threshold=dead_dye_threshold,
                )
            
        elif cell_type=="organoid":
            print(f"{get_current_time()} - Calculating organoid specific features...")
            df_tracks = calculate_organoid_specific_track_features(
                df_tracks,
                dead_dye_threshold=dead_dye_threshold,
                )
            #TODO Add organoid specific features
        else:
            print(f"{get_current_time()} - No cell type specified, skipping cell type specific features...")
        
        print(f"{get_current_time()} - Perform z-normalization on selected feature columns")
        df_tracks=normalize_track_features(
            df_tracks, 
            columns=[
                "mean_square_displacement",
                "speed",
                "mean_dead_dye"
            ]
        )
            
        tracks_out_path = Path(track_outdir, f"{sample_name}_{cell_type}_track_features.csv")
        print(f"{get_current_time()} - Writing output to {tracks_out_path}")
        df_tracks.to_csv(tracks_out_path, sep=",", index=False)
        
        # Adding a sample name for later combination of multiple track experiments
        df_tracks['sample_name']=sample_name
        df_tracks=df_tracks.sort_values(by=["sample_name", "TrackID", "relative_time"])
        df_all_tracks = pd.concat([df_all_tracks, df_tracks])
        
        end_time = time.time()
        h,m,s = format_time(start_time, end_time)
        print(f"###### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
        
    all_tracks_out_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_all_tracks.to_csv(all_tracks_out_path, index=False) 
    return(df_all_tracks)       

def calculate_image_based_track_features(
    df_tracks,
    cell_type,
    dead_channel,
    contact_threshold,
    element_size_x,
    element_size_y,
    element_size_z,
    #Paths
    raw_image_path,
    dead_mask_path,
    organoid_segments_path,
    tcell_segments_path,
    df_dead_mask_outpath,
    df_morphology_outpath,
    df_intensity_outpath,
    df_contacts_outpath,
    # Overwrite/redo df_intensity_outpath and df_contacts_outpath
    n_workers=1,
    overwrite=False,
    ):
    """
    This code calculates the various features for each timepoint in a track for each 
    separate experiment
    
    This codes works with either:
    - The generated segments and tracks from the BEHAV3D preprocessing modules
    - Imaris extracted statistics
    
    Output:
    - A .csv file containing all timepoints of all TrackIDs and their time-related 
      features
    """

    #TODO Add possibility to add multiple T cell types
    #TODO So both CD4 and CD8 segments, label the type for each track
    #TODO Then get distance and contact between all of them
    #TODO Perhaps allow for input of the track df that already has T cell type in there

    # Load in the images containing the organoid segments and T cell segments
    organoid_segments=load_image(organoid_segments_path)
    tcell_segments=load_image(tcell_segments_path)
    intensity_image = load_image(raw_image_path)
    dead_mask = load_image(dead_mask_path)
    
    if cell_type=="tcell":
        segments= tcell_segments
        segments_path = tcell_segments_path
    elif cell_type=="organoid":
        segments= organoid_segments
        segments_path = organoid_segments_path
        
    print(f"{get_current_time()} - Calculating morphology features...")
    morph_dtypes = {
        "TrackID": int,
        "position_t": int,
        "volume": float,
        "bbox_volume": float,
        "elongation": float,
        "extent": float,
        "equivalent_diameter": float,
        "major_axis_length": float,
        "minor_axis_length": float,
        "surface_area": float,
        "sphericity": float,
        "convex_volume": float,
        "orientation_vector": object,
    }
    
    if df_morphology_outpath.exists() and not overwrite:
        print("Morphology calculation .csv already exists. Loading in intensity information...")
        df_morphology = pd.read_csv(df_morphology_outpath, dtype=morph_dtypes)
        df_morphology["orientation_vector"] = df_morphology["orientation_vector"].apply(ast.literal_eval)
    else:
        df_morphology=calculate_morphology_features(
            segments_path=segments_path,
            n_workers=n_workers
        )
        df_morphology = df_morphology.astype(morph_dtypes)
        df_morphology.to_csv(df_morphology_outpath, sep=",", index=False)  
    df_tracks = pd.merge(df_tracks, df_morphology, how="left")
    
    print(f"{get_current_time()} - Calculating channel and especially death dye intensities...")
    if df_intensity_outpath.exists() and not overwrite:
        print("Intensity calculation .csv already exists. Loading in intensity information...")
        df_intensity = pd.read_csv(df_intensity_outpath)
    else:
        df_intensity=calculate_segment_intensity(
            segments=segments,
            intensity_image=intensity_image
        )
        if dead_channel is not None:
            df_intensity = df_intensity.rename(columns={f"mean_intensity_ch{dead_channel}":"mean_dead_dye"})
        df_intensity.to_csv(df_intensity_outpath, sep=",", index=False)
    df_tracks = pd.merge(df_tracks, df_intensity, how="left")
    
    print(f"{get_current_time()} - Calculating number of dead mask pixels...")
    if df_dead_mask_outpath.exists() and not overwrite:
        print("Dead mask calculation .csv already exists. Loading in intensity information...")
        df_dead_mask = pd.read_csv(df_dead_mask_outpath)
    else:
        df_dead_mask=calculate_dead_mask(
            segments=segments,
            dead_mask=dead_mask
        )
        df_dead_mask.to_csv(df_dead_mask_outpath, sep=",", index=False)
    df_tracks = pd.merge(df_tracks, df_dead_mask, how="left")
    
    print(f"{get_current_time()} - Calculating contact with organoids and other T cells...")
    print(f"Using a contact threshold of {contact_threshold} um")
    # Calculate the contact od each T cell with an organoid or a T cell
    # Explanation on how in the function itself
    contact_dtypes = {
                'TrackID': int,
                'position_t': float,
                'organoid_contact': bool,
                'organoid_contact_pixels': bool,
                'touching_organoids': str,
                'tcell_contact': bool,
                'tcell_contact_pixels': bool,
                'touching_tcells': str
                }
    if df_contacts_outpath.exists() and not overwrite:
        print("Contact .csv already exists. Loading in contact information...")
        df_contacts = pd.read_csv(
            df_contacts_outpath,
            dtype=contact_dtypes
            )

    else:
        df_contacts=calculate_organoid_and_tcell_contact(
            tcell_segments_path=tcell_segments_path,
            organoid_segments_path=organoid_segments_path,
            element_size_x=element_size_x,
            element_size_y=element_size_y,
            element_size_z=element_size_z,
            contact_threshold=contact_threshold,
            calculate_from=cell_type,
            n_workers=n_workers
        ) 
        df_contacts = df_contacts.astype(contact_dtypes)
        df_contacts.to_csv(df_contacts_outpath, sep=",", index=False)
    
    df_tracks = pd.merge(df_tracks, df_contacts, how="left")
    return(df_tracks)

def generalize_units_of_track_features(
    df_tracks,
    distance_unit,
    time_interval,
    time_unit
    ):
    # Converting the time and distance values to a default unit to allow comparison 
    # with differently provided units (defaults to µm and hours)
    print(f"{get_current_time()} - Converting distance and time unit to default um and hours...")   
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
    
    # Make default distance unit um over μm, as encoding of μm between python 
    # and e.g Excel can lead to formatting errors
    df_tracks["distance_unit"] = "um"
    df_tracks["time_unit"] = "h"
    return(df_tracks)
    
def calculate_tcell_specific_track_features(
    df_tracks,
    dead_dye_threshold,
    ):

    print(f"{get_current_time()} - Calculating cell death based on defined dead_dye_threshold {dead_dye_threshold}")
    df_tracks = calculate_death(df_tracks, threshold=dead_dye_threshold, threshold_column="mean_dead_dye")
    # Calculate various movement features such as speed and mean square displacement of the tracks  
    
    print(f"{get_current_time()} - Determining active contact of T cells")
    # Determining if a T cell is actively interacting with another T cell based on speed
    # More explanation at the top of this code
    df_tracks['list_touching_tcells'] = df_tracks['touching_tcells'].apply(
        lambda x: [] if pd.isna(x) or x=="" else list(map(int, x.split(',')))
        )
    
    active_interaction = []
    for _, row in df_tracks.iterrows():
        if row['tcell_contact']:
            max_mean_speed = max(row['mean_speed'], df_tracks.loc[df_tracks['SegmentID'].isin(row['list_touching_tcells']), 'mean_speed'].max())
            active_interaction.append(row['mean_speed'] == max_mean_speed)
        else:
            active_interaction.append(False)
    df_tracks=df_tracks.drop('list_touching_tcells', axis=1)      
    df_tracks['active_tcell_contact'] = active_interaction
    return(df_tracks)

def calculate_organoid_specific_track_features(
    df_tracks,
    dead_dye_threshold,
    ):

    print(f"{get_current_time()} - Calculating cell death based on nr dead pixels {dead_dye_threshold}")
    df_tracks = calculate_death(df_tracks, threshold=dead_dye_threshold, threshold_column="mean_dead_dye")
    # Calculate various movement features such as speed and mean square displacement of the tracks  
    
    return(df_tracks)


def calculate_imaris_track_features(
    df_tracks,
    organoid_contact_threshold,
    tcell_contact_threshold,
    distance_unit,
    ):
    """
    This code calculates the various features for each timepoint in a track for each 
    separate experiment
    
    This codes works with either:
    - The generated segments and tracks from the BEHAV3D preprocessing modules
    - Imaris extracted statistics
    
    Output:
    - A .csv file containing all timepoints of all TrackIDs and their time-related 
      features
    """
    
    print("Performing feature calculation from Imaris processing..")
    
    print(f"{get_current_time()} - Calculating contact with organoids... (From Imaris)")
    # Threshold the distance to organoid based on the supplied "organoid_contact_threshold"
    # This distance is calculated before in Imaris and supplied as a separate channel and
    # Extracted as a statistic (...Intensity_Min_Ch<#>_img=<#>.csv)
    df_tracks["organoid_contact"]=df_tracks["organoid_distance"]<=organoid_contact_threshold
    
    # Calculate the nearest Tcell based on the distance between the centroids of other T cells
    # Caution: This distance, unlike the BEHAV3D processing, is between centroids of cells, not borders
    # Thus the provided "tcell_contact_threshold" needs to reflect this
    print(f"{get_current_time()} - Calculating contact with T cells... (From Imaris)")
    print(f"Using a contact threshold of {tcell_contact_threshold}{distance_unit}")
    grouped = df_tracks.groupby('position_t')
    new_dfs = []

    ### The following code calculates the distances between all tracks of the same
    ### type (So between for example CD4 and other CD4 cells)
    touching_tcells_dict = {}
    # Calculate distances between a segment and all other segments with cdist
    for group_name, group_df in grouped:
        positions = group_df[['position_x', 'position_y', 'position_z']].values
        distances = cdist(positions, positions)
        np.fill_diagonal(distances, np.inf)
        distances_mask = distances <= tcell_contact_threshold

        tcell_contacts_list = []
        for i, row in enumerate(distances_mask):
            tcell_contacts = np.where(row)[0].tolist()
            tcell_contacts=[group_df.reset_index().loc[idx]["TrackID"] for idx in tcell_contacts]
            tcell_contacts_list.append(tcell_contacts)                 
        group_df['touching_tcells'] = tcell_contacts_list

        for track_id, tcell_contacts in zip(group_df['TrackID'], tcell_contacts_list):
            touching_tcells_dict[track_id] = tcell_contacts

        group_df['tcell_contact'] = group_df['touching_tcells'].apply(lambda x: len(x) > 0)

        ### As Imaris does not give interacting IDs, we can only add contact 
        ### but leave the actual interacting ID as unknown
        
        if "complementary_tcell_distance" in df_tracks.columns:
            # !! ADJUST THIS THRESHOLD OR MOVE TO AFTER ADDING ALL TRACKS
            group_df['tcell_contact'] = group_df['complementary_tcell_distance']<=tcell_contact_threshold
            group_df['touching_tcells'].append("unknown")
        group_df['touching_tcells'] = group_df['touching_tcells'].apply(lambda x: ",".join(map(str, x)) if isinstance(x, list) and len(x) > 0 else None)
        
        new_dfs.append(group_df)
    df_tracks = pd.concat(new_dfs, ignore_index=True)
    return(df_tracks)

def calculate_death(
    df_tracks,
    threshold,
    threshold_column="mean_dead_dye"
    ):
    # print(f"- Calculating cell death based on defined dead_dye_threshold {dead_dye_threshold}")
    df_tracks["dead"] = False
    
    # For any cell crossing the dead_dye_threshold, set the cell to dead. Any timepoint after this timepoint are
    # Also set to dead, even if the mean dead dye intensity goes under the threshold again
    for track_id in df_tracks["TrackID"].unique():
        track_df = df_tracks[df_tracks["TrackID"] == track_id]
        track_df_reset = track_df.reset_index(drop=True)
        threshold_indices = track_df_reset.reset_index(drop=True)[track_df_reset[threshold_column] >= threshold].index
        
        if not threshold_indices.empty:
            first_threshold_index = threshold_indices.min()
            df_tracks.loc[track_df.index[first_threshold_index:], "dead"] = True
    return df_tracks

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
           
def interpolate_missing_positions(
    df_tracks,
    cols_to_copy=[
        "TrackID",
        "organoid_contact",
        "organoid_contact_pixels",
        "touching_organoids",
        "tcell_contact",
        "tcell_contact_pixels",
        "touching_tcells",
        "distance_unit",
        "time_unit",
        "orientation_vector"
        ],
    cols_to_interpolate=[
        # "position_t", 
        # "position_z", 
        # "position_y", 
        # "position_x",
        # "mean_dead_dye"
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
    if  cols_to_interpolate is None or cols_to_interpolate == []:
        # Select all columns that are not in cols_to_copy
        cols_to_interpolate = df_tracks.columns.difference(cols_to_copy).tolist()
     
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
                    df_interpolated.at[idx, col] = df_interpolated.at[previous_idx, col]
        df_interpolated["interpolated"] = df_interpolated["interpolated"].astype("boolean").fillna(True)
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
    for t, tcell_stack in tqdm(enumerate(tcell_segments), total=len(tcell_segments)):
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


def _calculate_organoid_and_tcell_contact_single_timepoint(args):
    (
        t,
        tcell_segments_path,
        organoid_segments_path,
        element_size_x,
        element_size_y,
        element_size_z,
        contact_threshold,
        calculate_from
    ) = args
    
    tcell_segments = np.asarray(load_image(tcell_segments_path)[t])
    organoid_segments = np.asarray(load_image(organoid_segments_path)[t])
    
    if calculate_from == "tcell":
        segments_stack = tcell_segments
    elif calculate_from == "organoid":
        segments_stack = organoid_segments
    else:
        raise ValueError(f"calculate_from has to be either 'tcell' or 'organoid', got {calculate_from}")
    
    df_contacts = []
    segment_ids = np.unique(segments_stack)
    
    for segment_id in segment_ids:
        if segment_id == 0:
            continue
        
        stack_max_z, stack_max_y, stack_max_x = segments_stack.shape
        seg_locs = np.argwhere(segments_stack == segment_id)
        min_z, min_y, min_x = seg_locs.min(axis=0)
        max_z, max_y, max_x = seg_locs.max(axis=0)
        
        z_ext = 2 * math.ceil(contact_threshold / element_size_z)
        y_ext = 2 * math.ceil(contact_threshold / element_size_y)
        x_ext = 2 * math.ceil(contact_threshold / element_size_x)
        
        slicer = (
            slice(max(0, min_z - z_ext), min(stack_max_z, max_z + z_ext + 1)),
            slice(max(0, min_y - y_ext), min(stack_max_y, max_y + y_ext + 1)),
            slice(max(0, min_x - x_ext), min(stack_max_x, max_x + x_ext + 1))
        )
        
        tcell_cutout = tcell_segments[slicer]
        org_cutout = organoid_segments[slicer]
        seg_cutout = segments_stack[slicer]
        
        real_distances = distance_transform_edt(
            seg_cutout != segment_id,
            sampling=[element_size_z, element_size_y, element_size_x]
        )
        pix_distances = distance_transform_edt(seg_cutout != segment_id)

        organoid_contacts = [
            str(x) for x in np.unique(org_cutout[real_distances <= contact_threshold]) if x != 0
        ]
        if calculate_from == "organoid":
            organoid_contacts = [x for x in organoid_contacts if x != str(segment_id)]
        real_organoid_contact = len(organoid_contacts) > 0

        pix_organoid_contacts = [
            str(x) for x in np.unique(org_cutout[pix_distances <= 1.73]) if x != 0
        ]
        pix_organoid_contact = len(pix_organoid_contacts) > 0

        tcell_contacts = [
            str(x) for x in np.unique(tcell_cutout[real_distances <= contact_threshold]) if x != 0
        ]
        if calculate_from == "tcell":
            tcell_contacts = [x for x in tcell_contacts if x != str(segment_id)]
        real_tcell_contact = len(tcell_contacts) > 0

        pix_tcell_contacts = [
            str(x) for x in np.unique(tcell_cutout[pix_distances <= 1.73]) if x not in [0, segment_id]
        ]
        pix_tcell_contact = len(pix_tcell_contacts) > 0

        df_contacts.append(pd.DataFrame([{
            'TrackID': segment_id,
            'position_t': t,
            'organoid_contact': real_organoid_contact,
            'organoid_contact_pixels': pix_organoid_contact,
            'touching_organoids': ",".join(organoid_contacts) if real_organoid_contact else "",
            'tcell_contact': real_tcell_contact,
            'tcell_contact_pixels': pix_tcell_contact,
            'touching_tcells': ",".join(tcell_contacts) if real_tcell_contact else ""
        }]))

    return pd.concat(df_contacts)

def calculate_organoid_and_tcell_contact(
    tcell_segments_path,
    organoid_segments_path,
    contact_threshold,
    element_size_x,
    element_size_y,
    element_size_z,
    calculate_from="tcell",
    n_workers=1
    ):
    """
    Wrapper function to parallelize contact analysis.
    
    Returns:
        DataFrame of contact annotations.
    """
    tcell_segments = load_image(tcell_segments_path)
    timepoints = tcell_segments.shape[0]

    args_list = [
        (
            t,
            tcell_segments_path,
            organoid_segments_path,
            element_size_x,
            element_size_y,
            element_size_z,
            contact_threshold,
            calculate_from
        )
        for t in range(timepoints)
    ]

    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            results = list(tqdm(executor.map(_calculate_organoid_and_tcell_contact_single_timepoint, args_list), total=len(args_list)))
    else:
        results = [_calculate_organoid_and_tcell_contact_single_timepoint(args) for args in tqdm(args_list)]

    return pd.concat(results, ignore_index=True)

def calculate_segment_intensity(segments, intensity_image, calculation="mean"):
    """
    Calculates the intensity of a specific marker features for each segment.
    The calculation can be the minimum, maximum, mean or median
    """
    assert calculation in ["min", "max", "mean", "median"]
    intensity_image=np.transpose(intensity_image, axes=[0,2,3,4,1])
    df_intensity = []
    # for t, (tcell_stack, intensity_stack) in enumerate(zip(segments, intensity_image)):      
    for t, (tcell_stack, intensity_stack) in tqdm(enumerate(zip(segments, intensity_image)),total=len(segments)):      
        tcell_stack = np.asarray(tcell_stack)
        intensity_stack = np.asarray(intensity_stack)
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=intensity_stack, properties=['label', f'intensity_{calculation}']))
        properties["position_t"]=t
        df_intensity.append(properties)
    df_intensity = pd.concat(df_intensity)
    
    # Define a dictionary mapping old column names to new column names using regex
    column_mapping = {}
    for i in range(df_intensity.shape[1]):
        old_col_name = f'intensity_mean-{i}'
        new_col_name = f'mean_intensity_ch{i+1}'
        column_mapping[old_col_name] = new_col_name
    column_mapping["label"]="TrackID"
    df_intensity=df_intensity.rename(columns=column_mapping)
    return(df_intensity)

def calculate_dead_mask(segments, dead_mask):
    """
    Calculates the intensity of a specific marker features for each segment.
    The calculation can be the minimum, maximum, mean or median
    """
    df_intensity = []
    # for t, (tcell_stack, intensity_stack) in enumerate(zip(segments, intensity_image)):      
    for t, (tcell_stack, dead_mask_stack) in tqdm(enumerate(zip(segments, dead_mask)),total=len(segments)):      
        tcell_stack = np.asarray(tcell_stack)
        dead_mask_stack = np.asarray(dead_mask_stack)
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=dead_mask_stack, properties=['label', 'num_pixels', f'intensity_mean']))
        properties["position_t"]=t
        properties["nr_dead_mask_pixels"] = properties["num_pixels"] * properties["intensity_mean"]
        properties=properties.rename(columns={"label":"TrackID"})
        properties = properties[["TrackID", "position_t", "nr_dead_mask_pixels"]]
        df_intensity.append(properties)
    df_intensity = pd.concat(df_intensity)
    return(df_intensity)

def _calculate_morphology_single_timepoint(args):
    """Helper function to compute morphology features for a single timepoint."""
    t, segments_path, voxel_spacing = args
    segments = load_image(segments_path)
    stack = segments[t]
    stack = np.asarray(stack)
    properties = pd.DataFrame(
        regionprops_table(
            label_image=stack,
            properties=[
                "label", "area", "bbox_area", 
                "extent", "solidity", "equivalent_diameter",
                "major_axis_length", "minor_axis_length", "inertia_tensor", 
                "inertia_tensor_eigvals", "moments_central"
            ]
        )
    )
    properties.rename(columns={
        "label": "TrackID",
        "area": "volume",
        "bbox_area": "bbox_volume"
    }, inplace=True)
    properties["elongation"] = properties["major_axis_length"] / properties["minor_axis_length"]
    properties["position_t"] = t

    # Suppress known/harmless warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide")
        warnings.filterwarnings("ignore", message=".*convex hull image.*")
        

        surface_areas = []
        sphericities = []
        convex_volumes = []
        orientations = []

        for region_label in properties["TrackID"]:
            mask = (stack == region_label)
            coords = np.argwhere(mask)
            volume = properties.loc[properties["TrackID"] == region_label, "volume"].values[0]

            # Surface area + sphericity via marching cubes
            if np.count_nonzero(mask) >= 4:  # sanity check for marching_cubes
                try:
                    verts, faces, _, _ = marching_cubes(mask.astype(float), spacing=voxel_spacing)
                    surface_area = mesh_surface_area(verts, faces)
                    sphericity = (np.pi ** (1/3)) * ((6 * volume) ** (2/3)) / surface_area
                except Exception:
                    surface_area = np.nan
                    sphericity = np.nan
            else:
                surface_area = np.nan
                sphericity = np.nan
                
            surface_areas.append(surface_area)
            sphericities.append(sphericity)

            if coords.shape[0] >= 4:
                try:
                    scaled_coords = coords * np.array(voxel_spacing)
                    hull = ConvexHull(scaled_coords)
                    convex_volume = hull.volume
                except Exception:
                    convex_volume = np.nan
            else:
                convex_volume = np.nan
            convex_volumes.append(convex_volume)

            # Extract all the inertia tensor components from the DataFrame
            tensor_columns = [f"inertia_tensor-{i}-{j}" for i in range(3) for j in range(3)]

            orientations = []
            for idx, row in properties.iterrows():
                try:
                    tensor_values = np.array([row[col] for col in tensor_columns])
                    inertia_tensor = tensor_values.reshape((3, 3))
                    eigvals, eigvecs = np.linalg.eigh(inertia_tensor)
                    major_axis_vector = eigvecs[:, np.argmax(eigvals)]
                    orientations.append(major_axis_vector.tolist())
                except Exception:
                    orientations.append([np.nan, np.nan, np.nan])

        properties["surface_area"] = surface_areas
        properties["sphericity"] = sphericities
        properties["convex_volume"] = convex_volumes
        # Guard against divide-by-zero in solidity calculation
        with np.errstate(divide='ignore', invalid='ignore'):
            # Lower solidity = more protrusions or invaginations
            properties["solidity"] = properties["volume"] / properties["convex_volume"]
            
        # properties["solidity"] = properties["volume"] / properties["convex_volume"]
        properties["orientation_vector"] = orientations

    columns_to_drop = [col for col in properties.columns 
                   if col.startswith("inertia_tensor") 
                   or col.startswith("inertia_tensor_eigvals") 
                   or col.startswith("moments_central")]
    properties = properties.drop(columns=columns_to_drop)
    return properties

def calculate_morphology_features(segments_path, voxel_spacing=(1.0, 1.0, 1.0), n_workers=8):
    """
    Calculates morphological features (volume, shape, sphericity, etc.) for 3D segments.
    
    Parameters:
        segments_path: path to segments .zarr
        voxel_spacing: Tuple of (z, y, x) spacing in physical units (e.g., µm). Default is isotropic.
    
    Returns:
        pd.DataFrame with one row per object per time point.
    """
    segments=load_image(segments_path)
    timepoints = segments.shape[0]
    args_list = [(t, segments_path, voxel_spacing) for t in range(timepoints)]
    if n_workers > 1:
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            results = list(tqdm(executor.map(_calculate_morphology_single_timepoint, args_list), total=len(args_list)))
    else:
        results = [_calculate_morphology_single_timepoint(args) for args in tqdm(args_list) ]

    return pd.concat(results, ignore_index=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    metadata = pd.read_csv(config["metadata_csv_path"])
    run_behav3d_feature_extraction(config, metadata)