### TODO Perhaps set a "static" speed based on quantiles of mean_speed to define static and actively interacting cells
### However this excludes a lot of cells to actively interact, now we take the fastest of the contacting T cells as actively interacting

### TODO WHy is tcell_cntact TRUE everywhere?
"""
This script calculates the features of tracks for BEHAV3D analysis.

Flly flexible - works with ANY cell types defined in metadata using prefixes:
- or_ = Organoid types 
- im_ = Immune types 
- ot_ = Other types 

-------------------------------------
--------------- INPUT ---------------
-------------------------------------

.csv containing the following columns:
- TrackID       (The ID of the track a segments belongs to)
- SegmentID     (The unique ID of the segment of a specific timepoint)
- position_t    (The timepoint of the segment)
- position_z
- position_y
- position_x

-------------------------------------
--------------- OUTPUT --------------
-------------------------------------

# Features of a track at each timepoint per sample in the metadata csv (.csv)
- See "FEATURES TRACKS" below

# Combined summarized features for each track for all samples in metadata csv (.csv)
- Merged DataFrame with all features across samples

-------------------------------------
---------- FEATURES TRACKS ---------- 
-------------------------------------

DYNAMIC CONTACT FEATURES (generated for ALL cell types in metadata):
- {cell_type}_contact               (bool) - Real distance-based contact
- {cell_type}_contact_pixels        (bool) - Pixel-based contact (1.73 diagonal)
- touching_{cell_type}s             (str)  - Comma-separated TrackIDs
- active_{cell_type}_contact        (bool) - Active interaction (works for any cell type)

MORPHOLOGY FEATURES:
- nr_pixels, volume, bbox_volume, elongation, extent, equivalent_diameter
- major_axis_length, minor_axis_length, surface_area, sphericity
- convex_volume, orientation_vector

MOVEMENT FEATURES:
- displacement
- cumulative_displacement
- displacement_from_origin
- mean_square_displacement
- speed
- mean_speed
- interpolated
- time

INTENSITY FEATURES:
- mean_intensity_ch{N}  (for each channel in raw image)
- mean_dead_dye         (if dead_channel specified in metadata)

DEATH FEATURES:
- percentage_dead_mask
- nr_dead_mask_pixels
- increase_dead_mask
- dead                  (bool flag set when threshold crossed)

-------------------------------------
------ CONTACT FEATURE DETAILS ------
-------------------------------------

### {cell_type}_contact
- True/False
Per segment, creates a ZYX cutout with extended range around the segment border
and calculates a distance transform using physical spacing (µm). Any other cell
within "contact_threshold" µm counts as contacting.

### {cell_type}_contact_pixels  
- True/False
Same as above but uses pixel-based threshold (1.73 = diagonal distance) without
physical spacing. More lenient threshold for pixel-touching.

### touching_{cell_type}s
- String (comma-separated TrackIDs)
Contains TrackIDs of all contacting cells of this type. Empty string if none.
Self-contact excluded when calculating from the same cell type.

### active_{cell_type}_contact
- True/False
Identifies actively interacting cells vs passively contacted cells. When multiple
cells of the same type contact each other, ranks by mean_speed over rolling window
(default 10 timepoints). Only the cell with highest speed is marked as "active" interaction.
Works for any cell type (tcell, macro, nk, organoid, etc.).

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
calculating the "active_{cell_type}_contact"

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
from behav3d.core.utils import get_current_time, format_time, convert_time, convert_distance
from behav3d.io.images import load_image, convert_input_files_to_zarr
from tqdm import tqdm
from datetime import datetime

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

def run_feature_extraction(
    metadata, 
    config=None,
    output_dir=None,
    cell_type="tcell",
    features_choice=["movement", "intensity", "morphology", "contact", "death"],
    imaris=False,
    dead_mask_percentage_threshold=None,
    contact_threshold=None,
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
        # dead_dye_threshold=sample_metadata['dead_dye_threshold']
        
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
        
        # For Imaris: Collect all cell-type-specific contact thresholds from metadata
        contact_thresholds = {}
        if imaris:
            for col in sample_metadata.index:
                if col.endswith('_contact_threshold'):
                    # Extract cell type from column name (e.g., 'tcell_contact_threshold' -> 'tcell')
                    cell_type_name = col.replace('_contact_threshold', '')
                    if pd.notna(sample_metadata[col]):
                        contact_thresholds[cell_type_name] = sample_metadata[col]

        dead_channel=sample_metadata['dead_channel']
        
        print("###### Running track feature calculation")
        img_outdir = Path(output_dir, "images", sample_name)
        if not img_outdir.exists():
            img_outdir.mkdir(parents=True)

        analysis_outdir = Path(output_dir, "analysis", cell_type)
        track_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        track_intermediate_outdir= Path(track_outdir, "intermediate_results")
        feature_outdir = Path(analysis_outdir, "track_features")
        
        raw_image_path = sample_metadata['raw_image_path']
        
        # Dynamically find the current cell type's segments path
        current_cell_segments_path = None
        for prefix in ['or', 'im', 'ot']:
            col_name = f"{prefix}_{cell_type}_tracks_image_path"
            if col_name in sample_metadata.index and pd.notna(sample_metadata[col_name]):
                current_cell_segments_path = sample_metadata[col_name]
                break
        
        if current_cell_segments_path is None:
            raise ValueError(f"No tracks_image_path found for cell_type='{cell_type}' in sample {sample_name}")
        
        # Dynamically collect ALL organoid types' paths (for contact calculation)
        organoid_segments_paths = {}
        for col in sample_metadata.index:
            if col.startswith('or_') and col.endswith('_tracks_image_path'):
                parts = col.split('_')
                if len(parts) >= 4:
                    organoid_type = '_'.join(parts[1:-3])  # Extract organoid type name
                    if pd.notna(sample_metadata[col]):
                        organoid_segments_paths[organoid_type] = sample_metadata[col]
        
        # Dynamically collect ALL immune cell types' paths (for contact calculation)
        immune_segments_paths = {}
        for col in sample_metadata.index:
            if col.startswith('im_') and col.endswith('_tracks_image_path'):
                parts = col.split('_')
                if len(parts) >= 4:
                    immune_type = '_'.join(parts[1:-3])  # Extract immune type name
                    if pd.notna(sample_metadata[col]):
                        immune_segments_paths[immune_type] = sample_metadata[col]
        
        # Dynamically collect ALL other cell types' paths (for contact calculation)
        other_segments_paths = {}
        for col in sample_metadata.index:
            if col.startswith('ot_') and col.endswith('_tracks_image_path'):
                parts = col.split('_')
                if len(parts) >= 4:
                    other_type = '_'.join(parts[1:-3])  # Extract other type name
                    if pd.notna(sample_metadata[col]):
                        other_segments_paths[other_type] = sample_metadata[col]

        # Old: Construct dead_mask_path from output directory (kept for reference)
        # dead_mask_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
        
        # New: Read dead_mask_path from metadata (if it exists)
        dead_mask_path = None
        if 'dead_mask_path' in sample_metadata.index and pd.notna(sample_metadata['dead_mask_path']):
            dead_mask_path = Path(sample_metadata['dead_mask_path'])
            if not dead_mask_path.exists():
                print(f"⚠️ Warning: dead_mask_path in metadata does not exist: {dead_mask_path}")
                dead_mask_path = None
        
        print(f"{get_current_time()} - Converting all input files to .zarr for memory efficiency...")
        current_cell_segments_path, raw_image_path = convert_input_files_to_zarr(
            sample_name=sample_name,
            current_cell_segments_path=current_cell_segments_path,
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
        # Find the correct prefixed column (or_, im_, ot_)
        tracks_csv_col = None
        for prefix in ['or', 'im', 'ot']:
            col_name = f"{prefix}_{cell_type}_tracks_csv_path"
            if col_name in sample_metadata.index and pd.notna(sample_metadata[col_name]):
                tracks_csv_col = col_name
                break
        
        if tracks_csv_col is None:
            # Fallback to old non-prefixed format for backward compatibility
            tracks_csv_col = f"{cell_type}_tracks_csv_path"
        
        df_tracks_path = sample_metadata[tracks_csv_col]
        
        df_tracks=pd.read_csv(df_tracks_path, sep=",")
        # Adding a sample name for later combination of multiple track experiments
        df_tracks['sample_name']=sample_name
        
        print(f"{get_current_time()} - Generalizing the units of position and time to um and hours")
        df_tracks = generalize_units_of_track_features(
            df_tracks=df_tracks,
            distance_unit=distance_unit,
            time_interval=time_interval,
            time_unit=time_unit
            )
        time_interval=convert_time(time_interval, time_unit)
        time_unit = "h"
        
        required_features = {"morphology", "intensity", "contact", "death"}

        if any(feature in features_choice for feature in required_features):
            ### Calculate image based features (per timepoint)
            print(f"{get_current_time()} - Calculating single-timepoint image-based features...")
            df_intensity_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_intensity.csv")
            df_contacts_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_contact.csv")
            df_dead_mask_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_dead_mask.csv")
            df_morphology_outpath = Path(track_intermediate_outdir, f"{sample_name}_{cell_type}_morphology.csv")
            
            if imaris:
                df_tracks=calculate_imaris_track_features(
                    df_tracks=df_tracks,
                    cell_type=cell_type,
                    contact_threshold=contact_threshold,
                    contact_thresholds=contact_thresholds,
                    distance_unit=distance_unit,
                )
            else:
                df_tracks=calculate_image_based_track_features(
                    df_tracks=df_tracks,
                    cell_type=cell_type,
                    features_choice=features_choice,
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
                    current_cell_segments_path=current_cell_segments_path,
                    organoid_segments_paths=organoid_segments_paths,
                    immune_segments_paths=immune_segments_paths,
                    other_segments_paths=other_segments_paths,
                    overwrite=overwrite,
                    n_workers=n_workers
                )
        else:
            print(f"{get_current_time()} - Skipping image-based features as none requested in features_choice (morphology, intensity, contact, death)")
            
        
        print(f"{get_current_time()} - Interpolating missing timepoints based on time interval")
        # As sometimes 1 or several timepoints are missing in a track, interpolate these missing rows
        # Values are interpolated linearly, forward filled or left blank based on the column
        # More explanation within the function
        df_tracks= interpolate_missing_positions(df_tracks)

        if "movement" in features_choice:
            print(f"{get_current_time()} - Calculating movement features...")
            df_tracks=calculate_movement_features(
                df_tracks,
                time_interval = time_interval,
                rolling_meanspeed_window=rolling_meanspeed_window
                )
            df_tracks = df_tracks.sort_values(['TrackID', 'position_t'])
        else:
            print(f"{get_current_time()} - Skipping movement features as not requested in features_choice")
            
        if "contact" in features_choice:
            # Calculate death features if threshold specified and dead_channel exists
            if dead_mask_percentage_threshold is not None and dead_channel is not None and pd.notna(dead_channel):
                print(f"{get_current_time()} - Calculating cell death based on dead_mask_percentage_threshold {dead_mask_percentage_threshold}")
                df_tracks = calculate_death(df_tracks, threshold=dead_mask_percentage_threshold, threshold_column="percentage_dead_mask")
            
            # Calculate active contact for any cell type with same-type contacts
            touching_col = f'touching_{cell_type}s'
            if touching_col in df_tracks.columns:
                print(f"{get_current_time()} - Calculating active contact features for {cell_type}...")
                df_tracks = calculate_active_contact_features(
                    df_tracks,
                    cell_type=cell_type
                )
            else:
                print(f"{get_current_time()} - No same-type contacts found for {cell_type}, skipping active_contact calculation")
       
            
        tracks_out_path = Path(track_outdir, f"{sample_name}_{cell_type}_track_features.csv")
        print(f"{get_current_time()} - Writing output to {tracks_out_path}")
        df_tracks.to_csv(tracks_out_path, sep=",", index=False)
        
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
    
    features_choice,
    
    contact_threshold,
    element_size_x,
    element_size_y,
    element_size_z,
    #Paths
   
    dead_mask_path,
    current_cell_segments_path,
    organoid_segments_paths={},
    immune_segments_paths={},
    other_segments_paths={},
    raw_image_path = "",
    
    df_dead_mask_outpath="",
    df_morphology_outpath="",
    df_intensity_outpath="",
    df_contacts_outpath="",
    # Overwrite/redo df_intensity_outpath and df_contacts_outpath
    dead_mask_percentage_threshold=None,
    n_workers=1,
    overwrite=False,
    ):
    """
    This code calculates the various features for each timepoint in a track for each 
    separate experiment.
    
    Fully flexible - works with any cell type from metadata.
    
    Output:
    - A .csv file containing all timepoints of all TrackIDs and their time-related 
      features
    """

    # Load the current cell type's segments
    segments = load_image(current_cell_segments_path)
    segments_path = current_cell_segments_path
    
    # Load dead mask (only if path is provided)
    dead_mask = None
    if dead_mask_path is not None:
        dead_mask = load_image(dead_mask_path)
    
    # Load all organoid types' segments (for contact calculation)
    organoid_segments_dict = {}
    for org_type, org_path in organoid_segments_paths.items():
        organoid_segments_dict[org_type] = load_image(org_path)
    
    # Load all immune cell types' segments (for contact calculation)
    immune_segments_dict = {}
    for immune_type, immune_path in immune_segments_paths.items():
        immune_segments_dict[immune_type] = load_image(immune_path)
    
    # Load all other cell types' segments (for contact calculation)
    other_segments_dict = {}
    for other_type, other_path in other_segments_paths.items():
        other_segments_dict[other_type] = load_image(other_path)
    
    if "morphology" in features_choice:
        print(f"{get_current_time()} - Calculating morphology features...")
        morph_dtypes = {
            "TrackID": int,
            "position_t": int,
            "nr_pixels": int,
            "volume": float,
            "bbox_volume": float,
            "elongation": float,
            "extent": float,
            "equivalent_diameter": float,
            "major_axis_length": float,
            "minor_axis_length": float,
            "axis1_length": float,
            "axis2_length": float,
            "axis3_length": float,
            "oblateness": float,
            "prolateness": float,
            "surface_area": float,
            "sphericity": float,
            "convex_volume": float,
            "orientation_vector": object,
        }
        
        if df_morphology_outpath.exists() and not overwrite:
            print("Morphology calculation .csv already exists. Loading in morphology information...")
            df_morphology = pd.read_csv(df_morphology_outpath, dtype=morph_dtypes)
            df_morphology["orientation_vector"] = df_morphology["orientation_vector"].apply(ast.literal_eval)
        else:
            df_morphology=calculate_morphology_features(
                segments_path=segments_path,
                n_workers=n_workers,
                voxel_spacing=(element_size_z, element_size_y, element_size_x)
            )
            df_morphology = df_morphology.astype(morph_dtypes)
            
            if df_morphology_outpath != "":
                df_morphology.to_csv(df_morphology_outpath, sep=",", index=False)  
        df_tracks = pd.merge(df_tracks, df_morphology, how="left")
    
    else:
        print(f"{get_current_time()} - Skipping full morphology features, only computing fast nr_pixels & volume")
        df_basic = calculate_basic_morphology(
            segments_path=segments_path,
            voxel_spacing=(element_size_z, element_size_y, element_size_x),
            n_workers=n_workers,
        )

        basic_dtypes = {
            "TrackID": int,
            "position_t": int,
            "nr_pixels": int,
            "volume": float,
        }
        df_basic = df_basic.astype(basic_dtypes, copy=False)
        df_tracks = pd.merge(df_tracks, df_basic, how="left")

    if "intensity" in features_choice:
        print(f"{get_current_time()} - Calculating channel and especially death dye intensities...")
        
        intensity_image = load_image(raw_image_path)
        
        if df_intensity_outpath.exists() and not overwrite:
            print("Intensity calculation .csv already exists. Loading in intensity information...")
            df_intensity = pd.read_csv(df_intensity_outpath)
        else:
            df_intensity=calculate_segment_intensity(
                segments=segments,
                intensity_image=intensity_image
            )
            if dead_channel is not None and pd.notna(dead_channel):
                df_intensity = df_intensity.rename(columns={f"mean_intensity_ch{dead_channel}":"mean_dead_dye"})
            
            if df_intensity_outpath != "":
                df_intensity.to_csv(df_intensity_outpath, sep=",", index=False)
        df_tracks = pd.merge(df_tracks, df_intensity, how="left")
    else:
        print(f"{get_current_time()} - Skipping intensity features as not requested in features_choice")
    
    if "death" in features_choice:
        if dead_channel is not None and pd.notna(dead_channel) and dead_mask is not None:
            print(f"{get_current_time()} - Calculating number of dead mask pixels...")
            if df_dead_mask_outpath.exists() and not overwrite:
                print("Dead mask calculation .csv already exists. Loading in dead mask calculation information...")
                df_dead_mask = pd.read_csv(df_dead_mask_outpath)
            else:
                df_dead_mask=calculate_dead_mask(
                    segments=segments,
                    dead_mask=dead_mask
                )
                if df_dead_mask_outpath != "":
                    df_dead_mask.to_csv(df_dead_mask_outpath, sep=",", index=False)
            df_tracks = pd.merge(df_tracks, df_dead_mask, how="left")
        else:
            if dead_mask is None and dead_channel is not None:
                print(f"{get_current_time()} - Skipping death mask calculations: dead_mask_path not found in metadata or file doesn't exist")
            else:
                print(f"{get_current_time()} - Skipping death mask calculations: no dead_channel specified in metadata")
    else:
        print(f"{get_current_time()} - Skipping dead mask calculations as not requested in features_choice")
        
    if "contact" in features_choice:
        print(f"{get_current_time()} - Calculating contacts between {cell_type} and all other cell types...")
        print(f"Using a contact threshold of {contact_threshold} um")
        # Calculate contacts dynamically for all cell types found in metadata
        # The resulting DataFrame will have columns for ALL cell types (e.g., organoid1_contact, macro_contact, etc.)
        
        if df_contacts_outpath.exists() and not overwrite:
            print("Contact .csv already exists. Loading in contact information...")
            df_contacts = pd.read_csv(df_contacts_outpath)

        else:
            df_contacts=calculate_contact_features(
                current_cell_segments_path=current_cell_segments_path,
                organoid_segments_paths=organoid_segments_paths,
                immune_segments_paths=immune_segments_paths,
                other_segments_paths=other_segments_paths,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z,
                contact_threshold=contact_threshold,
                calculate_from=cell_type,
                n_workers=n_workers
            )
            if df_contacts_outpath != "":
                df_contacts.to_csv(df_contacts_outpath, sep=",", index=False)
        df_tracks = pd.merge(df_tracks, df_contacts, how="left")
    else:
        print(f"{get_current_time()} - Skipping contact calculations as not requested in features_choice")
        
    df_tracks = df_tracks.sort_values(
        by=["sample_name", "TrackID", "position_t"], 
        ascending=[True, True, True]
    )    
    new_order = ["sample_name"] + [c for c in df_tracks.columns.tolist() if c != "sample_name"]
    df_tracks = df_tracks[new_order]
    
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
    
def calculate_active_contact_features(df_tracks, cell_type):
    """
    Calculate active contact features for any cell type.
    
    Determines which cells are "actively" interacting vs "passively" contacted.
    When multiple cells of the same type are in contact, the one with the highest
    mean_speed is considered to be actively engaging.
    
    Args:
        df_tracks: DataFrame with track features including mean_speed and touching_{cell_type}s columns
        cell_type: The cell type to calculate active contact
    Returns:
        DataFrame with added active_{cell_type}_contact column
    """
    touching_col = f'touching_{cell_type}s'
    contact_col = f'{cell_type}_contact'
    active_contact_col = f'active_{cell_type}_contact'
    
    # Check if mean_speed column exists (requires movement features to be calculated first)
    if 'mean_speed' not in df_tracks.columns:
        print(f"{get_current_time()} - Warning: 'mean_speed' column not found. Skipping active contact calculation.")
        print(f"    Make sure 'movement' is included in features_choice to enable active contact detection.")
        # Set all active contacts to False since we can't determine activity without speed
        df_tracks[active_contact_col] = False
        return df_tracks
    
    print(f"{get_current_time()} - Determining active contact of {cell_type}")

    # --- Step 1: Explode touching cells column ---
    df_explode = (
        df_tracks[['TrackID', 'position_t', touching_col]]
        .dropna()
        .assign(**{touching_col: lambda d: d[touching_col].astype(str).str.split(',')})
        .explode(touching_col)
    )
    df_explode = df_explode[df_explode[touching_col].str.strip() != '']
    df_explode[touching_col] = df_explode[touching_col].astype(int)

    # --- Step 2: Get mean_speed of each touching cell ---
    speed_map = df_tracks[['TrackID', 'position_t', 'mean_speed']].rename(
        columns={'TrackID': touching_col, 'mean_speed': 'touching_speed'}
    )
    df_explode = df_explode.merge(speed_map, on=[touching_col, 'position_t'], how='left')

    # --- Step 3: Aggregate max touching speed for each cell ---
    max_touching_speed = df_explode.groupby(['TrackID', 'position_t'])['touching_speed'].max()

    # --- Step 4: Compare own speed with max of touching cells ---
    df_tracks = df_tracks.set_index(['TrackID', 'position_t'])
    df_tracks['max_touching_speed'] = max_touching_speed
    df_tracks['max_touching_speed'] = df_tracks['max_touching_speed'].fillna(-1)

    # Active contact = in contact AND moving faster than all touching cells
    df_tracks[active_contact_col] = (
        (df_tracks[contact_col]) & 
        (df_tracks['mean_speed'] >= df_tracks['max_touching_speed'])
    )

    df_tracks.reset_index(inplace=True)
    df_tracks.drop(columns=['max_touching_speed'], inplace=True)

    return df_tracks




def calculate_imaris_track_features(
    df_tracks,
    cell_type,
    contact_threshold,
    distance_unit,
    ):
    """
    Calculate contact features from Imaris-extracted statistics.
    Flexible - works with any cell type from metadata.
    
    Imaris provides pre-calculated distance columns in the CSV.
    This function thresholds those distances to determine contacts.
    
    Expected Imaris distance columns in df_tracks:
    - For organoid types: 'organoid_distance', 'organoid1_distance', 'organoid2_distance', etc.
    - For immune types: 'tcell_distance', 'macro_distance', etc.
    - For complementary contacts: 'complementary_{cell_type}_distance'
    
    Output:
    - DataFrame with contact columns for all detected cell types
    """
    
    print("Performing feature calculation from Imaris processing..")
    print(f"Using a contact threshold of {contact_threshold}{distance_unit}")
    
    # Detect all distance columns in the DataFrame
    distance_cols = [col for col in df_tracks.columns if col.endswith('_distance')]
    
    # Calculate contacts for each distance column
    for dist_col in distance_cols:
        # Extract cell type name from column (e.g., 'organoid_distance' -> 'organoid')
        other_cell_type = dist_col.replace('_distance', '').replace('complementary_', '')
        
        # Skip if this is the same cell type (will be handled by same-type contact below)
        if other_cell_type == cell_type:
            continue
        
        # Create contact column
        contact_col = f"{other_cell_type}_contact"
        df_tracks[contact_col] = df_tracks[dist_col] <= contact_threshold
        
        print(f"{get_current_time()} - Calculated contact with {other_cell_type} (From Imaris)")
    
    # Calculate same-type contacts (e.g., tcell-to-tcell, macro-to-macro)
    # This is based on centroid distances between cells of the same type
    print(f"{get_current_time()} - Calculating same-type contacts ({cell_type}-to-{cell_type})... (From Imaris)")
    
    grouped = df_tracks.groupby('position_t')
    new_dfs = []
    
    for group_name, group_df in grouped:
        positions = group_df[['position_x', 'position_y', 'position_z']].values
        distances = cdist(positions, positions)
        np.fill_diagonal(distances, np.inf)
        distances_mask = distances <= contact_threshold
        
        same_type_contacts_list = []
        for i, row in enumerate(distances_mask):
            contacts = np.where(row)[0].tolist()
            contact_ids = [group_df.reset_index().loc[idx]["TrackID"] for idx in contacts]
            same_type_contacts_list.append(contact_ids)
        
        group_df[f'touching_{cell_type}s'] = same_type_contacts_list
        group_df[f'{cell_type}_contact'] = group_df[f'touching_{cell_type}s'].apply(lambda x: len(x) > 0)
        
        # Check for complementary distance column (if Imaris calculated it)
        complementary_col = f'complementary_{cell_type}_distance'
        if complementary_col in df_tracks.columns:
            group_df[f'{cell_type}_contact'] = group_df[complementary_col] <= contact_threshold
            group_df[f'touching_{cell_type}s'].append("unknown")
        
        group_df[f'touching_{cell_type}s'] = group_df[f'touching_{cell_type}s'].apply(
            lambda x: ",".join(map(str, x)) if isinstance(x, list) and len(x) > 0 else None
        )
        
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
           
def interpolate_missing_positions(
    df_tracks,
    cols_to_copy=None,  # Will auto-detect contact columns + metadata
    cols_to_interpolate=None,
    col_to_none = [
        "SegmentID",
    ]
    ):
    """
    As not every track has a segment at every timepoint, interpolate the missing values of
    the missing timepoints.
    
    Fully dynamic - automatically detects ALL contact columns from dataframe.
    
    It interpolates various columns in different ways:
    -   Interpolates the numerical columns of [cols_to_interpolate] such as speed using linear
        interpolation
    -   Copies the columns of [cols_to_copy] using a forward fill from the last non-interpolated
        row of each TrackID (includes all *_contact, *_contact_pixels, touching_* columns)
    -   Puts None in any column not specified, such as SegmentID, as no actual segment exists
    """
    
    # Auto-detect contact columns if not specified
    if cols_to_copy is None:
        cols_to_copy = ["sample_name", "TrackID", "distance_unit", "time_unit", "orientation_vector", "principal_axes"]
        
        # Add metadata string columns that should be forward-filled, not interpolated
        for col in df_tracks.columns:
            # Add line_condition columns (organoid1_line_condition, tcell_line_condition, etc.)
            if col.endswith('_line_condition'):
                if col not in cols_to_copy:
                    cols_to_copy.append(col)
            # Add well and exp_nr columns
            if col in ['well', 'exp_nr']:
                if col not in cols_to_copy:
                    cols_to_copy.append(col)
        
        # Add all dynamically-generated contact columns
        for col in df_tracks.columns:
            if (col.endswith('_contact') or 
                col.endswith('_contact_pixels') or 
                col.startswith('touching_')):
                if col not in cols_to_copy:
                    cols_to_copy.append(col)
    
    # Interpolate missing timepoints so each calculation takes the same intervals
    if  cols_to_interpolate is None or cols_to_interpolate == []:
        # Select all columns that are not in cols_to_copy
        cols_to_interpolate = df_tracks.columns.difference(cols_to_copy).tolist()
        cols_to_interpolate = [col for col in cols_to_interpolate if col not in col_to_none]
    
    # Filter to only numeric columns - np.interp cannot handle object/string dtypes
    numeric_cols = df_tracks.select_dtypes(include=[np.number]).columns.tolist()
    cols_to_interpolate = [col for col in cols_to_interpolate if col in numeric_cols]
     
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
    rolling_meanspeed_window=5
    ):
    """
    Calculates various movement features for each timepoint of a track
    """

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
    
    def calculate_directional_persistence(track_coordinates, eps=1e-12):
        """
        Calculate directional persistence per timepoint compared to previous timepoint
        """
        steps = np.diff(track_coordinates, axis=0, prepend=track_coordinates[[0]])
        step_norms = np.apply_along_axis(np.linalg.norm, 1, steps)

        N = len(track_coordinates)
        persistence = np.zeros(N, dtype=float)

        # persistence needs previous step, so start at t=2
        for t in range(2, N):
            a = steps[t-1]
            b = steps[t]
            denom = step_norms[t-1] * step_norms[t]
            if denom > eps:
                persistence[t] = np.dot(a, b) / denom
            else:
                persistence[t] = 0.0
        persistence = np.clip(persistence, -1.0, 1.0)
        return persistence
    
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
        directional_persistence = calculate_directional_persistence(track_array_rel)

        # combine
        df_computed = pd.DataFrame({
            'displacement': displacement, 
            'cumulative_displacement': cumulative_displacement, 
            'displacement_from_origin': displacement_from_origin, 
            'mean_square_displacement':mean_square_displacement,
            'directional_persistence':directional_persistence
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

def calculate_cell_distance(
    current_cell_segments, 
    target_cell_segments, 
    element_size_x, 
    element_size_y, 
    element_size_z,
    current_cell_type="cell",
    target_cell_type="target"
    ):
    """
    Calculate distance from current cell type to target cell type.
    Fully flexible - works with any cell type combination.
    
    Args:
        current_cell_segments: Segments of the cell type being analyzed
        target_cell_segments: Segments of the target cell type
        element_size_x, element_size_y, element_size_z: Physical spacing in µm
        current_cell_type: Name of current cell type (for column naming)
        target_cell_type: Name of target cell type (for column naming)
    
    Returns:
        DataFrame with distance measurements for each cell at each timepoint
    """
    df_dist = []
    for t, current_stack in tqdm(enumerate(current_cell_segments), total=len(current_cell_segments)):
        target_stack = target_cell_segments[t,:,:,:]
        mask_target = np.ma.masked_where(target_stack==0, target_stack)
        dist_target = distance_transform_edt(mask_target.mask)
        real_dist_target = distance_transform_edt(
            mask_target.mask,
            sampling=[element_size_z, element_size_y, element_size_x]
        )
        properties_pix = pd.DataFrame(
            regionprops_table(
                label_image=current_stack, 
                intensity_image=dist_target, 
                properties=['label', 'intensity_min']
            )
        )
        properties_pix = properties_pix.rename(columns={"intensity_min": f"pix_distance_{target_cell_type}s"})
        properties_real = pd.DataFrame(
            regionprops_table(
                label_image=current_stack, 
                intensity_image=real_dist_target, 
                properties=['label', 'intensity_min']
            )
        )
        properties_real = properties_real.rename(columns={"intensity_min": f"real_distance_{target_cell_type}s"})
        properties = pd.merge(properties_pix, properties_real, how="left")
        properties["position_t"] = t
        df_dist.append(properties)
    df_dist = pd.concat(df_dist)
    df_dist = df_dist.rename(columns={"label": "TrackID"})
    df_dist[f"pix_{target_cell_type}_contact"] = df_dist[f"pix_distance_{target_cell_type}s"] <= 1.73
    return(df_dist)


def _calculate_contact_single_timepoint(args):
    """
    Calculate contacts between current cell type and ALL other cell types.
    Fully flexible - works with any combination of cell types.
    """
    (
        t,
        current_cell_segments_path,
        organoid_segments_paths,
        immune_segments_paths,
        other_segments_paths,
        element_size_x,
        element_size_y,
        element_size_z,
        contact_threshold,
        calculate_from
    ) = args
    
    # Load current cell type's segments for this timepoint
    current_segments = np.asarray(load_image(current_cell_segments_path)[t])
    
    # Load all organoid types' segments
    organoid_segments_dict = {}
    for org_type, org_path in organoid_segments_paths.items():
        organoid_segments_dict[org_type] = np.asarray(load_image(org_path)[t])
    
    # Load all immune cell types' segments
    immune_segments_dict = {}
    for immune_type, immune_path in immune_segments_paths.items():
        immune_segments_dict[immune_type] = np.asarray(load_image(immune_path)[t])
    
    # Load all other cell types' segments
    other_segments_dict = {}
    for other_type, other_path in other_segments_paths.items():
        other_segments_dict[other_type] = np.asarray(load_image(other_path)[t])
    
    df_contacts = []
    segment_ids = np.unique(current_segments)
    
    for segment_id in segment_ids:
        if segment_id == 0:
            continue
        
        stack_max_z, stack_max_y, stack_max_x = current_segments.shape
        seg_locs = np.argwhere(current_segments == segment_id)
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
        
        seg_cutout = current_segments[slicer]
        
        real_distances = distance_transform_edt(
            seg_cutout != segment_id,
            sampling=[element_size_z, element_size_y, element_size_x]
        )
        pix_distances = distance_transform_edt(seg_cutout != segment_id)
        
        contact_data = {
            'TrackID': segment_id,
            'position_t': t,
        }
        
        # Calculate contacts with ALL organoid types
        for org_type, org_segments in organoid_segments_dict.items():
            org_cutout = org_segments[slicer]
            
            org_contacts = [
                str(x) for x in np.unique(org_cutout[real_distances <= contact_threshold]) if x != 0
            ]
            # If calculating from this organoid type, exclude self-contact
            if calculate_from == org_type:
                org_contacts = [x for x in org_contacts if x != str(segment_id)]
            
            real_org_contact = len(org_contacts) > 0
            
            pix_org_contacts = [
                str(x) for x in np.unique(org_cutout[pix_distances <= 1.73]) if x != 0
            ]
            pix_org_contact = len(pix_org_contacts) > 0
            
            # Add columns for this organoid type
            contact_data[f'{org_type}_contact'] = real_org_contact
            contact_data[f'{org_type}_contact_pixels'] = pix_org_contact
            contact_data[f'touching_{org_type}s'] = ",".join(org_contacts) if real_org_contact else ""
        
        # Calculate contacts with ALL immune cell types
        for immune_type, immune_segments in immune_segments_dict.items():
            immune_cutout = immune_segments[slicer]
            
            immune_contacts = [
                str(x) for x in np.unique(immune_cutout[real_distances <= contact_threshold]) if x != 0
            ]
            # If calculating from this immune type, exclude self-contact
            if calculate_from == immune_type:
                immune_contacts = [x for x in immune_contacts if x != str(segment_id)]
            
            real_immune_contact = len(immune_contacts) > 0
            
            pix_immune_contacts = [
                str(x) for x in np.unique(immune_cutout[pix_distances <= 1.73]) if x != 0
            ]
            pix_immune_contact = len(pix_immune_contacts) > 0
            
            # Add columns for this immune type
            contact_data[f'{immune_type}_contact'] = real_immune_contact
            contact_data[f'{immune_type}_contact_pixels'] = pix_immune_contact
            contact_data[f'touching_{immune_type}s'] = ",".join(immune_contacts) if real_immune_contact else ""
        
        # Calculate contacts with ALL other cell types
        for other_type, other_segments in other_segments_dict.items():
            other_cutout = other_segments[slicer]
            
            other_contacts = [
                str(x) for x in np.unique(other_cutout[real_distances <= contact_threshold]) if x != 0
            ]
            # If calculating from this other type, exclude self-contact
            if calculate_from == other_type:
                other_contacts = [x for x in other_contacts if x != str(segment_id)]
            
            real_other_contact = len(other_contacts) > 0
            
            pix_other_contacts = [
                str(x) for x in np.unique(other_cutout[pix_distances <= 1.73]) if x != 0
            ]
            pix_other_contact = len(pix_other_contacts) > 0
            
            # Add columns for this other type
            contact_data[f'{other_type}_contact'] = real_other_contact
            contact_data[f'{other_type}_contact_pixels'] = pix_other_contact
            contact_data[f'touching_{other_type}s'] = ",".join(other_contacts) if real_other_contact else ""

        df_contacts.append(pd.DataFrame([contact_data]))




    return pd.concat(df_contacts)

def calculate_contact_features(
    current_cell_segments_path,
    organoid_segments_paths={},
    immune_segments_paths={},
    other_segments_paths={},
    contact_threshold=None,
    element_size_x=None,
    element_size_y=None,
    element_size_z=None,
    calculate_from=None,
    n_workers=1
    ):
    """
    Flexible contact calculation - works with any cell types from metadata.
    
    Args:
        current_cell_segments_path: Path to the current cell type being analyzed
        organoid_segments_paths: Dict of {organoid_type: path} for all organoid types
        immune_segments_paths: Dict of {immune_type: path} for all immune types
        other_segments_paths: Dict of {other_type: path} for all other types
        calculate_from: The cell type we're calculating features for (required for self-contact exclusion)
        
    Returns:
        DataFrame of contact annotations for the current cell type with ALL other cell types.
    """
    current_segments = load_image(current_cell_segments_path)
    timepoints = current_segments.shape[0]

    args_list = [
        (
            t,
            current_cell_segments_path,
            organoid_segments_paths,
            immune_segments_paths,
            other_segments_paths,
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
            results = list(tqdm(executor.map(_calculate_contact_single_timepoint, args_list), total=len(args_list)))
    else:
        results = [_calculate_contact_single_timepoint(args) for args in tqdm(args_list)]

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

def calculate_relative_increase(df, column, nr_timepoints_back, groupby="TrackID"):
    df = df.sort_values(by=[groupby, "position_t"]).copy()
    def relative_increase(group):
        values = group[column].values
        increases = []
        for i in range(len(group)):
            if i < nr_timepoints_back:
                prev = values[0]
            else:
                prev = values[i - nr_timepoints_back]
            increases.append(values[i] - prev)
        return pd.Series(increases, index=group.index)

    return df.groupby(groupby).apply(relative_increase, include_groups=False).reset_index(drop=True).values

# def calculate_relative_increase(df, column, nr_timepoints_back, groupby):
#     """
#     Calculate the percentage increase of a column compared to the value
#     """
#     df = df.sort_values(by=[groupby, "position_t"]).copy()

#     def compute_increase(group):
#         values = group[column].values
#         increases = []
#         for i in range(len(values)):
#             if i >= nr_timepoints_back:
#                 ref = values[i - nr_timepoints_back]
#             else:
#                 ref = values[0]
#             current = values[i]
#             if ref == 0:
#                 increase = float('inf') if current > 0 else 0.0
#             else:
#                 increase = (current - ref) / ref
#             increases.append(increase)
#         return pd.Series(increases, index=group.index)

#     # Use `transform`-like behavior by resetting the index before recombining
#     result = df.groupby(groupby).apply(compute_increase)
#     if isinstance(result.index, pd.MultiIndex):
#         result.index = result.index.droplevel(0)

#     return result

            
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
        properties.rename(columns={"intensity_mean":"percentage_dead_mask"}, inplace=True)
        properties["nr_dead_mask_pixels"] = properties["num_pixels"] * properties["percentage_dead_mask"]
        properties=properties.rename(columns={"label":"TrackID"})
        properties = properties[["TrackID", "position_t", "percentage_dead_mask", "nr_dead_mask_pixels"]]
        df_intensity.append(properties)
    df_intensity = pd.concat(df_intensity)
    df_intensity = df_intensity.sort_values(by=["TrackID", "position_t"]).reset_index(drop=True)
    df_intensity["increase_dead_mask"] = calculate_relative_increase(
            df = df_intensity,
            column="nr_dead_mask_pixels",
            nr_timepoints_back=10,
            groupby="TrackID"
        )
    return(df_intensity)

def _basic_counts_single_timepoint(args):
    """
    Worker: compute nr_pixels and volume for one timepoint using np.bincount.
    """
    t, segments_path, voxel_volume = args
    seg_t = load_image(segments_path)[t]  # expects a (Z, Y, X) label image at time t
    labels = np.asarray(seg_t, dtype=np.int64).ravel()

    if labels.size == 0:
        return pd.DataFrame(columns=["TrackID", "position_t", "nr_pixels", "volume"])

    # Count voxels per label (0 is background)
    counts = np.bincount(labels)
    if counts.size <= 1:
        # no foreground labels
        return pd.DataFrame(columns=["TrackID", "position_t", "nr_pixels", "volume"])

    # Build per-object table (skip background index 0)
    track_ids = np.nonzero(counts[1:] > 0)[0] + 1
    nr_pixels = counts[track_ids]
    volumes = nr_pixels.astype(np.float64) * float(voxel_volume)

    return pd.DataFrame(
        {
            "TrackID": track_ids.astype(np.int64),
            "position_t": np.full(track_ids.shape[0], t, dtype=np.int64),
            "nr_pixels": nr_pixels.astype(np.int64),
            "volume": volumes.astype(np.float64),
        }
    )

def calculate_basic_morphology(segments_path, voxel_spacing=(1.0, 1.0, 1.0), n_workers=8):
    """
    Fast calculation of nr_pixels and volume per object per timepoint.

    Parameters:
        segments_path: path to segments .zarr
        voxel_spacing: (z, y, x) physical spacing
        n_workers: number of processes

    Returns:
        pd.DataFrame with columns: TrackID, position_t, nr_pixels, volume
    """
    segments = load_image(segments_path)
    timepoints = int(segments.shape[0])
    vz, vy, vx = map(float, voxel_spacing)
    voxel_volume = vz * vy * vx

    args_list = [(t, segments_path, voxel_volume) for t in range(timepoints)]
    if n_workers and n_workers > 1:
        with ProcessPoolExecutor(max_workers=int(n_workers)) as ex:
            results = list(
                tqdm(
                    ex.map(_basic_counts_single_timepoint, args_list),
                    total=len(args_list),
                    # desc="Fast morphology (nr_pixels & volume)",
                )
            )
    else:
        results = [
            _basic_counts_single_timepoint(args)
            for args in tqdm(
                args_list, 
                # desc="Fast morphology (nr_pixels & volume)"
                )
        ]

    if not results:
        return pd.DataFrame(columns=["TrackID", "position_t", "nr_pixels", "volume"])

    return pd.concat(results, ignore_index=True)

def _calculate_morphology_single_timepoint(args):
    """
    Function to compute morphology features for a single timepoint 
    
    Features calculated:
    - volume
    - bbox_volume
    - extent
    - solidity
    - equivalent_diameter
    - major_axis_length
    - minor_axis_length
    - elongation
    - surface_area
    - sphericity
    - convex_volume
    - orientation_vector
    - oblateness
    - axis_length_a
    - axis_length_b
    - axis_length_c
    
    """
    t, segments_path, voxel_spacing = args
    segments = load_image(segments_path)
    stack = np.asarray(segments[t])

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*convex hull image.*")
        warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide")

        properties = pd.DataFrame(
            regionprops_table(
                label_image=stack,
                properties=[
                    "label", "num_pixels", "area", "bbox_area",
                    "extent", "solidity", "equivalent_diameter",
                    "major_axis_length", "minor_axis_length", "inertia_tensor",
                    "inertia_tensor_eigvals", "moments_central"
                ],
                spacing=voxel_spacing
            )
        )

    properties.rename(columns={
        "label": "TrackID",
        "area": "volume",
        "bbox_area": "bbox_volume",
        "num_pixels": "nr_pixels",
    }, inplace=True)

    # Derived scalar features
    properties["elongation"] = properties["major_axis_length"] / properties["minor_axis_length"]
    properties["position_t"] = t

    # Expensive per-region metrics
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="divide by zero encountered in scalar divide")
        warnings.filterwarnings("ignore", message=".*convex hull image.*")

        surface_areas = []
        sphericities = []
        convex_volumes = []

        # NEW: principal-axis lengths & ellipticity (3D)
        axis_length_a_list = []
        axis_length_b_list = []
        axis_length_c_list = []
        oblateness_list = []
        prolateness_list = []
        principal_axes_list = []  # each entry is 3x3, columns are unit vectors for a,b,c

        for region_label in properties["TrackID"]:
            mask = (stack == region_label)
            coords = np.argwhere(mask)
            volume = properties.loc[properties["TrackID"] == region_label, "volume"].values[0]

            # Surface area + sphericity via marching cubes (requires 3D support)
            z_coords = coords[:, 0] if coords.size else np.array([0])
            if (z_coords.max() - z_coords.min()) >= 1 and np.count_nonzero(mask) >= 4:
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

            # Convex hull volume (in physical units)
            if coords.shape[0] >= 4:
                try:
                    scaled_coords = coords * np.array(voxel_spacing)
                    hull = ConvexHull(scaled_coords, qhull_options='QJ')
                    convex_volume = hull.volume
                except Exception:
                    convex_volume = np.nan
            else:
                convex_volume = np.nan

            convex_volumes.append(convex_volume)

            if coords.shape[0] >= 3:
                try:
                    pts = coords.astype(float) * np.array(voxel_spacing, dtype=float)  # (N,3)
                    center = pts.mean(axis=0, keepdims=True)
                    X = pts - center

                    # SVD gives principal directions; Vt rows are unit vectors
                    U, S, Vt = np.linalg.svd(X, full_matrices=False)
                    V = Vt.T  # (3,3) columns are principal directions

                    # Project onto principal directions and get extent along each
                    proj = X @ V  # (N,3)
                    lengths = np.ptp(proj, axis=0)  # full lengths along each axis (a',b',c')

                    # Order by descending length: a >= b >= c
                    order = np.argsort(lengths)[::-1]
                    a, b, c = lengths[order]
                    V_sorted = V[:, order]  # columns aligned with (a,b,c)

                    # Oblate and prolate ellipticity
                    e_ob = 1 - (c / a) if a > 0 else np.nan
                    e_pro = 1 - (b / a) if a > 0 else np.nan
                    
                    axis_length_a_list.append(a)
                    axis_length_b_list.append(b)
                    axis_length_c_list.append(c)
                    oblateness_list.append(e_ob)
                    prolateness_list.append(e_pro)
                    principal_axes_list.append(V_sorted.tolist())
                except Exception:
                    axis_length_a_list.append(np.nan)
                    axis_length_b_list.append(np.nan)
                    axis_length_c_list.append(np.nan)
                    oblateness_list.append(np.nan)
                    prolateness_list.append(np.nan)
                    principal_axes_list.append([[np.nan, np.nan, np.nan]] * 3)
            else:
                axis_length_a_list.append(np.nan)
                axis_length_b_list.append(np.nan)
                axis_length_c_list.append(np.nan)
                oblateness_list.append(np.nan)
                principal_axes_list.append([[np.nan, np.nan, np.nan]] * 3)
            # ----------------------------------------------------------------

        properties["surface_area"] = surface_areas
        properties["sphericity"] = sphericities
        properties["convex_volume"] = convex_volumes

        # Guard against divide-by-zero in solidity calculation
        with np.errstate(divide='ignore', invalid='ignore'):
            properties["solidity"] = properties["volume"] / properties["convex_volume"]
            properties["surface_to_volume_ratio"] = properties["surface_area"] / properties["volume"]
        properties["axis1_length"] = axis_length_a_list
        properties["axis2_length"] = axis_length_b_list
        properties["axis3_length"] = axis_length_c_list
        properties["oblateness"] = oblateness_list
        properties["prolateness"] = prolateness_list
        # properties["principal_axes"] = principal_axes_list  # columns: a,b,c

    # ---- FIXED ORIENTATION COMPUTATION (O(R) instead of O(R^2)) ----
    tensor_columns = [f"inertia_tensor-{i}-{j}" for i in range(3) for j in range(3)]
    try:
        T = properties[tensor_columns].to_numpy(dtype=float).reshape(-1, 3, 3)  # (R,3,3)
        # Eigen-decomposition for all regions at once
        eigvals, eigvecs = np.linalg.eigh(T)  # eigvecs shape: (R,3,3), columns are eigenvectors
        max_idx = np.argmax(eigvals, axis=1)  # (R,)
        # Select major-axis eigenvector per region
        major_vecs = np.stack([eigvecs[i][:, max_idx[i]] for i in range(eigvecs.shape[0])], axis=0)
        properties["orientation_vector"] = major_vecs.tolist()
    except Exception:
        # Fallback if any tensors are missing/NaN
        properties["orientation_vector"] = [[np.nan, np.nan, np.nan]] * len(properties)
    # ----------------------------------------------------------------

    # Drop heavy intermediate columns
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

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser()
#     parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
#     parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
#     args = parser.parse_args()
#     with open(args.config, "r") as parameters:
#         config=yaml.load(parameters, Loader=yaml.SafeLoader)
#     metadata = pd.read_csv(config["metadata_csv_path"])
#     run_behav3d_feature_extraction(config, metadata)