import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt
from skimage.measure import regionprops_table
import argparse
import yaml
from tifffile import imread
from pathlib import Path
import h5py
import math

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
args = parser.parse_args()

def main(config):
    ### Process and extract features from the tracks
    # Get movement features (displacement, speed, etc.)
    print("### Running track feature calculation")
    with open("/Users/samdeblank/Library/CloudStorage/OneDrive-PrinsesMaximaCentrum/github/BEHAV3D-ilastik/scripts/data_processing/config.yml", "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    metadata=pd.read_csv("/Users/samdeblank/Library/CloudStorage/OneDrive-PrinsesMaximaCentrum/github/BEHAV3D-ilastik/configs/metadata_template_ilastik.csv")
    
    output_dir = config['output_dir']

    for _, sample in metadata.iterrows():
        sample_name = sample['sample_name']
        image_path = sample['image_path']
        image_internal_path = sample["image_internal_path"]
        red_lym_channel=sample['dead_dye_channel']

        element_size_x=sample['pixel_distance_xy']
        element_size_y=sample['pixel_distance_xy'] 
        element_size_z=sample['pixel_distance_z']
        distance_unit=sample['distance_unit']
        
        contact_threshold = sample["contact_threshold"]
        
        print(f"### Processing: {sample_name}")
        print("- Loading in tracks csv...")
        df_tracks_path = Path(output_dir, f"{sample_name}_tracks.csv")
        df_tracks=pd.read_csv(df_tracks_path, sep=",")
        
        print("- Calculating movement features...")
        df_tracks=calculate_movement_features(df_tracks)

        print("- Calculating contact with organoids and other T cells...")
        print(f"Using a contact threshold of {contact_threshold}{distance_unit}")
        # Get minimal distance of each segment to an organoid
        organoid_segments_path=Path(output_dir, f"{sample_name}_organoid_segments.tiff")
        organoid_segments=imread(organoid_segments_path)
        
        tcell_segments_path=Path(output_dir, f"{sample_name}_tcells_tracked.tiff")
        tcell_segments=imread(tcell_segments_path)
        
        df_contacts=calculate_organoid_and_tcell_contact(
            tcell_segments=tcell_segments,
            organoid_segments=organoid_segments,
            element_size_x=element_size_x,
            element_size_y=element_size_y,
            element_size_z=element_size_z,
            contact_threshold=contact_threshold
        )
        df_tracks = pd.merge(df_tracks, df_contacts)
        
        print("- Calculating death dye intensities...")
        
        # Split the path into file path and dataset path
        intensity_image = h5py.File(name=image_path, mode="r")[image_internal_path][:][red_lym_channel,:,:,:,:]
        df_red_lym_intensity=calculate_segment_intensity(
            tcell_segments=tcell_segments,
            intensity_image=intensity_image
        )
        df_tracks = pd.merge(df_tracks, df_red_lym_intensity)
        df_tracks=df_tracks.rename(columns={"intensity_mean":"dead_dye_mean"})

        df_tracks = df_tracks.sort_values(['TrackID', 'position_t'])
        
        tracks_out_path = Path(output_dir, f"{sample_name}_tracks_features.csv")
        print(f"- Writing output to {tracks_out_path}")
        df_tracks.to_csv(tracks_out_path, sep=",", index=False)
        
        print("### Done\n")
        
def calculate_movement_features(df_tracks):
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
        df_result = pd.concat([df_result, df_track[["position_t", "SegmentID"]]], axis=1)
        df_tracks_processed.append(df_result)
    df_tracks_processed = pd.concat(df_tracks_processed)
    df_tracks_processed=pd.merge(df_tracks, df_tracks_processed)
    return(df_tracks_processed)

def calculate_organoid_distance(
    tcell_segments, 
    organoid_segments, 
    element_size_x, 
    element_size_y, 
    element_size_z,
    contact_threshold
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
        properties=pd.merge(properties_pix,properties_real)
        properties["position_t"]=t
        df_dist_organoid.append(properties)
    df_dist_organoid = pd.concat(df_dist_organoid)
    df_dist_organoid=df_dist_organoid.rename(columns={"label":"TrackID"})
    df_dist_organoid["pix_organoid_contact"] =  df_dist_organoid["pix_distance_organoids"] <= 1.73
    return(df_dist_organoid)

def calculate_tcell_contact(
    tcell_segments, 
    element_size_x, 
    element_size_y, 
    element_size_z,
    contact_threshold
    ):
    """
    Calculates contact with other tcells by looping through
    the segments and cutting out a small area around it.
    It then calculates the euclidian distance of every pixel
    outside the segment and sees if any other segments are in this
    area.
    
    tcell_contact: 
    Based on the provided 'contact_threshold' value. Any other segment
    in this distance is seen as contacting
    
    tcell_contact_pixel:
    Based on direct pixel contact of one t cell and another, not influenced
    by the element_sizes
    """
    df_contact_tcell = []
    for t, tcell_stack in enumerate(tcell_segments):
        segment_ids = np.unique(tcell_stack)
  
        for segment_id in segment_ids:
            if segment_id==0:
                continue
            
            stack_max_z, stack_max_y, stack_max_x = tcell_stack.shape
            seg_locs = np.argwhere(tcell_stack==segment_id)
            min_z, min_y, min_x = seg_locs.min(axis=0)
            max_z, max_y, max_x = seg_locs.max(axis=0)
            
            z_ext = 2*math.ceil(contact_threshold / element_size_z)
            y_ext = 2*math.ceil(contact_threshold / element_size_y)
            x_ext = 2*math.ceil(contact_threshold / element_size_x)
            segment_cutout = tcell_stack[
                max(0, min_z-z_ext):min(stack_max_z, max_z+z_ext+1),
                max(0, min_y-y_ext):min(stack_max_y, max_y+y_ext+1),
                max(0, min_x-x_ext):min(stack_max_x, max_x+x_ext+1),
                ]
            
            real_dist_tcell=distance_transform_edt(
                segment_cutout!=segment_id,
                sampling=[element_size_z, element_size_y, element_size_x]
                )
            tcell_contacts = [str(x) for x in np.unique(segment_cutout[real_dist_tcell<=contact_threshold]) if x not in [0, segment_id]]
            real_tcell_contact = len(tcell_contacts)>0
            
            pix_dist_tcell=distance_transform_edt(
                segment_cutout!=segment_id
                )
            pix_tcell_contacts = [str(x) for x in np.unique(segment_cutout[pix_dist_tcell<= 1.73]) if x not in [0, segment_id]]
            pix_tcell_contact = len(pix_tcell_contacts)>0
            
            if real_tcell_contact:
                touching_tcells = ",".join(tcell_contacts.tolist())
                
            df_contact_tcell.append(pd.DataFrame([{
                'TrackID': segment_id, 
                'position_t': t,
                'tcell_contact': real_tcell_contact, 
                'tcell_contact_pixel': pix_tcell_contact,
                'touching_tcells': touching_tcells
            }]))
    
    df_contact_tcell=pd.concat(df_contact_tcell)
    return(df_contact_tcell)

def calculate_organoid_and_tcell_contact(
    tcell_segments,
    organoid_segments,
    element_size_x, 
    element_size_y, 
    element_size_z,
    contact_threshold
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
        segment_ids = np.unique(tcell_stack)
        org_stack = organoid_segments[t,:,:,:]
        for segment_id in segment_ids:
            if segment_id==0:
                continue
            
            stack_max_z, stack_max_y, stack_max_x = tcell_stack.shape
            seg_locs = np.argwhere(tcell_stack==segment_id)
            min_z, min_y, min_x = seg_locs.min(axis=0)
            max_z, max_y, max_x = seg_locs.max(axis=0)
            
            z_ext = 2*math.ceil(contact_threshold / element_size_z)
            y_ext = 2*math.ceil(contact_threshold / element_size_y)
            x_ext = 2*math.ceil(contact_threshold / element_size_x)
            segment_cutout = tcell_stack[
                max(0, min_z-z_ext):min(stack_max_z, max_z+z_ext+1),
                max(0, min_y-y_ext):min(stack_max_y, max_y+y_ext+1),
                max(0, min_x-x_ext):min(stack_max_x, max_x+x_ext+1),
                ]
            org_cutout = org_stack[
                max(0, min_z-z_ext):min(stack_max_z, max_z+z_ext+1),
                max(0, min_y-y_ext):min(stack_max_y, max_y+y_ext+1),
                max(0, min_x-x_ext):min(stack_max_x, max_x+x_ext+1),
                ]
            
            real_distances=distance_transform_edt(
                segment_cutout!=segment_id,
                sampling=[element_size_z, element_size_y, element_size_x]
                )
            pix_distances=distance_transform_edt(
                segment_cutout!=segment_id
                )
             
            organoid_contacts = [str(x) for x in np.unique(org_cutout[real_distances<=contact_threshold]) if x!=0]
            real_organoid_contact = len(organoid_contacts)>0
            
            pix_organoid_contacts = [str(x) for x in np.unique(org_cutout[pix_distances<= 1.73]) if x!=0]
            pix_organoid_contact = len(pix_organoid_contacts)>0
            
            tcell_contacts = [str(x) for x in np.unique(segment_cutout[real_distances<=contact_threshold]) if x not in [0, segment_id]]
            real_tcell_contact = len(tcell_contacts)>0

            pix_tcell_contacts = [str(x) for x in np.unique(segment_cutout[pix_distances<= 1.73]) if x not in [0, segment_id]]
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
    assert calculation in ["min", "max", "mean", "median"]
    
    df_intensities = []
    for t, tcell_stack in enumerate(tcell_segments):
        intensity_stack=intensity_image[t,:,:,:]
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=intensity_stack, properties=['label', f'intensity_{calculation}']))
        properties["position_t"]=t
        df_intensities.append(properties)
    df_intensities = pd.concat(df_intensities)
    df_intensities=df_intensities.rename(columns={"label":"TrackID"})
    return(df_intensities)

if __name__ == "__main__":
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    main(config)