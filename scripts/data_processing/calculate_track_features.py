import numpy as np
import pandas as pd
from scipy.ndimage.morphology import distance_transform_edt
from skimage.measure import regionprops_table
import argparse
import yaml

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
args = parser.parse_args()

with open(args.config, "r") as parameters:
    config=yaml.load(parameters, Loader=yaml.SafeLoader)

with open("/Users/samdeblank/Library/CloudStorage/OneDrive-PrinsesMaximaCentrum/github/BEHAV3D-ilastik/scripts/data_processing/config.yml", "r") as parameters:
    config=yaml.load(parameters, Loader=yaml.SafeLoader)
    
sample_name = config['sample_name']
data_path = config['data_path']
output_dir = config['output_dir']

red_lym_channel=config['red_lym_channel']

element_size_x=config['element_size_x']
element_size_y=config['element_size_y'] 
element_size_z=config['element_size_z']

def main(config):
### Process and extract features from the tracks
# Get movement features (displacement, speed, etc.)
    df_tracks=calculate_movement_features(
        df_tracks, 
        element_size_x=3.54, 
        element_size_y=3.54, 
        element_size_z=1.2
        )

    # Get minimal distance of each segment to an organoid
    organoid_segments=imread(seg_org_out_path)
    df_dist_org=calculate_organoid_distance(
        tcell_segments=im_track,
        organoid_segments=organoid_segments,
        element_size_x=3.54, 
        element_size_y=3.54, 
        element_size_z=1.2
    )
    df_tracks = pd.merge(df_tracks, df_dist_org)

    intensity_image=imread(data_path)[red_lym_channel,:,:,:,:]
    df_red_lym_intensity=calculate_segment_intensity(
        tcell_segments=im_track,
        intensity_image=intensity_image
    )
    df_tracks = pd.merge(df_tracks, df_red_lym_intensity)
    df_tracks=df_tracks.rename(columns={"intensity_mean":"dead_dye_mean"})

    df_tracks = df_tracks.sort_values(['track_id', 'position_t'])

def calculate_movement_features(df_tracks, element_size_x, element_size_y, element_size_z):
    # element_size_x=3.54
    # element_size_y=3.54
    # element_size_z=1.2

    df_tracks["real_x"]=df_tracks["position_x"]*element_size_x
    df_tracks["real_y"]=df_tracks["position_y"]*element_size_y
    df_tracks["real_z"]=df_tracks["position_z"]*element_size_z

    ## Convert the coordinates to time series
    def calculate_speed(a):
        """calculate the speed"""
        return (np.linalg.norm(a - [0,0,0]))
    def calculate_displacement(a):
        """calculate the displacement to the first timepoint"""
        return (np.linalg.norm(a - a[[0]]))
    def compute_MSD(path):
        totalsize=len(path)
        msd=[]
        for i in range(totalsize-1):
            j=i+1
            msd.append(np.sum((path[0:-j]-path[j::])**2)/float(totalsize-j))

        msd=np.array(msd)
        msd= np.insert(msd, 0, 0, axis=0)
        return msd

    ## split by unique trackID2 and process
    df_tracks_processed = []
    for track in df_tracks['track_id'].unique():
        df_track = df_tracks[df_tracks['track_id'] == track ].sort_values(by="position_t")
        df_track = df_track[['real_x', 'real_y', 'real_z']] ## select the data of interest
        
        # convert to array
        track_array = df_track.to_numpy()
        track_array_diff= np.diff(track_array,axis=0,prepend=track_array[[0]]) #normalized to previous raw coordinates
        arr_norm= track_array - track_array[[0]] #normalized to first raw coordinates
        # Compute per timepoint time stats
        array_track_dist=np.apply_along_axis(calculate_speed, 1, track_array_diff)
        array_cum_dist = np.cumsum(array_track_dist, axis = 0) 
        array_start_dist=np.apply_along_axis(calculate_displacement, 1, arr_norm)
        array_msd=compute_MSD(track_array)
        # combine
        d = {'speed': array_track_dist, 'cumulative_distance': array_cum_dist, 'displacement_from_start': array_start_dist,'mean_square_displacement':array_msd}
        df_computed = pd.DataFrame(data=d)
        df_result= pd.concat([df_track.reset_index(drop=True),df_computed.reset_index(drop=True)], axis=1)
        df_tracks_processed.append(df_result)
    df_tracks_processed = pd.concat(df_tracks_processed)
    df_tracks_processed=pd.merge(df_tracks, df_tracks_processed)
    return(df_tracks_processed)

def calculate_organoid_distance(tcell_segments, organoid_segments, element_size_x, element_size_y, element_size_z):
    df_dist_organoid = []
    for t, tcell_stack in enumerate(tcell_segments):
        org_stack = organoid_segments[t,:,:,:]
        mask_org= np.ma.masked_where(org_stack==0, org_stack)
        dist_org=distance_transform_edt(mask_org.mask)
        real_dist_org=distance_transform_edt(
            mask_org.mask,
            sampling=[element_size_z, element_size_y, element_size_z]
            )
        properties_pix=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=dist_org, properties=['label', 'intensity_min']))
        properties_pix=properties_pix.rename(columns={"intensity_min":"pix_distance_organoids"})
        properties_real=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=real_dist_org, properties=['label', 'intensity_min']))
        properties_real=properties_real.rename(columns={"intensity_min":"real_distance_organoids"})
        properties=pd.merge(properties_pix,properties_real)
        properties["position_t"]=t
        df_dist_organoid.append(properties)
    df_dist_organoid = pd.concat(df_dist_organoid)
    df_dist_organoid=df_dist_organoid.rename(columns={"label":"track_id"})
    df_dist_organoid["pix_organoid_contact"] =  df_dist_organoid["pix_distance_organoids"] <= 1.73
    return(df_dist_organoid)

def calculate_segment_intensity(tcell_segments, intensity_image, calculation="mean"):
    assert calculation in ["min", "max", "mean", "median"]
    
    df_intensities = []
    for t, tcell_stack in enumerate(tcell_segments):
        intensity_stack=intensity_image[t,:,:,:]
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=intensity_stack, properties=['label', f'intensity_{calculation}']))
        properties["position_t"]=t
        df_intensities.append(properties)
    df_intensities = pd.concat(df_intensities)
    df_intensities=df_intensities.rename(columns={"label":"track_id"})
    return(df_intensities)

if __name__ == "__main__":
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    main(config)