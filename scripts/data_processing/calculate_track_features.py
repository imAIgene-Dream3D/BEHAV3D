import numpy as np
import pandas as pd
from scipy.ndimage.morphology import distance_transform_edt
from skimage.measure import regionprops_table, regionprops


def calculate_movement_features(df_tracks, element_size_x, element_size_y, element_size_z):
    # element_size_x=3.54
    # element_size_y=3.54
    # element_size_z=1.2

    df_tracks["real_x"]=df_tracks["position_x"]*element_size_x
    df_tracks["real_y"]=df_tracks["position_y"]*element_size_y
    df_tracks["real_z"]=df_tracks["position_z"]*element_size_z

    ## Convert the coordinates to time series
    def my_func_speed(a):
        """calculate the speed"""
        return (np.linalg.norm(a - [0,0,0]))
    def my_func_displacement(a):
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
        array_track_dist=np.apply_along_axis(my_func_speed, 1, track_array_diff)
        array_cum_dist = np.cumsum(array_track_dist, axis = 0) 
        array_start_dist=np.apply_along_axis(my_func_displacement, 1, arr_norm)
        array_msd=compute_MSD(track_array)
        # combine
        d = {'track_dist': array_track_dist, 'cum_dist': array_cum_dist, 'displacement': array_start_dist,'msd':array_msd}
        df_computed = pd.DataFrame(data=d)
        df_result= pd.concat([df_track.reset_index(drop=True),df_computed.reset_index(drop=True)], axis=1)
        df_tracks_processed.append(df_result)
    df_tracks_processed = pd.concat(df_tracks_processed)
    return(df_tracks_processed)

def calculate_organoid_distance(tcell_segments, organoid_segments, element_size_x, element_size_y, element_size_z):
    df_dist_organoid = []
    for t, tcell_stack in enumerate(tcell_segments):
        org_stack = organoid_segments[t,:,:,:]
        mask_org= np.ma.masked_where(org_stack==0, org_stack)
        dist_org=distance_transform_edt(mask_org.mask)
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, intensity_image=dist_org, properties=['label', 'intensity_min']))
        properties["position_t"]=t
        df_dist_organoid.append(properties)
    df_dist_organoid = pd.concat(df_dist_organoid)
    df_dist_organoid=df_dist_organoid.rename(columns={"label":"track_id"})
    return(df_dist_organoid)

    