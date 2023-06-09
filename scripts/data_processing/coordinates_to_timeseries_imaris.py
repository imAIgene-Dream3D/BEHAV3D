# -*- coding: utf-8 -*-
"""
Created on Mon May 23 15:26:03 2022
This script processes the new data that is tracked. It takes as input the coordinates of each cell at each timepoint. 
@author: akrashen
"""
import numpy as np
import random
import os
import pandas as pd
# Import the data:

path='E:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BHVD_BEHAV3D/3.Analysis/Imaris/20210803_WT1_n3(1)_dist_trans_t_cell_imaris_Statistics' 
filenames =  os.listdir(path)
#Import all csv and concatenate
print(filenames)
random.seed(123)

dfs = []
for f in filenames:
    newdf = pd.read_csv(path +"/" +f,skiprows=3)
    newdf['basename']= f.split(".")[0]
    dfs.append(newdf)
df = pd.concat(dfs, axis=1)
df=df.iloc[:, [0,5,6,7,9,10,20,21,22]]
##df=df[df['dataset']=='output29']
# Add new column that would be unique to tracks
df['TrackID2'] = df['basename'].astype('str') + '_' + df['TrackID'].astype('str')  # this is the calculated field


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
def coord_to_series(data):  
    dfs_2 = []
    for track in data['TrackID2'].unique():
        df_sliced = data[data['TrackID2'] == track ]
        df_sliced2 = df_sliced[['Position X', 'Position Y', 'Position Z']] ## select the data of interest
        # convert to array
        array = df_sliced2.to_numpy()
        array1= np.diff(array,axis=0,prepend=array[[0]]) #normalized to previous raw coordinates
        arr_norm= array - array[[0]] #normalized to first raw coordinates
        # Compute per timepoint time stats
        array_track_dist=np.apply_along_axis(my_func_speed, 1, array1)
        array_cum_dist = np.cumsum(array_track_dist, axis = 0) 
        array_start_dist=np.apply_along_axis(my_func_displacement, 1, arr_norm)
        array_msd=compute_MSD(array)
        # combine
        d = {'track_dist': array_track_dist, 'cum_dist': array_cum_dist, 'displacement': array_start_dist,'msd':array_msd}
        df_computed = pd.DataFrame(data=d)
        df_result= pd.concat([df_sliced.reset_index(drop=True),df_computed.reset_index(drop=True)], axis=1)
        dfs_2.append(df_result)
    df_processed = pd.concat(dfs_2)
    return df_processed


df_processed = coord_to_series(df)

## Save output to process
os.chdir('E:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BHVD_BEHAV3D/3.Analysis/Imaris/processed_T_cells_stats_Imaris')
df_processed.to_csv("processed_imaris_wt1.csv", index=False)

