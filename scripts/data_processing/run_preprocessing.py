import os
import sys
import argparse
import subprocess
import yaml

from pathlib import Path
from tifffile import imread, imwrite
from ilastik import (
    run_ilastik_pixel_classifier, 
    run_ilastik_object_classifier,
    run_ilastik_object_splitter
)
from trackmate import run_trackmate
from calculate_track_features import calculate_movement_features, calculate_organoid_distance

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
args = parser.parse_args()

with open(args.config, "r") as parameters:
    config=yaml.load(parameters, Loader=yaml.SafeLoader)

with open("/Users/samdeblank/Library/CloudStorage/OneDrive-PrinsesMaximaCentrum/github/BEHAV3D-ilastik/scripts/data_processing/config.yml", "r") as parameters:
    config=yaml.load(parameters, Loader=yaml.SafeLoader)
    
data_path = config['data_path']
output_dir = config['output_dir']

ilastik_path = config['ilastik_path']
ilastik_pix_clas_model = config['ilastik_pixel_classifier_model']
ilastik_obj_clas_model = config['ilastik_object_classifier_model']
ilastik_postproc_model = config['ilastik_object_postprocessing_model']

trackmate_path = config['trackmate_path']

pix_prob_path=run_ilastik_pixel_classifier(
    raw_data=data_path,
    model=ilastik_pix_clas_model,
    output_dir=output_dir,
    ilastik_path=ilastik_path
)

preobj_out_path=run_ilastik_object_classifier(
    raw_data=data_path,
    pixel_probabilities=pix_prob_path,
    model=ilastik_obj_clas_model,
    output_dir=output_dir,
    ilastik_path=ilastik_path
)

finobj_out_path=run_ilastik_object_splitter(
    raw_data=data_path,
    segments=preobj_out_path,
    model=ilastik_postproc_model,
    output_dir=output_dir,
    ilastik_path=ilastik_path
)

### Add metadata to the image so TrackMate correctly reads the dimensions
im = imread(finobj_out_path)
imwrite(
    finobj_out_path,
    im.astype('uint16'),
    imagej=True,
    metadata={'axes':'TZYX'}
)

### Track the data using TrackMate
tracks=run_trackmate(str(finobj_out_path))
tracks.to_csv("/Users/samdeblank/surfdrive/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/testdata_medium_tracks.csv", sep=",", index=False)

### Assign the tracks to existing segments
# Add 1 to every track_id so 0 is not a track in the image (should be background)
tracks["track_id"]=tracks["track_id"]+1

# Loop through spots, link to segments in the image and replace label with track_id
im_track = np.zeros_like(im)
for idx, row in tracks.iterrows():
    t,z,y,x = row["position_t"],row["position_z"], row["position_y"], row["position_x"]
    t,z,y,x = round(t), round(z), round(y), round(x)
    corr_seg = im[t,z,y,x]
    im_track[t,:,:,:][im[t,:,:,:]==corr_seg]=row["track_id"]
    # im_track = im_track[im==corr_seg]=row["track_id"]
imwrite(
    finobj_out_path,
    im_track
) 

### Process and extract features from the tracks
# TODO THIS CHANGES THE TRACKS VARIABLE WHILE RUNNING A FUNCTION.. WEIRD
ttracks=calculate_movement_features(
    tracks, 
    element_size_x=3.54, 
    element_size_y=3.54, 
    element_size_z=1.2
    )

organoid_segments=imread("/Users/samdeblank/surfdrive/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/testdata_medium_organoids.tiff")
df_dist_org=calculate_organoid_distance(
    tcell_segments=im_track,
    organoid_segments=organoid_segments,
    element_size_x=3.54, 
    element_size_y=3.54, 
    element_size_z=1.2
)

pd.merge(tracks, df_dist_org, by=["track_id", "position_y"]
### Combine into single table with all features per track