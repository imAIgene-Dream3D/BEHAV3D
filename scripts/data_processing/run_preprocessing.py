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

### Process and extract features from the tracks

### Combine into single table with all features per track