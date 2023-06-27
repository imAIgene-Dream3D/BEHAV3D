import os
import sys
import argparse
import subprocess
import yaml
from pathlib import Path
from tifffile import imread, imwrite
from trackmate import run_trackmate

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


### Segment the data using Ilastik
def run_bash(script, path_bash="bash", stdin=None):
    """Run a bash-script. Returns (stdout, stderr), raises error on non-zero return code"""
    if sys.version_info[0]<3 or (sys.version_info[0]==3 and sys.version_info[1]<7):
        proc = subprocess.run(
            [path_bash, '-c', script],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
            )
    else:
        proc = subprocess.run(
            [path_bash, '-c', script],
            capture_output=True
            )
    stdout, stderr = proc.stdout, proc.stderr
    return(stdout.decode("utf-8") , stderr.decode("utf-8") , proc.returncode)

## Perform pixel classification
pix_prob_path = Path(output_dir, f"{Path(data_path).stem}_probabilities.h5")
command = (
    f"{ilastik_path} --headless "
    f"--project='{ilastik_pix_clas_model}' "
    # f"--output_format='multipage tiff' "
    # f"--output_filename_format={pix_prob_path} "
    f"--output_format='hdf5' "
    f"--output_internal_path /probabilities "
    f"--output_filename_format={pix_prob_path} "
    f"--output_axis_order=ctzyx "
    f"--export_source='Probabilities' "
    f"{data_path}"
)
out, error, returncode = run_bash(command)

## Perform object classification
preobj_out_path = Path(output_dir, f"{Path(data_path).stem}_preobjects.h5")
preobj_feature_csv_path = Path(output_dir, f"{Path(data_path).stem}_preobject_features.h5")
command = (
    f"{ilastik_path} --headless "
    f"--project='{ilastik_obj_clas_model}' "
    # f"--output_format='multipage tiff' "
    # f"--output_filename_format={preobj_out_path} "
    f"--output_format='hdf5' "
    f"--output_internal_path /pre_objects "
    f"--output_filename_format={preobj_out_path} "
    f"--table_filename={preobj_feature_csv_path} "
    f"--raw_data {data_path} " 
    f"--prediction_maps {pix_prob_path}/probabilities "

    f"--export_source='object identities' "
)
out, error, returncode = run_bash(command)

## Perform pixel classification
# finobj_out_path = Path(output_dir, f"{Path(data_path).stem}_finalobjects.h5")
finobj_out_path = Path(output_dir, f"{Path(data_path).stem}_finalobjects.tiff")
command = (
    f"{ilastik_path} --headless "
    f"--project='{ilastik_postproc_model}' "
    f"--output_format='multipage tiff' "
    #  f"--output_format='hdf5' "
    f"--output_filename_format={finobj_out_path} "
    # f"--output_internal_path /segments "
    # f"--output_axis_order=zyx "
    f"--raw_data {data_path} " 
    f"--segmentation_image {preobj_out_path}/pre_objects "
    f"--export_source='Object-Identities' "
)
out, error, returncode = run_bash(command)

im = imread(finobj_out_path)
imwrite(
   finobj_out_path,
    im.astype('uint16'),
    imagej=True,
    metadata={'axes':'TZYX'}
)

### Track the data using TrackMate
 # https://imagej.net/plugins/trackmate/scripting/scripting
 # https://forum.image.sc/t/trackmate-scripting-command-to-generate-full-track-statistics-table-from-analysis/10166/2
tracks=run_trackmate(str(finobj_out_path))
 
### Process and extract features from the tracks

### Combine into single table with all features per track