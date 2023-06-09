import os
import sys
import argparse
import subprocess
import yaml
from pathlib import Path

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
args = parser.parse_args()

with open(args.config, "r") as parameters:
    config=yaml.load(parameters, Loader=yaml.SafeLoader)

with open("/Users/samdeblank/Library/CloudStorage/OneDrive-PrinsesMaximaCentrum/github/BEHAV3D/scripts/data_processing/config.yml", "r") as parameters:
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
pix_out_path = Path(output_dir, f"{Path(data_path).stem}_probabilities.tiff")
command = (
    f"{ilastik_path} --headless "
    f"--project='{ilastik_pix_clas_model}' "
    f"--output_format='multipage tiff' "
    f"--output_filename_format={pix_out_path} "
    f"--output_axis_order=ctzyx "
    f"--export_source='Probabilities' "
    f"{data_path}"
)
out, error, returncode = run_bash(command)

## Perform object classification
preobj_out_path = Path(output_dir, f"{Path(data_path).stem}_preobjects.tiff")
command = (
    f"{ilastik_path} --headless "
    f"--project='{ilastik_obj_clas_model}' "
    f"--output_format='multipage tiff' "
    f"--output_filename_format={preobj_out_path} "
    # f"--output_axis_order=zyx "
    f"--export_source='simple segmentation' "
    f"{pix_out_path}"
)
out, error, returncode = run_bash(command)

## Perform pixel classification
finobj_out_path = Path(output_dir, f"{Path(data_path).stem}_finalobjects.tiff")
command = (
    f"{ilastik_path} --headless "
    f"--project='{ilastik_postproc_model}' "
    f"--output_format='multipage tiff' "
    f"--output_filename_format={finobj_out_path} "
    f"--output_axis_order=zyx "
    f"--export_source='object identities' "
    f"{preobj_out_path}"
)


print(command)
out, error, returncode = run_bash(command)
### Track the data using TrackMate

### Process and extract features from the tracks

### Combine into single table with all features per track