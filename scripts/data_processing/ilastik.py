import os
import sys
import subprocess

from pathlib import Path
from tifffile import imread, imwrite

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

def run_ilastik_pixel_classifier(raw_data, model, output_dir=None, ilastik_path="ilastik"):
    ## Perform pixel classification
    if output_dir is None:
        out_path = Path(Path(raw_data).parent, f"{Path(raw_data).stem}_probabilities.h5")
    else:
        out_path = Path(output_dir, f"{Path(raw_data).stem}_probabilities.h5")
    command = (
        f"{ilastik_path} --headless "
        f"--project='{model}' "
        # f"--output_format='multipage tiff' "
        # f"--output_filename_format={pix_prob_path} "
        f"--output_format='hdf5' "
        f"--output_internal_path /probabilities "
        f"--output_filename_format={out_path} "
        f"--output_axis_order=ctzyx "
        f"--export_source='Probabilities' "
        f"{raw_data}"
    )
    out, error, returncode = run_bash(command)
    print(out)
    print(error)
    return(out_path)
    
def run_ilastik_object_classifier(raw_data, pixel_probabilities, model, output_dir=None, ilastik_path="ilastik"):
    if output_dir is None:
        out_path = Path(Path(raw_data).parent, f"{Path(raw_data).stem}_preobjects.h5")
        out_csv_path = Path(Path(raw_data).parent, f"{Path(raw_data).stem}_preobject_features.h5")
    else:
        out_path = Path(output_dir, f"{Path(raw_data).stem}_preobjects.h5")
        out_csv_path = Path(output_dir, f"{Path(raw_data).stem}_preobject_features.h5")

    command = (
        f"{ilastik_path} --headless "
        f"--project='{model}' "
        # f"--output_format='multipage tiff' "
        # f"--output_filename_format={preobj_out_path} "
        f"--output_format='hdf5' "
        f"--output_internal_path /pre_objects "
        f"--output_filename_format={out_path} "
        f"--table_filename={out_csv_path} "
        f"--raw_data {raw_data} " 
        f"--prediction_maps {pixel_probabilities}/probabilities "

        f"--export_source='object identities' "
    )
    out, error, returncode = run_bash(command)
    print(out)
    print(error)
    return(out_path)

def run_ilastik_object_splitter(raw_data, segments, model, output_dir=None, ilastik_path="ilastik"):
    if output_dir is None:
        out_path = Path(Path(raw_data).parent, f"{Path(raw_data).stem}_finalobjects.tiff")
    else:
        out_path = Path(output_dir, f"{Path(raw_data).stem}_finalobjects.tiff")

    command = (
        f"{ilastik_path} --headless "
        f"--project='{model}' "
        f"--output_format='multipage tiff' "
        #  f"--output_format='hdf5' "
        f"--output_filename_format={out_path} "
        # f"--output_internal_path /segments "
        # f"--output_axis_order=zyx "
        f"--raw_data {raw_data} " 
        f"--segmentation_image {segments}/pre_objects "
        f"--export_source='Object-Identities' "
    )
    out, error, returncode = run_bash(command)
    print(out)
    print(error)
    return(out_path)