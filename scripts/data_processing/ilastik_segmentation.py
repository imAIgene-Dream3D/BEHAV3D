"""
Perform Ilastik segmentation of 4D datasets (czyx) with the following channels:
1 -
2 - 
3 - 
4 - 


Known errors:

If running with verbose=True, there will be an error code for object segmentation:
IndexError: index 1 is out of bounds for axis 0 with size 0

Segmentation is still performed however, and this error can be ignored
"""
import os
import sys
import argparse
import subprocess
import yaml
import h5py

from pathlib import Path
from tifffile import imread, imwrite

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
parser.add_argument('-m', '--metadata', type=str, help='path to a metadata.csv file that stores metadata per image', required=False)
parser.add_argument('-k', '--keep_all', action='store_true', help='Keep original files after checks if compression has been succesful', required=False)
parser.add_argument('-v', '--verbose', action='store_true', help='Verbose', required=False)

args = parser.parse_args()

with open(args.config, "r") as parameters:
    config=yaml.load(parameters, Loader=yaml.SafeLoader)

def main(config, metadata, keep_all=False):
    ilastik_path = config['ilastik_path']
    ilastik_pix_clas_model = config['ilastik_pixel_classifier_model']
    ilastik_org_seg_model = config['ilastik_organoid_segmentation_model']
    ilastik_tcell_seg_model = config['ilastik_tcell_segmentation_model']
    ilastik_postproc_model = config['ilastik_object_postprocessing_model']
    
    for _, sample in metadata.iterrows():
        sample_name = sample['sample_name']
        image_path = sample['image_path']
        image_internal_path = sample["image_internal_path"]
        full_image_path = f"{image_path}/{image_internal_path}"
        output_dir = config['output_dir']

        print(f"### Processing: {sample_name}")
         
        ### General pixel classifier
        print("- Running the ilastik pixel classifier...")
        pix_prob_path = Path(output_dir, f"{sample_name}_probabilities.h5")
        pix_prob_path=run_ilastik_pixel_classifier(
            raw_data=full_image_path,
            model=ilastik_pix_clas_model,
            out_path=pix_prob_path,
            ilastik_path=ilastik_path,
            verbose=verbose
        )

        ### Organoid segmentation
        print("- Performing initial organoid segmentation...")
        preseg_org_h5_path = Path(output_dir, f"{sample_name}_organoid_presegments.h5")
        preseg_org_h5_path=run_ilastik_object_segmentation(
            raw_data=full_image_path,
            pixel_probabilities=pix_prob_path,
            model=ilastik_org_seg_model,
            out_path=preseg_org_h5_path,
            ilastik_path=ilastik_path,
            verbose=verbose
        )

        print("- Performing organoid segment postprocessing (splitting)...")
        seg_org_h5_path = Path(output_dir, f"{sample_name}_organoid_segments.h5")
        seg_org_h5_path=run_ilastik_object_splitter(
            raw_data=full_image_path,
            segments=preseg_org_h5_path,
            model=ilastik_postproc_model,
            out_path=seg_org_h5_path,
            ilastik_path=ilastik_path,
            verbose=verbose
        )

        im = h5py.File(name=seg_org_h5_path, mode="r")["segments"][:].squeeze()
        # im = imread(seg_tcell_out_path)
        seg_org_tiff_path = Path(output_dir, f"{sample_name}_organoid_segments.tiff")
        imwrite(
            seg_org_tiff_path,
            im.astype('uint16'),
            imagej=True,
            metadata={'axes':'TZYX'}
        )

        ### Tcell segmentation
        print("- Performing initial T cell segmentation...")
        preseg_tcell_h5_path = Path(output_dir, f"{sample_name}_tcell_presegments.h5")
        preseg_tcell_h5_path=run_ilastik_object_segmentation(
            raw_data=full_image_path,
            pixel_probabilities=pix_prob_path,
            model=ilastik_tcell_seg_model,
            out_path=preseg_tcell_h5_path,
            ilastik_path=ilastik_path
        )

        print("- Performing T cell segment postprocessing (splitting)...")
        seg_tcell_h5_path = Path(output_dir, f"{sample_name}_tcell_segments.h5")
        seg_tcell_h5_path=run_ilastik_object_splitter(
            raw_data=full_image_path,
            segments=preseg_tcell_h5_path,
            model=ilastik_postproc_model,
            out_path=seg_tcell_h5_path,
            ilastik_path=ilastik_path
        )

        ### Add metadata to the image so TrackMate correctly reads the dimensions
        seg_tcell_tiff_path = Path(output_dir, f"{sample_name}_tcell_segments.tiff")
        im = h5py.File(name=seg_tcell_h5_path, mode="r")["segments"][:].squeeze()
        # im = imread(seg_tcell_out_path)
        imwrite(
            seg_tcell_tiff_path,
            im.astype('uint16'),
            imagej=True,
            metadata={'axes':'TZYX'}
        )
        
        if not keep_all:
            print("- Removing all intermediate files and keeping only final segment tiffs as '--keep_all' is not supplied...")
            os.remove(pix_prob_path)
            os.remove(preseg_org_h5_path)
            os.remove(seg_org_h5_path)
            os.remove(preseg_tcell_h5_path)
            os.remove(seg_tcell_h5_path)
        print("### Done\n")
        
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

def run_ilastik_pixel_classifier(
    raw_data, 
    model, 
    out_path=None, 
    output_dir=None, 
    ilastik_path="ilastik",
    verbose=False
    ):
    ## Perform pixel classification
    if out_path is None:
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
    if verbose:
        print("Running command:")
        print(command)
    
    out, error, returncode = run_bash(command)
    assert returncode==0, f"Pixel classification did not complete correctly\n\nOut:{out}\n\nError:{error}"

    return(out_path)
    
def run_ilastik_object_segmentation(
    raw_data, 
    pixel_probabilities, 
    model, 
    out_path=None, 
    output_dir=None, 
    ilastik_path="ilastik",
    verbose=False
    ):
    if out_path is None:
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
        # f"--table_filename={out_csv_path} "
        f"--raw_data {raw_data} " 
        f"--prediction_maps {pixel_probabilities}/probabilities "

        f"--export_source='object identities' "
    )
    if verbose:
        print("Running command:")
        print(command)
    out, error, returncode = run_bash(command)
    
    assert os.path.isfile(out_path), f"Object segmentation did not complete correctly\n\nOut:\n{out}\n\nError:\n{error}"

    return(out_path)

def run_ilastik_object_splitter(
    raw_data, 
    segments, 
    model, 
    out_path=None, 
    output_dir=None, 
    ilastik_path="ilastik",
    verbose=False
    ):
    if out_path is None:
        if output_dir is None:
            out_path = Path(Path(raw_data).parent, f"{Path(raw_data).stem}_finalobjects.tiff")
        else:
            out_path = Path(output_dir, f"{Path(raw_data).stem}_finalobjects.tiff")

    command = (
        f"{ilastik_path} --headless "
        f"--project='{model}' "
        # f"--output_format='multipage tiff' "
         f"--output_format='hdf5' "
        f"--output_filename_format={out_path} "
        f"--output_internal_path /segments "
        # f"--output_axis_order=zyx "
        f"--raw_data {raw_data} " 
        f"--segmentation_image {segments}/pre_objects "
        f"--export_source='Object-Identities' "
    )
    
    if verbose:
        print("Running command:")
        print(command)
    
    out, error, returncode = run_bash(command)
    assert returncode==0, f"Object splitting did not complete correctly\n\nOut:{out}\n\nError:{error}"
    
    return(out_path)

if __name__ == "__main__":
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    keep_all=args.keep_all
    metadata = pd.read_csv(args.metadata)
    verbose=args.verbose
    main(config, metadata, keep_all, verbose)