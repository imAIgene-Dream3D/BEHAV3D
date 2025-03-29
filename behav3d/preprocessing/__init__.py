from tifffile import imread, imwrite
from skimage.measure import regionprops_table
import pandas as pd
import numpy as np
from pathlib import Path
import math
from behav3d.utils.fileio import load_image,save_as_zarr

def convert_input_files_to_zarr(
    tcell_segments_path,
    organoid_segments_path,
    raw_image_path,
    outfolder = None,
    chunks=None,
    overwrite=False
    ):
    
    tcell_segments_path = Path(tcell_segments_path)
    organoid_segments_path = Path(organoid_segments_path)
    raw_image_path = Path(raw_image_path)
    
    if outfolder is None:
        tcell_zarr_out_path = tcell_segments_path.with_suffix(".zarr.zip")
        organoid_zarr_out_path = organoid_segments_path.with_suffix(".zarr.zip")
        raw_image_zarr_out_path = raw_image_path.with_suffix(".zarr.zip")
    else:   
        tcell_zarr_out_path = Path(outfolder, tcell_segments_path.stem + ".zarr.zip")
        organoid_zarr_out_path = Path(outfolder, organoid_segments_path.stem + ".zarr.zip")
        raw_image_zarr_out_path = Path(outfolder, raw_image_path.stem + ".zarr.zip")
    
    if tcell_zarr_out_path.exists() and not overwrite:
        print("Skipping conversion of tcell segments to zarr, as file already exists")
    else:
        tcell_segments = load_image(tcell_segments_path)
        chunksize = (1,) + tcell_segments.shape[1:]
        save_as_zarr(
            img=tcell_segments, 
            path=tcell_zarr_out_path, 
            chunks=chunks
            )
    
    if organoid_zarr_out_path.exists():
       print("Skipping conversion of organoid segments to zarr, as file already exists")
    else:
        organoid_segments = load_image(organoid_segments_path)
        chunksize = (1,) + organoid_segments.shape[1:]
        save_as_zarr(
            img=organoid_segments, 
            path=organoid_zarr_out_path, 
            chunks=chunks
            )

    if raw_image_zarr_out_path.exists():
       print("Skipping conversion of raw_image to zarr, as file already exists")
    else:
        img = load_image(raw_image_path)
        chunksize = (1,) + img.shape[1:]
        save_as_zarr(
            img=img, 
            path=raw_image_zarr_out_path, 
            chunks=chunks
            )
    
    return(
        tcell_zarr_out_path,
        organoid_zarr_out_path,
        raw_image_zarr_out_path
    )

def convert_segments_to_tracks(
        tracks_csv_path,
        segments_path,
        outpath,
        element_size_z,
        element_size_y,
        element_size_x
    ):
    
    ### Assign the tracks to existing segments
    # Loop through spots, link to segments in the image and replace label with TrackID
    print("- Assigning track ID to segmented image to create tracked image...")
    df_tracks = pd.read_csv(tracks_csv_path)
    tcell_segments = imread(segments_path)
    
    df_centroids = []
    # TODO 
    # Check if trackmate removes objects and still keep them in image
    for t, tcell_stack in enumerate(tcell_segments):
        properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, properties=['label', f'centroid']))
        properties["position_t"]=t
        df_centroids.append(properties)
    df_centroids = pd.concat(df_centroids)
    df_centroids["position_z"]=df_centroids["centroid-0"]*element_size_z
    df_centroids["position_y"]=df_centroids["centroid-1"]*element_size_y
    df_centroids["position_x"]=df_centroids["centroid-2"]*element_size_x
    
    tcells_tracked = np.zeros_like(tcell_segments)
    for _, row in df_tracks.iterrows():
        t,z,y,x = int(row["position_t"]),row["position_z"], row["position_y"], row["position_x"]
        corr_seg=None
        # There seems to be an issue where the tracking output gives 
        # e.g. 25.99355 and regionprops * elemen_size gives 25.99354999999999 (Floating point arithmetics issue)
        # To solve we use the math.isclose() with an abs_tol of 0.0001 to still match these values
        corr_seg = df_centroids[
            (df_centroids["position_x"].apply(lambda val: math.isclose(val, x, abs_tol=0.0001))) &
            (df_centroids["position_y"].apply(lambda val: math.isclose(val, y, abs_tol=0.0001))) &
            (df_centroids["position_z"].apply(lambda val: math.isclose(val, z, abs_tol=0.0001))) &
            (df_centroids["position_t"].apply(lambda val: math.isclose(val, t, abs_tol=0.0001)))
        ]
        assert len(corr_seg) > 0, f"Position of center segment corresponds to no tracked center, which is an error"
        assert len(corr_seg) <=1, f"Position of center segment corresponds to multiple tracked centers, which is an error"
        
        corr_seg=corr_seg["label"].values[0]
        assert corr_seg!=0, f"Position of center segment corresponds to background (0), which is an error"
        
        tcells_tracked[t,:,:,:][tcell_segments[t,:,:,:]==corr_seg]=row["TrackID"]
        # im_track = im_track[im==corr_seg]=row["TrackID"]
        
    imwrite(
        outpath,
        tcells_tracked
    ) 
    # return(tcells_tracked)