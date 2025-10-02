from tifffile import imread, imwrite
from skimage.measure import regionprops_table
import pandas as pd
import numpy as np
from pathlib import Path
import math
from behav3d.utils.fileio import load_image, save_as_zarr, append_to_zarr
import shutil
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
import multiprocessing

def convert_segments_to_tracks(
    df_tracks,
    segments,
    outpath,
    n_workers=None
    ):
    # tracked_img = np.zeros_like(segments)
    if n_workers is None:
        n_workers = multiprocessing.cpu_count()
        
    outpath = Path(outpath)
    if outpath.exists():
        if outpath.is_dir():
            # Remove the directory if it already exists
            shutil.rmtree(outpath)
        else:
            # Remove the file if it already exists
            outpath.unlink()
        # Remove the file if it already exists
        
    if outpath.suffix == ".zip":
        outpath = Path(outpath.parent, outpath.stem)
    
    assert outpath.suffix == ".zarr", "Supplied outpath is not .zarr or .zarr.zip"
    
    def _apply_segment_to_track(args):
        t_segment_img, t_df_tracks, t, track_out_path = args
        
        t_segment_img = np.asarray(t_segment_img)
        t_tracked_img = np.zeros_like(t_segment_img)
        
        for _, row in t_df_tracks.iterrows():
            # print(row["SegmentID"], row["TrackID"], (tracked_img[t]==row["SegmentID"]).any())
            t_tracked_img[t_segment_img==row["SegmentID"]] = row["TrackID"]
        return t_tracked_img
        
    args_list = [(t_segments, df_tracks[df_tracks["position_t"]==t] , t, outpath) for t, t_segments in enumerate(segments)]
    result=[]
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        result+=list(tqdm(executor.map(_apply_segment_to_track, args_list), total=len(args_list)))
    
    tracked_img = np.stack(result, axis=0)
    save_as_zarr(tracked_img, path=outpath)
    # for t, t_seg in tqdm(enumerate(segments), total=len(segments)):
    #     t_seg = np.asarray(t_seg)
    #     tracked_img = np.zeros_like(t_seg)
    #     t_df_tracks = df_tracks[df_tracks["position_t"]==t]
    #     for _, row in t_df_tracks.iterrows():
    #         # print(row["SegmentID"], row["TrackID"], (tracked_img[t]==row["SegmentID"]).any())
    #         tracked_img[t_seg==row["SegmentID"]] = row["TrackID"]

    #     tracked_img = np.expand_dims(tracked_img, axis=0)
    #     append_to_zarr(
    #         img=tracked_img, 
    #         outpath=outpath
    #         )
    return(outpath)

def convert_all_tracked_images_to_csv(
    metadata,
    output_dir,
    cell_type="organoid",
    overwrite=False,
    **kwargs
    ):
    for idx, sample in metadata.iterrows():
        sample_name=sample['sample_name']
        print(f"Tracking sample: {sample_name}")
        
        tracked_img_outdir = Path(output_dir, "images", sample_name)
        tracked_csv_outdir = Path(output_dir, "trackdata", sample_name, cell_type)
        
        # segmented_img_path = Path(tracked_img_outdir, f"{sample_name}_{cell_type}_segments.zarr")
        tracked_img_outpath = Path(tracked_img_outdir, f"{sample_name}_{cell_type}_tracked.zarr")
        segmented_img_path = tracked_img_outpath
        tracked_csv_outpath = Path(tracked_csv_outdir, f"{sample_name}_{cell_type}_tracks.csv")
        
        if not tracked_img_outdir.exists():
            tracked_img_outdir.mkdir(parents=True)
        if not tracked_csv_outdir.exists():
            tracked_csv_outdir.mkdir(parents=True)
        
        if (
            (
                not tracked_csv_outpath.exists() or 
                not tracked_img_outpath.exists()
            ) or overwrite
            ):
            element_size_x = sample["pixel_distance_xy"]
            element_size_y = sample["pixel_distance_xy"]
            element_size_z = sample["pixel_distance_z"]
            
            df_tracks = convert_tracked_image_to_csv(
                img_path=segmented_img_path,
                outpath=tracked_csv_outpath,
                element_size_x=element_size_x,
                element_size_y=element_size_y,
                element_size_z=element_size_z
            )
        else:
            print("Tracking already exists... Provide overwrite=True to overwrite... Loading existing tracking data")
            
        metadata.at[idx, f"{cell_type}_tracks_image_path"] = str(tracked_img_outpath)
        metadata.at[idx, f"{cell_type}_tracks_csv_path"] = str(tracked_csv_outpath)
    return(metadata)

def convert_tracked_image_to_csv(
    img_path,
    outpath=None,
    element_size_z=1,
    element_size_y=1,
    element_size_x=1
    ):
    segments = load_image(img_path)
    df_tracks = []
    for t, t_seg in tqdm(enumerate(segments), total=len(segments)):
        t_seg = np.asarray(t_seg)
        properties=pd.DataFrame(regionprops_table(label_image=t_seg, properties=['label', f'centroid']))
        properties["position_t"]=t
        df_tracks.append(properties)
    df_tracks = pd.concat(df_tracks)
    df_tracks.rename(
            columns={
                    'label': 'SegmentID',
                    'centroid-0': "pixel_position_z",
                    'centroid-1': "pixel_position_y",
                    'centroid-2': "pixel_position_x",
                }, 
            inplace=True
        )
    df_tracks["TrackID"]=df_tracks["SegmentID"]
    df_tracks["position_z"]=df_tracks["pixel_position_z"]*element_size_z
    df_tracks["position_y"]=df_tracks["pixel_position_y"]*element_size_y
    df_tracks["position_x"]=df_tracks["pixel_position_x"]*element_size_x
    df_tracks = df_tracks[["TrackID", "SegmentID", "position_t", "position_x", "position_y", "position_z", "pixel_position_x", "pixel_position_y", "pixel_position_z"]]
    if outpath is not None:
        outpath = Path(outpath)
        df_tracks.to_csv(outpath, sep=",", index=False)
        
    return(df_tracks)
        
# def convert_segments_to_tracks(
#         tracks_csv_path,
#         segments_path,
#         outpath,
#         element_size_z,
#         element_size_y,
#         element_size_x
#     ):
    
#     ### Assign the tracks to existing segments
#     # Loop through spots, link to segments in the image and replace label with TrackID
#     print("- Assigning track ID to segmented image to create tracked image...")
#     df_tracks = pd.read_csv(tracks_csv_path)
#     tcell_segments = imread(segments_path)
    
#     df_centroids = []
#     # TODO 
#     # Check if trackmate removes objects and still keep them in image
#     for t, tcell_stack in enumerate(tcell_segments):
#         properties=pd.DataFrame(regionprops_table(label_image=tcell_stack, properties=['label', f'centroid']))
#         properties["position_t"]=t
#         df_centroids.append(properties)
#     df_centroids = pd.concat(df_centroids)
#     df_centroids["position_z"]=df_centroids["centroid-0"]*element_size_z
#     df_centroids["position_y"]=df_centroids["centroid-1"]*element_size_y
#     df_centroids["position_x"]=df_centroids["centroid-2"]*element_size_x
    
#     tcells_tracked = np.zeros_like(tcell_segments)
#     for _, row in df_tracks.iterrows():
#         t,z,y,x = int(row["position_t"]),row["position_z"], row["position_y"], row["position_x"]
#         corr_seg=None
#         # There seems to be an issue where the tracking output gives 
#         # e.g. 25.99355 and regionprops * elemen_size gives 25.99354999999999 (Floating point arithmetics issue)
#         # To solve we use the math.isclose() with an abs_tol of 0.0001 to still match these values
#         corr_seg = df_centroids[
#             (df_centroids["position_x"].apply(lambda val: math.isclose(val, x, abs_tol=0.0001))) &
#             (df_centroids["position_y"].apply(lambda val: math.isclose(val, y, abs_tol=0.0001))) &
#             (df_centroids["position_z"].apply(lambda val: math.isclose(val, z, abs_tol=0.0001))) &
#             (df_centroids["position_t"].apply(lambda val: math.isclose(val, t, abs_tol=0.0001)))
#         ]
#         assert len(corr_seg) > 0, f"Position of center segment corresponds to no tracked center, which is an error"
#         assert len(corr_seg) <=1, f"Position of center segment corresponds to multiple tracked centers, which is an error"
        
#         corr_seg=corr_seg["label"].values[0]
#         assert corr_seg!=0, f"Position of center segment corresponds to background (0), which is an error"
        
#         tcells_tracked[t,:,:,:][tcell_segments[t,:,:,:]==corr_seg]=row["TrackID"]
#         # im_track = im_track[im==corr_seg]=row["TrackID"]
        
#     imwrite(
#         outpath,
#         tcells_tracked
#     ) 
#     # return(tcells_tracked)