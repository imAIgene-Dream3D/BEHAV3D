import napari
import numpy as np
import pandas as pd
from tifffile import imread, imwrite
from pathlib import Path
import h5py
import yaml
# df_tracks=df_tracks[df_tracks["relative_time"]<=30]

def backproject_data(
    config,
    sample_name,
    cell_type = "tcells",
    backproject_columns=[
        "TrackID",
        "cluster"
        ]
    ):
    
    sample_name="AIM_MB2_Exp58_Img003_donor899"
    
    output_dir = config['output_dir']
    
    metadata = pd.read_csv(config["metadata_csv_path"])
    df_sample = metadata[metadata["sample_name"]==sample_name]
    assert(sample_name in df_sample["sample_name"].values), f"Supplied sample name {sample_name} not in metadata"
 
    track_img_path = Path(output_dir, f"{sample_name}_{cell_type}_tracked.tiff")
    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
        ]
   
    track_img = imread(track_img_path)
    df_tracks_clustered=pd.read_csv(Path(output_dir, "BEHAV3D_UMAP_clusters.csv"))
    
    backproject_data = backproject_columns(
        track_img,
        df_tracks_clustered,
        columns=[
            "cluster",
            "mean_speed"
            ]
        )
    backproject_data = {**{"TrackID": track_img}, **backproject_data}
    
    # Output to .h5 format
    backproj_out_path = track_img_path.with_name(f"{track_img_path.stem}_backprojected.h5")
    out_h5 = h5py.File(name=backproj_out_path, mode="w")
    
    for k,v in backproject_data.items():
        out_h5.create_dataset(name=k, data=v)
        out_h5[k].attrs.update({  
                "type":"image",
                "elsize":elsize
                }
        )    
    out_h5.close()
    
    backproject_data
    view_napari(images=image_list, labels=["label"], names=["TrackID"])
    
# def backproject_data(
#     track_img_path,
#     elsize=[1,1,1,1]
#     df_tracks_clustered=None,
#     df_tracks_filtered=None,
#     backproj_out_path=None,
#     backproject_columns=[
#         "TrackID",
#         "cluster"
#         ]
#     ):
#     track_img_path=Path("/Users/samdeblank/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/BEHAV3D_run/ilastik_output/AIM_MB2_Exp58_Img003_donor899_tcells_tracked.tiff")
#     track_img_path=Path(track_img_path)
#     backproj_out_path = track_img_path.with_name(f"{track_img_path.stem}_backprojected.h5")

#     image_list = []
    
#     track_img = imread(track_img_path)
#     # if isinstance(segment_img, str):
#     #     segment_img = imread(segment_img)
#     #     segment_img = imread("/Users/samdeblank/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/BEHAV3D_run/ilastik_output/AIM_MB2_Exp58_Img003_donor899_tcells_tracked.tiff")
#     assert isinstance(track_img, np.ndarray), "segment_img should be type np.array or a path to a .tiff file"
    
#     df_tracks_clustered=pd.read_csv("/Users/samdeblank/Documents/1.projects/BHVD_BEHAV3D/BEHAV3D-ilastik/test/BEHAV3D_run/ilastik_output/BEHAV3D_UMAP_clusters.csv")
    
#     # Overwrite segment ID with cluster ID and add to image list
#     backproject_data = backproject_columns(
#         track_img,
#         df_tracks_clustered,
#         columns=[
#             "cluster",
#             "mean_speed"
#             ]
#         )
#     backproject_data = {**{"TrackID": track_img}, **backproject_data}
    
#     # Output to .h5 format
#     out_h5 = h5py.File(name=backproj_out_path, mode="w")
    
#     for k,v in backproject_data.items():
#         out_h5.create_dataset(name=k, data=v)
#         out_h5[k].attrs.update({  
#                 "type":"image",
#                 "elsize":
#                 }
#         )    
#     out_h5.close()
    
#     backproject_data
#     view_napari(images=image_list, labels=["label"], names=["TrackID"])

def backproject_columns(
    track_img,
    df_tracks_clustered,
    columns=["cluster", "mean_speed"]
    ):
    
    mapped_imgs = {col:np.zeros_like(track_img) for col in columns}
    for idx, row in df_tracks_clustered.iterrows():
        print(row["TrackID"])
        track_mask = track_img==row["TrackID"]
        for col in columns:
            mapped_imgs[col][track_mask]=row[col]
    return(mapped_imgs)
        
def view_napari(
    images=[],
    labels=[],
    names=[],
    elsize=[]
    ):
    
    viewer=napari.Viewer()
    
    if elsize is None or elsize==[]:
            elsize= [1] * images[0].ndim
            
    for idx, img in enumerate(images):
        if labels is None or labels==[]:
            label = None
        else:
            label=labels[idx]
        if names is None or labels==[]:
            name = None
        else:
            name=names[idx]
        if label=="label" or label=="segment":
            viewer.add_labels(img, name=name, scale=elsize)
        else:
            viewer.add_image(img, name=name, scale=elsize)
    for lay in viewer.layers: lay.visible = False
    napari.run()  
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    # metadata = pd.read_csv(config["metadata_csv_path"])
    run_tcell_analysis(config)