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
    columns=[],
    save=False
    ):
    
    sample_name="AIM_MB2_Exp58_Img003_donor899"
    
    output_dir = config['output_dir']
    img_dir = Path(output_dir, "images", sample_name)
    analysis_outdir = Path(output_dir, "analysis", "tcells")
    results_outdir = Path(analysis_outdir, "results")
    backproj_outdir = Path(analysis_outdir, "backprojection")
    
    if not backproj_outdir.exists():
        backproj_outdir.mkdir(parents=True)
        
    metadata = pd.read_csv(config["metadata_csv_path"])
    df_sample = metadata[metadata["sample_name"]==sample_name]
    assert(sample_name in df_sample["sample_name"].values), f"Supplied sample name {sample_name} not in metadata"
 
    track_img_path = Path(img_dir, f"{sample_name}_{cell_type}_tracked.tiff")
    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
        ]
   
    track_img = imread(track_img_path)
    
    # raw_h5_path = df_sample['image_path'].values[0]
    # raw_internal_path = df_sample["image_internal_path"].values[0]
        
    # raw_data = np.array(h5py.File(name=raw_h5_path, mode="r")[raw_internal_path])
    
    df_tracks_clustered=pd.read_csv(Path(results_outdir, "BEHAV3D_UMAP_clusters.csv"))
    track_img = np.where(np.isin(track_img, df_tracks_clustered["TrackID"].unique()), track_img, 0)

    trackid_data = {
        # "raw_data":{
        #     "img":raw_data, 
        #     "type":"image"
        #     },
        "TrackID":{
            "img":track_img, 
            "type":"label"
            }
        }
    if columns!=[]:
        columns=[x for x in df_tracks_clustered.columns if x not in metadata.columns.tolist()+["TrackID", "UMAP1", "UMAP2"]]

    backprojected_cols = backproject_columns(
        track_img,
        df_tracks_clustered,
        columns=columns
        )
    backproject_data = {**trackid_data, **backprojected_cols}
    
    view_napari(backproject_data, elsize)
    
    if save:
        # Output to .h5 format
        backproj_out_path = Path(backproj_outdir, f"{track_img_path.stem}_backprojected.h5")
        out_h5 = h5py.File(name=backproj_out_path, mode="w")
        

        for k,v in backproject_data.items():
            out_h5.create_dataset(name=k, data=v["img"], compression="gzip", compression_opts=1)
            out_h5[k].attrs.update({  
                    "type":v["type"],
                    "elsize":elsize
                    }
            )    
        out_h5.close()
    
def backproject_columns(
    track_img,
    df_tracks_clustered,
    columns=["ClusterID", "mean_speed"]
    ):
    
    mapped_imgs = {col:{"img":np.zeros_like(track_img),"type":"image"} for col in columns}
    for idx, row in df_tracks_clustered.iterrows():
        print(row["TrackID"])
        track_mask = track_img==row["TrackID"]
        for col in columns:
            mapped_imgs[col]["img"][track_mask]=row[col]
            if "ID" in col:
                mapped_imgs[col]["type"]="label"
    return(mapped_imgs)
        
def view_napari(
    backproject_data,
    elsize
    ):
    
    viewer=napari.Viewer()
    
    for k,v in backproject_data.items():
        if v["type"]=="label" or v["type"]=="segment":
            viewer.add_labels(v["img"], name=k, scale=elsize)
        else:
            viewer.add_image(v["img"], name=k, scale=elsize)
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