import napari
import numpy as np
import pandas as pd
from pathlib import Path
import h5py
import yaml
import time
import argparse
from behav3d.utils import format_time
from behav3d.utils.fileio import get_filepath_stem, load_image, load_zarr, save_as_zarr
from tqdm import tqdm
import dask.array as da
# df_tracks=df_tracks[df_tracks["relative_time"]<=30]

def backproject_behav3d(
    metadata,
    sample_name,
    config=None,
    output_dir=None,
    cell_type = "tcell",
    columns=[],
    save=False
    ):
    
    print(f"--------------- Backprojecting for {sample_name} ---------------")
    start_time = time.time()

    assert config is not None or all(
        [output_dir, metadata is not None]
    ), "Either 'config' or 'output_dir and metadata' parameters must be supplied"
    
    if not all([output_dir, metadata is not None]):
        output_dir = config['output_dir']
        metadata = pd.read_csv(config["metadata_csv_path"])
        
    analysis_outdir = Path(output_dir, "analysis", "tcell")
    results_outdir = Path(analysis_outdir, "results")
    backproj_outdir = Path(analysis_outdir, "backprojection")
    feature_outdir = Path(analysis_outdir, "track_features")
    img_outdir = Path(output_dir, "images", sample_name)
    
    if not backproj_outdir.exists():
        backproj_outdir.mkdir(parents=True)
        
    df_sample = metadata[metadata["sample_name"]==sample_name]
    assert(sample_name in df_sample["sample_name"].values), f"Supplied sample name {sample_name} not in metadata"
 
    raw_img_path = Path(df_sample["raw_image_path"].values[0])
    
    track_img_path = Path(df_sample["tcell_tracks_image_path"].values[0])
    
    org_track_img_path = Path(df_sample["organoid_tracks_image_path"].values[0])
    df_tracks_path = Path(df_sample["tcell_tracks_csv_path"].values[0])
    
    df_tracks_full_clustered_path = Path(results_outdir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv")
    df_tracks_full_clustered = pd.read_csv(df_tracks_full_clustered_path)
    df_tracks_full_clustered = df_tracks_full_clustered[df_tracks_full_clustered["sample_name"]==sample_name]
    
    df_tracks_all_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_tracks_all = pd.read_csv(df_tracks_all_path)
    df_tracks_all = df_tracks_all[df_tracks_all["sample_name"]==sample_name]
    
    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
        ]
    
    backproj_out_path = Path(backproj_outdir, f"{get_filepath_stem(track_img_path)}_backprojected.zarr")
    raw_img = load_image(raw_img_path)
    
    raw_img_data = {
            "raw_data":{
                "img":raw_img,
                "type":"image"
                }
            }
    
    print("- Loading in tracked segments")
    track_img = load_image(track_img_path)
    org_track_img = load_image(org_track_img_path)
    
    track_img = np.expand_dims(track_img, axis=1)
    org_track_img = np.expand_dims(org_track_img, axis=1)
    org_track_img = da.tile(org_track_img, (1,raw_img.shape[-4],1,1,1))
    
    dead_mask_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
    dead_mask = load_zarr(dead_mask_path)
    dead_mask = np.expand_dims(dead_mask, axis=1)
    dead_mask = da.tile(dead_mask, (1,raw_img.shape[-4],1,1,1))
    
    df_tracks_clustered=pd.read_csv(Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv"))
    filt_track_img = np.where(np.isin(track_img, df_tracks_clustered["TrackID"].unique()), track_img, 0)
    
    track_img = da.tile(track_img, (1,raw_img.shape[-4],1,1,1))
    
    trackid_data = {
        "Tcell_TrackID":{
            "img":track_img, 
            "type":"label"
            },
        "filtered_Tcell_TrackID":{
            "img":filt_track_img, 
            "type":"label"
            },
        "Organoid_TrackID":{
            "img":org_track_img, 
            "type":"label"
            },
        "Dead_Mask":{  
            "img":dead_mask, 
            "type":"label"
            }
        }
    
    if columns==[]:
        columns=[x for x in df_tracks_clustered.columns if x not in metadata.columns.tolist()+["TrackID", "UMAP1", "UMAP2"]]
        
    print("- Backprojecting all features onto each segment")
    backprojected_cols = backproject_columns(
        track_img=filt_track_img,
        df_tracks_clustered=df_tracks_clustered,
        zarr_outpath=backproj_out_path,
        columns=columns
        )
    
    label_columns = ["ClusterID"]
    for col in backprojected_cols.keys():
        ### Repeat the image so it fits with the raw image dimensionsw
        backprojected_cols[col]["img"] = da.tile(backprojected_cols[col]["img"], (1,raw_img.shape[-4],1,1,1))
        if col in label_columns:
            backprojected_cols[col]["type"] = "label"
    
    trackid_data["filtered_Tcell_TrackID"]["img"] = da.tile(trackid_data["filtered_Tcell_TrackID"]["img"], (1,raw_img.shape[-4],1,1,1))
    
    backproject_data = {**trackid_data, **backprojected_cols}
    visualize_data = {**raw_img_data, **backproject_data}
    
    print("- Visualizing backprojection in napari")
    view_napari(
        visualize_data, 
        df_tracks_full=df_tracks_all,
        df_tracks_clustered=df_tracks_full_clustered, 
        elsize=elsize)
    
    # if save:
    #     # Output to .h5 format
    #     save_backprojection(
    #         backproj_out_path=backproj_out_path,
    #         backprojection_data=backproject_data,
    #         elsize=elsize
    #     )
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return({
        "path": backproj_out_path,
        "data": visualize_data,
        "elsize": elsize
        })

def load_backprojection_h5(backprojection_h5_path):
    data_dict = {}
    with h5py.File(backprojection_h5_path, 'r') as backproj_h5:
        for dataset_name, dataset in backproj_h5.items():
            if 'type' in dataset.attrs:
                data_type=dataset.attrs['type']
            else:
                data_type="image"
            data_dict[dataset_name] = {
                'img': dataset[:],  # Assuming dataset is your actual data
                'type': data_type
            }
    return(data_dict)

def backproject_columns(
    track_img,
    zarr_outpath,
    df_tracks_clustered,
    columns=["ClusterID", "mean_speed"]
    ):
    zarr_outpath = Path(zarr_outpath)
    if Path(f"{zarr_outpath}.zip").exists():
        zarr_outpath = Path(f"{zarr_outpath}.zip")
        
    mapped_imgs = {col:{"img":None,"type":"image"} for col in columns}
    if zarr_outpath is not None:
        if zarr_outpath.exists():
            print("Backprojection .zarr already exists. Loading backprojection data...")
            for col in columns:
                mapped_imgs[col]["img"]=load_image(zarr_outpath, group=col)
        else:     
            if zarr_outpath.suffix==".zip":
                zarr_outpath = zarr_outpath.stem        
            for col in tqdm(columns, total=len(columns)):
                mapped_img = np.asarray(track_img.copy())
                col_dict = dict(zip(df_tracks_clustered['TrackID'], df_tracks_clustered[col]))
                mask = np.isin(mapped_img, list(col_dict.keys()))
                mapped_img[~mask]=0
                mapped_img[mask] = np.vectorize(col_dict.get)(mapped_img[mask])
                mapped_imgs[col]["img"]=mapped_img
                save_as_zarr(
                    img=mapped_img, 
                    path=zarr_outpath, 
                    chunks=(1,) + mapped_img.shape[1:],
                    group=col
                    )
        return(mapped_imgs)
        
def view_napari(
    backproject_data,
    df_tracks_full,
    df_tracks_clustered,
    elsize
    ):
    
    viewer=napari.Viewer()
    
    for idx, (k,v) in enumerate(backproject_data.items()):
        v["img"] = np.transpose(v["img"], (1,0,2,3,4))
        if v["type"]=="label" or v["type"]=="segment":
            viewer.add_labels(v["img"], name=k, scale=elsize)
        else:
            viewer.add_image(v["img"], name=k, scale=elsize, colormap='inferno')
    for lay in viewer.layers: lay.visible = False
    
    napari_tracks_clustered = df_tracks_clustered[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()
    features = {
        'ClusterID': df_tracks_clustered["ClusterID"].to_numpy()
    }
    napari_tracks_full = df_tracks_full[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()

    viewer.add_tracks(napari_tracks_clustered, name='Filtered T-cell Tracks', properties=features, features=features, color_by="ClusterID", tail_length=75)
    viewer.add_tracks(napari_tracks_full, name='All T-cell Tracks',  tail_length=75)
    napari.run()  
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    # metadata = pd.read_csv(config["metadata_csv_path"])
    backproject_behav3d(config)