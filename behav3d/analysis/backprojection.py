import napari
import numpy as np
import pandas as pd
from tifffile import imread, imwrite
from pathlib import Path
import h5py
import yaml
import time
import argparse
from behav3d import format_time
# df_tracks=df_tracks[df_tracks["relative_time"]<=30]

def backproject_behav3d(
    metadata,
    sample_name,
    config=None,
    output_dir=None,
    cell_type = "tcells",
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
        
    analysis_outdir = Path(output_dir, "analysis", "tcells")
    results_outdir = Path(analysis_outdir, "results")
    backproj_outdir = Path(analysis_outdir, "backprojection")
    
    if not backproj_outdir.exists():
        backproj_outdir.mkdir(parents=True)
        
    df_sample = metadata[metadata["sample_name"]==sample_name]
    assert(sample_name in df_sample["sample_name"].values), f"Supplied sample name {sample_name} not in metadata"
 
    raw_img_path = Path(df_sample["raw_image_path"].values[0])
    track_img_path = Path(df_sample["tcell_segments_path"].values[0])

    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
        ]
    
    backproj_out_path = Path(backproj_outdir, f"{track_img_path.stem}_backprojected.h5")
    raw_img = imread(raw_img_path)
    
    raw_img_data = {
            "raw_data":{
                "img":raw_img,
                "type":"image"
                }
            }
    
    print("- Loading in tracked segments")
    track_img = imread(track_img_path)
    
    df_tracks_clustered=pd.read_csv(Path(results_outdir, "BEHAV3D_UMAP_clusters.csv"))
    track_img = np.where(np.isin(track_img, df_tracks_clustered["TrackID"].unique()), track_img, 0)

    trackid_data = {
        "TrackID":{
            "img":track_img, 
            "type":"label"
            }
        }
    
    if columns==[]:
        columns=[x for x in df_tracks_clustered.columns if x not in metadata.columns.tolist()+["TrackID", "UMAP1", "UMAP2"]]

    print("- Backprojecting all features onto each segment")
    backprojected_cols = backproject_columns(
        track_img,
        df_tracks_clustered,
        columns=columns
        )
    backproject_data = {**trackid_data, **backprojected_cols}
    visualize_data = {**raw_img_data, **backproject_data}
    
    print("- Visualizing backprojection in napari")
    view_napari(visualize_data, elsize)
    
    if save:
        # Output to .h5 format
        save_backprojection(
            backproj_out_path=backproj_out_path,
            backprojection_data=backproject_data,
            elsize=elsize
        )
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
    df_tracks_clustered,
    columns=["ClusterID", "mean_speed"]
    ):
    
    # mapped_imgs = {col:{"img":np.zeros_like(track_img),"type":"image"} for col in columns}
    # for idx, row in df_tracks_clustered.iterrows():
    #     track_mask = track_img==row["TrackID"]
    #     for col in columns:
    #         mapped_imgs[col]["img"][track_mask]=row[col]
    #         if "ID" in col:
    #             mapped_imgs[col]["type"]="label"
    mapped_imgs = {col:{"img":np.copy(track_img),"type":"image"} for col in columns}
    for col in columns:
        print(f"Backprojecting: {col}")
        col_dict = dict(zip(df_tracks_clustered['TrackID'], df_tracks_clustered[col]))
        mask = np.isin(mapped_imgs[col]["img"], list(col_dict.keys()))
        mapped_imgs[col]["img"][~mask]=0
        mapped_imgs[col]["img"][mask] = np.vectorize(col_dict.get)(mapped_imgs[col]["img"][mask])

    return(mapped_imgs)
        
def view_napari(
    backproject_data,
    elsize
    ):
    
    viewer=napari.Viewer()
    
    for idx, (k,v) in enumerate(backproject_data.items()):
        if v["type"]=="label" or v["type"]=="segment":
            viewer.add_labels(v["img"], name=k, scale=elsize)
        else:
            viewer.add_image(v["img"], name=k, scale=elsize)
    for lay in viewer.layers: lay.visible = False
    napari.run()  

def save_backprojection(
    backproj_out_path,
    backprojection_data,
    elsize
    ):
    print (f"--------------- Saving backprojection to {backproj_out_path} ---------------")
    start_time = time.time()
    out_h5 = h5py.File(name=backproj_out_path, mode="w")
        
    for k,v in backprojection_data.items():
        out_h5.create_dataset(name=k, data=v["img"], compression="gzip", compression_opts=1)
        out_h5[k].attrs.update({  
                "type":v["type"],
                "elsize":elsize
                }
        )    
    out_h5.close()
    end_time = time.time()
    h,m,s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config=yaml.load(parameters, Loader=yaml.SafeLoader)
    # metadata = pd.read_csv(config["metadata_csv_path"])
    backproject_behav3d(config)