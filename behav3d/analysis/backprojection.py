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

def backproject_mean_features_behav3d(
    metadata,
    sample_name,
    config=None,
    output_dir=None,
    cell_type="tcell",
    columns=[],
    save=False
    ):
    """
    Wrapper for mean mode backprojection. Calls generic backprojection function.
    """
    return backproject_generic(
        metadata=metadata,
        sample_name=sample_name,
        config=config,
        output_dir=output_dir,
        cell_type=cell_type,
        columns=columns,
        save=save,
        mode="mean"
    )


def backproject_time_features_behav3d(
    metadata,
    sample_name,
    config=None,
    output_dir=None,
    cell_type="tcell",
    columns=[],
    save=False
    ):
    """
    Wrapper for time mode backprojection. Calls generic backprojection function.
    """
    return backproject_generic(
        metadata=metadata,
        sample_name=sample_name,
        config=config,
        output_dir=output_dir,
        cell_type=cell_type,
        columns=columns,
        save=save,
        mode="time"
    )


def backproject_generic(
    metadata,
    sample_name,
    config=None,
    output_dir=None,
    cell_type="tcell",
    columns=[],
    save=False,
    mode="mean"
    ):
    """
    Generic backprojection function that works with ANY cell type.
    
    Args:
        metadata: Metadata DataFrame
        sample_name: Name of sample to backproject
        config: Config dict (optional)
        output_dir: Output directory
        cell_type: Cell type name (e.g., 'tcell', 'organoid1', 'macro')
        columns: List of columns to backproject
        save: Whether to save backprojection to zarr
        mode: 'mean' or 'time' backprojection mode
    """
    
    print(f"--------------- Backprojecting {cell_type} for {sample_name} ({mode} mode) ---------------")
    start_time = time.time()

    assert config is not None or all(
        [output_dir, metadata is not None]
    ), "Either 'config' or 'output_dir and metadata' parameters must be supplied"
    
    if not all([output_dir, metadata is not None]):
        output_dir = config['output_dir']
        metadata = pd.read_csv(config["metadata_csv_path"])
        
    analysis_outdir = Path(output_dir, "analysis", cell_type)
    results_outdir = Path(analysis_outdir, "results")
    backproj_outdir = Path(analysis_outdir, "backprojection")
    feature_outdir = Path(analysis_outdir, "track_features")
    img_outdir = Path(output_dir, "images", sample_name)
    
    if not backproj_outdir.exists():
        backproj_outdir.mkdir(parents=True)
        
    df_sample = metadata[metadata["sample_name"]==sample_name]
    assert(sample_name in df_sample["sample_name"].values), f"Supplied sample name {sample_name} not in metadata"
 
    raw_img_path = Path(df_sample["raw_image_path"].values[0])
    
    # ✨ Generic: Look for {cell_type}_tracks_image_path column (with or without category prefix)
    track_img_col = f"{cell_type}_tracks_image_path"
    if track_img_col not in df_sample.columns:
        # Try with prefix patterns (im_, or_, ot_)
        for prefix in ['im_', 'or_', 'ot_']:
            alt_col = f"{prefix}{cell_type}_tracks_image_path"
            if alt_col in df_sample.columns:
                track_img_col = alt_col
                break
        else:
            raise ValueError(f"Cannot find tracks image path column for {cell_type}. Expected: {cell_type}_tracks_image_path or {prefix}{cell_type}_tracks_image_path")
    
    track_img_path = Path(df_sample[track_img_col].values[0])
    
    # Load UMAP clusters to get ClusterID
    df_umap_clusters_path = Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv")
    df_umap_clusters = pd.read_csv(df_umap_clusters_path)
    
    if mode == "mean":
        # For mean mode: load summarized features (track-averaged)
        df_tracks_summarized_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features_summarized.csv")
        df_tracks_full = pd.read_csv(df_tracks_summarized_path)
        df_tracks_full = df_tracks_full[df_tracks_full["sample_name"]==sample_name]
        
        # Merge with ClusterID for backprojection
        df_tracks_clustered = pd.merge(df_tracks_full, df_umap_clusters[["TrackID", "ClusterID"]], on='TrackID', how='left')
        
        # Also load time-varying data for napari track positions
        df_tracks_all_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
        df_tracks_all = pd.read_csv(df_tracks_all_path)
        df_tracks_all = df_tracks_all[df_tracks_all["sample_name"]==sample_name]
    else:
        # For time mode: load time-varying features
        df_tracks_all_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
        df_tracks_all = pd.read_csv(df_tracks_all_path)
        df_tracks_all = df_tracks_all[df_tracks_all["sample_name"]==sample_name]
        
        # Merge with ClusterID for backprojection
        df_tracks_clustered = pd.merge(df_tracks_all, df_umap_clusters[["TrackID", "ClusterID"]], on='TrackID', how='left')
        
        # For time mode, df_tracks_full is same as df_tracks_all
        df_tracks_full = df_tracks_all
    
    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
    ]
    
    backproj_out_path = Path(backproj_outdir, f"{get_filepath_stem(track_img_path)}_backprojected.zarr")
    raw_img = load_image(raw_img_path)
    
    raw_img_data = {
        "raw_data": {
            "img": raw_img,
            "type": "image"
        }
    }
    
    print(f"- Loading in tracked segments for {cell_type}")
    track_img = load_image(track_img_path)
    
    # Expand dims and tile
    track_img = np.expand_dims(track_img, axis=1)
    filt_track_img = np.where(np.isin(track_img, df_tracks_clustered["TrackID"].unique()), track_img, 0)
    track_img_tiled = da.tile(track_img, (1, raw_img.shape[-4], 1, 1, 1))
    
    # Build trackid_data dict
    trackid_data = {
        f"{cell_type}_TrackID": {
            "img": track_img_tiled,
            "type": "label"
        },
        f"filtered_{cell_type}_TrackID": {
            "img": filt_track_img,
            "type": "label"
        }
    }
    
    # ✨ Dynamically load ALL other cell types' tracks from metadata
    from behav3d.utils import (
        detect_immune_cell_types_from_metadata,
        detect_organoid_types_from_metadata,
        detect_other_cell_types_from_metadata
    )
    
    all_cell_types = (detect_immune_cell_types_from_metadata(metadata) + 
                      detect_organoid_types_from_metadata(metadata) + 
                      detect_other_cell_types_from_metadata(metadata))
    
    # Load tracks for other cell types (for visualization context)
    for other_ct in all_cell_types:
        if other_ct == cell_type:
            continue  # Skip the main cell type we're backprojecting
        
        other_track_col = f"{other_ct}_tracks_image_path"
        # Try with prefixes
        if other_track_col not in df_sample.columns:
            for prefix in ['im_', 'or_', 'ot_']:
                alt_col = f"{prefix}{other_ct}_tracks_image_path"
                if alt_col in df_sample.columns:
                    other_track_col = alt_col
                    break
        
        if other_track_col in df_sample.columns and pd.notna(df_sample[other_track_col].values[0]):
            try:
                other_track_path = Path(df_sample[other_track_col].values[0])
                if other_track_path.exists():
                    other_track_img = load_image(other_track_path)
                    other_track_img = np.expand_dims(other_track_img, axis=1)
                    other_track_img = da.tile(other_track_img, (1, raw_img.shape[-4], 1, 1, 1))
                    
                    trackid_data[f"{other_ct}_TrackID"] = {
                        "img": other_track_img,
                        "type": "label"
                    }
            except Exception as e:
                print(f"  ⚠️  Could not load {other_ct} tracks: {e}")
    
    # Load dead mask if available
    dead_mask_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
    if dead_mask_path.exists():
        try:
            dead_mask = load_zarr(dead_mask_path)
            dead_mask = np.expand_dims(dead_mask, axis=1)
            dead_mask = da.tile(dead_mask, (1, raw_img.shape[-4], 1, 1, 1))
            trackid_data["Dead_Mask"] = {
                "img": dead_mask,
                "type": "label"
            }
        except Exception as e:
            print(f"  ⚠️  Could not load dead mask: {e}")
    
    # Select columns to backproject
    if columns == []:
        columns = [x for x in df_tracks_clustered.columns 
                   if x not in metadata.columns.tolist() + ["TrackID", "UMAP1", "UMAP2"]]
    
    # DON'T strip 'mean_' prefix - both time and mean CSVs have mean_intensity_ch* columns
    # (not intensity_ch* as we initially thought)
    
    # Check if dead channel exists
    has_dead_channel = "dead_channel" in df_sample.columns and pd.notna(df_sample["dead_channel"].values[0])
    
    # ALWAYS filter columns (even if pre-selected in UI)
    valid_columns = []
    for col in columns:
        # Skip if column doesn't exist in dataframe
        if col not in df_tracks_clustered.columns:
            print(f"  Skipping missing column: {col}")
            continue
        
        # Skip string columns
        if col.startswith('touching_') or col == 'orientation_vector':
            print(f"  Skipping string column: {col}")
            continue
        
        # Skip dead features if no dead channel
        if not has_dead_channel and ('dead' in col.lower() or col == 'is_dead'):
            print(f"  Skipping dead feature (no dead channel): {col}")
            continue
        
        # Only include numeric columns
        if pd.api.types.is_numeric_dtype(df_tracks_clustered[col]):
            valid_columns.append(col)
        else:
            print(f"  Skipping non-numeric column: {col}")
    
    columns = valid_columns
    print(f"- Backprojecting {len(columns)} features onto each segment")
    backprojected_cols = backproject_columns(
        track_img=filt_track_img,
        df_tracks_clustered=df_tracks_clustered,
        zarr_outpath=backproj_out_path,
        columns=columns
    )
    
    # Mark ClusterID as label type
    label_columns = ["ClusterID"]
    for col in backprojected_cols.keys():
        backprojected_cols[col]["img"] = da.tile(backprojected_cols[col]["img"], (1, raw_img.shape[-4], 1, 1, 1))
        if col in label_columns:
            backprojected_cols[col]["type"] = "label"
    
    trackid_data[f"filtered_{cell_type}_TrackID"]["img"] = da.tile(
        trackid_data[f"filtered_{cell_type}_TrackID"]["img"], (1, raw_img.shape[-4], 1, 1, 1))
    
    backproject_data = {**trackid_data, **backprojected_cols}
    visualize_data = {**raw_img_data, **backproject_data}
    
    print("- Visualizing backprojection in napari")
    viewer = view_napari_generic(
        visualize_data,
        df_tracks_full=df_tracks_all,
        df_tracks_clustered=df_tracks_clustered,
        cell_type=cell_type,
        elsize=elsize
    )
    
    end_time = time.time()
    h, m, s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    
    print(f"\n✅ Backprojection finished (napari was launched inside the function).")
    return viewer


def load_backprojection_h5(backprojection_h5_path):
    data_dict = {}
    with h5py.File(backprojection_h5_path, 'r') as backproj_h5:
        for dataset_name, dataset in backproj_h5.items():
            if 'type' in dataset.attrs:
                data_type=dataset.attrs['type']
            else:
                data_type="image"
            data_dict[dataset_name] = {
                'img': dataset[:],
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
    """Legacy function for backward compatibility"""
    return view_napari_generic(
        backproject_data,
        df_tracks_full,
        df_tracks_clustered,
        cell_type="tcell",
        elsize=elsize
    )


def view_napari_generic(
    backproject_data,
    df_tracks_full,
    df_tracks_clustered,
    cell_type,
    elsize
    ):
    """
    Generic napari viewer for any cell type.
    
    Args:
        backproject_data: Dict of data layers to visualize
        df_tracks_full: DataFrame with all tracks (time-varying, for positions)
        df_tracks_clustered: DataFrame with clustered tracks (may be summarized or time-varying)
        cell_type: Cell type being visualized
        elsize: Element size for scaling
    """
    
    viewer = napari.Viewer()
    
    for idx, (k, v) in enumerate(backproject_data.items()):
        v["img"] = np.transpose(v["img"], (1, 0, 2, 3, 4))
        if v["type"] == "label" or v["type"] == "segment":
            viewer.add_labels(v["img"], name=k, scale=elsize)
        else:
            viewer.add_image(v["img"], name=k, scale=elsize, colormap='inferno')
    
    for lay in viewer.layers:
        lay.visible = False
    
    # For napari tracks: need time-varying positions with ClusterID
    # df_tracks_full always has positions, df_tracks_clustered has ClusterID
    df_for_napari = pd.merge(df_tracks_full, df_tracks_clustered[["TrackID", "ClusterID"]], on='TrackID', how='left')
    
    napari_tracks_clustered = df_for_napari[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()
    features = {
        'ClusterID': df_for_napari["ClusterID"].to_numpy()
    }
    napari_tracks_full = df_tracks_full[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()

    viewer.add_tracks(napari_tracks_clustered, name=f'Filtered {cell_type} Tracks', 
                      properties=features, features=features, color_by="ClusterID", tail_length=75)
    viewer.add_tracks(napari_tracks_full, name=f'All {cell_type} Tracks', tail_length=75)
    napari.run()
    return viewer
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Input parameters for automatic data transfer.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file that stores all required paths', required=False)
    args = parser.parse_args()
    with open(args.config, "r") as parameters:
        config = yaml.load(parameters, Loader=yaml.SafeLoader)
    # Use generic backprojection
    backproject_generic(config=config)
