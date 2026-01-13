import napari
import numpy as np
import pandas as pd
from pathlib import Path
import yaml
import time
import argparse
from behav3d.core.utils import format_time
from behav3d.io.images import get_filepath_stem, load_image, load_zarr, save_as_zarr
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
    Backproject mean (summarized) features onto tracked segments.
    Features are averaged per track, so each track has one value per feature.
    """
    
    print(f"--------------- Backprojecting {cell_type} for {sample_name} (mean mode) ---------------")
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
    
    # Find tracks image path column dynamically
    track_img_col = _find_track_image_column(df_sample, cell_type)
    track_img_path = Path(df_sample[track_img_col].values[0])
    
    # Load clustered data (summarized - one row per track)
    df_tracks_clustered = pd.read_csv(Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv"))
    
    # Load full track data (time-varying - for napari track positions)
    df_tracks_all_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_tracks_all = pd.read_csv(df_tracks_all_path)
    df_tracks_all = df_tracks_all[df_tracks_all["sample_name"]==sample_name]
    
    # Load clustered full data (time-varying with ClusterID)
    df_tracks_full_clustered_path = Path(results_outdir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv")
    df_tracks_full_clustered = pd.read_csv(df_tracks_full_clustered_path)
    df_tracks_full_clustered = df_tracks_full_clustered[df_tracks_full_clustered["sample_name"]==sample_name]
    
    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
    ]
    
    # Use separate zarr file for mean mode
    backproj_out_path = Path(backproj_outdir, f"{get_filepath_stem(track_img_path)}_backprojected_mean.zarr")
    raw_img = load_image(raw_img_path)
    
    raw_img_data = {
        "raw_data": {
            "img": raw_img,
            "type": "image"
        }
    }
    
    print("- Loading in tracked segments")
    track_img = load_image(track_img_path)
    
    track_img = np.expand_dims(track_img, axis=1)
    
    # Filter to only tracks with ClusterID
    filt_track_img = np.where(np.isin(track_img, df_tracks_clustered["TrackID"].unique()), track_img, 0)
    
    # Check filtered result (handle dask arrays)
    if hasattr(filt_track_img, 'compute'):
        unique_filtered = da.unique(filt_track_img).compute()
    else:
        unique_filtered = np.unique(filt_track_img)
    unique_filtered = unique_filtered[unique_filtered != 0]
    print(f"  Found {len(unique_filtered)} matching TrackIDs for filtered view")
    
    track_img = da.tile(track_img, (1, raw_img.shape[-4], 1, 1, 1))
    
    # Build trackid_data
    trackid_data = {
        f"{cell_type}_TrackID": {
            "img": track_img,
            "type": "label"
        },
        f"filtered_{cell_type}_TrackID": {
            "img": filt_track_img,
            "type": "label"
        }
    }
    
    # Load other cell type tracks for context
    trackid_data = _load_other_cell_tracks(trackid_data, metadata, df_sample, cell_type, raw_img)
    
    # Load dead mask if available
    trackid_data = _load_dead_mask(trackid_data, df_sample, img_outdir, sample_name, raw_img)
    
    # Select columns to backproject
    if columns == []:
        columns = [x for x in df_tracks_clustered.columns 
                   if x not in metadata.columns.tolist() + ["TrackID", "UMAP1", "UMAP2"]]
    
    # Always ensure ClusterID is included
    if "ClusterID" not in columns and "ClusterID" in df_tracks_clustered.columns:
        columns = ["ClusterID"] + columns
    
    # Filter to valid numeric columns
    columns = _filter_valid_columns(columns, df_tracks_clustered, df_sample)
    
    print(f"- Backprojecting {len(columns)} features onto each segment")
    backprojected_cols = backproject_columns(
        track_img=filt_track_img,
        df_tracks_clustered=df_tracks_clustered,
        zarr_outpath=backproj_out_path,
        columns=columns
    )
    
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
    viewer = view_napari(
        visualize_data,
        df_tracks_full=df_tracks_all,
        df_tracks_clustered=df_tracks_full_clustered,
        cell_type=cell_type,
        elsize=elsize
    )
    
    end_time = time.time()
    h, m, s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return viewer


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
    Backproject time-varying features onto tracked segments.
    Features vary per timepoint, so each track point has its own value.
    """
    
    print(f"--------------- Backprojecting {cell_type} for {sample_name} (time mode) ---------------")
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
    
    # Find tracks image path column dynamically
    track_img_col = _find_track_image_column(df_sample, cell_type)
    track_img_path = Path(df_sample[track_img_col].values[0])
    
    # Load UMAP clusters (to get ClusterID mapping)
    df_umap_clusters = pd.read_csv(Path(results_outdir, f"BEHAV3D_{cell_type}_UMAP_clusters.csv"))
    
    # Load full track data (time-varying)
    df_tracks_all_path = Path(feature_outdir, f"BEHAV3D_{cell_type}_combined_track_features.csv")
    df_tracks_all = pd.read_csv(df_tracks_all_path)
    df_tracks_all = df_tracks_all[df_tracks_all["sample_name"]==sample_name]
    
    # Merge ClusterID into time-varying data
    df_tracks_clustered = pd.merge(
        df_tracks_all, 
        df_umap_clusters[["TrackID", "ClusterID"]], 
        on='TrackID', 
        how='inner'  # Only keep tracks that have ClusterID
    )
    
    # Load clustered full data for napari
    df_tracks_full_clustered_path = Path(results_outdir, f"BEHAV3D_{cell_type}_combined_track_features_clustered.csv")
    df_tracks_full_clustered = pd.read_csv(df_tracks_full_clustered_path)
    df_tracks_full_clustered = df_tracks_full_clustered[df_tracks_full_clustered["sample_name"]==sample_name]
    
    elsize = [
        df_sample["pixel_distance_z"].values[0],
        df_sample["pixel_distance_xy"].values[0],
        df_sample["pixel_distance_xy"].values[0]
    ]
    
    # Use separate zarr file for time mode
    backproj_out_path = Path(backproj_outdir, f"{get_filepath_stem(track_img_path)}_backprojected_time.zarr")
    raw_img = load_image(raw_img_path)
    
    raw_img_data = {
        "raw_data": {
            "img": raw_img,
            "type": "image"
        }
    }
    
    print("- Loading in tracked segments")
    track_img = load_image(track_img_path)
    
    track_img = np.expand_dims(track_img, axis=1)
    
    # Filter to only tracks with ClusterID
    filt_track_img = np.where(np.isin(track_img, df_tracks_clustered["TrackID"].unique()), track_img, 0)
    
    # Check filtered result (handle dask arrays)
    if hasattr(filt_track_img, 'compute'):
        unique_filtered = da.unique(filt_track_img).compute()
    else:
        unique_filtered = np.unique(filt_track_img)
    unique_filtered = unique_filtered[unique_filtered != 0]
    print(f"  Found {len(unique_filtered)} matching TrackIDs for filtered view")
    
    track_img = da.tile(track_img, (1, raw_img.shape[-4], 1, 1, 1))
    
    # Build trackid_data
    trackid_data = {
        f"{cell_type}_TrackID": {
            "img": track_img,
            "type": "label"
        },
        f"filtered_{cell_type}_TrackID": {
            "img": filt_track_img,
            "type": "label"
        }
    }
    
    # Load other cell type tracks for context
    trackid_data = _load_other_cell_tracks(trackid_data, metadata, df_sample, cell_type, raw_img)
    
    # Load dead mask if available
    trackid_data = _load_dead_mask(trackid_data, df_sample, img_outdir, sample_name, raw_img)
    
    # Select columns to backproject
    if columns == []:
        columns = [x for x in df_tracks_clustered.columns 
                   if x not in metadata.columns.tolist() + ["TrackID", "UMAP1", "UMAP2"]]
    
    # Always ensure ClusterID is included
    if "ClusterID" not in columns and "ClusterID" in df_tracks_clustered.columns:
        columns = ["ClusterID"] + columns
    
    # Filter to valid numeric columns
    columns = _filter_valid_columns(columns, df_tracks_clustered, df_sample)
    
    print(f"- Backprojecting {len(columns)} features onto each segment")
    backprojected_cols = backproject_columns(
        track_img=filt_track_img,
        df_tracks_clustered=df_tracks_clustered,
        zarr_outpath=backproj_out_path,
        columns=columns
    )
    
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
    viewer = view_napari(
        visualize_data,
        df_tracks_full=df_tracks_all,
        df_tracks_clustered=df_tracks_full_clustered,
        cell_type=cell_type,
        elsize=elsize
    )
    
    end_time = time.time()
    h, m, s = format_time(start_time, end_time)
    print(f"### DONE - elapsed time: {h}:{m:02}:{s:02}\n")
    return viewer


def _find_track_image_column(df_sample, cell_type):
    """Find the correct column name for track image path."""
    track_img_col = f"{cell_type}_tracks_image_path"
    if track_img_col not in df_sample.columns:
        for prefix in ['im_', 'or_', 'ot_']:
            alt_col = f"{prefix}{cell_type}_tracks_image_path"
            if alt_col in df_sample.columns:
                return alt_col
        raise ValueError(f"Cannot find tracks image path column for {cell_type}")
    return track_img_col


def _load_other_cell_tracks(trackid_data, metadata, df_sample, cell_type, raw_img):
    """Load tracks from other cell types for visualization context."""
    from behav3d.core.metadata import (
        detect_immune_cell_types_from_metadata,
        detect_organoid_types_from_metadata,
        detect_other_cell_types_from_metadata
    )
    
    all_cell_types = (detect_immune_cell_types_from_metadata(metadata) + 
                      detect_organoid_types_from_metadata(metadata) + 
                      detect_other_cell_types_from_metadata(metadata))
    
    for other_ct in all_cell_types:
        if other_ct == cell_type:
            continue
        
        try:
            other_track_col = _find_track_image_column(df_sample, other_ct)
            if other_track_col in df_sample.columns and pd.notna(df_sample[other_track_col].values[0]):
                other_track_path = Path(df_sample[other_track_col].values[0])
                if other_track_path.exists():
                    other_track_img = load_image(other_track_path)
                    other_track_img = np.expand_dims(other_track_img, axis=1)
                    other_track_img = da.tile(other_track_img, (1, raw_img.shape[-4], 1, 1, 1))
                    trackid_data[f"{other_ct}_TrackID"] = {
                        "img": other_track_img,
                        "type": "label"
                    }
        except (ValueError, Exception):
            pass
    
    return trackid_data


def _load_dead_mask(trackid_data, df_sample, img_outdir, sample_name, raw_img):
    """Load dead mask if dead channel is enabled."""
    has_dead_channel = "dead_channel" in df_sample.columns and pd.notna(df_sample["dead_channel"].values[0])
    
    if has_dead_channel:
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
            except Exception:
                pass
    
    return trackid_data


def _filter_valid_columns(columns, df_tracks_clustered, df_sample):
    """Filter columns to only valid numeric ones."""
    has_dead_channel = "dead_channel" in df_sample.columns and pd.notna(df_sample["dead_channel"].values[0])
    
    valid_columns = []
    for col in columns:
        if col not in df_tracks_clustered.columns:
            print(f"  Skipping missing column: {col}")
            continue
        
        if col.startswith('touching_') or col == 'orientation_vector':
            print(f"  Skipping string column: {col}")
            continue
        
        if not has_dead_channel and ('dead' in col.lower() or col == 'is_dead'):
            continue
        
        if pd.api.types.is_numeric_dtype(df_tracks_clustered[col]):
            valid_columns.append(col)
    
    return valid_columns


def backproject_columns(
    track_img,
    zarr_outpath,
    df_tracks_clustered,
    columns=["ClusterID", "mean_speed"]
    ):
    """
    Map column values onto track image pixels.
    Each pixel value (TrackID) is replaced with the corresponding column value.
    """
    import zarr
    
    zarr_outpath = Path(zarr_outpath)
    if Path(f"{zarr_outpath}.zip").exists():
        zarr_outpath = Path(f"{zarr_outpath}.zip")
        
    mapped_imgs = {col: {"img": None, "type": "image"} for col in columns}
    
    if zarr_outpath is not None:
        if zarr_outpath.exists():
            # Check which columns already exist in zarr
            try:
                zarr_store = zarr.open(str(zarr_outpath), mode='r')
                existing_groups = list(zarr_store.keys()) if hasattr(zarr_store, 'keys') else []
            except Exception:
                existing_groups = []
            
            cols_to_load = [c for c in columns if c in existing_groups]
            cols_to_compute = [c for c in columns if c not in existing_groups]
            
            if cols_to_load:
                print(f"  Loading {len(cols_to_load)} features from existing .zarr")
                for col in cols_to_load:
                    mapped_imgs[col]["img"] = load_image(zarr_outpath, group=col)
            
            if cols_to_compute:
                print(f"  Computing {len(cols_to_compute)} new features")
                zarr_path_for_save = zarr_outpath.stem if zarr_outpath.suffix == ".zip" else zarr_outpath
                for col in tqdm(cols_to_compute, total=len(cols_to_compute)):
                    mapped_img = np.asarray(track_img.copy())
                    col_dict = dict(zip(df_tracks_clustered['TrackID'], df_tracks_clustered[col]))
                    mask = np.isin(mapped_img, list(col_dict.keys()))
                    mapped_img[~mask] = 0
                    col_values = np.array([col_dict.get(k, 0) for k in mapped_img[mask]])
                    col_values = np.nan_to_num(col_values, nan=0.0)
                    mapped_img[mask] = col_values
                    mapped_imgs[col]["img"] = mapped_img
                    save_as_zarr(
                        img=mapped_img,
                        path=zarr_path_for_save,
                        chunks=(1,) + mapped_img.shape[1:],
                        group=col
                    )
        else:
            if zarr_outpath.suffix == ".zip":
                zarr_outpath = zarr_outpath.stem
            for col in tqdm(columns, total=len(columns)):
                mapped_img = np.asarray(track_img.copy())
                col_dict = dict(zip(df_tracks_clustered['TrackID'], df_tracks_clustered[col]))
                mask = np.isin(mapped_img, list(col_dict.keys()))
                mapped_img[~mask] = 0
                col_values = np.array([col_dict.get(k, 0) for k in mapped_img[mask]])
                col_values = np.nan_to_num(col_values, nan=0.0)
                mapped_img[mask] = col_values
                mapped_imgs[col]["img"] = mapped_img
                save_as_zarr(
                    img=mapped_img,
                    path=zarr_outpath,
                    chunks=(1,) + mapped_img.shape[1:],
                    group=col
                )
        return mapped_imgs


def view_napari(
    backproject_data,
    df_tracks_full,
    df_tracks_clustered,
    cell_type,
    elsize
    ):
    """
    Visualize backprojection in napari.
    Like original: single tracks layer with color_by='ClusterID'.
    """
    
    viewer = napari.Viewer()
    
    # Add all image/label layers
    for k, v in backproject_data.items():
        v["img"] = np.transpose(v["img"], (1, 0, 2, 3, 4))
        if v["type"] == "label" or v["type"] == "segment":
            # Labels layer (like original) - hover shows the value directly
            viewer.add_labels(v["img"], name=k, scale=elsize)
        else:
            # Feature images - use inferno with auto contrast
            img_data = v["img"]
            if hasattr(img_data, 'compute'):
                img_np = np.asarray(img_data)
            else:
                img_np = img_data
            valid_vals = img_np[(img_np != 0) & np.isfinite(img_np)]
            if valid_vals.size > 0:
                vmin, vmax = np.percentile(valid_vals, [2, 98])
                if vmin >= vmax:
                    vmax = vmin + 1
            else:
                vmin, vmax = 0, 1
            viewer.add_image(v["img"], name=k, scale=elsize, colormap='inferno',
                           contrast_limits=[float(vmin), float(vmax)])
    
    # Hide all layers initially
    for lay in viewer.layers:
        lay.visible = False
    
    # Prepare tracks data for napari
    # df_tracks_clustered has ClusterID column (from _clustered.csv)
    napari_tracks_clustered = df_tracks_clustered[
        ["TrackID", "position_t", "position_z", "position_y", "position_x"]
    ].to_numpy()
    
    # Features dict - one value per point in tracks array
    features = {
        'ClusterID': df_tracks_clustered["ClusterID"].to_numpy().astype(float)
    }
    
    napari_tracks_full = df_tracks_full[
        ["TrackID", "position_t", "position_z", "position_y", "position_x"]
    ].to_numpy()
    
    # Debug: show cluster distribution
    cluster_values = features['ClusterID']
    unique_clusters = np.unique(cluster_values[~np.isnan(cluster_values)]).astype(int)
    valid_clusters = [c for c in unique_clusters if c > 0]
    print(f"  ClusterID distribution: {len(valid_clusters)} clusters: {list(valid_clusters)}")
    
    # Add filtered tracks with ClusterID coloring (like original)
    # Key fix: create layer first, then set color_by
    tracks_layer = viewer.add_tracks(
        napari_tracks_clustered,
        name=f'Filtered {cell_type} Tracks',
        features=features,
        tail_length=75
    )
    tracks_layer.color_by = 'ClusterID'
    tracks_layer.colormap = 'turbo'
    
    print(f"\n  To see ClusterID of a cell:")
    print(f"      - Enable the 'ClusterID' labels layer")
    print(f"      - Hover over a colored pixel to see the cluster number")
    
    # Add all tracks layer
    viewer.add_tracks(
        napari_tracks_full,
        name=f'All {cell_type} Tracks',
        tail_length=75,
        visible=False
    )
    
    # Show raw_data by default
    if 'raw_data' in [lay.name for lay in viewer.layers]:
        viewer.layers['raw_data'].visible = True
    
    napari.run()
    return viewer


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Backproject features onto tracked segments.')
    parser.add_argument('-c', '--config', type=str, help='path to a config.yml file', required=False)
    parser.add_argument('--cell_type', type=str, default='tcell', help='cell type to backproject')
    parser.add_argument('--mode', type=str, default='mean', choices=['mean', 'time'], help='backprojection mode')
    args = parser.parse_args()
    
    with open(args.config, "r") as parameters:
        config = yaml.load(parameters, Loader=yaml.SafeLoader)
    
    metadata = pd.read_csv(config["metadata_csv_path"])
    
    if args.mode == 'mean':
        backproject_mean_features_behav3d(config=config, cell_type=args.cell_type)
    else:
        backproject_time_features_behav3d(config=config, cell_type=args.cell_type)
