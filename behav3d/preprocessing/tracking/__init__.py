from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil

import napari
import numpy as np
import pandas as pd
from tqdm import tqdm

from skimage.measure import regionprops_table

from behav3d.io.images import append_to_zarr, save_as_zarr, load_image, load_zarr, write_zarr_parallel
from .multicolor_tracking_processing import (
    combine_multicolor_tracked_outputs,
    combine_multicolor_tracked_outputs_for_base,
)

def convert_segments_to_tracks(
    df_tracks,
    segments,
    outpath,
    n_workers=None
    ):
    n_workers = max(1, int(n_workers or 1))
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

    def _build_tracked_img(t, t_segments, t_pairs):
        t_segment_img = np.asarray(t_segments)
        t_tracked_img = np.zeros_like(t_segment_img)

        for segment_id, track_id in t_pairs:
            t_tracked_img[t_segment_img == segment_id] = track_id

        return t_tracked_img

    tracks_by_time = {
        int(t): group[["SegmentID", "TrackID"]].to_numpy(copy=False)
        for t, group in df_tracks.groupby("position_t", sort=False)
    }

    if n_workers == 1:
        for t, t_segments in tqdm(enumerate(segments), total=len(segments)):
            t_tracked_img = _build_tracked_img(t, t_segments, tracks_by_time.get(t, ()))
            append_to_zarr(
                img=np.expand_dims(t_tracked_img, axis=0),
                outpath=outpath,
            )
    else:
        if len(segments) == 0:
            raise ValueError("segments must contain at least one timepoint for parallel writing")
        first_seg = np.asarray(segments[0])
        write_zarr_parallel(
            outpath=outpath,
            shape=(len(segments),) + tuple(first_seg.shape),
            dtype=first_seg.dtype,
            overwrite=True,
        )

        def _write_timepoint(t):
            t_tracked_img = _build_tracked_img(t, segments[t], tracks_by_time.get(t, ()))
            write_zarr_parallel(outpath=outpath, index=t, data=t_tracked_img)
            return t

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            for _ in tqdm(executor.map(_write_timepoint, range(len(segments))), total=len(segments)):
                pass
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
    prop_cols = ['label', 'centroid-0', 'centroid-1', 'centroid-2']
    for t, t_seg in tqdm(enumerate(segments), total=len(segments)):
        t_seg = np.asarray(t_seg)
        properties = pd.DataFrame(regionprops_table(label_image=t_seg, properties=['label', 'centroid']))
        properties = properties.reindex(columns=prop_cols)
        properties["position_t"]=t
        df_tracks.append(properties)
    if df_tracks:
        df_tracks = pd.concat(df_tracks, ignore_index=True)
    else:
        df_tracks = pd.DataFrame(columns=prop_cols + ["position_t"])
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
        
def visualize_tracks(
    metadata_row,
    timepoint_range = None,
    channel_colors =("cyan", "yellow", "red", "green", "magenta", "blue", "gray", "turbo", "viridis", "plasma", "inferno", "twilight"),
) -> None:
    from behav3d.core.metadata import is_multicolor_celltype, is_combined_multicolor_celltype

    sample_name = metadata_row['sample_name']
    print(f"Sample selected: {sample_name}")
    
    raw_image_zarr = Path(metadata_row['raw_image_path'])
    print(f"Loading raw image: {raw_image_zarr}")
    
    # Collect ALL cell types with tracks:
    # 1. Prefixed types: or_/im_/ot_ columns
    # 2. Merged/grouped types: standalone columns like tcells_merged_tracks_image_path
    all_cell_tracks = {}
    suffixes_image = ("_tracks_image_path",)
    suffixes_csv = ("_tracks_csv_path",)

    for prefix in ['or', 'im', 'ot']:
        for col in metadata_row.index:
            if col.startswith(f"{prefix}_") and col.endswith("_tracks_image_path"):
                parts = col.split('_')
                if len(parts) >= 4:
                    cell_type = '_'.join(parts[1:-3])
                    if pd.notna(metadata_row[col]) and Path(metadata_row[col]).exists():
                        csv_col = f"{prefix}_{cell_type}_tracks_csv_path"
                        if csv_col in metadata_row.index and pd.notna(metadata_row[csv_col]):
                            all_cell_tracks[cell_type] = {
                                'image': Path(metadata_row[col]),
                                'csv': Path(metadata_row[csv_col]),
                                'prefix': prefix,
                            }

    # Merged / grouped outputs (no or_/im_/ot_ prefix)
    for col in metadata_row.index:
        if col.startswith(("or_", "im_", "ot_")):
            continue
        if not col.endswith("_tracks_image_path"):
            continue
        base = col[: -len("_tracks_image_path")]
        if not is_combined_multicolor_celltype(base):
            continue
        val = metadata_row[col]
        if pd.isna(val) or not str(val).strip():
            continue
        p = Path(str(val).strip())
        if not p.exists():
            continue
        csv_col = f"{base}_tracks_csv_path"
        csv_val = metadata_row.get(csv_col)
        all_cell_tracks[base] = {
            'image': p,
            'csv': Path(str(csv_val).strip()) if pd.notna(csv_val) and str(csv_val).strip() else None,
            'prefix': '',  # no prefix — merged/grouped
        }

    if not all_cell_tracks:
        raise ValueError(f"No tracked cell types found in sample {sample_name}. Expected prefixed columns (or_/im_/ot_)_{{cell_type}}_tracks_image_path")

    elsize_xy = metadata_row['pixel_distance_xy']
    elsize_z = metadata_row['pixel_distance_z']
    elsizes = (1, elsize_z, elsize_xy, elsize_xy)
    
    # Load raw image
    raw_image = load_zarr(raw_image_zarr)
    
    # Load all cell type tracks
    cell_data = {}
    for cell_type, paths in all_cell_tracks.items():
        print(f"Loading tracked labels for {cell_type}: {paths['image']}")
        csv_path = paths['csv']
        df = pd.read_csv(csv_path) if (csv_path is not None and csv_path.exists()) else pd.DataFrame()
        cell_data[cell_type] = {
            'tracks': load_zarr(paths['image']),
            'df': df,
            'prefix': paths['prefix']
        }
    
    print(f"Original image shape: {raw_image.shape}")

    # Slice timepoints
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, raw_image.shape[0] - 1))
        end_t = max(start_t, min(end_t, raw_image.shape[0] - 1))
        print(f"Slicing timepoints from {start_t} to {end_t}")
        
        raw_image = raw_image[start_t:end_t + 1]
        for cell_type in cell_data:
            cell_data[cell_type]['tracks'] = cell_data[cell_type]['tracks'][start_t:end_t + 1]
            df = cell_data[cell_type]['df']
            if not df.empty and 'position_t' in df.columns:
                cell_data[cell_type]['df'] = df[(df['position_t'] >= start_t) & (df['position_t'] <= end_t)].copy()
                cell_data[cell_type]['df']['position_t'] = cell_data[cell_type]['df']['position_t'] - start_t

    # Launch viewer
    viewer = napari.Viewer()

    assert raw_image.ndim == 5, f"raw_image should have 5 dimensions, but {raw_image.ndim} found."
        
    T, C, Z, Y, X = raw_image.shape
    for ch in range(C):
        channel_data = raw_image[:, ch]  # (T, Z, Y, X)
        viewer.add_image(
            channel_data,
            name=f"channel_{ch}",
            colormap=channel_colors[ch % len(channel_colors)],
            scale=elsizes,  # T, Z, Y, X
            blending="additive",
            channel_axis=None,
        )

    # Add all cell type tracks dynamically
    category_base_colors = {
        'or': ['magenta', 'red', 'orange'],      # organoids
        'im': ['cyan', 'blue', 'darkblue'],     # immune cells
        'ot': ['yellow', 'green', 'lime'],       # other cells
        '':   ['white', 'gray'],                 # merged/grouped (no prefix)
    }
    category_counts = {'or': 0, 'im': 0, 'ot': 0, '': 0}
    
    for cell_type, data in cell_data.items():
        prefix = data['prefix']

        # Determine visibility:
        # - individual per-channel multicolor layers → hidden (keep in viewer but not cluttering)
        # - everything else (including merged/grouped) → visible
        visible = not is_multicolor_celltype(cell_type)

        color_idx = category_counts.get(prefix, 0)
        colors_list = category_base_colors.get(prefix, ['white'])
        color = colors_list[color_idx % len(colors_list)]
        category_counts[prefix] = category_counts.get(prefix, 0) + 1
        
        try:
            tracked_labels = np.asarray(data['tracks'])
            viewer.add_labels(
                tracked_labels,
                name=f"{cell_type} segments (tracked)",
                scale=elsizes,
                visible=visible,
            )
        except Exception as e:
            print(f"  Skipping {cell_type} segments layer: {e}")

        df = data['df']
        if df.empty:
            print(f"  Skipping {cell_type} tracks: no CSV data")
            continue
        required = {"TrackID", "position_t"}
        if not required.issubset(set(df.columns)):
            print(f"  Skipping {cell_type} tracks: missing required CSV columns")
            continue
        pos_cols = ["TrackID", "position_t"]
        for col in ("pixel_position_z", "pixel_position_y", "pixel_position_x",
                    "position_z", "position_y", "position_x"):
            if col in df.columns:
                pos_cols.append(col)
            if len(pos_cols) == 5:  # TrackID + t + z + y + x
                break
        # Fall back to 4-column format if z is missing
        if len(pos_cols) < 5:
            for col in ("pixel_position_y", "pixel_position_x", "position_y", "position_x"):
                if col in df.columns and col not in pos_cols:
                    pos_cols.append(col)
                if len(pos_cols) == 4:
                    break
        track_coords = df[pos_cols].to_numpy()
        try:
            viewer.add_tracks(
                track_coords,
                name=f'{cell_type} Tracks',
                colormap=color,
                tail_length=20,
                visible=visible,
            )
        except Exception as e:
            print(f"  Skipping {cell_type} tracks layer: {e}")
    

    print("Launching Napari viewer...")
    napari.run()
    return viewer
