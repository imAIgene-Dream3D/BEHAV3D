from pathlib import Path
import napari
import numpy as np
from behav3d.utils.fileio import load_zarr
import pandas as pd

def visualize_tracks(
    metadata_row,
    timepoint_range = None,
    channel_colors =("cyan", "yellow", "red", "green", "magenta", "blue", "gray", "turbo", "viridis", "plasma", "inferno", "twilight"),
) -> None:
   
    sample_name = metadata_row['sample_name']
    print(f"Sample selected: {sample_name}")
    
    raw_image_zarr = Path(metadata_row['raw_image_path'])
    
    # Dynamically collect ALL cell types with tracks
    all_cell_tracks = {}
    for prefix in ['or', 'im', 'ot']:
        for col in metadata_row.index:
            if col.startswith(f"{prefix}_") and col.endswith("_tracks_image_path"):
                # Extract cell type name (e.g., or_organoid1_tracks_image_path -> organoid1)
                parts = col.split('_')
                if len(parts) >= 4:  # prefix_celltype_tracks_image_path
                    cell_type = '_'.join(parts[1:-3])  # everything between prefix and _tracks_image_path
                    if pd.notna(metadata_row[col]) and Path(metadata_row[col]).exists():
                        csv_col = f"{prefix}_{cell_type}_tracks_csv_path"
                        if csv_col in metadata_row.index and pd.notna(metadata_row[csv_col]):
                            all_cell_tracks[cell_type] = {
                                'image': Path(metadata_row[col]),
                                'csv': Path(metadata_row[csv_col]),
                                'prefix': prefix
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
        cell_data[cell_type] = {
            'tracks': load_zarr(paths['image']),
            'df': pd.read_csv(paths['csv']),
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
            # Slice the image data
            cell_data[cell_type]['tracks'] = cell_data[cell_type]['tracks'][start_t:end_t + 1]
            # Filter the CSV data to only include tracks within the timepoint range
            df = cell_data[cell_type]['df']
            cell_data[cell_type]['df'] = df[(df['position_t'] >= start_t) & (df['position_t'] <= end_t)].copy()
            # Adjust position_t to be relative to the new start (so it matches the sliced arrays)
            cell_data[cell_type]['df']['position_t'] = cell_data[cell_type]['df']['position_t'] - start_t

    # Launch viewer
    viewer = napari.Viewer()

    assert raw_image.ndim == 5, f"raw_image should have 5 dimensions, but {raw_image.ndim} found."
        
    T, C, Z, Y, X = raw_image.shape
    for ch in range(C):
        channel_data = raw_image[:, ch]  # (T, Z, Y, X)
        viewer.add_image(
            channel_data,
            name=f"channel_{ch+1}",
            colormap=channel_colors[ch % len(channel_colors)],
            scale=elsizes,  # T, Z, Y, X
            blending="additive",
            channel_axis=None,
        )

    # Add all cell type tracks dynamically
    # Define base colors for different categories
    category_base_colors = {
        'or': ['magenta', 'red', 'orange'],      # organoids
        'im': ['cyan', 'blue', 'turquoise'],     # immune cells
        'ot': ['yellow', 'green', 'lime']        # other cells
    }
    
    # Track how many of each category we've seen
    category_counts = {'or': 0, 'im': 0, 'ot': 0}
    
    for cell_type, data in cell_data.items():
        prefix = data['prefix']
        # Get color based on how many of this category we've already added
        color_idx = category_counts[prefix]
        if prefix not in category_base_colors:
            raise ValueError(f"Unknown prefix '{prefix}' for cell type '{cell_type}'. Expected 'or', 'im', or 'ot'.")
        colors_list = category_base_colors[prefix]
        color = colors_list[color_idx % len(colors_list)]  # Cycle through colors if more types than colors
        category_counts[prefix] += 1
        
        viewer.add_labels(data['tracks'], name=f"{cell_type} segments (tracked)", scale=elsizes)
        track_coords = data['df'][["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()
        viewer.add_tracks(track_coords, name=f'{cell_type} Tracks', colormap=color, tail_length=20)
    

    print("Launching Napari viewer...")
    napari.run()
    return viewer