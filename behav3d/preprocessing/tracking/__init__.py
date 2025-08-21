from pathlib import Path
import napari
import numpy as np
from behav3d.utils.fileio import load_zarr
import pandas as pd

def visualize_tracks(
    metadata_row,
    timepoint_range = None,
    channel_colors =("cyan", "yellow", "red", "green", "magenta", "blue"),
) -> None:
   
    sample_name = metadata_row['sample_name']
    print(f"Sample selected: {sample_name}")
    
    raw_image_zarr = Path(metadata_row['raw_image_path'])
    tcell_tracks_zarr = Path(metadata_row['tcell_tracks_image_path'])
    organoid_tracks_zarr = Path(metadata_row['organoid_tracks_image_path'])
    tcell_tracks_csv_path = Path(metadata_row['tcell_tracks_csv_path'])
    organoid_tracks_csv_path = Path(metadata_row['organoid_tracks_csv_path'])

    elsize_xy = metadata_row['pixel_distance_xy']
    elsize_z = metadata_row['pixel_distance_z']
    elsizes = (1, elsize_z, elsize_xy, elsize_xy)
    raw_image = load_zarr(raw_image_zarr)
    tcell_tracks = load_zarr(tcell_tracks_zarr)
    organoid_tracks = load_zarr(organoid_tracks_zarr)
    df_tcell_tracks = pd.read_csv(tcell_tracks_csv_path)
    df_organoid_tracks = pd.read_csv(organoid_tracks_csv_path)
    
    
    print(f"Original image shape: {raw_image.shape}")

    # Slice timepoints
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, img_arr.shape[0] - 1))
        end_t = max(start_t, min(end_t, img_arr.shape[0] - 1))
        print(f"Slicing timepoints from {start_t} to {end_t}")
        img_arr = img_arr[start_t:end_t + 1]
        
        raw_image = raw_image[start_t:end_t + 1]
        tcell_tracks = tcell_tracks[start_t:end_t + 1]
        organoid_tracks = organoid_tracks[start_t:end_t + 1]

    # for name in loaded_masks:
    #     loaded_masks[name] = np.asarray(loaded_masks[name])

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

    viewer.add_labels(tcell_tracks, name="T cell segments (tracked)", scale=elsizes)
    viewer.add_labels(organoid_tracks, name="Organoid segments (tracked)", scale=elsizes)
    
    tcell_track_coords = df_tcell_tracks[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()
    organoid_track_coords = df_organoid_tracks[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()

    viewer.add_tracks(tcell_track_coords, name='T-cell Tracks', colormap="cyan", tail_length=20)
    viewer.add_tracks(organoid_track_coords, name='Organoid Tracks',  colormap="cyan", tail_length=20)
    

    print("Launching Napari viewer...")
    napari.run()
    return viewer