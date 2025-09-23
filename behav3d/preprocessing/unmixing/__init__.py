from pathlib import Path
import napari
import numpy as np
from behav3d.utils.fileio import load_zarr
import pandas as pd

def visualize_unmix( 
                    metadata_row: pd.Series, 
                    timepoint_range=None, 
                    channel_colors =("cyan", "yellow", "red", "green", "magenta", "blue"),
                    ) -> None:
    """Launch napari viewer for the given metadata row (pandas Series)."""
    sample_name = metadata_row['sample_name']
    print(f"Sample selected: {sample_name}")
    
    unmix_img_zarr = Path(metadata_row['signal_unmixing_image_path']).expanduser()

    elsize_xy = metadata_row['pixel_distance_xy']
    elsize_z = metadata_row['pixel_distance_z']
    elsizes = (1, elsize_z, elsize_xy, elsize_xy)
    unmix_img = load_zarr(unmix_img_zarr)
    
    print(f"Original image shape: {unmix_img.shape}")

    # Slice timepoints
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, unmix_img.shape[0] - 1))
        end_t = max(start_t, min(end_t, unmix_img.shape[0] - 1))
        print(f"Slicing timepoints from {start_t} to {end_t}")        
        unmix_img = unmix_img[start_t:end_t + 1]

    # for name in loaded_masks:
    #     loaded_masks[name] = np.asarray(loaded_masks[name])

    # Launch viewer
    viewer = napari.Viewer()

    assert unmix_img.ndim == 5, f"raw_image should have 5 dimensions, but {unmix_img.ndim} found."
        
    T, C, Z, Y, X = unmix_img.shape
    for ch in range(C):
        channel_data = unmix_img[:, ch]  # (T, Z, Y, X)
        viewer.add_image(
            channel_data,
            name=f"channel_{ch+1}",
            colormap=channel_colors[ch % len(channel_colors)],
            scale=elsizes,  # T, Z, Y, X
            blending="additive",
            channel_axis=None,
        )

    print("Launching Napari viewer...")
    napari.run()
    return viewer