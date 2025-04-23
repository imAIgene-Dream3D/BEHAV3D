import numpy as np
import pandas as pd
from behav3d.utils.fileio import load_image
import napari

cell_type = 'tcell'
raw_path = '/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/test_windows/images/SMI_JM1_Exp004_Img14/SMI_JM1_Exp004_Img14.zarr'
seg_path = f'/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/test_windows/images/SMI_JM1_Exp004_Img14/SMI_JM1_Exp004_Img14_{cell_type}_tracked.zarr'
mask_dead_path = '/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/test_windows/images/SMI_JM1_Exp004_Img14/SMI_JM1_Exp004_Img14_mask_dead.zarr'
csv_path = f'/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/test_windows/analysis/{cell_type}/track_features/BEHAV3D_{cell_type}_combined_track_features.csv'


# Step 2: Load the CSV table
df = pd.read_csv(csv_path)  # or use pd.read_excel if Excel
seg = np.asarray(load_image(seg_path))
raw = np.asarray(load_image(raw_path))
dead = np.asarray(load_image(mask_dead_path))
# Step 3: Loop through each timepoint

# dead, raw, seg = dead[300:], raw[300:], seg[300:]
perc_dead = np.zeros_like(seg, dtype=np.float32)

for t in df['position_t'].unique():
    t=int(t)
    time_df = df[df['position_t'] == t]
    for _, row in time_df.iterrows():
        track_id = int(row['TrackID'])
        dead_pixels = row['percentage_dead_mask']  # or whatever column you need
        
        # Create a mask of where this TrackID appears in the image at this time
        mask = seg[t] == track_id
        
        # Replace with nr_dead_mask_pixels value
        perc_dead[t][mask] = dead_pixels

perc_dead = np.expand_dims(perc_dead, axis=1)  # Add channel dimension for visualization
perc_dead = np.repeat(perc_dead, 3, axis=1)  # Repeat for RGB visualization

# dead, raw, seg = dead[300:], raw[300:], seg[300:]
pixels_dead = np.zeros_like(seg, dtype=np.int32)

for t in df['position_t'].unique():
    t=int(t)
    time_df = df[df['position_t'] == t]
    for _, row in time_df.iterrows():
        track_id = int(row['TrackID'])
        dead_pixels = int(row['nr_dead_mask_pixels'])  # or whatever column you need
        
        # Create a mask of where this TrackID appears in the image at this time
        mask = seg[t] == track_id
        
        # Replace with nr_dead_mask_pixels value
        pixels_dead[t][mask] = dead_pixels
pixels_dead = np.expand_dims(pixels_dead, axis=1)  # Add channel dimension for visualization
pixels_dead = np.repeat(pixels_dead, 3, axis=1)  # R


mean_dead = np.zeros_like(seg, dtype=np.int32)

for t in df['position_t'].unique():
    t=int(t)
    time_df = df[df['position_t'] == t]
    for _, row in time_df.iterrows():
        track_id = int(row['TrackID'])
        dead_dye = row['mean_dead_dye']  # or whatever column you need
        
        # Create a mask of where this TrackID appears in the image at this time
        mask = seg[t] == track_id
        
        # Replace with nr_dead_mask_pixels value
        mean_dead[t][mask] = dead_dye
        
mean_dead = np.expand_dims(mean_dead, axis=1)  # Add channel dimension for visualization
mean_dead = np.repeat(mean_dead, 3, axis=1)  # R




dead = np.expand_dims(dead, axis=1)  # Add channel dimension for visualization
dead = np.repeat(dead, 3, axis=1)

viewer=napari.Viewer()
viewer.add_image(raw, name='raw')
viewer.add_image(dead, name='dead_pixels', colormap='gray')
viewer.add_image(perc_dead, name='percentage_dead_mask', colormap='inferno')
viewer.add_image(pixels_dead, name='nr_dead_mask_pixels', colormap='inferno')
viewer.add_image(mean_dead, name='mean_dead_dye', colormap='inferno')

# viewer.add_labels(seg, name='segmentation')
napari.run()

# Step 4: Save or visualize new_image
# e.g., np.save('relabelled_image.npy', new_image)