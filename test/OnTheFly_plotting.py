
from behav3d.utils import load_behav3d_metadata, check_behav3d_metadata
from behav3d.utils.fileio import load_image
from behav3d.utils.preprocessing import calc_z_projection
from scipy.ndimage import center_of_mass
import numpy as np
import pandas as pd
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import torch
import zarr
import dask.array as da
import napari
from napari_animation import Animation
from tqdm import tqdm
from qtpy.QtWidgets import QApplication

from napari.utils import nbscreenshot
from PIL import Image
import imageio.v3 as iio
from PIL import ImageChops
from skimage import measure
from scipy.ndimage import binary_erosion

output_dir = r"/Volumes/T7_Sam//BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE"
metadata_csv_path = r"/Volumes/T7_Sam//BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv"
metadata = load_behav3d_metadata(metadata_csv_path)


        
# sample_name = "ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)"
sample_name = "ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)"

# sample_name = "ROCHE_JM1_Exp011_Img09"

sample_metadata = metadata[metadata["sample_name"] == sample_name].iloc[0]

img_outdir = Path(output_dir, "images", sample_name)
raw_image_path = Path(img_outdir, f"{sample_name}.zarr")
org_tracks_path = Path(sample_metadata["organoid_tracks_image_path"])
tcell_tracks_path = Path(sample_metadata["tcell_tracks_image_path"])
tcell_tracks_csv_path = Path(sample_metadata["tcell_tracks_csv_path"])
dead_mask_path = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
organoid_track_csv_path = f"/Volumes/T7_Sam//BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/analysis/organoid/results/per_sample/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)/organoid/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)_organoid_track_analysis.csv"
tcell_segments_path = Path(sample_metadata["tcell_segments_image_path"])

raw_image = load_image(raw_image_path)
raw_image_tcell = raw_image[:,0]
raw_image_org = raw_image[:,1]
raw_image_death = raw_image[:,2]

org_tracks = load_image(org_tracks_path)
tcell_tracks = load_image(tcell_tracks_path)
tcell_segments = load_image(tcell_segments_path)
dead_mask = load_image(dead_mask_path)
tcell_tracks_csv = pd.read_csv(tcell_tracks_csv_path)
org_tracks_csv = pd.read_csv(organoid_track_csv_path)
# raw_image = calc_z_projection(raw_image, projection="max", z_axis=-3).compute()
# org_tracks = calc_z_projection(org_tracks, projection="max", z_axis=-3).compute()
# tcell_tracks = calc_z_projection(tcell_tracks, projection="max", z_axis=-3).compute()
# dead_mask = calc_z_projection(dead_mask, projection="max", z_axis=-3).compute()
dying_orgs = np.asarray(org_tracks.copy()).astype(np.float32)
for t, t_dying_orgs in tqdm(enumerate(dying_orgs), total = dying_orgs.shape[0]):
    t_df = org_tracks_csv[org_tracks_csv["position_t"] == t]
    # col_dict = dict(zip(t_df['TrackID'], t_df["dead"]))
    col_dict = dict(zip(t_df['TrackID'], t_df["dead"]))
    mask = np.isin(dying_orgs[t], list(col_dict.keys()))
    dying_orgs[t][~mask]=0
    dying_orgs[t][mask] = np.vectorize(col_dict.get)(dying_orgs[t][mask])


def compute_2d_borders_from_3d_labels(label_image_4d):
    # label_image_4d: (T, Z, Y, X)
    T, Z, Y, X = label_image_4d.shape

    # Max project along Z → shape (T, Y, X)
    projected = label_image_4d.max(axis=1)

    # Prepare border array with shape (T, 1, Y, X)
    borders = np.zeros((T, 1, Y, X), dtype=label_image_4d.dtype)

    for t in range(T):
        labels_2d = projected[t]
        for label in np.unique(labels_2d):
            if label == 0:
                continue
            mask = labels_2d == label
            eroded = binary_erosion(mask)
            border = mask & ~eroded
            borders[t, 0][border] = label

    return borders
org_borders = compute_2d_borders_from_3d_labels(org_tracks.compute())

#Visualized for behav3d(4): 5, 15
segment_id = 46
mask = org_tracks == segment_id

raw_image_tcell[~mask] = 0
raw_image_org[~mask] = 0
raw_image_death[~mask] = 0

# raw_image_org = raw_image[:,1]
# raw_image_death = raw_image[:,2]

napari_tracks_full = tcell_tracks_csv[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()

elsize = (sample_metadata["pixel_distance_z"], sample_metadata["pixel_distance_xy"], sample_metadata["pixel_distance_xy"])
viewer = napari.Viewer()

viewer.add_image(raw_image_tcell, name='Raw Image (T-cell)', colormap= "cyan", scale=elsize, blending='additive')
viewer.add_image(raw_image_org, name='Raw Image (Organoid)', colormap= "yellow", scale=elsize, blending='additive')
viewer.add_image(raw_image_death, name='Raw Image (Death)', colormap= "red", scale=elsize, blending='additive')
viewer.add_tracks(napari_tracks_full, name='All T-cell Tracks',  tail_length=20)
viewer.add_image(dead_mask, name='Dead Mask', scale=elsize, blending='additive')
viewer.add_labels(tcell_segments, name='T-cell Segments', scale=elsize, blending='additive')
viewer.add_labels(tcell_tracks, name='T-cell Tracks', opacity=0.5, scale=elsize)
viewer.add_labels(org_tracks, name='Organoid Tracks', opacity=0.5, scale=elsize)
viewer.add_labels(org_borders, name='Organoid Tracks (borders)', opacity=0.5, scale=elsize)
viewer.add_image(dying_orgs, name='Dying Organoids', colormap='yellow', scale=elsize, blending='additive')
# viewer.add_shapes(contours, shape_type='polygon', edge_color='red', name="Contours")
# viewer.close()


images = []

def crop_to_content(img_np, threshold=10):
    """
    Crop black borders from an RGB numpy image.
    Args:
        img_np: (H, W, 3) image
        threshold: pixel intensity threshold
    Returns:
        cropped image (np.ndarray)
    """
    gray = np.mean(img_np, axis=2)
    mask = gray > threshold
    if not np.any(mask):
        return img_np  # Return original if all black

    y_min, y_max = np.where(mask.any(axis=1))[0][[0, -1]]
    x_min, x_max = np.where(mask.any(axis=0))[0][[0, -1]]

    return img_np[y_min:y_max+1, x_min:x_max+1]

process_events = QApplication.processEvents
for i in tqdm(range(tcell_tracks.shape[0])):
    viewer.dims.set_current_step(0, i)
    process_events()
    img = viewer.screenshot(canvas_only=True)
    
    # Crop the black borders from the screenshot
    # pil_img = Image.fromarray(img)
    cropped = crop_to_content(img)
    images.append(cropped)

# Convert PIL images to numpy arrays
frames_np = [np.array(img) for img in images]

# Save as .mp4 with h.264 encoding
iio.imwrite(
    "/Users/s.deblank-3/Downloads/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)_DeadMask.mp4",
    frames_np,
    fps=20,
    codec="libx264",
    quality=8,  # 0–10 (10 = best)
    pixelformat="yuv420p",  # ensure compatibility
)







df = pd.read_csv(r"/Volumes/T7_Sam/\BHVD_BEHAV3D\BEHAV3D_python\runs\ROCHE\analysis\organoid\results\per_sample\ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)\organoid\ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)_organoid_general_analysis.csv")

time = df['position_t']
alive = df['nr_alive']
percentage_alive = df['percentage_alive']

# Set up figure and axes
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # second y-axis for percentage reference

# One orange line
line, = ax1.plot([], [], lw=2, color='blue', label='Alive')

# Axis limits and labels
ax1.set_xlim(0, time.max())
ax1.set_ylim(0, alive.max() + 1)
ax2.set_ylim(0, 110)

ax1.set_xlabel('Time')
ax1.set_ylabel('Alive Count')
ax2.set_ylabel('% Alive')
plt.title('Alive Organoids Over Time')

# Init
def init():
    line.set_data([], [])
    return line,

# Animate: 1 timepoint per frame
def animate(i):
    x = time[:i+1]
    y = alive[:i+1]
    line.set_data(x, y)
    return line,

# Animate: 450 timepoints, 20 FPS, 1 timepoint per frame
ani = animation.FuncAnimation(
    fig, animate, init_func=init,
    frames=len(time), interval=1000 / 20, blit=True  # 20 FPS
)
# Save the animation as MP4
outpath = r"C:\Users\Samde\Downloads\alive_organoids_animation.mp4"
ani.save(outpath, fps=20, extra_args=['-vcodec', 'libx264'])
