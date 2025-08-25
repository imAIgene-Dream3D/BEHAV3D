
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
from behav3d.utils.preprocessing import filter_median

output_dir = r"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE"
metadata_csv_path = r"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/metadata.csv"
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
organoid_track_csv_path = f"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/analysis/organoid/results/per_sample/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)/organoid/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)_organoid_track_analysis.csv"
organoid_track_csv_path = f"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/analysis/organoid/results/per_sample/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)/organoid/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)_organoid_track_analysis.csv"

tcell_segments_path = Path(sample_metadata["tcell_segments_image_path"])

raw_image = load_image(raw_image_path).compute()
raw_image_tcell = raw_image[:,0]
raw_image_org = raw_image[:,1]
raw_image_death = raw_image[:,2]

# from dask_image.ndfilters import median_filter
# from dask.diagnostics import ProgressBar


# def k(r): return 2*r + 1

# # Per-channel filtering in 3D (Z,Y,X), preserving time (T)
# tcell  = raw_image[:, 0]                         # (T, Z, Y, X)
# org    = raw_image[:, 1]
# death  = raw_image[:, 2]

# tcell_f = median_filter(tcell, size=(1, k(1), k(1), k(1)))   # radius=2
# org_f   = median_filter(org,   size=(1, k(1), k(1), k(1)))   # radius=3
# death_f = median_filter(death, size=(1, k(1), k(1), k(1)))   # radius=1

# Choose scheduler:
# Threads (good with SciPy/Dask C-code that releases GIL)
# with ProgressBar():  # shows a live text bar in the cell/output
#     tcell_f_np, org_f_np, death_f_np = da.compute(tcell_f, org_f, death_f)
# raw_image_tcell = tcell_f_np
# raw_image_org = org_f_np
# raw_image_death = death_f_np
# #median filter raw_image
# raw_image_tcell = filter_median(raw_image_tcell, radius=2)
# raw_image_org = filter_median(raw_image_org, radius=3)
# raw_image_death = filter_median(raw_image_death, radius=1)

org_tracks = load_image(org_tracks_path).compute()
tcell_tracks = load_image(tcell_tracks_path).compute()
tcell_segments = load_image(tcell_segments_path).compute()
dead_mask = load_image(dead_mask_path).compute()
tcell_tracks_csv = pd.read_csv(tcell_tracks_csv_path)
org_tracks_csv = pd.read_csv(organoid_track_csv_path)


# t_min = 125
# raw_image = raw_image[t_min:]
# raw_image_tcell = raw_image_tcell[t_min:]
# raw_image_org = raw_image_org[t_min:]
# raw_image_death = raw_image_death[t_min:]
# org_tracks = org_tracks[t_min:]
# tcell_tracks = tcell_tracks[t_min:]
# dead_mask = dead_mask[t_min:]

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


# def compute_2d_borders_from_3d_labels(label_image_4d):
#     # label_image_4/Volumes/T7_Sam (T, Z, Y, X)
#     T, Z, Y, X = label_image_4d.shape

#     # Max project along Z → shape (T, Y, X)
#     projected = label_image_4d.max(axis=1)

#     # Prepare border array with shape (T, 1, Y, X)
#     borders = np.zeros((T, 1, Y, X), dtype=label_image_4d.dtype)

#     for t in range(T):
#         labels_2d = projected[t]
#         for label in np.unique(labels_2d):
#             if label == 0:
#                 continue
#             mask = labels_2d == label
#             eroded = binary_erosion(mask)
#             border = mask & ~eroded
#             borders[t, 0][border] = label

#     return borders
# org_borders = compute_2d_borders_from_3d_labels(org_tracks.compute())

#Visualized for behav3d(4): 5, 15
# # !!! MAke LIST
# segment_ids = [3,7,13,39]
# mask = np.isin(org_tracks, segment_ids)
# mask  = mask.any(axis=0)    
# mask = np.broadcast_to(mask, (org_tracks.shape[0], *mask.shape))
# mask.shape


# raw_image_tcell[~mask] = 0
# raw_image_org[~mask] = 0
# raw_image_death[~mask] = 0

# raw_image_org = raw_image[:,1]
# raw_image_death = raw_image[:,2]

napari_tracks_full = tcell_tracks_csv[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()
# org_tracks_full = org_tracks_csv[["TrackID", "position_t", "position_z", "position_y", "position_x"]].to_numpy()
elsize = (sample_metadata["pixel_distance_z"], sample_metadata["pixel_distance_xy"], sample_metadata["pixel_distance_xy"])
viewer = napari.Viewer()

viewer.add_labels(dying_orgs.astype(int), name='Dying Organoids', scale=elsize, )
viewer.add_image(raw_image_tcell, name='Raw Image (T-cell)', colormap= "cyan", scale=elsize, blending='additive')
viewer.add_image(raw_image_org, name='Raw Image (Organoid)', colormap= "yellow", scale=elsize, blending='additive')
viewer.add_image(raw_image_death, name='Raw Image (Death)', colormap= "red", scale=elsize, blending='additive')
viewer.add_tracks(napari_tracks_full, name='All T-cell Tracks',  tail_length=20)
viewer.add_image(dead_mask, name='Dead Mask', scale=elsize, blending='additive')
# viewer.add_labels(tcell_segments, name='T-cell Segments', scale=elsize, blending='additive')
viewer.add_labels(tcell_tracks, name='T-cell Tracks', opacity=0.5, scale=elsize)
viewer.add_labels(org_tracks, name='Organoid Tracks', opacity=0.5, scale=elsize)
# viewer.add_labels(org_borders, name='Organoid Tracks (borders)', opacity=0.5, scale=elsize)
# viewer.add_labels(dying_orgs, name='Dying Organoids', colormap='yellow', scale=elsize, blending='additive')
# viewer.add_shapes(contours, shape_type='polygon', edge_color='red', name="Contours")
# viewer.close()


images = []

def crop_to_content(img_np, threshold=10):
    """
    Crop black borders from an RGB numpy image.
    Args:
        img_np: (H, W, 3) image
        threshol/Volumes/T7_Sam pixel intensity threshold
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
    "/Users/s.deblank-3/surfdrive/Documents/Presentations/videos/20250825_BEHAV3DExample3D_DyingOrganoids.mp4",
    frames_np,
    fps=40,
    codec="libx264",
    quality=8,  # 0–10 (10 = best)
    pixelformat="yuv420p",  # ensure compatibility
)


# Load your CSV
df = pd.read_csv('/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/analysis/organoid/results/per_sample/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)/organoid/ROCHE_JM1_Exp016-1_Img02_MAGEA4_TCB_Behav3d(4)_organoid_general_analysis.csv')

time = df['position_t']
alive = df['nr_alive']
percentage_alive = df['percentage_alive']

# --- Figure & Axes ---
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()  # second y-axis for percentage reference

# Set black background for figure and axes
fig.patch.set_facecolor('black')
ax1.set_facecolor('black')
ax2.set_facecolor('black')

# One white line
line, = ax1.plot([], [], lw=2, color='red', label='Alive')

# Axis limits and labels
ax1.set_xlim(0, time.max())
ax1.set_ylim(0, alive.max() + 1)
ax2.set_ylim(0, 110)

ax1.set_xlabel('Time', color='white')
ax1.set_ylabel('Alive Count', color='white')
ax2.set_ylabel('% Alive', color='white')
plt.title('Alive Organoids Over Time', color='white')

# Make all ticks and spines white
for ax in [ax1, ax2]:
    ax.tick_params(axis='x', colors='white')
    ax.tick_params(axis='y', colors='white')
    for spine in ax.spines.values():
        spine.set_color('white')

# --- Init function ---
def init():
    line.set_data([], [])
    return line,

# --- Animate function ---
def animate(i):
    x = time[:i+1]
    y = alive[:i+1]
    line.set_data(x, y)
    return line,

# --- Animate: 20 FPS ---
ani = animation.FuncAnimation(
    fig, animate, init_func=init,
    frames=len(time), interval=1000 / 20, blit=True
)

# Save as MP4
outpath = r"/Users/s.deblank-3/Downloads/alive_organoids_animation_blackbg.mp4"
ani.save(outpath, fps=20, extra_args=['-vcodec', 'libx264'])