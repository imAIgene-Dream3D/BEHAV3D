from pathlib import Path
import numpy as np
import time
import shutil

import napari
from magicgui import magicgui
from magicgui.widgets import PushButton
from qtpy.QtWidgets import QPlainTextEdit, QWidget, QVBoxLayout, QApplication

from aicspylibczi import CziFile

from skimage import data, segmentation, feature, future
from skimage.measure import label
from skimage.segmentation import watershed, relabel_sequential

from sklearn.ensemble import RandomForestClassifier
# from napari_apoc import PixelClassifier  # Commented out - use scikit-learn instead
from scipy.ndimage import binary_fill_holes, find_objects

from behav3d.utils.preprocessing import open_mask, dilate_mask
from behav3d.utils.segmentation import segment_size_filter, get_border_segments, remove_boundary_segments, calculate_edt, segment_2d_filter
from behav3d.utils.fileio import save_as_zarr, load_zarr, load_image, append_to_zarr

from tqdm import tqdm

from functools import partial

import zarr
import dask.array as da
import gc

from picasso.nn_picasso import PICASSOnn
import os
import random
from scipy.ndimage import gaussian_filter, median_filter
import torch

## TODO create a function for BEHAV3D notebook

def signal_unmixing(
    metadata,
    output_dir,
    manual_dim_order=None,
    sink_channel=None,
    source_channels=None,
    sum_scale=None,
    timepoints=None,
    train_z=[4,5,6],  
    train_t=[3],
    ref_z=4,
    bg_percentile=1,
    median_size=3,
    gaussian_sigma=4,
    visualize_napari=True,
    training_log=True
):

    """
    Perform signal unmixing on microscopy image data and save the result in Zarr format.

    Parameters
    ----------
    metadata : object
        Metadata for the input image (contains acquisition details and paths).
    output_dir : str or Path
        Directory where the unmixed Zarr file will be saved.
    manual_dim_order : str, optional
        Explicit dimension order string (e.g. "TCZYX"). If None, will attempt to infer from metadata.
    sink_channel : int, default=1
        Index of the sink channel (0-based). This channel is the one to be unmixed. If not provided, defaults to tcell channel.
    source_channels : list of int, default=[0, 2]
        Indices of the source channels (0-based). These channels contribute to the sink channel. If not provided, defaults to live and dead channels.
    sum_scale : list of float, default=[20, 30]
        Scaling factors for each of the source channels. Must be the same length as `source_channels`.
    timepoints : int, optional
        Number of time frames to process. If None, process all available time frames.
    train_z : int or list of int, default=[5]
        Number of z-slices to use for training the unmixing model, or a list of specific z-slices. If 0 will use only `ref_z`.
    train_t : int or list of int, default=[3]
        Number of time frames to use for training, or a list of specific frames. If 0 will use only t=0.
    ref_z : int, default=4
        Reference z-slice index (0-based) to use if `train_z` is not provided.
    bg_percentile : float, default=1
        Percentile for estimating the background signal (e.g., 1 = 1st percentile).
    median_filter_size : int, default=3
        Size of the median filter to apply during preprocessing.
    gaussian_sigma : float, default=4
        Sigma value for Gaussian smoothing applied during preprocessing.
    visualize_napari : bool, default=True
        Whether to visualize the results in Napari for user inspection.
    training_log : bool, default=True
        Whether to print training logs during model fitting.

    Returns
    -------

    zarr files with unmixed data saved in the $OUTPUT_PATH/image/SignalUnmixing directory.

    Notes
    -----

    - The function will visualize the result in Napari (if enabled by workflow).
    - User can interactively decide whether to accept or rerun the unmixing with different parameters.

    """

    # Set random seeds for reproducibility
    np.random.seed(42)
    random.seed(42)
    torch.manual_seed(42)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(42)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

    # --- Background removal and signal enhancement for ch0 and ch1 ---
    def remove_background_and_enhance(img, bg_value, max_value, scale=1.0, median_size=3, gaussian_sigma=1.0, norm_max=65535):
        """
        Apply scaling first, then median filter, then subtract background, then normalize.
        """
        from scipy.ndimage import median_filter
        # Apply scaling first
        scaled = img.astype(np.float32) * scale
        # Apply median filter
        gauss = gaussian_filter(scaled, sigma=gaussian_sigma)    
        filtered = median_filter(gauss, size=median_size)


        # Remove background (subtract bg_value, clip at 0)
        bg_removed = np.clip(filtered - bg_value, 0, None)
        # Normalize to norm_max
        if max_value > 0:
            norm = np.clip(bg_removed / max_value * norm_max, 0, norm_max)
        else:
            norm = bg_removed
        return norm.astype(np.uint16)
    
    import sys, os
    from contextlib import contextmanager

    # Suppress stdout
    @contextmanager
    def suppress_stdout():
        with open(os.devnull, "w") as devnull:
            old_stdout = sys.stdout
            sys.stdout = devnull
            try:
                yield
            finally:
                sys.stdout = old_stdout

    
    # --- Main function body ---
    # if n_workers is None:
    #     n_workers = multiprocessing.cpu_count()
        
    signal_unmix_outdir = Path(output_dir, "images", "SignalUnmixing")
    signal_unmix_outdir.mkdir(exist_ok=True, parents=True)

    for idx, sample in metadata.iterrows():
        
        sample_name = sample['sample_name']
        
        #TODO: What should the user input in metadata?
        tcell_ch=sample['tcell_channel'] # Channel to unmix by default
        live_ch=sample['live_channel']
        dead_ch=sample['dead_channel']

        # Set the default channels and scales if not provided
        if sink_channel is None:
            sink_channel = live_ch
        if source_channels is None:
            source_channels = [0, 1]
        if sum_scale is None:
            sum_scale = [20, 30]
        if not train_z:
            train_z = [ref_z]
        if not train_t:
            train_t = [0]

        raw_image_path = Path(sample['raw_image_path'])
        raw_image_zarr =  Path(signal_unmix_outdir, sample_name, f"{sample_name}.zarr")

        if timepoints is not None:
            output_zarr_path = Path(signal_unmix_outdir, sample_name, f"{sample_name}_unmixed_t{timepoints}.zarr")
            all_outputs_path = Path(signal_unmix_outdir, sample_name, f"{sample_name}_all_outputs_t{timepoints}.zarr")
        else:
            output_zarr_path = os.path.join(signal_unmix_outdir, sample_name, f"{sample_name}_full_unmixed.zarr")
            all_outputs_path = os.path.join(signal_unmix_outdir, sample_name, f"{sample_name}_all_outputs.zarr")
        
        if not raw_image_zarr.exists():
            print(f"- Converting raw image to .zarr for memory efficiency...")
            images = load_image(raw_image_path)
            print(f"Original image shape: {images.shape}")
            default_order = ("T", "C", "Z", "Y", "X")  # T, C, Z, Y, X
            if manual_dim_order is not None:
                # Convert to tuple of characters if string
                if isinstance(manual_dim_order, str):
                    manual_dim_order = tuple(manual_dim_order)
                if manual_dim_order != default_order:
                    # Compute permutation: for each axis in default_order, find its index in manual_dim_order
                    perm = [manual_dim_order.index(ax) for ax in default_order]
                    print(f"Transposing image from {manual_dim_order} to {default_order} using permutation {perm}")
                    images = images.transpose(perm)
                    print(f"New image shape: {images.shape}")
                else:
                    print("Image is already in default order (T, C, Z, Y, X).")
            else:
                print("No manual_dim_order provided, assuming image is already in default order.")
            chunksize = (1,) + images.shape[1:]
            save_as_zarr(
                img=images, 
                path=raw_image_zarr, 
                chunks=chunksize
            )

        # --- Read from Zarr and process only the first N timepoints ---
        lazy_data = da.from_zarr(raw_image_zarr)
        T, C, Z, Y, X = lazy_data.shape
        N = timepoints if timepoints is not None else T
        lazy_data = lazy_data[:N]
        print("Zarr data shape:", lazy_data.shape)


        # --- Calculate global background and max for ch0 and ch1 (for normalization) ---
        # Use all data for global stats (could use only training set if desired)
        ch0_data = lazy_data[:, source_channels[0]].compute().flatten()
        ch1_data = lazy_data[:, source_channels[1]].compute().flatten()
        bg0 = np.percentile(ch0_data, bg_percentile)
        bg1 = np.percentile(ch1_data, bg_percentile)
        max0 = (ch0_data - bg0).clip(min=0).max()
        max1 = (ch1_data - bg1).clip(min=0).max()
        print(f"ch0: BG={bg0:.2f}, Max={max0:.2f}")
        print(f"ch1: BG={bg1:.2f}, Max={max1:.2f}")

        # --- TRAINING: concatenate user-specified (t, z) for each source, train on that image only ---
        mixing_matrix = np.array([[1], [-1]])

        # Integrate compatibility with both int and list inputs for train_t and train_z
        if train_z is None:
            train_z = [ref_z]
        elif isinstance(train_z, int):
            train_z = list(np.linspace(0, Z, train_z, endpoint=False, dtype=int))

        if train_t is None:
            train_t = [0]
        if isinstance(train_t, int):
            train_t = list(np.linspace(0, N, train_t, endpoint=False, dtype=int))

        # First unmix: sink=ch2, source=ch0 (bg removed, median filtered, normalized, scaled)
        sink_images = []
        source_images = []
        for t in train_t:
            t_data = lazy_data[t].compute()
            for z in train_z:
                sink_images.append(t_data[sink_channel][z])
                raw_source = t_data[source_channels[0]][z] #* sum_scale[0]
                bg_removed = remove_background_and_enhance(raw_source, bg0, max0,scale=sum_scale[0], median_size=median_size, gaussian_sigma=gaussian_sigma)
                source_images.append(bg_removed)
        sink_image = np.concatenate(sink_images, axis=0)
        source_image = np.concatenate(source_images, axis=0)
        train_images = [sink_image, source_image]
        print(f"sink_image shape: {sink_image.shape}")
        print(f"source_image shape: {source_image.shape}")
        print("Training first unmixing model...")
        if training_log:
            model1 = PICASSOnn(mixing_matrix)
            _ = list(model1.train_loop(train_images, max_iter=1000, lr=0.001))
        else:
            with suppress_stdout():
                model1 = PICASSOnn(mixing_matrix)
                _ = list(model1.train_loop(train_images, max_iter=1000, lr=0.001))        
        alpha1 = model1.mixing_parameters[0, :]
        bg1_model = model1.mixing_parameters[1, :].T

        # Second unmix: sink=unmixed1, source=ch1 (bg removed, median filtered, normalized, scaled)
        sink_images2 = []
        source_images2 = []
        for t in train_t:
            t_data = lazy_data[t].compute()
            for z in train_z:
                # Apply first unmix to get sink for second unmix
                sink_im = t_data[sink_channel][z]
                raw_source0 = t_data[source_channels[0]][z] #* args.sum_scale0
                src0 = remove_background_and_enhance(raw_source0, bg0, max0, scale=sum_scale[0],median_size=median_size, gaussian_sigma=gaussian_sigma)
                ims1 = np.array([sink_im.flatten(), src0.flatten()]).T
                unmixed1 = (ims1 - bg1_model).clip(min=0).dot(alpha1)
                unmixed1 = unmixed1.clip(min=0).astype('uint16').reshape(sink_im.shape)
                sink_images2.append(unmixed1)
                # Source for second unmix
                raw_source1 = t_data[source_channels[1]][z] #* args.sum_scale1
                src1 = remove_background_and_enhance(raw_source1, bg1, max1, scale=sum_scale[1],median_size=median_size, gaussian_sigma=gaussian_sigma)
                source_images2.append(src1)
        sink_image2 = np.concatenate(sink_images2, axis=0)
        source_image2 = np.concatenate(source_images2, axis=0)
        train_images2 = [sink_image2, source_image2]
        print(f"sink_image2 shape: {sink_image2.shape}")
        print(f"source_image2 shape: {source_image2.shape}")
        print("Training second unmixing model...")
        if training_log:
            model2 = PICASSOnn(mixing_matrix)
            _ = list(model2.train_loop(train_images2, max_iter=1000, lr=0.001))
        else:
            with suppress_stdout():
                model2 = PICASSOnn(mixing_matrix)
                _ = list(model2.train_loop(train_images2, max_iter=1000, lr=0.001))
        alpha2 = model2.mixing_parameters[0, :]
        bg2 = model2.mixing_parameters[1, :].T

        # --- INFERENCE: apply learned parameters to each (t, z) independently ---
        output_zarr = zarr.open(output_zarr_path, mode='w', shape=(N, Z, Y, X), dtype='uint16', chunks=(1, 1, Y, X))
        for t in tqdm(range(N), desc="Processing timepoints"):
            t_data = lazy_data[t].compute()
            sink_stack = t_data[sink_channel]
            source_stack0 = t_data[source_channels[0]] #* sum_scale[0]
            source_stack1 = t_data[source_channels[1]] #* sum_scale[1]
            unmixed1 = np.zeros_like(sink_stack, dtype=np.uint16)
            unmixed2 = np.zeros_like(sink_stack, dtype=np.uint16)
            for z in range(sink_stack.shape[0]):
                # First unmix
                sink_im = sink_stack[z]
                src0 = remove_background_and_enhance(source_stack0[z], bg0, max0, scale=sum_scale[0], median_size=median_size, gaussian_sigma=gaussian_sigma)
                ims1 = np.array([sink_im.flatten(), src0.flatten()]).T
                out1 = (ims1 - bg1_model).clip(min=0).dot(alpha1)
                out1 = out1.clip(min=0).astype('uint16').reshape(sink_im.shape)
                unmixed1[z] = out1
                # Second unmix
                src1 = remove_background_and_enhance(source_stack1[z], bg1, max1, scale=sum_scale[1], median_size=median_size, gaussian_sigma=gaussian_sigma)
                ims2 = np.array([out1.flatten(), src1.flatten()]).T
                out2 = (ims2 - bg2).clip(min=0).dot(alpha2)
                out2 = out2.clip(min=0).astype('uint16').reshape(sink_im.shape)
                unmixed2[z] = out2
            output_zarr[t] = unmixed2
        print(f"All timepoints processed and saved to {output_zarr_path}")

        # --- Visualize the result in napari ---
        if visualize_napari:
            try:
                original = da.from_zarr(raw_image_zarr)[:N].compute()  # (N, C, Z, Y, X)
                unmixed = da.from_zarr(output_zarr_path)[:N].compute()  # (N, Z, Y, X)
                unmixed = unmixed[:, None, ...]  # (N, 1, Z, Y, X)

                # Compute enhanced ch0 and ch1 for visualization
                enhanced_ch0 = np.zeros_like(original[:, source_channels[0]], dtype=np.uint16)  # (N, Z, Y, X)
                enhanced_ch1 = np.zeros_like(original[:, source_channels[1]], dtype=np.uint16)
                for t in range(N):
                    for z in range(original.shape[2]):
                        enhanced_ch0[t, z] = remove_background_and_enhance(original[t, source_channels[0], z], bg0, max0, scale=sum_scale[0], median_size=median_size, gaussian_sigma=gaussian_sigma)
                        enhanced_ch1[t, z] = remove_background_and_enhance(original[t, source_channels[1], z], bg1, max1, scale=sum_scale[1], median_size=median_size, gaussian_sigma=gaussian_sigma)

                available_colormaps = [
                    "magenta", "green", "gray", "yellow", "red", "cyan", "orange", "blue"
                ]
                viewer = napari.Viewer()
                # Visualize all original channels as before
                for c in range(original.shape[1]):
                    viewer.add_image(
                        original[:, c],
                        name=f'Original ch{c}',
                        contrast_limits=[np.min(original[:, c]), np.max(original[:, c])],
                        colormap=available_colormaps[c % len(available_colormaps)],
                        blending='additive',
                        metadata={'axes': 'TZYX'}
                    )
                # Add enhanced ch0 and ch1 as extra layers
                viewer.add_image(
                    enhanced_ch0,
                    name=f'Enhanced source ch{source_channels[0]}',
                    contrast_limits=[0, 65535],
                    colormap='magenta',
                    blending='additive',
                    metadata={'axes': 'TZYX'}
                )
                viewer.add_image(
                    enhanced_ch1,
                    name=f'Enhanced source ch{source_channels[1]}',
                    contrast_limits=[0, 65535],
                    colormap='green',
                    blending='additive',
                    metadata={'axes': 'TZYX'}
                )
                # Visualize final unmixed result
                viewer.add_image(
                    unmixed[:, 0],
                    name='Unmixed Final',
                    contrast_limits=[0, 65535],
                    colormap='yellow',
                    blending='additive',
                    metadata={'axes': 'TZYX'}
                )
                napari.run()
            except Exception as e:
                print(f"Napari visualization skipped: {e}")
        else:
            print("Napari visualization skipped as per user setting.")

        # At the end, save all outputs to a new Zarr group
        # Combine unmixed with original channels and save
        original = da.from_zarr(raw_image_zarr)[:N].compute()  # (N, C, Z, Y, X)
        unmixed = da.from_zarr(output_zarr_path)[:N].compute()  # (N, Z, Y, X)
        final_output = np.concatenate((original, unmixed[:, None, ...]), axis=1)  # (N, C+1, Z, Y, X)
        chunksize = (1,) + final_output.shape[1:]
        save_as_zarr(
            img=final_output, 
            path=all_outputs_path, 
            chunks=chunksize
            )
        print(f"Saved all outputs to {all_outputs_path}")

        ## Update metadata and return
        metadata.at[idx, 'tcell_channel'] = final_output.shape[1]-1  # Update tcell channel index
        metadata.at[idx, 'raw_image_path'] = str(all_outputs_path)

    return metadata
