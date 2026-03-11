import os

from pathlib import Path
import numpy as np
import pandas as pd
import time
import shutil

import napari
from magicgui import magicgui
from magicgui.widgets import PushButton
from qtpy.QtWidgets import QPlainTextEdit, QWidget, QVBoxLayout, QApplication, QLabel, QCheckBox, QDoubleSpinBox, QSpinBox, QPushButton as QtPushButton

from aicspylibczi import CziFile

from skimage import data, segmentation, feature, future
from skimage import filters as skimage_filters
from skimage.measure import label
from skimage.segmentation import watershed, relabel_sequential

from sklearn.ensemble import RandomForestClassifier
# from napari_apoc import PixelClassifier  # Commented out - use scikit-learn instead
from scipy import ndimage
from scipy.ndimage import binary_fill_holes, find_objects

from behav3d.preprocessing.segmentation import segment_size_filter, get_border_segments, remove_boundary_segments, segment_2d_filter
from behav3d.preprocessing import open_mask, dilate_mask, calculate_edt, zeropad_image_to_match_shape
from behav3d.io.images import load_image, get_image_shape, load_zarr, save_as_zarr, append_to_zarr

import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from tqdm import tqdm

import joblib
from functools import partial

import zarr
import dask.array as da
import gc

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False

try:
    import threadpoolctl
    _HAS_THREADPOOLCTL = True
except ImportError:
    _HAS_THREADPOOLCTL = False


def _log_mem(tag, log_func=print):
    """Log current process memory usage (RSS in MB)."""
    if _HAS_PSUTIL:
        proc = psutil.Process(os.getpid())
        rss_mb = proc.memory_info().rss / (1024 * 1024)
        log_func(f"[MEM] {tag}: {rss_mb:.0f} MB")
    else:
        log_func(f"[MEM] {tag}: psutil not installed, skipping memory log")

## TODO create a function for BEHAV3D notebook
'''
def calculate_image_features(
    image,
    intensity=True,
    edges=True,
    texture=False,
    # sigma_min=sigma_min,
    sigma_max=16,
    channel_axis=0,
    ):
    """
    Calculate multiscale basic features for a given image.
    Possible to save as a joblib and load later
    """
    img = np.asarray(image)
    orig_shape = img.shape
    ndim = img.ndim
    
    # Identify spatial axes (everything except channel_axis)
    if channel_axis is None:
        spatial_axes = list(range(ndim))
    else:
        spatial_axes = [ax for ax in range(ndim) if ax != channel_axis]
    
    # If image is 2D, squeeze pseudo-3D-axis
    squeeze_axes = [ax for ax in spatial_axes if orig_shape[ax] == 1]
    if squeeze_axes:
        img_squeezed = np.squeeze(img, axis=tuple(squeeze_axes))
    else:
        img_squeezed = img
        
    feats = features_func(
        img_squeezed,
        intensity=intensity,
        edges=edges,
        texture=texture,
        channel_axis=channel_axis,
        sigma_max=sigma_max,
        )
    
    n_features = feats.shape[-1]
    orig_spatial_shape = [orig_shape[ax] for ax in spatial_axes]
    feats_restored = feats.reshape((*orig_spatial_shape, n_features))
    return feats_restored
'''
# ---------------------------------------------------------------------------
# Sigma values in µm — optimised via random-forest feature importance.
# They are converted to pixel units at runtime using the sample's pixel size.
# ---------------------------------------------------------------------------
SIGMA_UM = [0.89, 1.77, 3.54, 14.16]


def make_features(pixel_size_xy_um=1.0, pixel_size_z_um=1.0, verbose=True):
    """
    Build a feature extractor whose Gaussian sigmas are expressed in **pixels**,
    derived from the physical SIGMA_UM values and the given pixel sizes.

    Parameters
    ----------
    pixel_size_xy_um : float
        XY pixel size in µm/px.
    pixel_size_z_um : float
        Z pixel size in µm/px.
    verbose : bool
        If True, print the sigma conversion table.

    Returns
    -------
    callable
        A function with signature ``f(image)`` that returns multiscale features.
    """
    sigmas_px = []
    for s_um in SIGMA_UM:
        sz = s_um / pixel_size_z_um
        sxy = s_um / pixel_size_xy_um
        sigmas_px.append((sz, sxy, sxy))

    if verbose:
        print("Sigma values (µm → px) used for feature extraction:")
        print(f"  Pixel size  :  xy = {pixel_size_xy_um:.4f} µm/px   z = {pixel_size_z_um:.4f} µm/px")
        print(f"  {'σ (µm)':>10s}  {'σ_z (px)':>10s}  {'σ_xy (px)':>10s}")
        print(f"  {'-'*36}")
        for s_um, (sz, sxy, _) in zip(SIGMA_UM, sigmas_px):
            print(f"  {s_um:10.2f}  {sz:10.2f}  {sxy:10.2f}")

    def _features(image):
        all_features = []
        for sigma in sigmas_px:
            for func in [skimage_filters.gaussian,
                         skimage_filters.sobel]:
                if image.ndim > 3 and image.shape[0] > 1:  # channel axis present
                    per_ch = []
                    for ch in range(image.shape[0]):
                        filtered = func(image[ch], sigma=sigma) if func is skimage_filters.gaussian else func(skimage_filters.gaussian(image[ch], sigma=sigma))
                        per_ch.append(filtered)
                    all_features.extend(per_ch)
                else:
                    filtered = func(image, sigma=sigma) if func is skimage_filters.gaussian else func(skimage_filters.gaussian(image, sigma=sigma))
                    all_features.append(filtered)
        return np.stack(all_features, axis=-1)

    return _features


def format_sigma_table(pixel_size_xy_um, pixel_size_z_um):
    """Return a formatted string showing σ (µm) → σ (px) conversion."""
    lines = []
    lines.append("Sigma values (µm → px) used for feature extraction:")
    lines.append(f"  Pixel size  :  xy = {pixel_size_xy_um:.4f} µm/px   z = {pixel_size_z_um:.4f} µm/px")
    lines.append(f"  {'σ (µm)':>10s}  {'σ_z (px)':>10s}  {'σ_xy (px)':>10s}")
    lines.append(f"  {'-'*36}")
    for s_um in SIGMA_UM:
        sz = s_um / pixel_size_z_um
        sxy = s_um / pixel_size_xy_um
        lines.append(f"  {s_um:10.2f}  {sz:10.2f}  {sxy:10.2f}")
    lines.append(f"  SIGMA_UM = {SIGMA_UM}")
    return "\n".join(lines)


# Module-level fallback: pixel_size=1.0 → sigmas in µm used as-is in pixels.
# Wraps the real function so it prints a warning the first time it's called.
_fallback_compute = make_features(1.0, 1.0, verbose=False)
_fallback_warned = False

def compute_features(image):
    """Fallback feature extractor — warns that pixel sizes were not provided."""
    global _fallback_warned
    if not _fallback_warned:
        print(
            "⚠️  compute_features called without pixel-size calibration.\n"
            "    Sigma values are being used in µm as-is (pixel_size = 1.0).\n"
            "    Please provide pixel_distance_xy / pixel_distance_z in the metadata\n"
            "    so that sigmas are correctly converted to pixels."
        )
        _fallback_warned = True
    return _fallback_compute(image)


def postprocess_mask(mask, fill_holes=True, opening_nr_pixels=1):
    """
    Postprocess a binary mask by filling holes and/or applying opening.
    
    Parameters
    ----------
    mask : np.ndarray
        Binary mask (bool or 0/1). Can be 2D (Y, X), 3D (Z, Y, X), or 4D (T, Z, Y, X)
    fill_holes : bool
        Whether to fill holes in the mask
    opening_nr_pixels : int
        Number of pixels for morphological opening (0 to disable)
    
    Returns
    -------
    np.ndarray
        Postprocessed binary mask (same shape and dtype as input)
    """
    mask = mask.astype(bool)  # Ensure binary
    
    if fill_holes:  
        # Handle different dimensionalities
        if mask.ndim == 2:  # (Y, X)
            mask = binary_fill_holes(mask)
        elif mask.ndim == 3:  # (Z, Y, X)
            filled_mask = np.zeros_like(mask)
            for i in range(mask.shape[0]):
                filled_mask[i] = binary_fill_holes(mask[i])
            mask = filled_mask
        elif mask.ndim == 4:  # (T, Z, Y, X)
            filled_mask = np.zeros_like(mask)
            for t in range(mask.shape[0]):
                for z in range(mask.shape[1]):
                    filled_mask[t, z] = binary_fill_holes(mask[t, z])
            mask = filled_mask
        else:
            raise ValueError(f"Unsupported mask dimensionality: {mask.ndim}D")
    
    if opening_nr_pixels > 0 and opening_nr_pixels is not None:
        mask = open_mask(mask, nr_pixels=opening_nr_pixels)
    
    return mask.astype(bool)

def refine_segments(
    segments: np.ndarray,
    edt: np.ndarray,
    mask: np.ndarray,
    thr: float,
    segment_size_min: int,
    n_workers: int | None = 1,
    out_dtype=np.uint16,
):
    """
    One refinement pass:
      - relabel segments sequentially (so find_objects indices match label ids)
      - for each label, run refine_segment(local_mask, local_edt, thr, ...)
      - stitch refined local labels into a global label image with unique IDs

    Returns:
      new_segments (same shape as mask)
    """
    # Make labels contiguous so: slices[i] corresponds to label (i+1)
    segments, _, _ = relabel_sequential(segments)

    slices = find_objects(segments)
    args_list = []

    for i, slc in enumerate(slices):
        if slc is None:
            continue

        label_id = i + 1
        local_mask = (segments[slc] == label_id)
        if not np.any(local_mask):
            continue

        local_edt = edt[slc]
        minc = [s.start for s in slc]
        args_list.append((local_mask, local_edt, thr, segment_size_min, minc))

    if n_workers is None:
        n_workers = multiprocessing.cpu_count()

    if len(args_list) == 0:
        return np.zeros_like(mask, dtype=out_dtype)

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        results = list(executor.map(_refine_segment, args_list))

    new_segments = np.zeros_like(mask, dtype=np.uint32)  # safer while assembling
    current_label = 0

    for local_segment, minc in results:
        z0, y0, x0 = minc

        nonzero = local_segment > 0
        if not np.any(nonzero):
            continue

        local_segment = local_segment.astype(np.uint32, copy=False)
        local_max = int(local_segment.max())

        local_segment[nonzero] += current_label

        zz, yy, xx = np.where(nonzero)
        new_segments[z0 + zz, y0 + yy, x0 + xx] = local_segment[nonzero]

        current_label += local_max

    if out_dtype is not None:
        if np.iinfo(out_dtype).max < new_segments.max():
            raise ValueError(
                f"Label overflow: max label {new_segments.max()} exceeds dtype {out_dtype}"
            )
        new_segments = new_segments.astype(out_dtype, copy=False)

    return new_segments

def _refine_segment(args):
    """Refine a single segment by reapplying EDT-based splitting."""
    # start_time = time.time()
    # label_id, segment_mask, full_edt, edt_threshold, segment_size_min = args
    local_mask, local_edt, edt_thr_refined, segment_size_min, minc = args
    # if np.sum(segment_mask) == 0:
    #     return np.zeros_like(segment_mask)
   
    local_seeds = label(local_edt >= edt_thr_refined)
    # print("###", label_id, "refine_segment time elapsed: ", time.time() - start_time)
    if np.max(local_seeds) < 2:
        return (local_mask.astype(np.int32), tuple(minc))
    
    new_seg = watershed(-local_edt, markers=local_seeds, mask=local_mask)
    # print("###", label_id, "refine_segment time elapsed: ", time.time() - start_time)
    new_seg = segment_size_filter(new_seg, size_min=segment_size_min)
    
    # If EDT threshold is too high and all segments were removed, use original segment
    if np.max(new_seg) == 0:
        return (local_mask.astype(np.int32), tuple(minc))
    
    # print("###", label_id, "refine_segment time elapsed: ", time.time() - start_time)
    new_seg = watershed(-local_edt, markers=new_seg, mask=local_mask)
    # print("###", label_id, "refine_segment time elapsed: ", time.time() - start_time)
    new_seg, _, _ = relabel_sequential(new_seg)
    # print("###", label_id, "refine_segment time elapsed: ", time.time() - start_time)
    return (new_seg, tuple(minc))

def segment_mask(mask, edt_thr=1.5, edt_thr_refined=None, segment_size_min=15, use_dims=2, n_workers=1):
    edt = calculate_edt(mask, use_dims=use_dims)

    # Normalize refined thresholds
    if edt_thr_refined is not None and not isinstance(edt_thr_refined, list):
        edt_thr_refined = [edt_thr_refined]

    # Start from connected components as initial segments
    segments = label(mask)

    # Pass 1: use edt_thr
    segments = refine_segments(
        segments=segments,
        edt=edt,
        mask=mask,
        thr=edt_thr,
        segment_size_min=segment_size_min,
        n_workers=n_workers,
        out_dtype=np.uint16,
    )

    # Pass 2..N: refined thresholds
    if edt_thr_refined is not None:
        for thr in edt_thr_refined:
            segments = refine_segments(
                segments=segments,
                edt=edt,
                mask=mask,
                thr=thr,
                segment_size_min=segment_size_min,
                n_workers=n_workers,
                out_dtype=np.uint16,
            )
    segments = segment_size_filter(segments, size_min=segment_size_min)
    segments = segment_2d_filter(segments)
    return segments


def _apply_classifier(args):
    clf, path, outpath, idx = args
    features = np.asarray(load_image(path, mode="r")[idx])
    prediction = zarr.open(outpath, mode="r+")
    prediction[idx] = future.predict_segmenter(features, clf)

# def apply_classifier_in_memory(args):
#     clf, path, outpath, idx = args
#     features = np.asarray(load_image(path, mode="r")[idx])
#     prediction = zarr.open(outpath, mode="r+")
#     prediction[idx] = future.predict_segmenter(features, clf)


def apply_classifier(classifier, features_outpath, pred_labels_outpath, n_workers=4):
    shape = load_image(features_outpath).shape
    pred_labels = da.zeros(shape[:-1], chunks=(1,) + shape[1:-1], dtype='int16')
    save_as_zarr(pred_labels, pred_labels_outpath)

    # Disable internal RF parallelism to prevent nested thread explosion
    # (ThreadPoolExecutor already parallelizes across timepoints)
    original_n_jobs = getattr(classifier, 'n_jobs', 1)
    classifier.n_jobs = 1

    args_list = [(classifier, str(features_outpath), str(pred_labels_outpath), idx) for idx in range(pred_labels.shape[0])]
    test=[]
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        test+=list(tqdm(executor.map(_apply_classifier, args_list), total=len(args_list)))

    classifier.n_jobs = original_n_jobs
    
    pixel_class = np.asarray(load_image(pred_labels_outpath))
    
    return pixel_class
        
def train_pixel_classifier(
    output_dir,
    metadata=None,
    examples_per_sample = 3,
    sample_specific_classifier=False,
    n_workers=None,
    manual_dim_order=None,
    images=None,
    #two_org_types=False,  # Deprecated
    organoid_types=None,
    immune_types=None,
    other_types=None,
    initial_params=None,
    on_params_changed=None,
    ):
    """
    Train a pixel classifier for segmentation with support for multiple cell types.
    Args:
        output_dir: Output directory
        metadata: Metadata DataFrame
        examples_per_sample: Number of timepoints per sample
        sample_specific_classifier: Whether to use sample-specific classifier
        n_workers: Number of workers
        manual_dim_order: Optional. Tuple/list specifying the order for transpose, e.g. (1,0,2,3,4) for (C,T,Z,Y,X)
        two_org_types: Deprecated, use organoid_types instead
        organoid_types: List of organoid type names (e.g., ['organoid1', 'organoid2'])
        immune_types: List of immune cell type names (e.g., ['tcell', 'macro'])
        other_types: List of other cell type names (e.g., ['tum'])
    """
    assert images is not None or metadata is not None, "Either supply a BEHAV3D metadata table or a list of images"
    if n_workers is None:
        n_workers = multiprocessing.cpu_count()
    
    # Validate and cap n_workers to available CPU count
    max_workers = multiprocessing.cpu_count()
    if n_workers > max_workers:
        safe_workers = max(1, max_workers // 2)
        print(f"⚠️ Requested {n_workers} workers but only {max_workers} CPUs available. Using {safe_workers} workers (half of available) for system stability.")
        n_workers = safe_workers
    elif n_workers < 1:
        print(f"⚠️ Invalid n_workers={n_workers}. Using 1 worker (sequential).")
        n_workers = 1
    
    # Detect cell types from metadata if not provided
    from behav3d.core.metadata import (
        detect_organoid_types_from_metadata,
        detect_immune_cell_types_from_metadata,
        detect_other_cell_types_from_metadata,
        has_dead_channel
    )
    
    if metadata is not None:
        if organoid_types is None:
            organoid_types = detect_organoid_types_from_metadata(metadata)
        if immune_types is None:
            immune_types = detect_immune_cell_types_from_metadata(metadata)
        if other_types is None:
            other_types = detect_other_cell_types_from_metadata(metadata)
        has_death = has_dead_channel(metadata)
    else:
        has_death = False  # Default to False if no metadata
    
    # Default to empty lists if still None
    organoid_types = organoid_types or []
    immune_types = immune_types or []
    other_types = other_types or []
    
    # Combine all cell types
    all_cell_types = organoid_types + immune_types + other_types
    
    print(f"Training classifier for cell types:")
    print(f"   Organoids: {organoid_types}")
    print(f"   Immune: {immune_types}")
    print(f"   Other: {other_types}")
    print(f"Dead channel: {has_death}")
        
    pixel_class_outdir = Path(output_dir, "images", "PixelClassification")
    pixel_class_outdir.mkdir(exist_ok=True, parents=True)
    
    features_outpath = Path(pixel_class_outdir, 'PixelClassifier_Features.zarr')
    image_outpath = Path(pixel_class_outdir, 'PixelClassifier_Images.zarr')
    
    
    if images is None:
        if not features_outpath.exists() or not image_outpath.exists():
            if image_outpath.exists():
                shutil.rmtree(image_outpath)
            if features_outpath.exists():
                shutil.rmtree(features_outpath)
            
            all_images = []
            all_features = []
            
            channels = []
            for _, sample in metadata.iterrows():
                raw_image_path = Path(sample['raw_image_path'])
                raw_image_shape = get_image_shape(raw_image_path)
                nr_channels = raw_image_shape[-4]
                channels.append(nr_channels)
                
            shapes_info = "\n".join(
                [f"{sample['sample_name']}: {get_image_shape(Path(sample['raw_image_path']))}"
                for _, sample in metadata.iterrows()]
            )

            assert len(set(channels)) == 1, (
                "All samples must have the same number of channels.\n"
                f"Image shapes:\n{shapes_info}"
            )
                
            for idx, sample in metadata.iterrows():
                
                sample_name = sample['sample_name']
                        
                raw_image_path = Path(sample['raw_image_path'])
                raw_image_zarr =  Path(output_dir, "images", sample_name, f"{sample_name}.zarr")
                
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
                    # else:
                    #     print("No manual_dim_order provided, assuming image is already in default order.")
                    chunksize = (1,) + images.shape[1:]
                    save_as_zarr(
                        img=images, 
                        path=raw_image_zarr, 
                        chunks=chunksize
                    )
                            
                images = load_image(raw_image_zarr)
                max_t = images.shape[0]-1
                print(images.shape)
                idc = np.linspace(0, max_t, examples_per_sample, dtype=int)
                print(f"Taking timepoints: {idc}")
                
                # sample_images = [images[t, [tcell_ch, live_ch, dead_ch]] for t in idc]
                sample_images = [images[t] for t in idc]
                all_images+=sample_images
            
            # Calculate the maximum shape that the images need to be in each dimension
            max_shape = tuple(max(img.shape[i] for img in all_images) for i in range(images[0].ndim))

            # --- Build feature extractor from metadata pixel sizes ---
            pixel_sizes_xy = metadata['pixel_distance_xy'].unique()
            pixel_sizes_z  = metadata['pixel_distance_z'].unique()
            assert len(pixel_sizes_xy) == 1 and len(pixel_sizes_z) == 1, (
                "All samples must have the same pixel size for consistent feature extraction.\n"
                f"  pixel_distance_xy values: {pixel_sizes_xy}\n"
                f"  pixel_distance_z  values: {pixel_sizes_z}"
            )
            pixel_size_xy_um = float(pixel_sizes_xy[0])
            pixel_size_z_um  = float(pixel_sizes_z[0])
            extract_features = make_features(pixel_size_xy_um, pixel_size_z_um, verbose=True)

            print(f"Calculating features...")
            for img in tqdm(all_images):
                append_to_zarr(
                    np.expand_dims(
                        zeropad_image_to_match_shape(
                            img = extract_features(img),
                            target_shape=max_shape[-3:],
                            axes=[-4, -3, -2]
                        ), axis=0), 
                    features_outpath
                    )
            
            
            all_images = [zeropad_image_to_match_shape(img, max_shape) for img in all_images]
            
            print(f"Max shape found for padding input images: {max_shape}")
            # Allow manual override of dimension order for transpose
            all_images = da.stack(all_images)
            all_images = all_images.transpose(1, 0, 2, 3, 4) ## this order corresponds to (C, T, Z, Y, X)
            save_as_zarr(all_images, image_outpath)
            del all_images
            gc.collect()
            
            all_images = load_zarr(image_outpath)
            all_features = load_zarr(features_outpath)
        else:
            all_images = load_image(image_outpath)
            all_features = load_image(features_outpath) 
        all_images = np.asarray(all_images)
    else:     
        all_images = []
        all_features = []
        for img in images:
            feature_img = compute_features(img)
            all_images.append(img)
            all_features.append(feature_img)
        
        all_images = np.stack(all_images, axis=0)
        all_images = all_images.transpose(1, 0, 2, 3, 4)
        save_as_zarr(all_images, image_outpath)
        all_features = np.stack(all_features, axis=0)
        save_as_zarr(all_features, features_outpath)
    
    def segment_and_update(
        pixel_class_outdir,
        only_segment=False,
        n_workers: int = 16,
        log=print,
        **kwargs  # Dynamic parameters: edt_thresholds, segment_size_mins, opening_nr_pixels, fill_holes
        ):
        """
        Dynamic segmentation function that handles different cell types.
        Parameters are passed as kwargs with prefixes:
        - {cell_type}_edt_threshold: EDT threshold for segmentation
        - {cell_type}_segment_size_min: Minimum segment size
        - {cell_type}_opening_nr_pixels: Mask opening pixels
        - {cell_type}_fill_holes: Whether to fill holes (bool)
        """
        print("segment_and_update() called!")
        
        # Parse kwargs into separate dictionaries
        edt_thresholds = {}
        segment_size_mins = {}
        opening_nr_pixels_dict = {}
        fill_holes_dict = {}
        
        # Define suffixes for parameter parsing
        SUFFIX_EDT = '_edt_threshold'
        SUFFIX_SIZE = '_segment_size_min'
        SUFFIX_OPENING = '_opening_nr_pixels'
        SUFFIX_FILLHOLES = '_fill_holes'
        
        for key, value in kwargs.items():
            if key.endswith(SUFFIX_EDT):
                cell_type = key[:-len(SUFFIX_EDT)]
                edt_thresholds[cell_type] = value
            elif key.endswith(SUFFIX_SIZE):
                cell_type = key[:-len(SUFFIX_SIZE)]
                segment_size_mins[cell_type] = value
            elif key.endswith(SUFFIX_OPENING):
                cell_type = key[:-len(SUFFIX_OPENING)]
                opening_nr_pixels_dict[cell_type] = value
            elif key.endswith(SUFFIX_FILLHOLES):
                cell_type = key[:-len(SUFFIX_FILLHOLES)]
                fill_holes_dict[cell_type] = value
        
        print(f"Received EDT thresholds: {edt_thresholds}")
        print(f"Received segment size mins: {segment_size_mins}")
        print(f"Received opening pixels: {opening_nr_pixels_dict}")
        print(f"Received fill holes: {fill_holes_dict}")
        
        start_time = time.time()
        _log_mem("segment_and_update START", log)
        log("###### Running Segmentation\n")
        log(f"EDT Thresholds: {edt_thresholds}\n")
        log(f"Segment Size Mins: {segment_size_mins}\n")
        log(f"Opening Pixels: {opening_nr_pixels_dict}\n")
        log(f"Fill Holes: {fill_holes_dict}\n")
        log(f"only_segment = {only_segment}\n")
        log(f"all_cell_types = {all_cell_types}\n")
        QApplication.processEvents()
        
        # Access the label layer and feature image
        image_data = all_images  # Use the original all_images data
        log(f"image_data shape: {image_data.shape if hasattr(image_data, 'shape') else 'N/A'}\n")
        log(f"all_features shape: {all_features.shape if hasattr(all_features, 'shape') else 'N/A'}\n")
        QApplication.processEvents()
        
        # Save user labels for all cell types
        cell_type_labels = {}  # Store {cell_type: label_data}
        
        # Save death labels (only if present)
        if has_death:
            log("Saving death labels...\n")
            QApplication.processEvents()
            dead_label_layer = viewer.layers['User Provided Labels (Dead)']
            dead_label_data = dead_label_layer.data
            dead_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
            # Remove existing zarr to avoid ContainsArrayError
            if dead_labels_outpath.exists():
                shutil.rmtree(dead_labels_outpath)
            save_as_zarr(dead_label_data, dead_labels_outpath)
            log("Saved death labels\n")
            QApplication.processEvents()
        
        # Save dynamic cell type labels and apply postprocessing for training
        for cell_type in all_cell_types:
            layer_name = f'User Provided Labels ({cell_type.capitalize()})'
            label_layer = viewer.layers[layer_name]
            label_data = label_layer.data.copy()
            
            # Apply postprocessing to user labels before training
            # Extract foreground mask (label 2) for postprocessing
            # Get postprocessing parameters for this cell type
            # opening_nr_pixels = opening_nr_pixels_dict.get(cell_type, 
            #     (3 if cell_type in organoid_types else 0))
            # fill_holes = fill_holes_dict.get(cell_type, True)
            
            # # Extract foreground mask and postprocess
            # fg_mask = (label_data == 2).astype(bool)
            # if np.any(fg_mask):
            #     processed_fg = postprocess_mask(
            #         fg_mask, 
            #         fill_holes=fill_holes, 
            #         opening_nr_pixels=int(opening_nr_pixels)
            #     )
            #     # Reconstruct labels: background=1, foreground=2
            #     label_data = np.where(processed_fg, 2, 
            #                          np.where(label_data == 1, 1, 0)).astype(label_data.dtype)
            
            labels_outpath = Path(pixel_class_outdir, f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr')
            # Remove existing zarr to avoid ContainsArrayError
            if labels_outpath.exists():
                shutil.rmtree(labels_outpath)
            save_as_zarr(label_data, labels_outpath)
            cell_type_labels[cell_type] = label_data
            log(f"Saved {cell_type} labels (postprocessed for training)\n")
            QApplication.processEvents()


        def train_classifier(user_labels, features):
            flat_label_data = user_labels.ravel()
            flat_features = features.reshape(-1, features.shape[-1])  # shape: (N_total, 90)

            # Get 1D indices where labels exist
            label_indices = np.flatnonzero(flat_label_data > 0)

            selected_features = np.asarray(flat_features[label_indices]) # (N_selected, 90)
            selected_labels = flat_label_data[label_indices]   
            
            nr_bg_pix = int(np.sum(selected_labels==1))
            nr_fg_pix = int(np.sum(selected_labels==2))
            total_pix = nr_bg_pix + nr_fg_pix
            
            log(f"Found {nr_bg_pix} background pixels")
            log(f"Found {nr_fg_pix} foreground pixels")
            
            class_weights = {
                1: nr_bg_pix / total_pix,
                2: nr_fg_pix / total_pix,
            }
            # Use scikit-learn RandomForestClassifier (fallback if napari_apoc not available)
            clf = RandomForestClassifier(
                n_estimators=50,
                n_jobs=-1, 
                max_depth=20, 
                # max_samples=0.05,
                class_weight=class_weights
            )
            
            clf = future.fit_segmenter(selected_labels, selected_features, clf)
            return clf
        
        # Train classifiers for all cell types
        classifiers = {}  # Store {cell_type: classifier}
        clf_death = None  # Initialize to None
        
        _log_mem("before classifier training", log)
        log(f"\n--- Starting classifier training (only_segment={only_segment}) ---\n")
        QApplication.processEvents()
        
        if not only_segment:
            # Train death classifier (only if death channel is present)
            if has_death:
                log("\n### Training Random Forest Classifier (Cell Death)")
                clf_death = train_classifier(dead_label_data, all_features)
                death_random_forest_outpath = Path(pixel_class_outdir, 'PixelClassifier_Death.joblib')
                log(f"Saving to {death_random_forest_outpath}")
                joblib.dump(clf_death, death_random_forest_outpath)
                QApplication.processEvents()
            else:
                log("\n### Skipping death classifier (no dead channel)")
                QApplication.processEvents()
        
            # Train classifiers for all detected cell types
            for cell_type in all_cell_types:
                log(f"\n### Training Random Forest Classifier ({cell_type.capitalize()})")
                clf = train_classifier(cell_type_labels[cell_type], all_features)
                clf_outpath = Path(pixel_class_outdir, f'PixelClassifier_{cell_type.capitalize()}.joblib')
                log(f"Saving to {clf_outpath}")
                joblib.dump(clf, clf_outpath)
                classifiers[cell_type] = clf
                QApplication.processEvents()
        
        
        # Apply classifiers to predict pixels
        pred_masks = {}  # Store {cell_type: predicted_mask}
        
        # Define prediction paths
        pred_death_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_Death_PredictedLabels.zarr')
        pred_labels_paths = {
            cell_type: Path(pixel_class_outdir, f'PixelClassifier_{cell_type.capitalize()}_PredictedLabels.zarr')
            for cell_type in all_cell_types
        }
        
        if not only_segment:
            # Clear old predictions
            if has_death and pred_death_labels_outpath.exists():
                shutil.rmtree(pred_death_labels_outpath)
            for path in pred_labels_paths.values():
                if path.exists():
                    shutil.rmtree(path)
            
            # Predict death pixels (only if death channel is present)
            if has_death:
                log("\n### Predicting Death Pixels")
                _log_mem("before death prediction", log)
                QApplication.processEvents()
                log("   Starting apply_classifier for death...")
                QApplication.processEvents()
                pred_death_mask = apply_classifier(clf_death, features_outpath, pred_death_labels_outpath, n_workers=n_workers)
                _log_mem("after death prediction", log)
                log("    Death prediction finished!")
                QApplication.processEvents()
                viewer.layers["Pixel Classification (Dead)"].data = pred_death_mask
                log("    Death layer updated!")
                QApplication.processEvents()
            else:
                log("\n### Skipping death prediction (no dead channel)")
                QApplication.processEvents()

            # Predict pixels for all cell types
            for cell_type in all_cell_types:
                log(f"\n### Predicting {cell_type.capitalize()} Pixels")
                _log_mem(f"before {cell_type} prediction", log)
                QApplication.processEvents()
                pred_mask = apply_classifier(classifiers[cell_type], features_outpath, pred_labels_paths[cell_type], n_workers=n_workers)
                _log_mem(f"after {cell_type} prediction", log)

                pred_masks[cell_type] = pred_mask
                viewer.layers[f"Pixel Classification ({cell_type.capitalize()})"].data = pred_mask
        
        else:
            # Load existing predictions
            if has_death:
                log("\n### Loading Death Prediction Mask")
                QApplication.processEvents()
                pred_death_mask = viewer.layers["Pixel Classification (Dead)"].data
            
            # Load predictions for all cell types
            for cell_type in all_cell_types:
                log(f"\n### Loading {cell_type.capitalize()} Prediction Mask")
                QApplication.processEvents()
                pred_mask = viewer.layers[f"Pixel Classification ({cell_type.capitalize()})"].data
                pred_masks[cell_type] = pred_mask
            
        # Segment instances for all cell types
        _log_mem("before segmentation", log)
        log("\n### Segment Cell Instances")
        QApplication.processEvents()
        
        # Segment each cell type using EDT watershed
        segmented_cells = {}  # Store {cell_type: segmented_mask}
        for cell_type in all_cell_types:
            edt_threshold = edt_thresholds.get(cell_type, 
                (12.0 if cell_type in organoid_types else (2.5 if cell_type in immune_types else 1.0)))
            segment_size_min = segment_size_mins.get(cell_type,
                (1000 if cell_type in organoid_types else (10 if cell_type in immune_types else 10)))
              # Get postprocessing parameters for this cell type
            opening_nr_pixels = opening_nr_pixels_dict.get(cell_type, 
                (3 if cell_type in organoid_types else 0))
            fill_holes = fill_holes_dict.get(cell_type, True)
            
            log(f"\n### Segmenting {cell_type.capitalize()} (EDT threshold={edt_threshold}, min_size={segment_size_min})")
            QApplication.processEvents()

            pred_mask = pred_masks[cell_type]
            n_timepoints = pred_mask.shape[0]
            log(f"   Processing {n_timepoints} timepoints...")
            log(f"   Mask shape: {pred_mask.shape}")
            QApplication.processEvents()

            segmented_timepoints = []
            for t_idx in range(n_timepoints):
                log(f"   [T{t_idx+1}] Loading mask...")
                QApplication.processEvents()

                mask_t = pred_mask[t_idx]
                # Remove background (label 1), keep only foreground (label 2)
                mask_t = (mask_t == 2)  # bool mask is fine
                fg_pixels = int(mask_t.sum())
                log(f"   [T{t_idx+1}] Foreground pixels: {fg_pixels}")
                QApplication.processEvents()

                if fg_pixels == 0:
                    segmented = np.zeros_like(mask_t, dtype=np.uint16)
                    log(f"   [T{t_idx+1}] Empty mask, skipped")
                    QApplication.processEvents()
                else:
                    # Use your unified pipeline (EDT + refine passes)
                    log(f"   [T{t_idx+1}] Segmenting via segment_mask()...")
                    QApplication.processEvents()

                    # log(f"\n### Preprocessing {cell_type.capitalize()} classifier mask (opening_nr_pixels={opening_nr_pixels}, fill_holes={fill_holes})")
                    # Convert label mask to binary for postprocessing
                    # pred_mask has labels: 0=background, 1=background_label, 2=foreground
                    # Extract foreground mask (label 2) and postprocess
                    # fg_mask = (pred_mask == 2).astype(bool)
                    mask_t = postprocess_mask(mask_t, fill_holes=fill_holes, opening_nr_pixels=int(opening_nr_pixels))
                
                    segmented = segment_mask(
                        mask=mask_t,
                        edt_thr=edt_threshold,          # first pass threshold
                        edt_thr_refined=None,           # or e.g. [edt_threshold + 0.5, edt_threshold + 1.0]
                        segment_size_min=segment_size_min,
                        use_dims=2,                     # IMPORTANT: 2D frames (change to 3 for 3D)
                        n_workers=1,                    # or whatever you want
                    )

                    segmented = segmented.astype(np.uint16, copy=False)
                    log(f"   [T{t_idx+1}] Done! labels={int(segmented.max())}")
                    QApplication.processEvents()

                segmented_timepoints.append(segmented)

            full_seg = np.stack(segmented_timepoints, axis=0)
            segmented_cells[cell_type] = full_seg
            _log_mem(f"after {cell_type} segmentation (all timepoints)", log)
            log(f"   {cell_type.capitalize()} segmentation complete!")
            QApplication.processEvents()

            viewer.layers[f"{cell_type.capitalize()} Segments"].data = full_seg

        # Save images
        image_outpath = Path(pixel_class_outdir, 'PixelClassifier_Images.zarr')
        # Remove existing zarr to avoid ContainsArrayError
        if image_outpath.exists():
            shutil.rmtree(image_outpath)
        save_as_zarr(image_data, image_outpath)

        _log_mem("segment_and_update END", log)
        log(f"\n###### DONE time elapsed: {time.time() - start_time:.2f} s")
    
    def save_pixel_classification(log=print):
        """Save all user-provided labels for all cell types"""
        # Save death labels (only if present)
        if has_death:
            dead_label_layer = viewer.layers['User Provided Labels (Dead)']
            dead_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
            # Remove existing zarr to avoid ContainsArrayError
            if dead_labels_outpath.exists():
                shutil.rmtree(dead_labels_outpath)
            save_as_zarr(dead_label_layer.data, dead_labels_outpath)
            log(f"Saved death labels to: {dead_labels_outpath}")
        
        # Save labels for all cell types
        for cell_type in all_cell_types:
            layer_name = f'User Provided Labels ({cell_type.capitalize()})'
            label_layer = viewer.layers[layer_name]
            labels_outpath = Path(pixel_class_outdir, f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr')
            # Remove existing zarr to avoid ContainsArrayError
            if labels_outpath.exists():
                shutil.rmtree(labels_outpath)
            save_as_zarr(label_layer.data, labels_outpath)
            log(f"Saved {cell_type} labels to: {labels_outpath}")
        
        log("All user labels saved!")
            
    # Create Napari viewer
    viewer = napari.Viewer()
    
    # Split channels and add them as separate colored layers
    # all_images shape is (channels, time, z, y, x)
    n_channels = all_images.shape[0]
    
    # Define colors for different channels (extend if more channels exist)
    channel_colors = ['cyan', 'yellow', 'red', 'green', 'magenta', 'blue', 
                   'gray', 'turbo', 'viridis', 'plasma', 'inferno', 'twilight']

    for ch in range(n_channels):
        channel_data = all_images[ch]  # Shape: (time, z, y, x)
        
        # Flatten and filter out zero pixels
        flat_vals = np.asarray(channel_data).reshape(-1)
        nonzero_vals = flat_vals[flat_vals > 0]

        if nonzero_vals.size > 0:
            channel_percentile = float(np.percentile(nonzero_vals, 99.8))
            contrast_limits = (0, channel_percentile)
        else:
            print(f"⚠️ Channel {ch} appears empty. Using fallback contrast limits.")
            contrast_limits = (0, 1e-3)  # or any small dummy range
        
        # Add channel as separate layer with color
        channel_name = f"Channel {ch}"
        
        img_layer = viewer.add_image(
            channel_data, 
            name=channel_name, 
            contrast_limits=contrast_limits,
            colormap=channel_colors[ch] if ch < len(channel_colors) else 'gray',
            blending='additive', # This allows channels to blend together
            opacity=0.8
        )
        img_layer.contrast_limits_range = (0, float(channel_data.max()))

    # DYNAMIC LAYER CREATION FOR ALL CELL TYPES
    # Create user label layers and predicted label layers for each detected cell type
    user_layers = {}
    pixelclass_layers = {}
    segment_layers = {}
    
    # Add death labels (only if dead channel is present)
    if has_death:
        dead_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
        if dead_labels_outpath.exists():
            print("Loading existing user labelled dead data")
            dead_user_labels = np.asarray(load_zarr(dead_labels_outpath))
        else:
            dead_user_labels = np.zeros(all_images.shape[1:]).astype(np.int16)
        user_layers["User Provided Labels (Dead)"] = dead_user_labels
        
        # Predicted death labels
        pred_death_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_Death_PredictedLabels.zarr')
        if pred_death_labels_outpath.exists():
            pixelclass_layers["Pixel Classification (Dead)"] = np.asarray(load_zarr(pred_death_labels_outpath))
        else:
            pixelclass_layers["Pixel Classification (Dead)"] = np.zeros(all_images.shape[1:]).astype(np.int16)
    
    # Dynamic cell type layers
    for cell_type in all_cell_types:
        # User labels
        user_label_path = Path(pixel_class_outdir, f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr')
        if user_label_path.exists():
            print(f"Loading existing user labelled {cell_type} data")
            user_labels = np.asarray(load_zarr(user_label_path))
        else:
            user_labels = np.zeros(all_images.shape[1:]).astype(np.int16)
        user_layers[f"User Provided Labels ({cell_type.capitalize()})"] = user_labels
        
        # Predicted labels
        pred_label_path = Path(pixel_class_outdir, f'PixelClassifier_{cell_type.capitalize()}_PredictedLabels.zarr')
        if pred_label_path.exists():
            pixelclass_layers[f"Pixel Classification ({cell_type.capitalize()})"] = np.asarray(load_zarr(pred_label_path))
        else:
            pixelclass_layers[f"Pixel Classification ({cell_type.capitalize()})"] = np.zeros(all_images.shape[1:]).astype(np.int16)
        
        # Segment layers (for live visualization)
        segment_layers[f"{cell_type.capitalize()} Segments"] = np.zeros(all_images.shape[1:]).astype(np.int16)
    
    # Add all user label layers to viewer
    for name, data in user_layers.items():
        layer = viewer.add_labels(data, name=name, opacity=0.8)
        layer.blending = 'additive'
    
    # Add all pixel classification layers to viewer
    for name, data in pixelclass_layers.items():
        layer = viewer.add_labels(data, name=name, opacity=0.8, visible=False)
        layer.blending = 'additive'
   
    # Add all segment layers to viewer
    for name, data in segment_layers.items():
        layer = viewer.add_labels(data, name=name, opacity=0.8, visible=False)
        layer.blending = 'additive'

    log_output = QPlainTextEdit()
    log_output.setReadOnly(True)
    log_widget = QWidget()
    layout = QVBoxLayout()
    layout.addWidget(log_output)
    log_widget.setLayout(layout)
    viewer.window.add_dock_widget(log_widget, area="right", name="Log Output")
    
    # Create base GUI with only segment checkbox
    
    update_function = partial(
        segment_and_update, 
        pixel_class_outdir=pixel_class_outdir,
        n_workers=n_workers,
        log=log_output.appendPlainText
    )
    gui = magicgui(update_function, 
                only_segment={"widget_type": "Checkbox", "text": "Only Segment"},
                auto_call=False
    )
    
    # Create dynamic spinboxes for all cell types (organoids, immune, other)
    # Store spinboxes in dictionaries organized by parameter type
    edt_spinboxes = {}
    segment_size_spinboxes = {}
    opening_spinboxes = {}
    fill_holes_checkboxes = {}
    
    # Helper function to add spinboxes for a cell type
    def add_cell_type_spinboxes(cell_type, category):
        """Add all parameter spinboxes for a given cell type"""
        # Default values based on category
        if category == 'organoid':
            default_edt = 12.0
            default_segment_size = 1000
            default_opening = 3
            default_fill_holes = True
        elif category == 'immune':
            default_edt = 2.5
            default_segment_size = 10
            default_opening = 0
            default_fill_holes = True
        else:  # other
            default_edt = 1.0
            default_segment_size = 10
            default_opening = 0
            default_fill_holes = True
        
        # Override defaults with initial_params if provided
        if initial_params is not None:
            default_edt = initial_params.get(f"{cell_type}_edt_threshold", default_edt)
            default_segment_size = initial_params.get(f"{cell_type}_segment_size_min", default_segment_size)
            default_opening = initial_params.get(f"{cell_type}_opening_nr_pixels", default_opening)
            default_fill_holes = initial_params.get(f"{cell_type}_fill_holes", default_fill_holes)
        
        # Read spinbox limits from initial_params (fall back to defaults per category)
        edt_min = float(initial_params.get(f"{cell_type}_edt_min", 0)) if initial_params else 0
        edt_max = float(initial_params.get(f"{cell_type}_edt_max", 50.0)) if initial_params else 50.0
        edt_step = float(initial_params.get(f"{cell_type}_edt_step", 0.5)) if initial_params else 0.5
        size_min = int(initial_params.get(f"{cell_type}_segment_size_min_limit", 0)) if initial_params else 0
        size_max = int(initial_params.get(f"{cell_type}_segment_size_max_limit", 100000)) if initial_params else 100000
        size_step = int(initial_params.get(f"{cell_type}_segment_size_step", 10)) if initial_params else 10
        opening_min = int(initial_params.get(f"{cell_type}_opening_min", 0)) if initial_params else 0
        opening_max = int(initial_params.get(f"{cell_type}_opening_max", 10)) if initial_params else 10
        opening_step = int(initial_params.get(f"{cell_type}_opening_step", 1)) if initial_params else 1
        
        # EDT threshold spinbox
        label = QLabel(f"{cell_type} EDT threshold")
        gui.native.layout().addWidget(label)
        edt_spinbox = QDoubleSpinBox()
        edt_spinbox.setMinimum(edt_min)
        edt_spinbox.setMaximum(edt_max)
        edt_spinbox.setSingleStep(edt_step)
        edt_spinbox.setValue(default_edt)
        edt_spinbox.setObjectName(f"{cell_type}_edt_threshold")
        edt_spinboxes[cell_type] = edt_spinbox
        gui.native.layout().addWidget(edt_spinbox)
        
        # Segment size min spinbox
        label = QLabel(f"{cell_type} min segment size")
        gui.native.layout().addWidget(label)
        segment_size_spinbox = QSpinBox()
        segment_size_spinbox.setMinimum(size_min)
        segment_size_spinbox.setMaximum(size_max)
        segment_size_spinbox.setSingleStep(size_step)
        segment_size_spinbox.setValue(int(default_segment_size))
        segment_size_spinbox.setObjectName(f"{cell_type}_segment_size_min")
        segment_size_spinboxes[cell_type] = segment_size_spinbox
        gui.native.layout().addWidget(segment_size_spinbox)
        
        # Opening pixels spinbox
        label = QLabel(f"{cell_type} opening pixels")
        gui.native.layout().addWidget(label)
        opening_spinbox = QSpinBox()
        opening_spinbox.setMinimum(opening_min)
        opening_spinbox.setMaximum(opening_max)
        opening_spinbox.setSingleStep(opening_step)
        opening_spinbox.setValue(int(default_opening))
        opening_spinbox.setObjectName(f"{cell_type}_opening_nr_pixels")
        opening_spinboxes[cell_type] = opening_spinbox
        gui.native.layout().addWidget(opening_spinbox)
        
        # Fill holes checkbox
        fill_holes_cb = QCheckBox(f"{cell_type} fill holes")
        fill_holes_cb.setChecked(default_fill_holes)
        fill_holes_cb.setObjectName(f"{cell_type}_fill_holes")
        fill_holes_checkboxes[cell_type] = fill_holes_cb
        gui.native.layout().addWidget(fill_holes_cb)
    
    # Add spinboxes for all cell types
    for organoid_type in organoid_types:
        add_cell_type_spinboxes(organoid_type, 'organoid')
    
    for immune_type in immune_types:
        add_cell_type_spinboxes(immune_type, 'immune')
    
    for other_type in other_types:
        add_cell_type_spinboxes(other_type, 'other')
    
    # Override the call behavior to include all parameters
    def custom_call(*args, **kwargs):
        """Custom call function that collects all parameters from spinboxes"""
        print("Custom call function triggered")  # Debug
        
        # Collect all parameters
        all_params = {}
        
        # EDT thresholds
        for cell_type, spinbox in edt_spinboxes.items():
            all_params[f"{cell_type}_edt_threshold"] = spinbox.value()
        
        # Segment size mins
        for cell_type, spinbox in segment_size_spinboxes.items():
            all_params[f"{cell_type}_segment_size_min"] = int(spinbox.value())
        
        # Opening pixels
        for cell_type, spinbox in opening_spinboxes.items():
            all_params[f"{cell_type}_opening_nr_pixels"] = int(spinbox.value())
        
        # Fill holes
        for cell_type, checkbox in fill_holes_checkboxes.items():
            all_params[f"{cell_type}_fill_holes"] = checkbox.isChecked()
        
        print(f"Collected parameters: {all_params}")  # Debug
        
        # Save parameters via callback if provided
        if on_params_changed is not None:
            try:
                on_params_changed(all_params)
                print("Parameters saved via callback")
            except Exception as e:
                print(f"Warning: Could not save parameters via callback: {e}")
        
        # Call with all parameters as kwargs, with error handling
        try:
            return segment_and_update(
                pixel_class_outdir=pixel_class_outdir,
                only_segment=gui.only_segment.value,
                n_workers=n_workers,
                log=log_output.appendPlainText,
                **all_params
            )
        except Exception as e:
            import traceback
            error_msg = f"\nERROR: {str(e)}\n{traceback.format_exc()}"
            print(error_msg)
            log_output.appendPlainText(error_msg)
            raise
    
    # Replace default call button behavior with custom function
    # Disconnect existing connections to avoid double-calling
    try:
        gui.call_button.clicked.disconnect()
    except (TypeError, RuntimeError):
        pass  # No connections exist yet
    
    # Connect ONLY the button click to our custom function (not gui.called to avoid double execution)
    gui.call_button.clicked.connect(lambda: custom_call())
    
    viewer.window.add_dock_widget(gui)
        
    save_button = PushButton(label="Save User Labels")
    save_function = partial(
        save_pixel_classification,
        log=log_output.appendPlainText
    )
    save_button.clicked.connect(save_function)
    gui.native.layout().addWidget(save_button.native)

    # Add "Clear Labels" buttons per user label layer
    gui.native.layout().addWidget(QLabel("--- Clear Labels ---"))
    
    if has_death:
        clear_dead_btn = QtPushButton("Clear Labels: Dead")
        def _make_clear_fn(layer_name):
            def _clear():
                lyr = viewer.layers[layer_name]
                lyr.data = np.zeros_like(lyr.data)
            return _clear
        clear_dead_btn.clicked.connect(_make_clear_fn('User Provided Labels (Dead)'))
        gui.native.layout().addWidget(clear_dead_btn)
    
    for cell_type in all_cell_types:
        layer_name = f'User Provided Labels ({cell_type.capitalize()})'
        clear_btn = QtPushButton(f"Clear Labels: {cell_type.capitalize()}")
        clear_btn.clicked.connect(_make_clear_fn(layer_name))
        gui.native.layout().addWidget(clear_btn)

    #napari.run()
    return viewer

'''def _process_single_timepoint(args):
    """
    Worker function to process a single timepoint in parallel.
    
    Parameters
    ----------
    args : tuple
        (t, t_img, classifiers, clf_death, all_cell_types, organoid_types, immune_types, 
         other_types, all_edt_thresholds, all_segment_size_mins, all_opening_nr_pixels,
         all_fill_holes, only_segment, loaded_masks, has_death, sample_extract_features)
    
    Returns
    -------
    dict
        Contains masks and segments for all cell types at this timepoint
    """
    (t, t_img, classifiers, clf_death, all_cell_types, organoid_types, immune_types, 
     other_types, all_edt_thresholds, all_segment_size_mins, all_opening_nr_pixels,
     all_fill_holes, only_segment, loaded_masks, has_death, sample_extract_features) = args
    
    result = {
        'timepoint': t,
        'masks': {},
        'segments': {},
        'death_mask': None
    }
    
    # Extract features once for all classifiers
    if not only_segment:
        features = sample_extract_features(t_img)
    
    # Apply classifiers and get predictions for each cell type
    if not only_segment:
        # Apply pixel classifiers for each cell type
        for cell_type in all_cell_types:
            # Predict pixels using classifier
            pred_mask = future.predict_segmenter(features, classifiers[cell_type])
            # Convert to binary mask: label 1 = background -> 0, label 2 = foreground -> 1
            pred_mask[pred_mask > 0] -= 1  # Convert 1->0 (bg), 2->1 (fg)
            
            # Apply postprocessing
            opening_nr_pixels = all_opening_nr_pixels.get(cell_type, 0)
            fill_holes = all_fill_holes.get(cell_type, True)
            
            # Extract foreground mask (now label 1 after subtraction) for postprocessing
            fg_mask = (pred_mask == 1).astype(bool)
            processed_fg = postprocess_mask(fg_mask, fill_holes=fill_holes, opening_nr_pixels=int(opening_nr_pixels))
            # Reconstruct binary mask: 0=background, 1=foreground
            pred_mask = processed_fg.astype(pred_mask.dtype)
            
            result['masks'][cell_type] = pred_mask
        
        # Apply death classifier (only if present)
        if has_death and clf_death is not None:
            death_mask = future.predict_segmenter(features, clf_death)
            death_mask[death_mask > 0] -= 1
            result['death_mask'] = death_mask
    
    # Segment each cell type
    for cell_type in all_cell_types:
        edt_threshold = all_edt_thresholds.get(cell_type, 1.0)
        segment_size_min = all_segment_size_mins.get(cell_type, 10)
        
        if only_segment and cell_type in loaded_masks:
            mask = np.asarray(loaded_masks[cell_type][t])
        elif cell_type in result['masks']:
            mask = result['masks'][cell_type]
        else:
            continue
        
        # Segment the mask
        seg = segment_mask(
            mask=mask,
            edt_thr=edt_threshold,
            edt_thr_refined=None,
            segment_size_min=segment_size_min,
            use_dims=2
        )
        
        result['segments'][cell_type] = seg
    
    return result'''

def run_pixel_classifier_segmentation(
    output_dir,
    metadata,
    organoid_edt_thresholds=None,  # Dict: {organoid_type: threshold}
    immune_edt_thresholds=None,    # Dict: {immune_type: threshold}
    other_edt_thresholds=None,      # Dict: {other_type: threshold}
    organoid_segment_size_mins=None,  # Dict: {organoid_type: min_size}
    immune_segment_size_mins=None,    # Dict: {immune_type: min_size}
    other_segment_size_mins=None,      # Dict: {other_type: min_size}
    organoid_opening_nr_pixels=None,   # Dict: {organoid_type: pixels}
    immune_opening_nr_pixels=None,     # Dict: {immune_type: pixels}
    other_opening_nr_pixels=None,      # Dict: {other_type: pixels}
    organoid_fill_holes=None,          # Dict: {organoid_type: bool}
    immune_fill_holes=None,            # Dict: {immune_type: bool}
    other_fill_holes=None,             # Dict: {other_type: bool}
    timepoint_range=None,
    clf_organoid_paths=None,        # Dict: {organoid_type: path}
    clf_immune_paths=None,          # Dict: {immune_type: path}
    clf_other_paths=None,           # Dict: {other_type: path}
    clf_death_path=None,
    only_segment=False,
    overwrite_existing=False,
    n_workers=4,
    log_callback=print,
    ):

    """
    Run pixel classifier segmentation with support for multiple cell types.
    
    Parameters:
    -----------
    organoid_edt_thresholds : dict
        EDT thresholds for organoid types
    immune_edt_thresholds : dict
        EDT thresholds for immune types
    other_edt_thresholds : dict
        EDT thresholds for other cell types
    clf_organoid_paths, clf_immune_paths, clf_other_paths : dict
        Paths to classifier files for each cell type
    n_workers : int
        Number of parallel workers for timepoint processing (default: 4)
        Will be capped at available CPU count to prevent overloading
    log_callback : callable
        Function to handle log messages (default: print)
    """
    
    # Use log_callback instead of print
    log = log_callback

    # Validate and cap n_workers to available CPU count
    max_workers = multiprocessing.cpu_count()
    if n_workers > max_workers:
        safe_workers = max_workers
        print(f"⚠️ Requested {n_workers} workers but only {max_workers} CPUs available. Using {safe_workers} workers.")
        n_workers = safe_workers
    elif n_workers < 1:
        log(f"⚠️ Invalid n_workers={n_workers}. Using 1 worker (sequential).")
        n_workers = 1
    
    # Detect cell types from metadata
    from behav3d.core.metadata import (
        detect_organoid_types_from_metadata,
        detect_immune_cell_types_from_metadata,
        detect_other_cell_types_from_metadata,
        has_dead_channel
    )
    
    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)
    has_death = has_dead_channel(metadata)
    
    all_cell_types = organoid_types + immune_types + other_types
    
    print(f"Detected cell types: organoids={organoid_types}, immune={immune_types}, other={other_types}")
    print(f"Dead channel: {has_death}")
    _log_mem("run_pixel_classifier_segmentation START")
    
    # Initialize parameter dictionaries with defaults if not provided
    if organoid_edt_thresholds is None:
        organoid_edt_thresholds = {ct: 12.0 for ct in organoid_types}
    if immune_edt_thresholds is None:
        immune_edt_thresholds = {ct: 2.5 for ct in immune_types}
    if other_edt_thresholds is None:
        other_edt_thresholds = {ct: 1.0 for ct in other_types}
    
    if organoid_segment_size_mins is None:
        organoid_segment_size_mins = {ct: 1000 for ct in organoid_types}
    if immune_segment_size_mins is None:
        immune_segment_size_mins = {ct: 10 for ct in immune_types}
    if other_segment_size_mins is None:
        other_segment_size_mins = {ct: 10 for ct in other_types}
    
    if organoid_opening_nr_pixels is None:
        organoid_opening_nr_pixels = {ct: 3 for ct in organoid_types}
    if immune_opening_nr_pixels is None:
        immune_opening_nr_pixels = {ct: 0 for ct in immune_types}
    if other_opening_nr_pixels is None:
        other_opening_nr_pixels = {ct: 0 for ct in other_types}
    
    if organoid_fill_holes is None:
        organoid_fill_holes = {ct: True for ct in organoid_types}
    if immune_fill_holes is None:
        immune_fill_holes = {ct: True for ct in immune_types}
    if other_fill_holes is None:
        other_fill_holes = {ct: True for ct in other_types}
    
    # Combine all parameters
    all_edt_thresholds = {**organoid_edt_thresholds, **immune_edt_thresholds, **other_edt_thresholds}
    all_segment_size_mins = {**organoid_segment_size_mins, **immune_segment_size_mins, **other_segment_size_mins}
    all_opening_nr_pixels = {**organoid_opening_nr_pixels, **immune_opening_nr_pixels, **other_opening_nr_pixels}
    all_fill_holes = {**organoid_fill_holes, **immune_fill_holes, **other_fill_holes}
    
    # Setup classifier directory
    pixelclass_dir = Path(output_dir, "images", "PixelClassification")
    pixelclass_dir.mkdir(exist_ok=True, parents=True)
    
    # Load classifiers for each cell type
    classifiers = {}
    
    # Helper function to load/copy classifier
    def _load_classifier(cell_type, provided_path=None):
        default_name = f'PixelClassifier_{cell_type.capitalize()}.joblib'
        default_path = pixelclass_dir / default_name
        
        if provided_path:
            provided_path = Path(provided_path)
            if not provided_path.samefile(default_path):
                shutil.copy(provided_path, default_path)
            clf_path = default_path
        else:
            clf_path = default_path
        
        if clf_path.exists():
            return joblib.load(clf_path)
        else:
            raise FileNotFoundError(f"Classifier not found for {cell_type}: {clf_path}")
    
    # Load organoid classifiers
    organoid_clf_paths = clf_organoid_paths or {}
    for org_type in organoid_types:
        classifiers[org_type] = _load_classifier(org_type, organoid_clf_paths.get(org_type))
    
    # Load immune classifiers
    immune_clf_paths = clf_immune_paths or {}
    for immune_type in immune_types:
        classifiers[immune_type] = _load_classifier(immune_type, immune_clf_paths.get(immune_type))
    
    # Load other classifiers
    other_clf_paths = clf_other_paths or {}
    for other_type in other_types:
        classifiers[other_type] = _load_classifier(other_type, other_clf_paths.get(other_type))
    
    # MEJORA 1: clf.n_jobs=1 — prevent RF from spawning all cores per prediction
    # The RF saved from napari has n_jobs=-1 (all system cores).
    # Without this, each predict call launches N threads.
    # With n_workers=4 workers × N threads = thread explosion → crash.
    # With n_jobs=1: each worker uses exactly 1 core for RF. You control the total.
    for clf in classifiers.values():
        clf.n_jobs = 1
    print("✓ Classifiers set to n_jobs=1 for stable parallel prediction")
    
    # Load death classifier (only if death channel is present)
    clf_death = None
    if has_death:
        death_clf_default = pixelclass_dir / 'PixelClassifier_Death.joblib'
        if clf_death_path and not Path(clf_death_path).samefile(death_clf_default):
            shutil.copy(clf_death_path, death_clf_default)
        if death_clf_default.exists():
            clf_death = joblib.load(death_clf_default)
            clf_death.n_jobs = 1  # Same fix for death classifier
            print(f"Loaded death classifier from: {death_clf_default}")
        else:
            log(f"⚠️ Warning: Death channel detected but classifier not found at: {death_clf_default}")
    else:
        log("Skipping death classifier (no dead channel)")
    
    _log_mem("after loading classifiers")
    
    # Ensure all segment path columns exist with correct dtype (object/string)
    for cell_type in all_cell_types:
        if cell_type in organoid_types:
            prefix = 'or'
        elif cell_type in immune_types:
            prefix = 'im'
        elif cell_type in other_types:
            prefix = 'ot'
        else:
            continue
        
        path_col = f'{prefix}_{cell_type}_segments_image_path'
        if path_col not in metadata.columns:
            metadata[path_col] = pd.NA
        # Ensure column has object dtype to prevent FutureWarning
        metadata[path_col] = metadata[path_col].astype('object')
    
    # Process each sample
    for idx, sample in metadata.iterrows():
        print(f"Processing sample: {sample['sample_name']}")
        _log_mem(f"sample {sample['sample_name']} START")
        start_time = time.time()
        sample_name = sample['sample_name']
        
        # Setup paths
        raw_image_path = Path(sample['raw_image_path'])
        raw_image_zarr = Path(output_dir, "images", sample_name, f"{sample_name}.zarr")
        img_outdir = Path(output_dir, "images", sample_name)
        img_outdir.mkdir(parents=True, exist_ok=True)
            
        # Create output paths for all cell types
        segments_outpaths = {}
        mask_outpaths = {}
        
        for cell_type in all_cell_types:
            segments_outpaths[cell_type] = Path(img_outdir, f"{sample_name}_{cell_type}_segments.zarr")
            mask_outpaths[cell_type] = Path(img_outdir, f"{sample_name}_{cell_type}_mask.zarr")
        
        death_mask_outpath = Path(img_outdir, f"{sample_name}_mask_dead.zarr")

        # Load or convert raw image
        if not raw_image_zarr.exists():
            img = load_image(raw_image_path)
            save_as_zarr(img, raw_image_zarr)
        img = load_image(raw_image_zarr)
        print(f"  Image shape: {img.shape}")
        _log_mem(f"after loading image zarr handle")
        
        # Check if already segmented
        all_segments_exist = all(p.exists() for p in segments_outpaths.values())
        
        # -----------------------------------------------------------------
        # Resume / skip / overwrite decision
        # -----------------------------------------------------------------
        processing_marker = img_outdir / ".seg_processing"
        resuming = False     # will be True when we continue an interrupted run

        if overwrite_existing:
            # Overwrite: start fresh
            for seg_path in segments_outpaths.values():
                if seg_path.exists():
                    shutil.rmtree(seg_path)
            for mask_path in mask_outpaths.values():
                if mask_path.exists():
                    shutil.rmtree(mask_path)
            if death_mask_outpath.exists():
                shutil.rmtree(death_mask_outpath)
            # Remove any stale markers
            for f in img_outdir.glob(".seg_done_*"):
                f.unlink()
            if processing_marker.exists():
                processing_marker.unlink()

        elif processing_marker.exists():
            # Interrupted run detected --> RESUME
            # zarrs exist but may be incomplete,open them in r+ mode
            # only reprocess timepoints whose .seg_done marker is missing.
            resuming = True
            log(f"  ⚡ Interrupted run detected for {sample_name} — resuming")

        elif not only_segment:
            # Check if sample is already fully done
            skip_sample = False
            for cell_type in all_cell_types:
                prefix = 'or' if cell_type in organoid_types else ('im' if cell_type in immune_types else 'ot')
                path_col = f'{prefix}_{cell_type}_segments_image_path'
                existing_path = sample.get(path_col)
                if pd.notna(existing_path) and str(existing_path).strip():
                    if Path(str(existing_path)).exists():
                        skip_sample = True
                        break
            if not skip_sample and all(p.exists() for p in segments_outpaths.values()):
                skip_sample = True

            if skip_sample:
                log(f"  Sample {sample_name} already segmented. Skipping (Overwrite existing unchecked).")
                continue

        # If we are NOT resuming and NOT overwriting, clean for masks too
        if not resuming and not overwrite_existing:
            for seg_path in segments_outpaths.values():
                if seg_path.exists():
                    shutil.rmtree(seg_path)
            if not only_segment:
                if death_mask_outpath.exists():
                    shutil.rmtree(death_mask_outpath)
                for mask_path in mask_outpaths.values():
                    if mask_path.exists():
                        shutil.rmtree(mask_path)

        # Load existing masks if only_segment mode
        loaded_masks = {}
        if only_segment:
            # If a previous classification left the processing flag, masks may be incomplete
            if processing_marker.exists() and not resuming:
                raise RuntimeError(
                    f"only_segment mode: {sample_name} has a .seg_processing marker, "
                    f"meaning a previous full classification was interrupted. "
                    f"Masks may be incomplete — re-run full classification first."
                )
            missing_masks = []
            for cell_type, mask_path in mask_outpaths.items():
                if mask_path.exists():
                    loaded_masks[cell_type] = load_image(mask_path)
                else:
                    missing_masks.append(cell_type)
            if missing_masks:
                raise FileNotFoundError(
                    f"only_segment mode: mask(s) not found for {missing_masks} "
                    f"in {img_outdir}. Run full classification first."
                )
            # Validate mask timepoint dimensions match the raw image
            for cell_type, mask_arr in loaded_masks.items():
                if mask_arr.shape[0] < img.shape[0]:
                    raise ValueError(
                        f"only_segment mode: {cell_type} mask has {mask_arr.shape[0]} "
                        f"timepoints but raw image has {img.shape[0]}. "
                        f"Mask appears truncated — re-run full classification for this sample."
                    )

        # Determine which timepoints to process
        if timepoint_range is not None:
            start_t, end_t = timepoint_range
            start_t = max(0, min(start_t, img.shape[0]-1))
            end_t = max(start_t, min(end_t, img.shape[0]-1))
            timepoint_indices = list(range(start_t, end_t + 1))
            log(f"  Processing timepoints: {start_t} to {end_t} (total: {len(timepoint_indices)})")
        else:
            timepoint_indices = list(range(img.shape[0]))
            log(f"  Processing all timepoints: 0 to {img.shape[0]-1}")
        
        n_total = len(timepoint_indices)
        # Mask/segment spatial shape: skip T and C dimensions
        # predict_segmenter output is (Z, Y, X)
        spatial_shape = tuple(img.shape[2:])  # (Z, Y, X), NOT (C, Z, Y, X)
        
        # Build per-sample feature extractor from pixel size metadata
        sample_pixel_xy = float(sample.get('pixel_distance_xy', 1.0))
        sample_pixel_z  = float(sample.get('pixel_distance_z', 1.0))
        sample_extract_features = make_features(sample_pixel_xy, sample_pixel_z, verbose=True)
        
        # MEJORA 2: Dynamic RAM check — reduce workers if images are too large
        # Compute real feature count from a tiny dummy array (microseconds, ~4 KB).
        # This is just a measurement — no side effects on real data.
        dummy = np.zeros((img.shape[1],) + (8, 8, 8), dtype=np.float32)
        n_features = sample_extract_features(dummy).shape[-1]
        del dummy
        
        # Accurate memory estimate per worker:
        # - features array: spatial_voxels × n_features × 8 bytes (float64)
        # - raw timepoint held briefly: timepoint_bytes
        # - intermediates (masks, segments): ~50% overhead
        spatial_voxels = int(np.prod(spatial_shape))
        features_bytes = spatial_voxels * n_features * 8  # float64
        timepoint_bytes = int(np.prod(img.shape[1:])) * 2  # int16 = 2 bytes
        estimated_memory_per_worker = int(features_bytes + timepoint_bytes)
        print(f"  Feature probe: {n_features} features → {features_bytes/1e9:.2f} GB features array, "
              f"{estimated_memory_per_worker/1e9:.2f} GB/worker estimated")
        
        if _HAS_PSUTIL:
            available_ram = psutil.virtual_memory().available
            safe_workers = max(1, int(available_ram * 0.8 // estimated_memory_per_worker))
            actual_workers = min(n_workers, safe_workers)
            
            if actual_workers < n_workers:
                print(
                    f"⚠️ Large images ({timepoint_bytes/1e9:.2f} GB/timepoint, "
                    f"~{estimated_memory_per_worker/1e9:.2f} GB/worker estimated). "
                    f"Available RAM: {available_ram/1e9:.1f} GB. "
                    f"Reducing workers {n_workers} → {actual_workers} to avoid crash."
                )
            else:
                print(
                    f"  Memory OK: {available_ram/1e9:.1f} GB available, "
                    f"~{estimated_memory_per_worker/1e9:.2f} GB/worker → using {actual_workers} workers"
                )
        else:
            actual_workers = n_workers
            print(f"  psutil not installed, skipping RAM check. Using {actual_workers} workers.")
        
        _log_mem(f"before processing loop ({n_total} timepoints, {actual_workers} workers)")
        
        # ---------------------------------------------------------
        # Create the sample-level processing marker
        # ---------------------------------------------------------
        processing_marker.touch()

        # ---------------------------------------------------------
        # Determine which timepoints are already done (resume mode)
        # ---------------------------------------------------------
        done_indices = set()
        if resuming:
            for p in img_outdir.glob(".seg_done_*"):
                try:
                    done_indices.add(int(p.name.split("_")[-1]))
                except ValueError:
                    pass

        # Filter out completed timepoints
        if done_indices:
            args_list_filtered = [(t_i, t) for t_i, t in enumerate(timepoint_indices) if t_i not in done_indices]
            log(f"  Resume: {len(done_indices)}/{n_total} timepoints already done, {len(args_list_filtered)} remaining")
        else:
            args_list_filtered = None  # will build later

        # MEJORA 3: Pre-allocated zarrs with fixed size — thread-safe index writes
        # Instead of append_to_zarr (which requires strict sequential order),
        # create zarrs with final size and each worker writes zarr[t_i] = data.
        # Two workers never touch the same index.
        # Result is always zarr[0]=T0, zarr[1]=T1... regardless of completion order.
        zarr_mode = "r+" if resuming else "w"
        zarr_kw = dict(shape=(n_total,) + spatial_shape, chunks=(1,) + spatial_shape)

        if not only_segment:
            zarr_masks = {}
            zarr_segs = {}
            for cell_type in all_cell_types:
                zarr_masks[cell_type] = zarr.open(
                    str(mask_outpaths[cell_type]),
                    mode=zarr_mode, dtype="int16", **zarr_kw,
                )
                zarr_segs[cell_type] = zarr.open(
                    str(segments_outpaths[cell_type]),
                    mode=zarr_mode, dtype="uint16", **zarr_kw,
                )
            zarr_death = None
            if has_death and clf_death is not None:
                zarr_death = zarr.open(
                    str(death_mask_outpath),
                    mode=zarr_mode, dtype="int16", **zarr_kw,
                )
        else:
            zarr_segs = {}
            for cell_type in all_cell_types:
                zarr_segs[cell_type] = zarr.open(
                    str(segments_outpaths[cell_type]),
                    mode=zarr_mode, dtype="uint16", **zarr_kw,
                )
            zarr_masks = None
            zarr_death = None
        
        # -----------------------------------------------------------------
        # Per-timepoint worker function
        # -----------------------------------------------------------------
        def process_timepoint(args, _ff=sample_extract_features):
            t_i, t = args
            t_start = time.time()
            
            # Read 1 timepoint from zarr (lazy read, materialized here)
            t_img = np.asarray(img[t])
            
            if not only_segment:
                # MEJORA 4: threadpoolctl limits feature extraction to 1 internal thread
                # multiscale_basic_features (scikit-image) internally uses all cores.
                # Without control: 4 workers × 16 cores = 64 competing threads.
                # With limits=1: 4 workers × 1 thread = actual_workers cores total.
                if _HAS_THREADPOOLCTL:
                    with threadpoolctl.threadpool_limits(limits=1):
                        features = _ff(t_img)
                else:
                    features = _ff(t_img)
                del t_img  # free immediately — features is all we need from here
                
                for cell_type in all_cell_types:
                    # RF already has n_jobs=1 (Mejora 1) → 1 core per worker
                    pred_mask = future.predict_segmenter(features, classifiers[cell_type])
                    pred_mask[pred_mask > 0] -= 1
                    
                    opening_nr_pixels_val = all_opening_nr_pixels.get(cell_type, 0)
                    fill_holes_val = all_fill_holes.get(cell_type, True)
                    fg_mask = (pred_mask == 1).astype(bool)
                    processed_fg = postprocess_mask(fg_mask, fill_holes=fill_holes_val, opening_nr_pixels=int(opening_nr_pixels_val))
                    pred_mask = processed_fg.astype(np.int16)
                    
                    # Write mask by exact index (Mejora 3)
                    # t_i=0 is always T0, t_i=1 is always T1, etc.
                    zarr_masks[cell_type][t_i] = pred_mask
                    
                    edt_threshold = all_edt_thresholds.get(cell_type, 1.0)
                    segment_size_min = all_segment_size_mins.get(cell_type, 10)
                    seg = segment_mask(
                        mask=pred_mask.astype(bool),
                        edt_thr=edt_threshold,
                        edt_thr_refined=None,
                        segment_size_min=segment_size_min,
                        use_dims=2
                    )
                    
                    zarr_segs[cell_type][t_i] = seg
                    del pred_mask, fg_mask, processed_fg, seg
                
                # Death classifier
                if has_death and clf_death is not None and zarr_death is not None:
                    death_mask = future.predict_segmenter(features, clf_death)
                    death_mask[death_mask > 0] -= 1
                    zarr_death[t_i] = death_mask
                    del death_mask
                
                del features
            else:
                # only_segment mode: read existing masks, segment, write
                for cell_type in all_cell_types:
                    if cell_type not in loaded_masks:
                        continue
                    mask = np.asarray(loaded_masks[cell_type][t])
                    
                    edt_threshold = all_edt_thresholds.get(cell_type, 1.0)
                    segment_size_min = all_segment_size_mins.get(cell_type, 10)
                    seg = segment_mask(
                        mask=mask,
                        edt_thr=edt_threshold,
                        edt_thr_refined=None,
                        segment_size_min=segment_size_min,
                        use_dims=2
                    )
                    
                    zarr_segs[cell_type][t_i] = seg
                    del mask, seg
                del t_img  # free in only_segment branch too
            
            gc.collect()
            
            # Mark this timepoint as fully written
            (img_outdir / f".seg_done_{t_i:06d}").touch()
            
            t_elapsed = time.time() - t_start
            if (t_i + 1) % 10 == 0 or t_i == 0 or (t_i + 1) == n_total:
                _log_mem(f"T{t} ({t_i+1}/{n_total}, {t_elapsed:.1f}s)")
            
            return t_i, t
        
        # MEJORA 5: ThreadPoolExecutor at timepoint level
        # Instead of 1 sequential timepoint with N chaotic internal threads,
        # actual_workers timepoints in parallel with 1 thread each.
        # Speedup ≈ actual_workers×, stable on any machine.
        if args_list_filtered is None:
            args_list = list(enumerate(timepoint_indices))
        else:
            args_list = args_list_filtered

        n_to_process = len(args_list)
        print(f"  Starting parallel processing with {actual_workers} workers ({n_to_process} timepoints)...")
        
        with ThreadPoolExecutor(max_workers=actual_workers) as executor:
            futures_map = {
                executor.submit(process_timepoint, args): args
                for args in args_list
            }
            completed = 0
            for fut in tqdm(as_completed(futures_map), total=n_to_process, desc=f"  {sample_name}"):
                try:
                    t_i, t = fut.result()
                    completed += 1
                    if completed % 10 == 0 or completed == 1 or completed == n_to_process:
                        print(f"  Progress: {completed}/{n_to_process} timepoints completed")
                except Exception as e:
                    args = futures_map[fut]
                    print(f"  ❌ Error in timepoint {args}: {e}")
                    raise
        
        # ---------------------------------------------------------
        # All timepoints done — clean up markers
        # ---------------------------------------------------------
        for f in img_outdir.glob(".seg_done_*"):
            f.unlink()
        if processing_marker.exists():
            processing_marker.unlink()
        print(f"  ✓ All {n_to_process} timepoints written successfully (markers cleaned)")
                
        # Update metadata with segment paths (with prefixes)
        for cell_type in all_cell_types:
            # Determine prefix based on category
            if cell_type in organoid_types:
                prefix = 'or'
            elif cell_type in immune_types:
                prefix = 'im'
            elif cell_type in other_types:
                prefix = 'ot'
            else:
                continue  # Skip unknown types
            
            path_col = f'{prefix}_{cell_type}_segments_image_path'
            metadata.at[idx, path_col] = str(segments_outpaths[cell_type])
            log(f"  Saved {cell_type} segments: {segments_outpaths[cell_type]}")
        
        # Update dead_mask_path in metadata if death channel is present AND column exists
        if has_death and 'dead_mask_path' in metadata.columns and death_mask_outpath.exists():
            metadata['dead_mask_path'] = metadata['dead_mask_path'].astype('object')
            metadata.at[idx, 'dead_mask_path'] = str(death_mask_outpath)
            log(f"  Updated dead_mask_path: {death_mask_outpath}")
        
        _log_mem(f"sample {sample_name} END")
        print(f"  Sample {sample_name} completed in {time.time() - start_time:.1f}s")
    
    _log_mem("run_pixel_classifier_segmentation END")
    return metadata