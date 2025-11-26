import os

from pathlib import Path
import numpy as np
import pandas as pd
import time
import shutil

import napari
from magicgui import magicgui
from magicgui.widgets import PushButton, FloatSlider
from qtpy.QtWidgets import QPlainTextEdit, QWidget, QVBoxLayout, QApplication, QLabel

from aicspylibczi import CziFile

from skimage import data, segmentation, feature, future
from skimage.measure import label
from skimage.segmentation import watershed, relabel_sequential

from sklearn.ensemble import RandomForestClassifier
# from napari_apoc import PixelClassifier  # Commented out - use scikit-learn instead
from scipy import ndimage
from scipy.ndimage import binary_fill_holes, find_objects

from behav3d.utils.preprocessing import open_mask, dilate_mask, zeropad_image_to_match_shape
from behav3d.utils.segmentation import segment_size_filter, get_border_segments, remove_boundary_segments, calculate_edt, segment_2d_filter
from behav3d.utils.fileio import save_as_zarr, load_zarr, load_image, append_to_zarr, get_image_shape

import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from tqdm import tqdm

import joblib
from functools import partial

import zarr
import dask.array as da
import gc

## TODO create a function for BEHAV3D notebook

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

sigma_min = 1
sigma_max = 16
features_func = partial(
        feature.multiscale_basic_features,
        intensity=True,
        edges=True,
        texture=False,
        # sigma_min=sigma_min,
        sigma_max=sigma_max,
        channel_axis=0,
    )

def postprocess_mask(mask, fill_holes=True, opening_nr_pixels=1):
    if fill_holes:  
        # mask = binary_fill_holes(mask)
        # Assume `mask` is a 3D binary array: (Z, Y, X)
        filled_mask = np.zeros_like(mask)
        for i in range(mask.shape[0]):
            filled_mask[i] = binary_fill_holes(mask[i])
        mask = filled_mask
    if opening_nr_pixels>0 and opening_nr_pixels is not None:
        mask = open_mask(mask)
    return(mask)


def refine_segment(args):
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
    # print("###", label_id, "refine_segment time elapsed: ", time.time() - start_time)
    new_seg = watershed(-local_edt, markers=new_seg, mask=local_mask)
    # print("###", label_id, "refine_segment time elapsed: ", time.time() - start_time)
    new_seg, _, _ = relabel_sequential(new_seg)
    # print("###", label_id, "refine_segment time elapsed: ", time.time() - start_time)
    return (new_seg, tuple(minc))

def segment_mask(mask, edt_thr=1.5, edt_thr_refined=[2, 2.5, 3], segment_size_min=15, use_dims=3, n_workers=1):
    offset = 1
    start_time = time.time()
    # Step 1: Initial segmentation
    edt = calculate_edt(mask, use_dims=use_dims)
    seeds = label(edt >= edt_thr)
    segments = watershed(-edt, markers = seeds, mask=mask)
    seeds2 = label(mask * (segments==0))
    seeds2[seeds2!=0] += seeds.max()
    # Relabel last segments to keep unique labels
    segments[segments==0]=seeds2[segments==0]
    segments = segment_size_filter(segments, size_min=segment_size_min)
    segments = watershed(-edt, markers = segments, mask=mask)
    # segments, _, _ = relabel_sequential(segments, offset)
    
    # If edt_thr_refined ois not list, turn into list
    if not isinstance(edt_thr_refined, list):
        if edt_thr_refined is not None:
            edt_thr_refined = [edt_thr_refined]
        
    if edt_thr_refined is not None:
        for thr in edt_thr_refined:
            # Step 2: Prepare for per-label refinement
            start_time = time.time()
            unique_labels = np.unique(segments)
            unique_labels = unique_labels[unique_labels != 0]  # Exclude background
            args_list = []
    
            slices = find_objects(segments)
            for i, slc in enumerate(slices):
                label_id = i + 1  # Because label 0 is background and ignored

                if slc is None:
                    continue  # This label is not present

                local_mask = segments[slc] == label_id
                local_edt = edt[slc]
                minc = [s.start for s in slc]

                args_list.append((local_mask, local_edt, thr, segment_size_min, minc))

            if n_workers is None:
                n_workers = multiprocessing.cpu_count()

            results = []
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                results += list(executor.map(refine_segment, args_list))
           
            current_label = 0
            segments = np.zeros_like(mask, dtype=np.uint16)
            for result in results:
                local_segment, minc = result
                z0, y0, x0 = minc
                z1, y1, x1 = z0 + local_segment.shape[0], y0 + local_segment.shape[1], x0 + local_segment.shape[2]

                # Assign to global array with label offset
                nonzero_mask = local_segment > 0
                local_segment[nonzero_mask] += current_label
                segments[z0:z1, y0:y1, x0:x1][nonzero_mask] = local_segment[nonzero_mask]

                current_label = segments.max() + 1  # Update for next segment
    else:
        segments, _, _ = relabel_sequential(segments, offset)
        
    # FIlter out 2D segments
    segments = segment_2d_filter(segments)
    return(segments)


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
    args_list = [(classifier, str(features_outpath), str(pred_labels_outpath), idx) for idx in range(pred_labels.shape[0])]
    test=[]
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        test+=list(tqdm(executor.map(_apply_classifier, args_list), total=len(args_list)))
    
    pixel_class = load_image(pred_labels_outpath)
    pixel_class[pixel_class==1] = 0
    # pixel_class[pixel_class>0] -= 1
    # pixel_class = pixel_class.compute()
    pixel_class = np.asarray(pixel_class)

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
    from behav3d.utils import (
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

            print(f"Calculating features...")
            for img in tqdm(all_images):
                append_to_zarr(
                    np.expand_dims(
                        zeropad_image_to_match_shape(
                            img = features_func(img),
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
            feature_img = features_func(img)
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
        **edt_thresholds  # Dynamic EDT thresholds for each cell type
        ):
        """
        Dynamic segmentation function that handles different cell types.
        EDT thresholds are passed as kwargs
        """
        print("✅ segment_and_update() called!")  # Debug
        print(f"Received EDT thresholds: {edt_thresholds}")  # Debug
        start_time = time.time()
        log("###### Running Segmentation\n")
        log(f"EDT Thresholds: {edt_thresholds}\n")
        QApplication.processEvents()
        
        # Access the label layer and feature image
        image_data = all_images  # Use the original all_images data
        
        # === Save user labels for all cell types ===
        cell_type_labels = {}  # Store {cell_type: label_data}
        
        # Save death labels (only if present)
        if has_death:
            dead_label_layer = viewer.layers['User Provided Labels (Dead)']
            dead_label_data = dead_label_layer.data
            dead_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
            save_as_zarr(dead_label_data, dead_labels_outpath)
            log("Saved death labels")
        
        # Save dynamic cell type labels
        for cell_type in all_cell_types:
            layer_name = f'User Provided Labels ({cell_type.capitalize()})'
            label_layer = viewer.layers[layer_name]
            label_data = label_layer.data
            labels_outpath = Path(pixel_class_outdir, f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr')
            save_as_zarr(label_data, labels_outpath)
            cell_type_labels[cell_type] = label_data
            log(f"Saved {cell_type} labels")


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
        
        # === Train classifiers for all cell types ===
        classifiers = {}  # Store {cell_type: classifier}
        clf_death = None  # Initialize to None
        
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
        
        
        # === Apply classifiers to predict pixels ===
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
                QApplication.processEvents()
                pred_death_mask = apply_classifier(clf_death, features_outpath, pred_death_labels_outpath, n_workers=n_workers)
                viewer.layers["Pixel Classification (Dead)"].data = pred_death_mask
            else:
                log("\n### Skipping death prediction (no dead channel)")
                QApplication.processEvents()

            # Predict pixels for all cell types
            for cell_type in all_cell_types:
                log(f"\n### Predicting {cell_type.capitalize()} Pixels")
                QApplication.processEvents()
                pred_mask = apply_classifier(classifiers[cell_type], features_outpath, pred_labels_paths[cell_type], n_workers=n_workers)
                pred_masks[cell_type] = pred_mask
                viewer.layers[f"Pixel Classification ({cell_type.capitalize()})"].data = pred_mask
        
        else:
            # Load existing predictions
            if has_death:
                log("\n### Loading Death Prediction Mask")
                QApplication.processEvents()
                pred_death_mask = viewer.layers["Pixel Classification (Dead)"].data
                pred_death_mask[pred_death_mask==1] = 0
                viewer.layers["Pixel Classification (Dead)"].data = pred_death_mask
            
            # Load predictions for all cell types
            for cell_type in all_cell_types:
                log(f"\n### Loading {cell_type.capitalize()} Prediction Mask")
                QApplication.processEvents()
                pred_mask = viewer.layers[f"Pixel Classification ({cell_type.capitalize()})"].data
                pred_mask[pred_mask==1] = 0
                viewer.layers[f"Pixel Classification ({cell_type.capitalize()})"].data = pred_mask
                pred_masks[cell_type] = pred_mask
            
        # === Segment instances for all cell types ===
        log("\n### Segment Cell Instances")
        QApplication.processEvents()
        
        # Segment each cell type using EDT watershed
        segmented_cells = {}  # Store {cell_type: segmented_mask}
        
        for cell_type in all_cell_types:
            edt_threshold = edt_thresholds.get(cell_type, 1.0)  # Default if not specified
            log(f"\n### Segmenting {cell_type.capitalize()} (EDT threshold={edt_threshold})")
            QApplication.processEvents()
            
            pred_mask = pred_masks[cell_type]
            
            # Segment each timepoint
            segmented_timepoints = []
            for t_idx in range(pred_mask.shape[0]):
                mask_t = pred_mask[t_idx]
                # Remove background (label 1), keep only foreground (label 2)
                mask_t = (mask_t == 2).astype(np.uint8)
                
                # Apply watershed segmentation with EDT
                if mask_t.sum() > 0:
                    # Distance transform
                    distance = ndimage.distance_transform_edt(mask_t)
                    # Find peaks above threshold
                    local_max = distance > edt_threshold
                    markers = label(local_max)
                    # Watershed
                    segmented = watershed(-distance, markers, mask=mask_t)
                else:
                    segmented = np.zeros_like(mask_t, dtype=np.int32)
                
                segmented_timepoints.append(segmented)
            
            full_seg = np.stack(segmented_timepoints, axis=0)
            segmented_cells[cell_type] = full_seg
            viewer.layers[f"{cell_type.capitalize()} Segments"].data = full_seg
            
        # Save images
            image_outpath = Path(pixel_class_outdir, 'PixelClassifier_Images.zarr')
            save_as_zarr(image_data, image_outpath)
            
            log(f"\n###### DONE time elapsed: {time.time() - start_time:.2f} s")
    
    def save_pixel_classification(log=print):
        """Save all user-provided labels for all cell types"""
        # Save death labels (only if present)
        if has_death:
            dead_label_layer = viewer.layers['User Provided Labels (Dead)']
            dead_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
            save_as_zarr(dead_label_layer.data, dead_labels_outpath)
            log(f"Saved death labels to: {dead_labels_outpath}")
        
        # Save labels for all cell types
        for cell_type in all_cell_types:
            layer_name = f'User Provided Labels ({cell_type.capitalize()})'
            label_layer = viewer.layers[layer_name]
            labels_outpath = Path(pixel_class_outdir, f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr')
            save_as_zarr(label_layer.data, labels_outpath)
            log(f"Saved {cell_type} labels to: {labels_outpath}")
        
        log("✅ All user labels saved!")
            
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
        channel_name = f"Channel {ch+1}"
        if ch < len(channel_colors):
            channel_name = f"Channel {ch+1} ({channel_colors[ch]})"
        
        img_layer = viewer.add_image(
            channel_data, 
            name=channel_name, 
            contrast_limits=contrast_limits,
            colormap=channel_colors[ch] if ch < len(channel_colors) else 'gray',
            blending='additive'  # This allows channels to blend together
        )
        img_layer.contrast_limits_range = (0, float(channel_data.max()))

    # ==== DYNAMIC LAYER CREATION FOR ALL CELL TYPES ====
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
        layer = viewer.add_labels(data, name=name, opacity=0.3)
    
    # Add all pixel classification layers to viewer
    for name, data in pixelclass_layers.items():
        layer = viewer.add_labels(data, name=name, opacity=0.3, visible=False)
   
    # Add all segment layers to viewer
    for name, data in segment_layers.items():
        viewer.add_labels(data, name=name, opacity=0.7, visible=False)

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
    
    # Create dynamic sliders for organoid types
    organoid_sliders = {}
    for organoid_type in organoid_types:
        # Create a label for the slider
        label = QLabel(f"{organoid_type} edt threshold")
        gui.native.layout().addWidget(label)
        
        # Create the slider (use organoid parameters)
        slider = FloatSlider(
            value=12.0,  # Default organoid threshold
            min=0.5,
            max=20.0,
            step=0.5,
            name=f"{organoid_type}_edt_threshold"
        )
        organoid_sliders[organoid_type] = slider
        # Add slider to gui layout
        gui.native.layout().addWidget(slider.native)
    
    # Create dynamic sliders for immune cell types
    immune_sliders = {}
    for immune_type in immune_types:
        # Create a label for the slider
        label = QLabel(f"{immune_type} edt threshold")
        gui.native.layout().addWidget(label)
        
        # Create the slider (use immune cell parameters)
        slider = FloatSlider(
            value=2.5,  # Default immune cell threshold
            min=0.5,
            max=15.0,
            step=0.5,
            name=f"{immune_type}_edt_threshold"
        )
        immune_sliders[immune_type] = slider
        # Add slider to gui layout
        gui.native.layout().addWidget(slider.native)
    
    # Create dynamic sliders for other cell types
    other_sliders = {}
    for other_type in other_types:
        # Create a label for the slider
        label = QLabel(f"{other_type} edt threshold")
        gui.native.layout().addWidget(label)
        
        # Create the slider (use default parameters for other types)
        slider = FloatSlider(
            value=1.0,  # Default "other" threshold
            min=0.5,
            max=20.0,
            step=0.5,
            name=f"{other_type}_edt_threshold"
        )
        other_sliders[other_type] = slider
        # Add slider to gui layout
        gui.native.layout().addWidget(slider.native)
    
    # Override the call behavior to include all thresholds
    def custom_call(*args, **kwargs):
        """Custom call function that collects EDT thresholds from all sliders"""
        print("Custom call function triggered")  # Debug
        
        # Get all EDT thresholds from sliders
        all_edt_thresholds = {}
        all_edt_thresholds.update({org_type: slider.value for org_type, slider in organoid_sliders.items()})
        all_edt_thresholds.update({immune_type: slider.value for immune_type, slider in immune_sliders.items()})
        all_edt_thresholds.update({other_type: slider.value for other_type, slider in other_sliders.items()})
        
        print(f"Collected EDT thresholds: {all_edt_thresholds}")  # Debug
        
        # Call with all thresholds as kwargs
        return segment_and_update(
            pixel_class_outdir=pixel_class_outdir,
            only_segment=gui.only_segment.value,
            n_workers=n_workers,
            log=log_output.appendPlainText,
            **all_edt_thresholds
        )
    
    # Replace default call button behavior with custom function
    # First, disconnect all existing connections
    try:
        gui.called.disconnect()
    except (TypeError, RuntimeError):
        # No connections exist yet
        pass
    
    # Connect our custom call function
    gui.called.connect(custom_call)
    
    # Also ensure the button directly triggers our function
    try:
        gui.call_button.clicked.disconnect()
    except (TypeError, RuntimeError):
        pass
    gui.call_button.clicked.connect(lambda: custom_call())
    
    viewer.window.add_dock_widget(gui)
        
    save_button = PushButton(label="Save User Labels")
    save_function = partial(
        save_pixel_classification,
        log=log_output.appendPlainText
    )
    save_button.clicked.connect(save_function)
    gui.native.layout().addWidget(save_button.native)

    napari.run()

# def apply_pixel_classifier(
#     classifier,
#     img,
#     t
#     ):
#     clf_path, path, outpath, idx = args
#     clf = joblib.load(clf_path)
#     features = np.asarray(load_image(path, mode="r")[idx])
#     prediction = zarr.open(outpath, mode="r+")
#     prediction[idx] = future.predict_segmenter(features, clf)
'''def _run_single_timepoint_segmentation(
    t_img,
    clf_org,
    clf_org_2,
    clf_tcell,
    clf_death,
    organoid_edt_threshold=6,
    organoid_min_size=1000,
    only_segment=False,
    tcell_mask=None,
    organoid_mask=None,
    organoid_mask_2=None,
    ):
    if only_segment:
        assert tcell_mask is not None, "tcell_segments must be provided when only_segment is True"
        assert organoid_mask is not None, "organoid_segments must be provided when only_segment is True"
        if clf_org_2 is not None:
            assert organoid_mask_2 is not None, "organoid_segments_2 must be provided when only_segment is True and clf_org_2 is not None"
        death_mask=None
    else:
        features = features_func(t_img)
        # result = future.predict_segmenter(features, clf)
            
        # print("\n### Predicting Organoid Pixels")
        organoid_mask = future.predict_segmenter(features, clf_org)
        organoid_mask[organoid_mask>0] -= 1

        if clf_org_2 is not None:
            organoid_mask_2 = future.predict_segmenter(features, clf_org_2)
            organoid_mask_2[organoid_mask_2>0] -= 1
        
        # print("\n### Predicting T-cell Pixels")
        tcell_mask = future.predict_segmenter(features, clf_tcell)
        tcell_mask[tcell_mask>0] -= 1
        
        # print("\n### Predicting Death Pixels")
        death_mask = future.predict_segmenter(features, clf_death)
        death_mask[death_mask>0] -= 1
    
    
    # print("\n### Segmenting Organoids and T-cells")
    if clf_org_2 is not None:
        two_org_types = True
        seg_organoid, seg_organoid_2, seg_tcell = segment_tcell_and_organoid(
            args = (organoid_mask, organoid_mask_2, tcell_mask, organoid_edt_threshold, two_org_types)
        )
        # return(seg_organoid, seg_organoid_2, seg_tcell, pred_death_mask)
    else:
        two_org_types = False
        seg_organoid, seg_tcell = segment_tcell_and_organoid(
            args = (organoid_mask, tcell_mask, organoid_edt_threshold, two_org_types),  
        )
        
        # return(seg_organoid, seg_tcell, pred_death_mask, pred_org_mask, pred_tcell_mask)
    return_dict = {
        "segments_organoid": seg_organoid,
        "segments_tcell": seg_tcell,
        "death_mask": death_mask,
        "mask_tcell": tcell_mask,
        "mask_organoid": organoid_mask
    }
    
    if two_org_types:
        return_dict["segments_organoid_2"] = seg_organoid_2
    
    return return_dict'''

def _process_single_timepoint(args):
    """
    Worker function to process a single timepoint in parallel.
    
    Parameters
    ----------
    args : tuple
        (t, t_img, classifiers, clf_death, all_cell_types, organoid_types, immune_types, 
         other_types, all_edt_thresholds, only_segment, loaded_masks, has_death)
    
    Returns
    -------
    dict
        Contains masks and segments for all cell types at this timepoint
    """
    (t, t_img, classifiers, clf_death, all_cell_types, organoid_types, immune_types, 
     other_types, all_edt_thresholds, only_segment, loaded_masks, has_death) = args
    
    result = {
        'timepoint': t,
        'masks': {},
        'segments': {},
        'death_mask': None
    }
    
    # Extract features once for all classifiers
    if not only_segment:
        features = features_func(t_img)
    
    # Apply classifiers and get predictions for each cell type
    if not only_segment:
        # Apply pixel classifiers for each cell type
        for cell_type in all_cell_types:
            # Predict pixels using classifier
            pred_mask = future.predict_segmenter(features, classifiers[cell_type])
            # Convert to binary mask (label 2 = foreground, label 1 = background)
            pred_mask[pred_mask > 0] -= 1  # Convert 1->0 (bg), 2->1 (fg)
            result['masks'][cell_type] = pred_mask
        
        # Apply death classifier (only if present)
        if has_death and clf_death is not None:
            death_mask = future.predict_segmenter(features, clf_death)
            death_mask[death_mask > 0] -= 1
            result['death_mask'] = death_mask
    
    # Segment each cell type
    for cell_type in all_cell_types:
        edt_threshold = all_edt_thresholds.get(cell_type, 1.0)
        
        if only_segment and cell_type in loaded_masks:
            mask = np.asarray(loaded_masks[cell_type][t])
        elif cell_type in result['masks']:
            mask = result['masks'][cell_type]
        else:
            continue
        
        # Determine segment size based on cell type category
        if cell_type in organoid_types:
            segment_size_min = 1000
        elif cell_type in immune_types:
            segment_size_min = 10
        else:
            segment_size_min = 10
        
        # Segment the mask
        seg = segment_mask(
            mask=mask,
            edt_thr=edt_threshold,
            edt_thr_refined=[edt_threshold+0.5, edt_threshold+1.0, edt_threshold+1.5] if (cell_type in immune_types or cell_type in other_types) else None,
            segment_size_min=segment_size_min,
            use_dims=3
        )
        
        result['segments'][cell_type] = seg
    
    return result

def run_pixel_classifier_segmentation(
    output_dir,
    metadata,
    organoid_edt_thresholds=None,  # Dict: {organoid_type: threshold}
    immune_edt_thresholds=None,    # Dict: {immune_type: threshold}
    other_edt_thresholds=None,      # Dict: {other_type: threshold}
    timepoint_range=None,
    clf_organoid_paths=None,        # Dict: {organoid_type: path}
    clf_immune_paths=None,          # Dict: {immune_type: path}
    clf_other_paths=None,           # Dict: {other_type: path}
    clf_death_path=None,
    only_segment=False,
    overwrite_existing=False,
    n_workers=4,
    # deprecated
    #corganoid_edt_threshold=None,
    #clf_org_path=None,
    #clf_org2_path=None,
    #clf_tcell_path=None,
    #two_org_types=False,
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
    """
    
    # Validate and cap n_workers to available CPU count
    max_workers = multiprocessing.cpu_count()
    if n_workers > max_workers:
        safe_workers = max(1, max_workers // 2)
        print(f"⚠️ Requested {n_workers} workers but only {max_workers} CPUs available. Using {safe_workers} workers (half of available) for system stability.")
        n_workers = safe_workers
    elif n_workers < 1:
        print(f"⚠️ Invalid n_workers={n_workers}. Using 1 worker (sequential).")
        n_workers = 1
    
    # Detect cell types from metadata
    from behav3d.utils import (
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
    
    # Initialize threshold dictionaries with defaults if not provided
    if organoid_edt_thresholds is None:
        organoid_edt_thresholds = {ct: 12.0 for ct in organoid_types}
    if immune_edt_thresholds is None:
        immune_edt_thresholds = {ct: 2.5 for ct in immune_types}
    if other_edt_thresholds is None:
        other_edt_thresholds = {ct: 1.0 for ct in other_types}
    
    # Combine all thresholds
    all_edt_thresholds = {**organoid_edt_thresholds, **immune_edt_thresholds, **other_edt_thresholds}
    
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
    
    # Load death classifier (only if death channel is present)
    clf_death = None
    if has_death:
        death_clf_default = pixelclass_dir / 'PixelClassifier_Death.joblib'
        if clf_death_path and not Path(clf_death_path).samefile(death_clf_default):
            shutil.copy(clf_death_path, death_clf_default)
        if death_clf_default.exists():
            clf_death = joblib.load(death_clf_default)
            print(f"Loaded death classifier from: {death_clf_default}")
        else:
            print(f"⚠️ Warning: Death channel detected but classifier not found at: {death_clf_default}")
    else:
        print("Skipping death classifier (no dead channel)")
    
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
        
        # Check if already segmented
        all_segments_exist = all(p.exists() for p in segments_outpaths.values())
        
        if all_segments_exist and not overwrite_existing and not only_segment:
            print("  Already segmented, skipping")
        else:   
            # Remove existing outputs if overwriting
            for cell_type, seg_path in segments_outpaths.items():
                if seg_path.exists():
                    shutil.rmtree(seg_path)
            
            # Load existing masks if only_segment mode
            loaded_masks = {}
            if only_segment:
                for cell_type, mask_path in mask_outpaths.items():
                    if mask_path.exists():
                        loaded_masks[cell_type] = load_image(mask_path)
                    else:
                        print(f"  Warning: Mask not found for {cell_type}: {mask_path}")
            else:
                # Remove mask files if not in only_segment mode
                if death_mask_outpath.exists():
                    shutil.rmtree(death_mask_outpath)
                for mask_path in mask_outpaths.values():
                    if mask_path.exists():
                        shutil.rmtree(mask_path)
                                             
            # Determine which timepoints to process
            if timepoint_range is not None:
                start_t, end_t = timepoint_range
                start_t = max(0, min(start_t, img.shape[0]-1))
                end_t = max(start_t, min(end_t, img.shape[0]-1))
                timepoint_indices = list(range(start_t, end_t + 1))
                print(f"  Processing timepoints: {start_t} to {end_t} (total: {len(timepoint_indices)})")
            else:
                timepoint_indices = list(range(img.shape[0]))
                print(f"  Processing all timepoints: 0 to {img.shape[0]-1}")
            
            # Prepare arguments for parallel processing
            args_list = [
                (t, img[t], classifiers, clf_death, all_cell_types, organoid_types, 
                 immune_types, other_types, all_edt_thresholds, only_segment, 
                 loaded_masks, has_death)
                for t in timepoint_indices
            ]
            
            # Process timepoints in parallel
            print(f"  Using {n_workers} parallel workers")
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                results = list(tqdm(
                    executor.map(_process_single_timepoint, args_list),
                    total=len(args_list),
                    desc=f"  {sample_name}"
                ))
            
            # Save results in timepoint order
            for result in sorted(results, key=lambda x: x['timepoint']):
                # Save masks
                if not only_segment:
                    for cell_type in all_cell_types:
                        if cell_type in result['masks']:
                            append_to_zarr(
                                np.expand_dims(result['masks'][cell_type], axis=0),
                                mask_outpaths[cell_type]
                            )
                    
                    # Save death mask
                    if result['death_mask'] is not None:
                        append_to_zarr(
                            np.expand_dims(result['death_mask'], axis=0),
                            death_mask_outpath
                        )
                
                # Save segments
                for cell_type in all_cell_types:
                    if cell_type in result['segments']:
                        append_to_zarr(
                            np.expand_dims(result['segments'][cell_type], axis=0),
                            segments_outpaths[cell_type]
                        )
                
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
        
        print(f"  Sample {sample_name} completed in {time.time() - start_time:.1f}s")
    
    return metadata