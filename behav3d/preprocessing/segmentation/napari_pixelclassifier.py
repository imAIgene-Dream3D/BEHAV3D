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

import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from tqdm import tqdm

import joblib
from functools import partial

import zarr
import dask.array as da
import gc

## TODO create a function for BEHAV3D notebook


sigma_max = 16
features_func = partial(
        feature.multiscale_basic_features,
        intensity=True,
        edges=True,
        texture=True,
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

def segment_tcell_and_organoid(
    # mask_organoid,
    # mask_tcell,
    args,
    tcell_edt_threshold=1.5,
    tcell_edt_threshold_refined=[2, 2.5, 3],
    tcell_segment_size_min=10,
    # organoid_edt_threshold=12,
    organoid_segment_size_min=1000,
    ):
    """
    Segment T-cells and organoids from the prediction mask.
    Parameters
    ----------
    prediction_mask : np.ndarray
        The prediction mask from the pixel classifier.
    tcell_val : int, optional
        The value of the T-cell pixels in the prediction mask. The default is 3.
    organoid_val : int, optional
        The value of the organoid pixels in the prediction mask. The default is 2.
    tcell_edt_threshold : float, optional
        The threshold for the EDT to split T-cell segments. The default is 2.5.
    organoid_edt_threshold : float, optional
        The threshold for the EDT to split organoid segments. The default is 10.
    tcell_segment_size_min : int, optional 
        The minimum size of T-cell segments. The default is 10.
    organoid_segment_size_min : int, optional
        The minimum size of organoid segments. The default is 500.
    use_dims : int, optional
        The number of dimensions to use for the EDT calculation. The default is 3.
    Returns
    -------
    segments : np.ndarray
        The segmented T-cells and organoids.
    """
    organoid, tcell, organoid_edt_threshold =  args
    tcell = postprocess_mask(tcell, opening_nr_pixels=0)
    organoid = postprocess_mask(organoid, opening_nr_pixels=3)

    start_time = time.time()
    organoid = segment_mask(
        mask=organoid,
        segment_size_min=organoid_segment_size_min,
        edt_thr=organoid_edt_threshold,
        edt_thr_refined=None,
        use_dims=3
    )
    
    start_time = time.time()
    tcell = segment_mask(
        mask=tcell,
        segment_size_min=tcell_segment_size_min,
        edt_thr=tcell_edt_threshold,
        edt_thr_refined=tcell_edt_threshold_refined,
        use_dims=3
    )
    return(organoid, tcell)

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
    pixel_class = pixel_class.compute()
    return pixel_class
        
def train_pixel_classifier(
    output_dir,
    metadata,
    examples_per_sample = 3,
    sample_specific_classifier=False,
    n_workers=None,
    manual_dim_order=None  # Optional: tuple/list for custom dimension order in transpose
    ):
    """
    Train a pixel classifier for segmentation.
    Args:
        output_dir: Output directory
        metadata: Metadata DataFrame
        examples_per_sample: Number of timepoints per sample
        sample_specific_classifier: Whether to use sample-specific classifier
        n_workers: Number of workers
        manual_dim_order: Optional. Tuple/list specifying the order for transpose, e.g. (1,0,2,3,4) for (C,T,Z,Y,X)
    """
    if n_workers is None:
        n_workers = multiprocessing.cpu_count()
        
    pixel_class_outdir = Path(output_dir, "images", "PixelClassification")
    pixel_class_outdir.mkdir(exist_ok=True, parents=True)
    
    features_outpath = Path(pixel_class_outdir, 'PixelClassifier_Features.zarr')
    image_outpath = Path(pixel_class_outdir, 'PixelClassifier_Images.zarr')
 
    if not features_outpath.exists() or not image_outpath.exists():
        if image_outpath.exists():
            shutil.rmtree(image_outpath)
        if features_outpath.exists():
            shutil.rmtree(features_outpath)
        
        all_images = []
        all_features = []
        for idx, sample in metadata.iterrows():
            
            sample_name = sample['sample_name']
            
            tcell_ch=sample['tcell_channel']
            live_ch=sample['live_channel']
            dead_ch=sample['dead_channel']
            
            print(f"Calculating features for: {sample_name}")
            
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
            # for idx in idc:
            #     sample_images.append(images[idx])      
            
            all_images+=sample_images
            
            for img in tqdm(sample_images):
                append_to_zarr(np.expand_dims(features_func(img), axis=0), features_outpath)
            
            # def calculate_features(args):
            #     path, t = args
            #     img = load_image(path)
            #     img = img[t]
            #     return features_func(img)
            
            # args_list = [(raw_image_zarr, t) for t in idc]
            # if n_workers > 1:
            #     with ProcessPoolExecutor(max_workers=n_workers) as executor:
            #         all_features+=list(tqdm(executor.map(calculate_features, args_list), total=len(idc)))
            # else:
            #     all_features+=[calculate_features(args) for args in tqdm(args_list) ]
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
    
    def segment_and_update(
        pixel_class_outdir,
        only_segment=False,
        tcell_edt_threshold=2.5,
        organoid_edt_threshold: float = 12,
        n_workers: int = 16,
        log=print
        ):
        start_time = time.time()
        log("###### Running Segmentation\n")
        QApplication.processEvents()
        # Access the label layer and feature image
        # Get the first channel layer for reference (all channels have same shape)
        image_layer = viewer.layers[f'Channel 1 (red)'] if 'Channel 1 (red)' in viewer.layers else viewer.layers[0]
        image_data = all_images  # Use the original all_images data
        
        org_label_layer = viewer.layers['User Provided Labels (Organoid)']
        org_label_data = org_label_layer.data
        
        tcell_label_layer = viewer.layers['User Provided Labels (Tcell)']
        tcell_label_data = tcell_label_layer.data
        
        dead_label_layer = viewer.layers['User Provided Labels (Dead)']
        dead_label_data = dead_label_layer.data
        
        org_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserOrganoidLabels.zarr')
        tcell_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserTcellLabels.zarr')
        dead_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
        save_as_zarr(org_label_data, org_labels_outpath)
        save_as_zarr(tcell_label_data, tcell_labels_outpath)
        save_as_zarr(dead_label_data, dead_labels_outpath)
    
        def train_classifier(user_labels, features):
            flat_label_data = user_labels.ravel()
            flat_features = features.reshape(-1, features.shape[-1])  # shape: (N_total, 90)

            # Get 1D indices where labels exist
            label_indices = np.flatnonzero(flat_label_data > 0)

            selected_features = flat_features[label_indices].compute()  # (N_selected, 90)
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
        
        if not only_segment:
            log("\n### Training Random Forest Classifier (Organoids)")
            clf_organoids = train_classifier(org_label_data, all_features)
            org_random_forest_outpath = Path(pixel_class_outdir, 'PixelClassifier_Organoid.joblib')
            log("Saving RandomForest, Sparse labels and input images to {org_random_forest_outpath}")
            joblib.dump(clf_organoids, org_random_forest_outpath)
            QApplication.processEvents()
            
            log("\n### Training Random Forest Classifier (T-cells)")
            clf_tcells = train_classifier(tcell_label_data, all_features)
            tcell_random_forest_outpath = Path(pixel_class_outdir, 'PixelClassifier_Tcell.joblib')
            log("Saving RandomForest, Sparse labels and input images to {tcell_random_forest_outpath}")
            joblib.dump(clf_tcells, tcell_random_forest_outpath)
            QApplication.processEvents()
            
            log("\n### Training Random Forest Classifier (Cell Death)")
            clf_death = train_classifier(dead_label_data, all_features)
            death_random_forest_outpath = Path(pixel_class_outdir, 'PixelClassifier_Death.joblib')
            log("Saving RandomForest, Sparse labels and input images to {death_random_forest_outpath}")
            joblib.dump(clf_death, death_random_forest_outpath)
            QApplication.processEvents()
        
        
        pred_org_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_Organoid_PredictedLabels.zarr')
        pred_tcell_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_Tcell_PredictedLabels.zarr')
        pred_death_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_Death_PredictedLabels.zarr')
        # parts_outpath = Path(pred_labels_outpath.parent, "pred_parts")
        
        if not only_segment:
            if pred_org_labels_outpath.exists():
                shutil.rmtree(pred_org_labels_outpath)
            if pred_tcell_labels_outpath.exists():
                shutil.rmtree(pred_tcell_labels_outpath)
            if pred_death_labels_outpath.exists():
                shutil.rmtree(pred_death_labels_outpath)

            log("\n### Predicting Organoid Pixels")
            QApplication.processEvents()
            pred_org_mask = apply_classifier(clf_organoids, features_outpath, pred_org_labels_outpath)
            
            log("\n### Predicting T-cell Pixels")
            QApplication.processEvents()
            pred_tcell_mask = apply_classifier(clf_tcells, features_outpath, pred_tcell_labels_outpath)
            
            log("\n### Predicting Death Pixels")
            QApplication.processEvents()
            pred_death_mask = apply_classifier(clf_death, features_outpath, pred_death_labels_outpath)
            
            viewer.layers["Pixel Classification (Organoid)"].data = pred_org_mask
            viewer.layers["Pixel Classification (Tcell)"].data = pred_tcell_mask
            viewer.layers["Pixel Classification (Dead)"].data = pred_death_mask

        else:
            log("\n### Loading Organoid Prediction Mask")
            QApplication.processEvents()
            pred_org_mask = viewer.layers["Pixel Classification (Organoid)"].data
            pred_org_mask[pred_org_mask==1] = 0
            viewer.layers["Pixel Classification (Organoid)"].data = pred_org_mask
            
            log("\n### Loading T-cell Prediction Mask")
            QApplication.processEvents()
            pred_tcell_mask = viewer.layers["Pixel Classification (Tcell)"].data
            pred_tcell_mask[pred_tcell_mask==1] = 0
            viewer.layers["Pixel Classification (Tcell)"].data = pred_tcell_mask
            
            log("\n### Loading Death Prediction Mask")
            QApplication.processEvents()
            pred_death_mask = viewer.layers["Pixel Classification (Dead)"].data
            pred_death_mask[pred_death_mask==1] = 0
            viewer.layers["Pixel Classification (Dead)"].data = pred_death_mask
            
        log("\n### Segment Cell and Organoid Instances")
        QApplication.processEvents()
        args_list = [(pred_org_mask[idx], pred_tcell_mask[idx], organoid_edt_threshold) for idx in range(org_label_data.shape[0])]
        results=[]

        with ThreadPoolExecutor(max_workers=n_workers) as executor:
            results+=list(tqdm(executor.map(segment_tcell_and_organoid, args_list), total=len(args_list)))

        full_seg_organoid, full_seg_tcell = zip(*results)

        full_seg_organoid = np.stack(full_seg_organoid, axis=0)
        full_seg_tcell = np.stack(full_seg_tcell, axis=0)
        
        log("\n### Updating Napari Data")
        
        viewer.layers["Organoid Segments"].data = full_seg_organoid
        viewer.layers["Tcell Segments"].data = full_seg_tcell
        
        image_outpath = Path(pixel_class_outdir, 'PixelClassifier_Images.zarr')
        save_as_zarr(image_data, image_outpath)
        
        log(f"\n###### DONE time elapsed: {time.time() - start_time:.2f} s")
    
    def save_pixel_classification(log=print):
        label_layer = viewer.layers['User Provided Labels']
        label_data = label_layer.data

        labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserLabels.zarr')
        save_as_zarr(label_data, labels_outpath)
        log(f"Saved Pixel Classification to: {labels_outpath}")
            
    # Create Napari viewer
    viewer = napari.Viewer()
    
    # Split channels and add them as separate colored layers
    # all_images shape is (channels, time, z, y, x)
    n_channels = all_images.shape[0]
    
    # Define colors for different channels
    channel_colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow']
    
    for ch in range(n_channels):
        channel_data = all_images[ch]  # Shape: (time, z, y, x)
        
        # Calculate contrast limits for this channel
        channel_percentile = float(np.percentile(channel_data.reshape(-1), 99))
        contrast_limits = (0, channel_percentile)
        
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
        img_layer.contrast_limits_range = (0, channel_data.max())

    org_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserOrganoidLabels.zarr')
    tcell_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserTcellLabels.zarr')
    dead_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
    
    
    if org_labels_outpath.exists():
        print("Loading existing user labelled organoid data")
        org_user_labels = np.asarray(load_zarr(org_labels_outpath))
    else:
        org_user_labels = np.zeros(all_images.shape[1:]).astype(np.int16)
        
    if tcell_labels_outpath.exists():
        print("Loading existing user labelled Tcell data")
        tcell_user_labels = np.asarray(load_zarr(tcell_labels_outpath))
    else:
        tcell_user_labels = np.zeros(all_images.shape[1:]).astype(np.int16)
        
    if dead_labels_outpath.exists():
        print("Loading existing user labelled dead data")
        dead_user_labels = np.asarray(load_zarr(dead_labels_outpath))
    else:
        dead_user_labels = np.zeros(all_images.shape[1:]).astype(np.int16)
    
    user_layers = {
        "User Provided Labels (Organoid)": org_user_labels,
        "User Provided Labels (Tcell)": tcell_user_labels,
        "User Provided Labels (Dead)": dead_user_labels,
    }

    for name, data in user_layers.items():
        layer = viewer.add_labels(data, name=name, opacity=0.3)
    
    pred_org_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_Organoid_PredictedLabels.zarr')
    pred_tcell_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_Tcell_PredictedLabels.zarr')
    pred_death_labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_Death_PredictedLabels.zarr')
    if not pred_org_labels_outpath.exists():
        pixelclass_layers = {
            "Pixel Classification (Organoid)": np.zeros(all_images.shape[1:]).astype(np.int16),
            "Pixel Classification (Tcell)": np.zeros(all_images.shape[1:]).astype(np.int16),
            "Pixel Classification (Dead)": np.zeros(all_images.shape[1:]).astype(np.int16),
        }
    else:
        pixelclass_layers = {
            "Pixel Classification (Organoid)": np.asarray(load_zarr(pred_org_labels_outpath)),
            "Pixel Classification (Tcell)": np.asarray(load_zarr(pred_tcell_labels_outpath)),
            "Pixel Classification (Dead)": np.asarray(load_zarr(pred_death_labels_outpath)),
        }
    
    for name, data in pixelclass_layers.items():
        layer = viewer.add_labels(data, name=name, opacity=0.3, visible=False)
   
    viewer.add_labels(np.zeros(all_images.shape[1:]).astype(np.int16), name="Organoid Segments", opacity=0.7, visible=False)
    viewer.add_labels(np.zeros(all_images.shape[1:]).astype(np.int16), name="Tcell Segments", opacity=0.7, visible=False)

    log_output = QPlainTextEdit()
    log_output.setReadOnly(True)
    log_widget = QWidget()
    layout = QVBoxLayout()
    layout.addWidget(log_output)
    log_widget.setLayout(layout)
    viewer.window.add_dock_widget(log_widget, area="right", name="Log Output")
    
    update_function = partial(
        segment_and_update, 
        pixel_class_outdir=pixel_class_outdir,
        n_workers=n_workers,
        log=log_output.appendPlainText
        )
    gui = magicgui(update_function, 
                tcell_edt_threshold={"widget_type": "FloatSlider", "min": 1.0, "max": 15.0, "step": 0.5},
                organoid_edt_threshold={"widget_type": "FloatSlider", "min": 1.0, "max": 20.0, "step": 0.5},
                only_segment={"widget_type": "Checkbox", "text": "Only Segment"}
                )
    viewer.window.add_dock_widget(gui)
        
    save_button = PushButton(label="Save User Labels")
    save_function = partial(
        save_pixel_classification,
        log =log_output.appendPlainText
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
def _run_single_timepoint_segmentation(
    t_img,
    clf_org,
    clf_tcell,
    clf_death,
    organoid_edt_threshold=6,
    ):
    features = features_func(t_img)
    # result = future.predict_segmenter(features, clf)
        
    # print("\n### Predicting Organoid Pixels")
    pred_org_mask = future.predict_segmenter(features, clf_org)
    pred_org_mask[pred_org_mask>0] -= 1
    
    # print("\n### Predicting T-cell Pixels")
    pred_tcell_mask = future.predict_segmenter(features, clf_tcell)
    pred_tcell_mask[pred_tcell_mask>0] -= 1
    
    # print("\n### Predicting Death Pixels")
    pred_death_mask = future.predict_segmenter(features, clf_death)
    pred_death_mask[pred_death_mask>0] -= 1
    
    # print("\n### Segmenting Organoids and T-cells")
    seg_organoid, seg_tcell = segment_tcell_and_organoid(
        args = (pred_org_mask, pred_tcell_mask, organoid_edt_threshold),  
    )
    
    return(seg_organoid, seg_tcell, pred_death_mask)
    
def run_pixel_classifier_segmentation(
    output_dir,
    metadata,
    organoid_edt_threshold=12,
    timepoint_range=None
    ):
    
    clf_org_path = Path(output_dir, "images", "PixelClassification", 'PixelClassifier_Organoid.joblib')
    clf_tcell_path = Path(output_dir, "images", "PixelClassification", 'PixelClassifier_Tcell.joblib')
    clf_death_path = Path(output_dir, "images", "PixelClassification", 'PixelClassifier_Death.joblib')
  
    clf_org = joblib.load(clf_org_path)
    clf_tcell = joblib.load(clf_tcell_path)
    clf_death = joblib.load(clf_death_path)
    
    for idx, sample in metadata.iterrows():
        print(f"Processing sample: {sample['sample_name']}") 
        start_time = time.time()
        sample_name = sample['sample_name']
        
        tcell_ch=sample['tcell_channel']
        live_ch=sample['live_channel']
        dead_ch=sample['dead_channel']
        
        raw_image_path = Path(sample['raw_image_path'])
        raw_image_zarr =  Path(output_dir, "images", sample_name, f"{sample_name}.zarr")
        img_outdir = Path(output_dir, "images", sample_name)
        if not img_outdir.exists():
            img_outdir.mkdir(parents=True)
            
        tcell_segments_outpath = Path(img_outdir, f"{sample_name}_tcell_segments.zarr")
        organoid_segments_outpath = Path(img_outdir, f"{sample_name}_organoid_segments.zarr")
        death_mask_outpath = Path(img_outdir, f"{sample_name}_mask_dead.zarr")
        
        if not raw_image_zarr.exists():
            img = load_image(raw_image_path)
            # img = img[:, [tcell_ch, live_ch, dead_ch]]
            save_as_zarr(img, raw_image_zarr)
        img = load_image(raw_image_zarr)
        print(img.shape)
        
        if (tcell_segments_outpath.exists() and 
            organoid_segments_outpath.exists()):
            print("Already segmented, skipping")
        else:   
            if tcell_segments_outpath.exists():
                shutil.rmtree(tcell_segments_outpath)
            if organoid_segments_outpath.exists():
                shutil.rmtree(organoid_segments_outpath)
            
            # Determine which timepoints to process
            if timepoint_range is not None:
                start_t, end_t = timepoint_range
                # Ensure range is within bounds
                start_t = max(0, min(start_t, img.shape[0]-1))
                end_t = max(start_t, min(end_t, img.shape[0]-1))
                timepoint_indices = list(range(start_t, end_t + 1))
                print(f"Processing timepoints: {start_t} to {end_t} (total: {len(timepoint_indices)})")
            else:
                timepoint_indices = list(range(img.shape[0]))
                print(f"Processing all timepoints: 0 to {img.shape[0]-1}")
            
            for t in tqdm(timepoint_indices, total=len(timepoint_indices)):
                t_img = img[t]
                seg_organoid, seg_tcell, pred_death_mask = _run_single_timepoint_segmentation(
                    t_img=t_img,
                    clf_org=clf_org,
                    clf_tcell=clf_tcell,
                    clf_death=clf_death,
                    organoid_edt_threshold=organoid_edt_threshold
                )
                append_to_zarr(np.expand_dims(seg_organoid, axis=0), organoid_segments_outpath)
                append_to_zarr(np.expand_dims(seg_tcell, axis=0), tcell_segments_outpath)
                append_to_zarr(np.expand_dims(pred_death_mask, axis=0), death_mask_outpath)
                
        metadata.at[idx, "tcell_segments_image_path"] = str(tcell_segments_outpath)
        metadata.at[idx, "organoid_segments_image_path"] = str(organoid_segments_outpath)
    return(metadata)