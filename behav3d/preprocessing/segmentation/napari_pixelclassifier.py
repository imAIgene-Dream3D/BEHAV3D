from pathlib import Path
import numpy as np
import time
import shutil

import napari
from magicgui import magicgui
from magicgui.widgets import PushButton

from aicspylibczi import CziFile

from skimage import data, segmentation, feature, future
from skimage.measure import label
from skimage.segmentation import watershed, relabel_sequential

from sklearn.ensemble import RandomForestClassifier
from scipy.ndimage import binary_fill_holes

from behav3d.utils.preprocessing import open_mask, dilate_mask
from behav3d.utils.segmentation import segment_size_filter, get_border_segments, remove_boundary_segments, calculate_edt
from behav3d.utils.fileio import save_as_zarr, load_zarr, load_image, append_to_zarr

from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

import joblib
from functools import partial

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
        
        mask = binary_fill_holes(mask)
    if opening_nr_pixels>0 and opening_nr_pixels is not None:
        mask = open_mask(mask)
    return(mask)

def segment_mask(mask, segment_splitting_edt, segment_size_min, use_dims):
    offset=1
    edt = calculate_edt(mask, use_dims=use_dims)
    seeds = label(edt > segment_splitting_edt)
    segments = watershed(mask, markers = seeds, mask=mask)
    seeds2 = label(mask * (segments==0))
    seeds2[seeds2!=0] += seeds.max()
    # Relabel last segments to keep unique labels
    segments[segments==0]=seeds2[segments==0]
    segments = segment_size_filter(segments, size_min=segment_size_min)
    segments, _, _ = relabel_sequential(segments, offset)
    return(segments)

def train_pixel_classifier(
    output_dir,
    metadata,
    examples_per_sample = 4,
    sample_specific_classifier=False,
    n_workers=16
    ):
    
    pixel_class_outdir = Path(output_dir, "images", "PixelClassification")
    pixel_class_outdir.mkdir(exist_ok=True)
    
    all_images = []
    all_features = []
    
    features_outpath = Path(pixel_class_outdir, 'PixelClassifier_Features.zarr')
    image_outpath = Path(pixel_class_outdir, 'PixelClassifier_Images.zarr')
    
    ### TEST
    if not features_outpath.exists() or not image_outpath.exists():
        if image_outpath.exists():
            shutil.rmtree(image_outpath)
        if features_outpath.exists():
            shutil.rmtree(features_outpath)
            
        for idx, sample in metadata.iterrows():
            
            sample_name = sample['sample_name']
            
            tcell_ch=sample['tcell_channel']
            live_ch=sample['live_channel']
            dead_ch=sample['dead_channel']
            
            print(f"Calculating features for: {sample_name}")
            
            raw_image_path = Path(sample['raw_image_path'])
            raw_image_zarr =  Path(output_dir, "images", sample_name, f"{sample_name}.zarr")
            
            images = load_image(raw_image_zarr)
            max_t = images.shape[0]-1
            
            idc = np.linspace(0, max_t, examples_per_sample, dtype=int)
            print(f"Taking timepoints: {idc}")
            
            sample_images = [images[t, [tcell_ch, live_ch, dead_ch]] for t in idc]
            # for idx in idc:
            #     sample_images.append(images[idx])      
            
            all_images+=sample_images
            
            for img in tqdm(sample_images):
                append_to_zarr(np.expand_dims(features_func(img), axis=0), features_outpath)
            # all_features+=[features_func(img) for img in tqdm(sample_images)]
            
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
        all_images = da.stack(all_images)
        all_images = all_images.transpose(1, 0, 2, 3, 4)
        save_as_zarr(all_images, image_outpath)
        del all_images
        gc.collect()
        all_images = load_zarr(image_outpath)
        
        # all_features = da.stack(all_features)
        # save_as_zarr(all_features, features_outpath)
        # del all_features
        # gc.collect()
        all_features = load_zarr(features_outpath)
    else:
        all_images = load_image(image_outpath)
        all_features = load_image(features_outpath)     

    def segment_and_update(
        pixel_class_outdir,
        tcell_edt_threshold: float = 2.5,
        ):
        start_time = time.time()
        print("###### Running Segmentation")
        # Access the label layer and feature image
        
        image_layer = viewer.layers['Image']
        image_data = image_layer.data
        
        label_layer = viewer.layers['User Provided Labels']
        label_data = label_layer.data
        
        print("Training Random Forest Classifier")
        clf = RandomForestClassifier(
            n_estimators=100,
            n_jobs=-1, 
            max_depth=20, 
            max_samples=0.05,
            class_weight="balanced"
            )
        
        clf = future.fit_segmenter(label_data, all_features, clf)
        
        print("Predicting Background, T-cells and Organoid pixels")
        pixel_class = future.predict_segmenter(all_features, clf)

        print("Segment Cell and Organoid Instances")
        full_seg_organoid=[]
        full_seg_tcell=[]
        for idx, t_pixel_class in enumerate(pixel_class):
            mask_organoid = t_pixel_class==2
            mask_tcell = t_pixel_class==3
            
            mask_tcell = postprocess_mask(mask_tcell, opening_nr_pixels=1)
            mask_organoid = postprocess_mask(mask_organoid, opening_nr_pixels=3)

            seg_organoid = segment_mask(
                mask=mask_organoid,
                segment_size_min=500,
                segment_splitting_edt=10,
                use_dims=3
            )
            
            seg_tcell = segment_mask(
                mask=mask_tcell,
                segment_size_min=10,
                segment_splitting_edt=tcell_edt_threshold,
                use_dims=3
            )

            full_seg_organoid.append(seg_organoid)
            full_seg_tcell.append(seg_tcell)
        full_seg_organoid = np.stack(full_seg_organoid, axis=0)
        full_seg_tcell = np.stack(full_seg_tcell, axis=0)
        
        print("Updating Napari Data")
        viewer.layers["Pixel Classification"].data = pixel_class
        viewer.layers["Organoid Segments"].data = full_seg_organoid
        viewer.layers["Tcell Segments"].data = full_seg_tcell
        
        print("Saving RandomForest, Sparse labels and input images to {pixel_class_outdir}")
        ### SAVE 
        random_forest_outpath = Path(pixel_class_outdir, 'PixelClassifier_RandomForest.joblib')
        joblib.dump(clf, random_forest_outpath)
        
        labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserLabels.zarr')
        save_as_zarr(label_data, labels_outpath)
        
        image_outpath = Path(pixel_class_outdir, 'PixelClassifier_Images.zarr')
        save_as_zarr(image_data, image_outpath)
        
        print(f"###### DONE time elapsed: {time.time() - start_time:.2f} s")
        
    # Create Napari viewer
    viewer = napari.Viewer()
    img_layer = viewer.add_image(all_images, name="Image", contrast_limits=(0,float(da.percentile(all_images[1,0].reshape(-1), 99).compute()[0])))
    img_layer.contrast_limits_range = (0, int(all_images.max().compute()))

    labels_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserLabels.zarr')
    
    if labels_outpath.exists():
        print("Loading existing user labelled data")
        user_labels = np.asarray(load_zarr(labels_outpath))
    else:
        user_labels = np.zeros(all_images.shape[1:]).astype(np.int16)
        
    viewer.add_labels(user_labels, name="User Provided Labels", opacity=0.3)
    viewer.add_labels(np.zeros(all_images.shape[1:]).astype(np.int16), name="Pixel Classification", opacity=0.7, visible=False)
    viewer.add_labels(np.zeros(all_images.shape[1:]).astype(np.int16), name="Organoid Segments", opacity=0.7, visible=False)
    viewer.add_labels(np.zeros(all_images.shape[1:]).astype(np.int16), name="Tcell Segments", opacity=0.7, visible=False)

    update_function = partial(segment_and_update, pixel_class_outdir=pixel_class_outdir)
    gui = magicgui(update_function, 
                tcell_edt_threshold={"widget_type": "FloatSlider", "min": 1.0, "max": 15.0, "step": 0.5}
                )
    viewer.window.add_dock_widget(gui)
    napari.run()
    
def run_pixel_classifier(
    classifier_path,
    metadata,
    ):
    
    clf = joblib.load(classifier_path)
    
    for idx, sample in metadata.iterrows():
        print(f"Processing sample: {sample['sample_name']}")
        start_time = time.time()
        sample_name = sample['sample_name']
        raw_image_path = Path(sample['raw_image_path'])
        
        img_outdir = Path(output_dir, "images", sample_name)
        if not img_outdir.exists():
            img_outdir.mkdir(parents=True)

        img = load_image(raw_image_path)
        img = np.transpose(img, [1,0,2,3,4])
        
        metadata.at[idx, "tcell_segments_image_path"] = str(segmenter.tcell_segments_outpath)
        metadata.at[idx, "organoid_tracks_image_path"] = str(segmenter.organoid_segments_outpath)
    return(metadata)