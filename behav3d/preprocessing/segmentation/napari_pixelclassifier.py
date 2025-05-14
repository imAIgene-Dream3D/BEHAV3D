from aicspylibczi import CziFile
from pathlib import Path
import numpy as np
from magicgui import magicgui
from magicgui.widgets import PushButton
import napari
from skimage import data, segmentation, feature, future
from sklearn.ensemble import RandomForestClassifier
from functools import partial
from scipy.ndimage import binary_fill_holes
from behav3d.utils.preprocessing import open_mask, dilate_mask
from skimage.measure import label
from behav3d.utils.segmentation import segment_size_filter, get_border_segments, remove_boundary_segments, calculate_edt
from skimage.segmentation import watershed, relabel_sequential


## TODO create a function for BEHAV3D notebook


path = r"D:\BHVD_BEHAV3D\BEHAV3D_python\data\Jess_ROCHE\ROCHE_JM1_Exp011_Img04_169M_50KTcells_withTCB.czi"
path = Path(path)
czi = CziFile(path)

max_t = czi.get_dims_shape()[0]["T"][1] - 1

img_f, _ = czi.read_image(T=0)
img_f = np.squeeze(img_f)
img_m, _ = czi.read_image(T=int(max_t/2))
img_m = np.squeeze(img_m)
img_l, _ = czi.read_image(T=max_t)
img_l = np.squeeze(img_l)

img = np.stack([img_f, img_m, img_l])
feature_img = []
for idx, ch_img in enumerate(img):
    print(idx)
    sigma_min = 1
    sigma_max = 30
    features_func = partial(
            feature.multiscale_basic_features,
            intensity=True,
            edges=True,
            texture=True,
            sigma_min=sigma_min,
            # sigma_max=sigma_max,
            channel_axis=0,
        )
    feature_img.append(features_func(ch_img))
feature_img = np.stack(feature_img, axis=0)
    
img = np.transpose(img, [1,0,2,3,4])
    
    
# viewer=napari.Viewer()
# viewer.add_image(np.transpose(feature_img, [0,4,1,2,3]), name="image")
# napari.run()  
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
    
def segment_and_update(
    tcell_edt_threshold: float = 2.5,
    ):
    print("Running segmentation")
    current_timepoint = int(viewer.dims.current_step[0])
    # feature_img_t = feature_img[current_timepoint]
    
    # Access the label layer and feature image
    label_layer = viewer.layers['User Provided Labels']
    label_data = label_layer.data
    
    clf = RandomForestClassifier(
        n_estimators=100,
        n_jobs=-1, 
        max_depth=20, 
        max_samples=0.05,
        class_weight="balanced"
        )
    
    clf = future.fit_segmenter(label_data, feature_img, clf)
    pixel_class = future.predict_segmenter(feature_img, clf)

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
            mask=mask_organoid,
            segment_size_min=10,
            segment_splitting_edt=tcell_edt_threshold,
            use_dims=3
        )

        full_seg_organoid.append(seg_organoid)
        full_seg_tcell.append(seg_tcell)
    full_seg_organoid = np.stack(full_seg_organoid, axis=0)
    full_seg_tcell = np.stack(full_seg_tcell, axis=0)
        
    viewer.layers["Pixel Classification"].data = pixel_class
    viewer.layers["Organoid Segments"].data = full_seg_organoid
    viewer.layers["Tcell Segments"].data = full_seg_tcell
    print("DONE")
    
# Create Napari viewer
viewer = napari.Viewer()
img_layer = viewer.add_image(img, name="Image", contrast_limits=(0,np.percentile(img[1], 99)))
img_layer.contrast_limits_range = (img.min(), img.max())

viewer.add_labels(np.zeros(img.shape[1:]).astype(np.int16), name="User Provided Labels", opacity=0.3)
viewer.add_labels(np.zeros(img.shape[1:]).astype(np.int16), name="Pixel Classification", opacity=0.7, visible=False)
viewer.add_labels(np.zeros(img.shape[1:]).astype(np.int16), name="Organoid Segments", opacity=0.7, visible=False)
viewer.add_labels(np.zeros(img.shape[1:]).astype(np.int16), name="Tcell Segments", opacity=0.7, visible=False)
# Create interactive sliders
# Create interactive sliders

# Create interactive sliders
gui = magicgui(segment_and_update, 
               tcell_edt_threshold={"widget_type": "FloatSlider", "min": 1.0, "max": 15.0, "step": 0.5}
              )
viewer.window.add_dock_widget(gui)
napari.run()
# viewer = NapariViewer(
#     img = [img_f, mask_dead, segmentation_results["segments"]],
#     label = ["image", "image",  "label"],    
# )
# viewer.view()
