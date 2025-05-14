from aicspylibczi import CziFile
from pathlib import Path
import numpy as np
from magicgui import magicgui
import napari

path = input("Provide a path to .czi to test configuration on (it takes first, middle and last timepoint)")

path = "/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/data/Jess_ROCHE/ROCHE_JM1_Exp011_Img04_169M_50KTcells_withTCB.czi"
path = Path(path)
czi = CziFile(path)

max_t = czi.get_dims_shape()[0]["T"][1] - 1

img_f, _ = czi.read_image(T=0)
img_f = np.squeeze(img_f)
img_f = calc_z_projection(img_f, z_axis=-3, projection='max')


img_m, _ = czi.read_image(T=int(max_t/2))
img_m = np.squeeze(img_m)
img_m= calc_z_projection(img_m, z_axis=-3, projection='max')

img_l, _ = czi.read_image(T=max_t)
img_l = np.squeeze(img_l)
img_l = calc_z_projection(img_l, z_axis=-3, projection='max')

img = np.stack([img_f, img_m, img_l])

def segment_and_update(
    # img=img,
    # tcell_ch = None,
    # live_ch = 1,
    # dead_ch =  2,
    live_SNR: float = 6,
    dead_SNR: float = 6,
    dead_peaks_SNR: float = 6,
    tcell_SNR: float =15,
    segment_splitting_edt: int = 10,
    remove_border_segments = False,
    segment_size_min = 200,
    ):

    current_timepoint = int(viewer.dims.current_step[0])
    img_t = img[current_timepoint]
    segmentation_results = segment_organoids(
        img_t,
        use_dims=2,
        segments_prev_tp=None,
        border_segments_ids=None,
        tcell_ch = None,
        live_ch = 1,
        dead_ch =  2,
        live_SNR = live_SNR,
        dead_SNR = dead_SNR,
        # dead_peaks_SNR = 6,
        tcell_SNR = tcell_SNR,
        remove_border_segments = remove_border_segments,
        segment_size_min = segment_size_min,
        segment_splitting_edt = segment_splitting_edt,
        return_segments_only=False
    )

    # self.segments_t0=self.segments.copy()
    # self.border_segments_ids = segmentation_results["border_segments_ids"]
    # self.channel_identity = segmentation_results["channel_identity"]
    # Filter the dead signal based on per segment internsity values
    # Perform second round of thresholding to get death peak signal as the peaks of dead signal 
    # are often higher than the background signal within organoids
    # dead_peaks_SNR = 6
    mask_dead = np.zeros_like(segmentation_results["mask_dead"])

    for segment_id in np.unique(segmentation_results["unfiltered_segments"]):
        if segment_id==0:
            continue
        segment_mask = segmentation_results["unfiltered_segments"]==segment_id
        dead_foreground_thr = np.percentile(segmentation_results["smooth_dead"][segment_mask], 20)*dead_peaks_SNR
        # self.dead_foreground_segment_thr[int(segment_id)]=float(dead_foreground_thr)
        mask_dead[(segmentation_results["smooth_dead"] >= dead_foreground_thr) &  (segment_mask)]=1

    # Only keep dead signal inside of a segment
    # As dead mask filtering is dependent per organoid base level dead dye
    # Border segments wouldnt do this filtering and thus show a lot of noise in the dead mask
    # That's why this is filtered
    mask_dead[segmentation_results["unfiltered_segments"]==0]=0
    mask_dead[segmentation_results["mask_tcell"]==1]=0


    # Update layers
    if "mask_dead" in viewer.layers:
        viewer.layers["mask_dead"].data = mask_dead
    else:
        viewer.add_labels(mask_dead, name="mask_dead")

    if "segments" in viewer.layers:
        viewer.layers["segments"].data = segmentation_results["segments"]
    else:
        viewer.add_labels(segmentation_results["segments"], name="segments")

# Create Napari viewer
viewer = napari.Viewer()
viewer.add_image(img, name="image")

# Create interactive sliders
gui = magicgui(segment_and_update, 
               live_SNR={"widget_type": "FloatSlider", "min": 1.0, "max": 15.0, "step": 0.5},
               dead_SNR={"widget_type": "FloatSlider", "min": 1.0, "max": 15.0, "step": 0.5},
               dead_peaks_SNR={"widget_type": "FloatSlider", "min": 1.0, "max": 15.0, "step": 0.5},
               tcell_SNR={"widget_type": "FloatSlider", "min": 1, "max": 15, "step": 0.5},
               segment_splitting_edt={"min": 1, "max": 15, "step": 1},
               segment_size_min = {"min": 0, "max": 1000, "step": 50},
               remove_border_segments={"widget_type": "Checkbox"}
              )
gui.changed.connect(segment_and_update)  # Auto-update when sliders change

def update_on_time_change(event):
    segment_and_update(
        live_SNR=gui.live_SNR.value, 
        dead_SNR=gui.dead_SNR.value, 
        dead_peaks_SNR=gui.dead_peaks_SNR.value, 
        tcell_SNR=gui.tcell_SNR.value, 
        segment_size_min=gui.segment_size_min.value, 
        segment_splitting_edt=gui.segment_splitting_edt.value,
        remove_border_segments=gui.remove_border_segments.value
    )
    
viewer.dims.events.current_step.connect(update_on_time_change)
viewer.window.add_dock_widget(gui)

# Initial segmentation
segment_and_update()

napari.run()
# viewer = NapariViewer(
#     img = [img_f, mask_dead, segmentation_results["segments"]],
#     label = ["image", "image",  "label"],    
# )
# viewer.view()
