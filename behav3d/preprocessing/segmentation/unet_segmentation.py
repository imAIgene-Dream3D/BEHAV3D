import sys
import time
import logging

import numpy as np
from scipy.ndimage import binary_fill_holes
from skimage.segmentation import watershed, relabel_sequential
from skimage.measure import label
from behav3d.preprocessing import calculate_edt, filter_median, sauvola_thresholding, dilate_mask, open_mask
from behav3d.preprocessing.segmentation import segment_size_filter, get_border_segments, remove_boundary_segments
from pathlib import Path
import torch
from behav3d.preprocessing.segmentation import unet_segmentation
from behav3d.io.images import load_image, save_as_zarr, load_zarr, get_filepath_stem
import zarr
import dask.array as da
from tqdm import tqdm

class BEHAV3D_Unet_Segmenter():
    def __init__(
        self, 
        img=None, 
        img_path=None,
        output_dir=None,
        sample_name=None,
        
        model=None,
        model_path=None,
        
        tcell_ch=0, 
        live_ch=1, 
        dead_ch=2, 
        
        use_dims = 3,
        
        unet_prob_thr = 0.7,
        
        dead_opening_nr_pixels = 1,
        dead_sauvola_window_size=4,
        dead_smooth_radius=0,
        dead_SNR=8,
        dead_peaks_SNR=3,
        
        #Segmentation
        organoids_segments_prev_tp = None,
        organoid_opening_nr_pixels = 3,
        organoid_dilation_nr_pixels = 2,
        organoid_segment_splitting_edt=10,
        
        tcell_opening_nr_pixels = 1,
        tcell_segment_splitting_edt = 2.5,
        tcell_segment_size_min=50,
        
        #Segment filtering
        remove_border_segments=False,
        organoid_start_segment_size_min=1000,
        organoid_segment_size_min=100,
        
        logger=None,
        verbose=True,
        overwrite=False,
        **kwargs
        ):
        self.img = img
        self.img_path = img_path
        if sample_name is None:
            if img_path is not None:
                sample_name = get_filepath_stem(img_path)
            else:
                raise ValueError("No sample name provided and no image path provided to extract sample name from")
        self.sample_name = sample_name
        
        if self.img_path is not None:
            self.img_path = Path(self.img_path)
            
        if output_dir is None and self.img_path is not None:
            self.output_dir = self.img_path.parent
        elif output_dir is None:
            self.output_dir = Path.cwd()
        else:
            self.output_dir = Path(output_dir)
            
        self.model= model
        self.tcell_ch = tcell_ch
        self.live_ch = live_ch
        self.dead_ch = dead_ch
        
        self.logger = logger

        self.verbose = verbose
        if self.logger is None:
            self.logger = self.setup_logging(verbose=self.verbose)
        
        self.overwrite = overwrite
        
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")    
        self.model_path = model_path
        if model_path is not None:
            self.model = torch.load(model_path, map_location=self.device)
            
        self.mask_tcell = None
        self.mask_organoids = None
        self.mask_dead = None
        
        self.dead_SNR = dead_SNR
        self.dead_peaks_SNR = dead_peaks_SNR
        
        self.unet_prob_thr = unet_prob_thr
        
        self.use_dims                       = use_dims
        self.tcell_opening_nr_pixels        = tcell_opening_nr_pixels
        self.organoid_opening_nr_pixels     = organoid_opening_nr_pixels
        self.dead_opening_nr_pixels         = dead_opening_nr_pixels
        self.dead_sauvola_window_size       = dead_sauvola_window_size
        self.dead_smooth_radius             = dead_smooth_radius
        
        self.organoids_segments_prev_tp     = organoids_segments_prev_tp
        self.organoid_dilation_nr_pixels    = organoid_dilation_nr_pixels
        self.organoid_segment_splitting_edt = organoid_segment_splitting_edt
        self.tcell_segment_splitting_edt    = tcell_segment_splitting_edt
        
        #Segment filtering
        self.organoid_start_segment_size_min      = organoid_start_segment_size_min
        self.organoid_segment_size_min      = organoid_segment_size_min
        self.tcell_segment_size_min         = tcell_segment_size_min
        self.remove_border_segments         = remove_border_segments
        
    def log_if_verbose(self, message, level="info"):
        if self.verbose:
            log_func = getattr(self.logger, level, self.logger.info)
            log_func(message)
        
    def apply_unet(self, img, unet_prob_thr=None, model_path=None, normalize=True):
        if unet_prob_thr is not None:
            self.unet_prob_thr = unet_prob_thr
        
        if model_path is not None:
            self.model_path = model_path
            self.model = torch.load(model_path, map_location=self.device)
        self.model.eval()

        model_input = img[[self.tcell_ch,  self.live_ch, self.dead_ch]]
        if normalize:
            try:
                model_input=(model_input/np.iinfo(model_input.dtype).max).astype(np.float32)
            except:
                model_input=(model_input/np.finfo(model_input.dtype).max).astype(np.float32)
        model_input = np.expand_dims(model_input, axis=0)
        model_input = torch.from_numpy(model_input.astype(np.float32)).float().to(self.device)
        
        def multiclass_predict(preds):
            preds_softmax = torch.nn.functional.softmax(preds, dim=1)
            max_idx = torch.argmax(preds, dim=1, keepdim=True)
            preds_thr = torch.empty(preds.shape, device=preds.device)
            preds_thr.zero_()
            preds_thr.scatter_(1, max_idx, 1)
            return(preds_thr, preds_softmax, max_idx)

        with torch.no_grad():
            output = self.model(model_input)
            preds_thr, preds_prob, _ = multiclass_predict(output)
            # preds_thr = preds_thr.cpu().numpy()
            preds_prob = np.squeeze(preds_prob.cpu().numpy())
            preds_thr = preds_prob >= self.unet_prob_thr
        
        mask_organoid, mask_tcell = preds_thr[1], preds_thr[2]
        return(mask_organoid, mask_tcell)
        # self.mask_organoid = preds_thr[1]
        # self.mask_tcell = preds_thr[2]
    
    @staticmethod
    def segment_foreground(
        img,
        smooth_radius,
        SNR_threshold,
        sauvola_window_size,
        use_dims,
        SNR_percentage,
        ):
        if smooth_radius>0:
            smooth_img = filter_median(img, radius=smooth_radius, use_dimensions=use_dims)
        else:
            smooth_img = img
            
        SNR_thr_absmin = float(np.percentile(smooth_img, SNR_percentage)*SNR_threshold)
        # fg_mask = sauvola_thresholding(
        #     smooth_img, 
        #     use_dimensions=use_dims, 
        #     window_size=sauvola_window_size, 
        #     absmin=SNR_thr_absmin
        #     )
        fg_mask = smooth_img >= SNR_thr_absmin
        return(fg_mask)
    
    def mask_dead_fg(self, img):
         # Masking the dead channel of the image
        mask_dead = self.segment_foreground(
                img=img[self.dead_ch],
                smooth_radius=self.dead_smooth_radius,
                SNR_threshold=self.dead_SNR,
                sauvola_window_size=self.dead_sauvola_window_size,
                use_dims=self.use_dims,
                SNR_percentage=90
            )
        
        # Remove single separeted pixel values of dead mask to reduce noise, real signal is often several pixels
        if self.dead_opening_nr_pixels!=0 and self.dead_opening_nr_pixels is not None:
            mask_dead_opened = open_mask(mask_dead, use_dimensions=2, nr_pixels=self.dead_opening_nr_pixels)
            # Refill more pixel-precise original dead mask with opened dead mask to reconnect dead pixels belonging to real signal
            mask_dead = watershed(mask_dead, markers = mask_dead_opened, mask=mask_dead)

        return(mask_dead)
    
    def segment_filtering_mask_dead_peaks(self):
        # Filter the dead signal based on per segment internsity values
        # Perform second round of thresholding to get death peak signal as the peaks of dead signal 
        # are often higher than the background signal within organoids
        smooth_dead = filter_median(
            self.img[self.dead_ch], 
            radius=self.dead_smooth_radius, 
            use_dimensions=self.use_dims
            )
        
        self.mask_dead_peaks = np.zeros_like(self.mask_dead)
        self.dead_foreground_segment_thr = {}
        for segment_id in np.unique(self.unfiltered_organoid_segments):
            if segment_id==0:
                continue
            segment_mask = self.unfiltered_organoid_segments==segment_id
            dead_foreground_thr = np.percentile(smooth_dead[segment_mask], 20)*self.dead_peaks_SNR
            self.dead_foreground_segment_thr[int(segment_id)]=float(dead_foreground_thr)
            self.mask_dead_peaks[(smooth_dead >= dead_foreground_thr) &  (segment_mask)]=1
        
        # Only keep dead signal inside of a segment
        # As dead mask filtering is dependent per organoid base level dead dye
        # Border segments wouldnt do this filtering and thus show a lot of noise in the dead mask
        # That's why this is filtered
        self.mask_dead_peaks[self.unfiltered_organoid_segments==0]=0
        self.mask_dead_peaks[self.combined_mask_tcell==1]=0
        
    def postprocess_all_masks(self):
        self.mask_tcell = self.postprocess_mask(self.mask_tcell, self.tcell_opening_nr_pixels)
        self.mask_organoid = self.postprocess_mask(self.mask_organoid, self.organoid_opening_nr_pixels)
        self.mask_dead = self.postprocess_mask(self.mask_dead, self.dead_opening_nr_pixels)
        
    @staticmethod
    def postprocess_mask(mask, fill_holes=True, opening_nr_pixels=1):
        if fill_holes:
            mask = binary_fill_holes(mask)
        if opening_nr_pixels>0 and opening_nr_pixels is not None:
            mask = open_mask(mask)
        return(mask)
        # if closing_nr_pixels>0:
        #     mask = close_mask(mask)
    
    def segment_organoids(self, mask_organoid, mask_dead, mask_tcell):
        ### For organoids specifically, remove the T cell masks from dead
        ### If not, the organoid mask will expand through dead T cells, giving wrong segmentation
        ### Especially when photoswitching, this will lead to wrong mask photowitching near, but not interacting T cells
        filt_dead_mask = mask_dead.copy()
        filt_dead_mask[mask_tcell==1]=0
        if self.dead_opening_nr_pixels!=0:
           filt_dead_mask = open_mask(filt_dead_mask, use_dimensions=self.use_dims, nr_pixels=self.dead_opening_nr_pixels)
        combined_mask = np.logical_or(filt_dead_mask, mask_organoid)
        
        offset=1
        if self.organoids_segments_prev_tp is None:
            # For the first timepoint, segments are determined more details (see function description)
            segments = self.segment_mask(
                mask=mask_organoid,
                segment_size_min=self.organoid_start_segment_size_min,
                segment_splitting_edt=self.organoid_segment_splitting_edt,
                use_dims=self.use_dims
            )
            combined_mask_org_dilated = dilate_mask(combined_mask, nr_pixels=self.organoid_dilation_nr_pixels)
            unfiltered_organoid_segments = watershed(combined_mask_org_dilated, markers = segments, mask=combined_mask_org_dilated)
            unfiltered_organoid_segments[combined_mask==0]=0
            self.border_segments_ids = get_border_segments(unfiltered_organoid_segments)
            
        else:
            # self.log_if_verbose(f"- Performing segmentation based on segments of previous timepoint")

            # For the subsequent timepoint, segments are copied from previous timepoint and refinedc
            seeds = self.organoids_segments_prev_tp.copy()
            seeds[combined_mask==0]=0
            # seeds = keep_largest_connected_components(seeds)
            unfiltered_organoid_segments = watershed(combined_mask, markers = seeds, mask=combined_mask)
            # logger.info(f"{np.unique(unfiltered_organoid_segments)}")
            # segments = segment_size_filter(segments, size_min=segment_size_min)
            mask_org_dilated = dilate_mask(combined_mask, nr_pixels=self.organoid_dilation_nr_pixels)
            unfiltered_organoid_segments = watershed(mask_org_dilated, markers = unfiltered_organoid_segments, mask=mask_org_dilated)
            unfiltered_organoid_segments[combined_mask==0]=0
            unfiltered_organoid_segments = segment_size_filter(unfiltered_organoid_segments, size_min=self.organoid_segment_size_min)
        
        if self.remove_border_segments:
            organoid_segments = remove_boundary_segments(
                unfiltered_organoid_segments,
                border_segments=self.border_segments_ids
                ) 
        else:
            organoid_segments = unfiltered_organoid_segments.copy()
        
        self.organoids_segments_prev_tp = unfiltered_organoid_segments.copy()
        return(organoid_segments, unfiltered_organoid_segments)
    
    def segment_tcells(self, mask_tcell):
       tcell_segments = self.segment_mask(
                mask=mask_tcell,
                segment_size_min=self.tcell_segment_size_min,
                segment_splitting_edt=self.tcell_segment_splitting_edt,
                use_dims=self.use_dims
            )
       return(tcell_segments)
    
    @staticmethod
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
    
    def run_single_timepoint(self, img):
        # self.logger.info("Applying supplied Unet Model for organoid and T-cell segmentation")
        mask_organoid, mask_tcell = self.apply_unet(img)
        
        # self.logger.info("- Masking the dead channel")
        mask_dead = self.mask_dead_fg(img)
        
        # self.logger.info("- Postprocessing all the created masks")  
        mask_tcell = self.postprocess_mask(mask_tcell, opening_nr_pixels=self.tcell_opening_nr_pixels)
        mask_organoid = self.postprocess_mask(mask_organoid, opening_nr_pixels=self.organoid_opening_nr_pixels)
        mask_dead = self.postprocess_mask(mask_dead, opening_nr_pixels=0)
        
        # self.logger.info("- Segmenting the organoids (with dead cells)")
        organoid_segments, unfiltered_organoid_segments = self.segment_organoids(
            mask_organoid=mask_organoid,
            mask_dead=mask_dead,
            mask_tcell=mask_tcell
        )
        
        # self.logger.info("- Segmenting the T cells")
        tcell_segments = self.segment_tcells(mask_tcell)
        
        return(tcell_segments, organoid_segments, mask_dead)
        
    def run(self, img=None):
        
        if img is not None:
            self.img = img
            
        elif self.img_path is not None:
            image_zarr_outpath = Path(self.output_dir, f"{self.sample_name}.zarr")
            if not Path(image_zarr_outpath).exists():
                self.logger.info(f"Convert image to .zarr for dask processing")
                self.logger.info(f"Saving to {image_zarr_outpath}")
                self.img = load_image(self.img_path)                
                chunks = (1,) + self.img.shape[1:]
                save_as_zarr(
                    img=self.img, 
                    path=image_zarr_outpath, 
                    chunks=chunks
                    )
            self.img = load_zarr(image_zarr_outpath)
        else:
            self.logger.error("- No image provided... Nothing to run processing on")
            
        if self.img.ndim < 4:
            self.logger.info("Provided images does not have enough dimensions for a BEHAV3D experiment and segmentation... Exiting")
        elif self.img.ndim == 4:
            self.logger.info("Segmenting single timepoint BEHAV3D image")
            self.tcell_segments, self.organoid_segments, self.mask_dead = self.run_single_timepoint(self.img)
        elif self.img.ndim == 5:
            self.logger.info(f"Segmenting multiple timepoint BEHAV3D image: {self.img.shape[0]} timepoints")
            self.first_timepoint = True
            
            self.tcell_segments_outpath = Path(self.output_dir, f"{self.sample_name}_tcell_segments.zarr")
            self.organoid_segments_outpath = Path(self.output_dir, f"{self.sample_name}_organoid_tracked.zarr")
            self.mask_dead_outpath = Path(self.output_dir, f"{self.sample_name}_mask_dead.zarr")
            
            if (Path(f"{self.tcell_segments_outpath}").exists() and 
                Path(f"{self.organoid_segments_outpath}").exists() and 
                Path(f"{self.mask_dead_outpath}").exists() and 
                self.overwrite==False
                ):
                self.logger.info("Segmentation already exists... Provide overwrite=True to overwrite")
                
                self.tcell_segments = load_zarr(Path(self.tcell_segments_outpath))
                self.organoid_segments = load_zarr(Path(self.organoid_segments_outpath))
                self.mask_dead = load_zarr(Path(f"{self.mask_dead_outpath}"))
                # return(self.tcell_segments, self.organoid_segments, self.mask_dead)
            else:
                for t, t_img in tqdm(enumerate(self.img), total=self.img.shape[0]):
                # for t, t_img in enumerate(self.img):
                    t_img = np.asarray(t_img)
                    tcell_segments, organoid_segments, mask_dead = self.run_single_timepoint(t_img)
                    
                    tcell_segments = np.expand_dims(tcell_segments, axis=0)
                    organoid_segments = np.expand_dims(organoid_segments, axis=0)
                    mask_dead = np.expand_dims(mask_dead, axis=0)
                    
                    self._append_to_zarr(tcell_segments, self.tcell_segments_outpath)
                    self._append_to_zarr(organoid_segments, self.organoid_segments_outpath)
                    self._append_to_zarr(mask_dead, self.mask_dead_outpath)
                    self.first_timepoint = False
                
                self.tcell_segments = load_zarr(self.tcell_segments_outpath)
                self.organoid_segments = load_zarr(self.organoid_segments_outpath)
                self.mask_dead = load_zarr(self.mask_dead_outpath)
        
            
    def _append_to_zarr(self, img, outpath):
        """
        Append a timepoint to an existing .zarr array
        If non-existent, create the .zarr array
        """
        outpath = Path(outpath)
        
        if not outpath.exists() or self.first_timepoint:
            zarr_file = zarr.open(
                outpath, 
                mode='w', 
                shape=(0,) + img.shape[1:], 
                chunks=(1,) + img.shape[1:], 
                dtype=img.dtype
                )
        else:
            zarr_file = zarr.open(outpath, mode='a')
        zarr_file.append(img)
          
    @staticmethod
    def setup_logging(log_path=None, verbose=True, logname="logger"):
        """
        Set up the logger for writing output to a logfile
        
        If self.interactive==True, also output to the terminal
        """
        logger = logging.getLogger(logname)
        if verbose:
            logger.setLevel(logging.INFO)
        else:
            logger.setLevel(logging.WARNING)
        
        if verbose or log_path is None:
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter('%(asctime)s - [%(levelname)s] - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            console_handler.setFormatter(console_formatter)
            logger.addHandler(console_handler)
            # logger.info("Logging to console")
        
        if log_path is not None:
            log_path=Path(log_path)
            if log_path.exists():
                open(log_path, 'w').close()
            file_handler = logging.FileHandler(log_path)
            if verbose:
                console_handler.setLevel(logging.INFO)
            else:
                console_handler.setLevel(logging.WARNING)
            file_formatter = logging.Formatter('%(asctime)s - [%(levelname)s] - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
            logger.info(f"Logging to file: {log_path}")
        return(logger)

def run_behav3d_unet_segmentation(
    metadata, 
    output_dir,
    model_path,
    tcell_ch=0, 
    live_ch=1, 
    dead_ch=2,
    use_dims = 3,
    unet_prob_thr = 0.7,
    logger=None,
    **kwargs
    ):
    for idx, sample in metadata.iterrows():
        print(f"Processing sample: {sample['sample_name']}")
        start_time = time.time()
        sample_name = sample['sample_name']
        raw_image_path = Path(sample['raw_image_path'])
        
        img_outdir = Path(output_dir, "images", sample_name)
        if not img_outdir.exists():
            img_outdir.mkdir(parents=True)

        # img = load_image(raw_image_path)
        segmenter = BEHAV3D_Unet_Segmenter(
            sample_name=sample_name,
            img_path=raw_image_path,
            output_dir=img_outdir,
            model_path=model_path,
            tcell_ch=sample['tcell_channel'], 
            live_ch=sample['live_channel'], 
            dead_ch=sample['dead_channel'],
            use_dims=use_dims,
            unet_prob_thr=unet_prob_thr,
            logger=logger,
            **kwargs
            )
        segmenter.run()
        
        metadata.at[idx, "tcell_segments_image_path"] = str(segmenter.tcell_segments_outpath)
        metadata.at[idx, "organoid_tracks_image_path"] = str(segmenter.organoid_segments_outpath)
    return(metadata)

def run_only_death_segmentation():
    sample_name = "ROCHE_IC1_Exp033_Img08"
    output_dir = f"/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/runs/ROCHE/images/{sample_name}"
    raw_image_path = Path(output_dir, f"{sample_name}.zarr")
    unet_model_path = "/Volumes/T7_Sam/BHVD_BEHAV3D/BEHAV3D_python/Unet_3D_model.pth"

    # mask_test = load_image(mask_dead_outpath)
    img = load_image(raw_image_path)
    img.shape

    segmenter = BEHAV3D_Unet_Segmenter(
            sample_name=sample_name,
            img_path=raw_image_path,
            output_dir=output_dir,
            model_path=unet_model_path,
            tcell_ch=0, 
            live_ch=1, 
            dead_ch=2,
            use_dims=3,
            unet_prob_thr=0.5,
            logger=None,
            )
    segmenter.img = img

    segmenter.mask_dead_outpath = Path(segmenter.output_dir, f"{segmenter.sample_name}_mask_dead.zarr")
    segmenter.dead_ch

    segmenter.first_timepoint = True


    for t, t_img in tqdm(enumerate(segmenter.img), total=segmenter.img.shape[0]):
        t_img = np.asarray(t_img)
        mask_dead = segmenter.mask_dead_fg(t_img)
        mask_dead = segmenter.postprocess_mask(mask_dead, segmenter.dead_opening_nr_pixels)
        mask_dead = np.expand_dims(mask_dead, axis=0)
        segmenter._append_to_zarr(mask_dead, segmenter.mask_dead_outpath)
        segmenter.first_timepoint = False
            # organoid_segments = segmenter.organoid_segments
            # tcell_segments = segmenter.tcell_segments
            # dead_mask = segmenter.mask_dead
            # return(organoid_segments, tcell_segments, dead_mask)
            