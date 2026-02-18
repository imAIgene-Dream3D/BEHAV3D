
import napari
import yaml
from magicgui.widgets import create_widget
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox, 
    QStackedWidget, QPushButton, QGroupBox, QFormLayout, 
    QSpinBox, QDoubleSpinBox, QCheckBox, QFileDialog, QScrollArea,
    QPlainTextEdit, QTextEdit
)
from qtpy.QtCore import Qt
import pandas as pd
from pathlib import Path
import traceback
import numpy as np
import dask.array as da
import shutil
import gc
import os
from functools import partial
from skimage import feature
# Scikit-learn imports for training
from sklearn.ensemble import RandomForestClassifier
import joblib


from behav3d.io.images import load_image, get_image_shape, save_as_zarr, append_to_zarr, load_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape
from behav3d.preprocessing.segmentation.napari_pixelclassifier import (
    train_pixel_classifier, 
    run_pixel_classifier_segmentation,
    segment_mask,
    postprocess_mask
)


# Feature calculation function (same as in napari_pixelclassifier.py)
sigma_min = 1
sigma_max = 16
features_func = partial(
        feature.multiscale_basic_features,
        intensity=True,
        edges=True,
        texture=False,
        sigma_max=sigma_max,
        channel_axis=0,
    )

class SegmentationTab(QWidget):
    def __init__(self, viewer: napari.Viewer, metadata_loader):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)

        # Method Selection
        method_group = QGroupBox("Segmentation Method")
        method_layout = QHBoxLayout()
        self.method_combo = QComboBox()
        self.method_combo.addItems(["Pixel Classifier (Random Forest)", "Cellpose (Deep Learning)", "Import segmentation"])
        self.method_combo.currentIndexChanged.connect(self._on_method_changed)
        method_layout.addWidget(QLabel("Method:"))
        method_layout.addWidget(self.method_combo)
        method_group.setLayout(method_layout)
        layout.addWidget(method_group)

        # Stacked Widget for Method Parameters
        self.param_stack = QStackedWidget()
        
        # 1. Pixel Classifier Page
        self.pixel_classifier_page = PixelClassifierWidget(self.viewer, self.metadata_loader, log_callback=self._log)
        self.param_stack.addWidget(self.pixel_classifier_page)

        # 2. Cellpose Page
        self.cellpose_page = CellposeWidget(self.viewer, self.metadata_loader)
        self.param_stack.addWidget(self.cellpose_page)

        # 3. Import Page
        self.import_page = ImportWidget(self.viewer, self.metadata_loader)
        self.param_stack.addWidget(self.import_page)


        layout.addWidget(self.param_stack)
        
        # Add stretch to keep things at the top
        # layout.addStretch()


        # Log Window
        self.log = QTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumHeight(120)
        self.log.setStyleSheet("font-family: monospace; font-size: 11px;")
        layout.addWidget(QLabel("Log"))
        layout.addWidget(self.log)

    def _log(self, msg):
        import datetime
        timestamp = datetime.datetime.now().strftime("%H:%M:%S")
        self.log.append(f"[{timestamp}] {msg}")
        self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum())

    def _on_method_changed(self, index):
        self.param_stack.setCurrentIndex(index)

def _cfg_get(cfg: dict, dotted_key: str, default=None):
    """Get a value from a nested dict using dotted key notation."""
    cur = cfg
    for part in dotted_key.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return cur


class PixelClassifierWidget(QWidget):
    def __init__(self, viewer, metadata_loader, log_callback=None):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self.log = log_callback if log_callback else print
        self._detect_cell_types()
        self._init_ui()
        self.all_images = None
        self.all_features = None

    def _detect_cell_types(self):
        from behav3d.core.metadata import (
            detect_organoid_types_from_metadata,
            detect_immune_cell_types_from_metadata,
            detect_other_cell_types_from_metadata,
            has_dead_channel
        )
        metadata = self.metadata_loader.metadata
        
        if metadata is None:
            self.organoid_types = []
            self.immune_types = []
            self.other_types = []
            self.has_death = False
            self.all_cell_types = []
            return

        self.organoid_types = detect_organoid_types_from_metadata(metadata)
        self.immune_types = detect_immune_cell_types_from_metadata(metadata)
        self.other_types = detect_other_cell_types_from_metadata(metadata)
        self.has_death = has_dead_channel(metadata)
        self.all_cell_types = self.organoid_types + self.immune_types + self.other_types

    def _init_ui(self):
        # Connect to metadata updates
        if hasattr(self.metadata_loader, 'metadata_loaded'):
            self.metadata_loader.metadata_loaded.connect(self._on_metadata_updated)

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        content = QWidget()
        content_layout = QVBoxLayout()
        content_layout.setContentsMargins(4, 4, 4, 4)
        content_layout.setSpacing(4)
        content.setLayout(content_layout)
        scroll.setWidget(content)
        layout.addWidget(scroll)

        # Read saved config values
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})

        # Training Configuration
        train_group = QGroupBox("1. Labeling & Training")
        train_form = QFormLayout()
        train_form.setContentsMargins(6, 6, 6, 6)
        train_form.setSpacing(4)
        train_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        
        self.spin_examples = QSpinBox()
        self.spin_examples.setValue(int(pc.get("examples_per_sample", 3)))
        self.spin_examples.setRange(1, 10)
        self.spin_examples.setMaximumWidth(70)
        train_form.addRow("Examples/sample:", self.spin_examples)
        
        self.check_sample_specific = QCheckBox("Per-sample classifier")
        self.check_sample_specific.setChecked(bool(pc.get("sample_specific_classifier", False)))
        train_form.addRow(self.check_sample_specific)
        
        self.spin_workers = QSpinBox()
        self.spin_workers.setValue(int(pc.get("workers", os.cpu_count() or 4)))
        self.spin_workers.setRange(1, 128)
        self.spin_workers.setMaximumWidth(70)
        train_form.addRow("Workers:", self.spin_workers)
        
        self.btn_load_training = QPushButton("Load Training Data")
        self.btn_load_training.setToolTip("Clears viewer and loads selected timepoints for labeling")
        self.btn_load_training.setStyleSheet("background-color: #007bff; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_load_training.clicked.connect(self._on_load_training_clicked)
        
        self.btn_save_labels = QPushButton("Save Labels")
        self.btn_save_labels.setToolTip("Save all user-provided label layers to disk")
        self.btn_save_labels.setStyleSheet("background-color: #fd7e14; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_save_labels.clicked.connect(self._save_user_labels)
        
        self.btn_train_update = QPushButton("Train Classifier")
        self.btn_train_update.setStyleSheet("background-color: #28a745; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_train_update.clicked.connect(self._on_train_clicked)
        
        train_btn_row = QHBoxLayout()
        train_btn_row.setSpacing(4)
        train_btn_row.addWidget(self.btn_save_labels)
        train_btn_row.addWidget(self.btn_train_update)
        
        train_layout = QVBoxLayout()
        train_layout.addLayout(train_form)
        train_layout.addWidget(self.btn_load_training)
        train_layout.addLayout(train_btn_row)
        
        train_group.setLayout(train_layout)
        
        content_layout.addWidget(train_group)

        # Segmentation Parameters (Dynamic)
        self.param_group = QGroupBox("2. Segmentation Parameters")
        self.param_layout = QVBoxLayout()
        self.param_group.setLayout(self.param_layout)

        self.btn_resegment = QPushButton("Test Segmentation Parameters")
        self.btn_resegment.setToolTip("Re-runs segmentation on the current slice/time point using current parameters.")
        self.btn_resegment.setStyleSheet("background-color: #17a2b8; color: white; font-weight: bold; border-radius: 4px; padding: 6px;")
        self.btn_resegment.clicked.connect(self._on_resegment_clicked)
        
        self.param_widgets = {} # {cell_type: {param_name: widget}}
        self._refresh_params_ui() # Populate initially
        
        content_layout.addWidget(self.param_group)

        # Execution Configuration
        seg_group = QGroupBox("3. Batch Segmentation")
        seg_form = QFormLayout()
        seg_form.setContentsMargins(6, 6, 6, 6)
        seg_form.setSpacing(4)
        seg_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        
        self.check_overwrite = QCheckBox("Overwrite existing")
        self.check_overwrite.setChecked(bool(pc.get("overwrite_existing", False)))
        seg_form.addRow(self.check_overwrite)

        self.check_process_all = QCheckBox("All timepoints")
        self.check_process_all.setChecked(bool(pc.get("use_all_timepoints", True)))
        self.check_process_all.toggled.connect(self._on_process_all_toggled)
        seg_form.addRow(self.check_process_all)
        
        use_all = bool(pc.get("use_all_timepoints", True))
        self.spin_t_start = QSpinBox()
        self.spin_t_start.setRange(0, 9999)
        self.spin_t_start.setValue(int(pc.get("tp_start", 0)))
        self.spin_t_start.setEnabled(not use_all)
        self.spin_t_start.setMaximumWidth(60)
        self.spin_t_end = QSpinBox()
        self.spin_t_end.setRange(0, 9999)
        self.spin_t_end.setValue(int(pc.get("tp_end", 0)))
        self.spin_t_end.setEnabled(not use_all)
        self.spin_t_end.setMaximumWidth(60)
        
        t_range_layout = QHBoxLayout()
        t_range_layout.setSpacing(4)
        t_range_layout.addWidget(QLabel("T:"))
        t_range_layout.addWidget(self.spin_t_start)
        t_range_layout.addWidget(QLabel("–"))
        t_range_layout.addWidget(self.spin_t_end)
        t_range_layout.addStretch()
        seg_form.addRow("Range:", t_range_layout)
        
        self.btn_run_segmentation = QPushButton("Run Batch Segmentation")
        self.btn_run_segmentation.setStyleSheet("background-color: #28a745; color: white; font-weight: bold; border-radius: 4px; padding: 6px; font-size: 12px;")
        self.btn_run_segmentation.clicked.connect(self._on_run_segmentation_clicked)
        
        seg_layout = QVBoxLayout()
        seg_layout.addLayout(seg_form)
        seg_layout.addWidget(self.btn_run_segmentation)
        seg_group.setLayout(seg_layout)
        
        content_layout.addWidget(seg_group)
        content_layout.addStretch()

    def _on_metadata_updated(self):
        self.log("Metadata updated, refreshing Segmentation Tab...")
        self._detect_cell_types()
        self._refresh_params_ui()

    def _refresh_params_ui(self):
        # Clear existing
        while self.param_layout.count():
            item = self.param_layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
        
        self.param_widgets = {}

        if not self.all_cell_types:
            self.param_layout.addWidget(QLabel("No cell types detected or no metadata loaded."))
            return

        # Read saved per-cell-type parameters from config
        pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})

        for cell_type in self.all_cell_types:
            ct_group = QGroupBox(f"{cell_type.capitalize()}")
            ct_form = QFormLayout()
            ct_form.setContentsMargins(6, 4, 6, 4)
            ct_form.setSpacing(3)
            ct_form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)

            # Determine defaults based on category
            if cell_type in self.organoid_types:
                def_edt, def_size, def_open, def_fill = 12.0, 1000, 3, True
            elif cell_type in self.immune_types:
                def_edt, def_size, def_open, def_fill = 2.5, 10, 0, True
            else:
                def_edt, def_size, def_open, def_fill = 1.0, 10, 0, True
            
            # Override defaults with saved config values if present
            saved_edt = pc.get(f"{cell_type}_edt_threshold", def_edt)
            saved_size = pc.get(f"{cell_type}_segment_size_min", def_size)
            saved_open = pc.get(f"{cell_type}_opening_nr_pixels", def_open)
            saved_fill = pc.get(f"{cell_type}_fill_holes", def_fill)

            w_edt = QDoubleSpinBox()
            w_edt.setRange(0.1, 50.0)
            w_edt.setSingleStep(0.5)
            w_edt.setValue(float(saved_edt))
            w_edt.setMaximumWidth(70)
            ct_form.addRow("EDT:", w_edt)
            
            w_size = QSpinBox()
            w_size.setRange(1, 100000)
            w_size.setValue(int(saved_size))
            w_size.setMaximumWidth(80)
            ct_form.addRow("Min size:", w_size)
            
            w_open = QSpinBox()
            w_open.setRange(0, 10)
            w_open.setValue(int(saved_open))
            w_open.setMaximumWidth(60)
            ct_form.addRow("Opening:", w_open)
            
            w_fill = QCheckBox("Fill holes")
            w_fill.setChecked(bool(saved_fill))
            ct_form.addRow(w_fill)
            
            ct_group.setLayout(ct_form)
            self.param_layout.addWidget(ct_group)
            
            self.param_widgets[cell_type] = {
                'edt': w_edt,
                'min_size': w_size,
                'opening': w_open,
                'fill_holes': w_fill
            }
        
        # Add Test Button at the end of params list
        self.param_layout.addWidget(self.btn_resegment)

    def _persist_params(self):
        """Save all segmentation widget values to behav3d_parameters.yml."""
        pc = self.metadata_loader.behav3d_parameters.setdefault("pixel_classifier", {})
        pc["examples_per_sample"] = int(self.spin_examples.value())
        pc["sample_specific_classifier"] = bool(self.check_sample_specific.isChecked())
        pc["workers"] = int(self.spin_workers.value())
        pc["use_all_timepoints"] = bool(self.check_process_all.isChecked())
        pc["tp_start"] = int(self.spin_t_start.value())
        pc["tp_end"] = int(self.spin_t_end.value())
        pc["overwrite_existing"] = bool(self.check_overwrite.isChecked())

        # Save per-cell-type parameters
        for cell_type in self.all_cell_types:
            if cell_type in self.param_widgets:
                w = self.param_widgets[cell_type]
                pc[f"{cell_type}_edt_threshold"] = float(w['edt'].value())
                pc[f"{cell_type}_segment_size_min"] = int(w['min_size'].value())
                pc[f"{cell_type}_opening_nr_pixels"] = int(w['opening'].value())
                pc[f"{cell_type}_fill_holes"] = bool(w['fill_holes'].isChecked())

        # Write to YAML file
        out_dir = self.metadata_loader.output_dir
        if out_dir:
            params_path = Path(out_dir) / "behav3d_parameters.yml"
            try:
                with open(params_path, "w") as f:
                    yaml.safe_dump(self.metadata_loader.behav3d_parameters, f, sort_keys=False)
            except Exception as e:
                self.log(f"Warning: Could not save parameters: {e}")

    def _save_user_labels(self):
        """Save all user-provided label layers to zarr on disk."""
        if not self.metadata_loader.output_dir:
            self.log("Cannot save labels: no output directory set.")
            return

        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)

        saved_any = False

        if self.has_death:
            layer_name = 'User Provided Labels (Dead)'
            if layer_name in self.viewer.layers:
                outpath = pixel_class_outdir / 'PixelClassifier_UserDeadLabels.zarr'
                if outpath.exists():
                    shutil.rmtree(outpath)
                save_as_zarr(self.viewer.layers[layer_name].data, outpath)
                self.log(f"Saved Dead labels to {outpath.name}")
                saved_any = True

        for cell_type in self.all_cell_types:
            layer_name = f'User Provided Labels ({cell_type.capitalize()})'
            if layer_name in self.viewer.layers:
                outpath = pixel_class_outdir / f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr'
                if outpath.exists():
                    shutil.rmtree(outpath)
                save_as_zarr(self.viewer.layers[layer_name].data, outpath)
                self.log(f"Saved {cell_type} labels to {outpath.name}")
                saved_any = True

        if not saved_any:
            self.log("No user label layers found in viewer to save.")
        else:
            self.log("All user labels saved.")

    def _on_process_all_toggled(self, checked):
        self.spin_t_start.setEnabled(not checked)
        self.spin_t_end.setEnabled(not checked)

    def _on_run_segmentation_clicked(self):
        self.log("Starting batch segmentation...")
        self._persist_params()
        try:
            # Get parameters
            overwrite = self.check_overwrite.isChecked()
            
            # Timepoint range
            if self.check_process_all.isChecked():
                timepoint_range = None
                self.log("Processing ALL timepoints.")
            else:
                t_start = self.spin_t_start.value()
                t_end = self.spin_t_end.value()
                if t_start > t_end:
                    self.log("Error: Start timepoint must be <= End timepoint.")
                    return
                timepoint_range = (t_start, t_end)
                self.log(f"Processing timepoints: {t_start} - {t_end}")
            
            # Build per-category parameter dicts matching run_pixel_classifier_segmentation signature
            org_edt, org_size, org_open, org_fill = {}, {}, {}, {}
            imm_edt, imm_size, imm_open, imm_fill = {}, {}, {}, {}
            oth_edt, oth_size, oth_open, oth_fill = {}, {}, {}, {}

            for cell_type in self.all_cell_types:
                ct_key = cell_type.lower()
                if ct_key in self.param_widgets:
                    w = self.param_widgets[ct_key]
                    edt = float(w['edt'].value())
                    sz  = int(w['min_size'].value())
                    opn = int(w['opening'].value())
                    fh  = bool(w['fill_holes'].isChecked())
                else:
                    edt, sz, opn, fh = 2.5, 10, 0, True

                if cell_type in self.organoid_types:
                    org_edt[cell_type] = edt; org_size[cell_type] = sz
                    org_open[cell_type] = opn; org_fill[cell_type] = fh
                elif cell_type in self.immune_types:
                    imm_edt[cell_type] = edt; imm_size[cell_type] = sz
                    imm_open[cell_type] = opn; imm_fill[cell_type] = fh
                else:
                    oth_edt[cell_type] = edt; oth_size[cell_type] = sz
                    oth_open[cell_type] = opn; oth_fill[cell_type] = fh

            # Call the segmentation function
            run_pixel_classifier_segmentation(
                metadata=self.metadata_loader.metadata,
                output_dir=self.metadata_loader.output_dir,
                organoid_edt_thresholds=org_edt or None,
                immune_edt_thresholds=imm_edt or None,
                other_edt_thresholds=oth_edt or None,
                organoid_segment_size_mins=org_size or None,
                immune_segment_size_mins=imm_size or None,
                other_segment_size_mins=oth_size or None,
                organoid_opening_nr_pixels=org_open or None,
                immune_opening_nr_pixels=imm_open or None,
                other_opening_nr_pixels=oth_open or None,
                organoid_fill_holes=org_fill or None,
                immune_fill_holes=imm_fill or None,
                other_fill_holes=oth_fill or None,
                only_segment=True,
                n_workers=self.spin_workers.value(),
                log_callback=self.log,
                overwrite_existing=overwrite,
                timepoint_range=timepoint_range,
            )
            self.log("Batch segmentation finished successfully!")

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during batch segmentation: {e}")

    def _on_load_training_clicked(self):
        try:
            self.log("Loading training data...")
            self._persist_params()
            self.viewer.layers.clear()
            
            if not self.metadata_loader.metadata_loaded:
                self.log("Metadata not loaded!")
                return

            output_dir = Path(self.metadata_loader.output_dir)
            pixel_class_outdir = output_dir / "images" / "PixelClassification"
            pixel_class_outdir.mkdir(parents=True, exist_ok=True)
            
            features_outpath = pixel_class_outdir / 'PixelClassifier_Features.zarr'
            image_outpath = pixel_class_outdir / 'PixelClassifier_Images.zarr'
            
            metadata = self.metadata_loader.metadata
            examples_per_sample = self.spin_examples.value()
            
            # --- CACHING: Check if cached data can be reused ---
            pc = _cfg_get(self.metadata_loader.behav3d_parameters, "pixel_classifier", {})
            saved_examples = pc.get("examples_per_sample", None)
            can_reuse_cache = (
                features_outpath.exists()
                and image_outpath.exists()
                and saved_examples == examples_per_sample
            )
            
            if can_reuse_cache:
                # --- FAST PATH: load from cached zarr files ---
                self.log("Loading cached images and features from disk (same examples_per_sample)...")
                stacked_images = np.asarray(load_image(image_outpath))
                stacked_features = np.asarray(load_image(features_outpath))
                self.log(f"  Cached images shape: {stacked_images.shape}")
                self.log(f"  Cached features shape: {stacked_features.shape}")
            else:
                # --- SLOW PATH: regenerate images and features ---
                if can_reuse_cache is False and saved_examples is not None and saved_examples != examples_per_sample:
                    self.log(f"examples_per_sample changed ({saved_examples} → {examples_per_sample}), regenerating...")
                else:
                    self.log("No cached data found, generating images and features...")
                
                all_images = []
                
                self.log(f"Checking {len(metadata)} samples...")
                
                for idx, sample in metadata.iterrows():
                    sample_name = sample['sample_name']
                    raw_image_path = Path(sample['raw_image_path'])
                    
                    if not raw_image_path.exists():
                        self.log(f"Error: Raw image path not found: {raw_image_path}")
                        continue
                    
                    raw_image_zarr = output_dir / "images" / sample_name / f"{sample_name}.zarr"
                    
                    if not raw_image_zarr.exists():
                        self.log(f"- Converting {sample_name} to .zarr...")
                        try:
                            images = load_image(raw_image_path)
                            save_as_zarr(img=images, path=raw_image_zarr)
                        except Exception as e:
                            self.log(f"Error loading/converting {sample_name}: {e}")
                            traceback.print_exc()
                            continue
                    
                    try:
                        images = load_image(raw_image_zarr)
                        if len(images.shape) == 5:
                             max_t = images.shape[0]-1
                        elif len(images.shape) == 4:
                            self.log(f"Warning: Image {sample_name} has shape {images.shape}. Assuming T,C,Z,Y,X logic might fail.")
                            max_t = images.shape[0]-1
                        else:
                            self.log(f"Warning: unexpected shape {images.shape}")
                            max_t = 0
                            
                        idc = np.linspace(0, max_t, examples_per_sample, dtype=int)
                        idc = np.unique(idc)
                        
                        self.log(f"  {sample_name}: Taking timepoints {idc}")
                        
                        sample_images = [images[t] for t in idc]
                        all_images += sample_images
                        
                    except Exception as e:
                         self.log(f"Error reading {raw_image_zarr}: {e}")
                         traceback.print_exc()

                if not all_images:
                    self.log("No images loaded.")
                    return

                max_shape = [0] * all_images[0].ndim
                for img in all_images:
                    for i in range(img.ndim):
                        max_shape[i] = max(max_shape[i], img.shape[i])
                max_shape = tuple(max_shape)
                
                self.log(f"Max shape for padding: {max_shape}")
                
                # Remove old cached files before regenerating
                if image_outpath.exists(): shutil.rmtree(image_outpath)
                if features_outpath.exists(): shutil.rmtree(features_outpath)
                
                self.log("Calculating features (this may take a while)...")
                
                processed_images = []
                processed_features = []
                
                for i, img in enumerate(all_images):
                    self.log(f"Processing image {i+1}/{len(all_images)}")
                    try:
                        feats = features_func(img)
                        img_padded = zeropad_image_to_match_shape(img, max_shape)
                        processed_images.append(img_padded)
                        processed_features.append(feats)
                    except Exception as e:
                        self.log(f"Error calculating features: {e}")
                        traceback.print_exc()

                if not processed_images:
                     return

                stacked_images = np.stack(processed_images)
                
                feat_max_shape = [0] * processed_features[0].ndim
                for f in processed_features:
                    for i in range(f.ndim):
                         feat_max_shape[i] = max(feat_max_shape[i], f.shape[i])
                feat_max_shape = tuple(feat_max_shape)
                
                processed_features_padded = []
                for f in processed_features:
                    processed_features_padded.append(zeropad_image_to_match_shape(f, feat_max_shape))
                
                stacked_features = np.stack(processed_features_padded)
                
                save_as_zarr(stacked_images, image_outpath)
                save_as_zarr(stacked_features, features_outpath)
                
                self.log("Images and features saved to disk.")
            
            # --- Common path: store in memory and update config ---
            self.all_images = stacked_images
            self.all_features = stacked_features
            
            # Save the current examples_per_sample to config so we can detect changes next time
            self._persist_params()
            
            # --- Add Layers to Viewer ---
            self.log("Adding layers to viewer...")
            channel_colors = ['cyan', 'yellow', 'red', 'green', 'magenta', 'blue', 'gray']
            
            n_channels = stacked_images.shape[1]
            for ch in range(n_channels):
                channel_data = stacked_images[:, ch, :, :, :]
                
                layer = self.viewer.add_image(
                    channel_data,
                    name=f"Channel {ch+1}",
                    colormap=channel_colors[ch % len(channel_colors)],
                    blending='additive'
                )
                try:
                    flat = channel_data.flatten()
                    flat = flat[flat > 0]
                    if flat.size > 0:
                        clim = np.percentile(flat, 99.8)
                        layer.contrast_limits = (0, clim)
                except:
                    pass
            
            # Create label layer dimensions
            dims = stacked_images.shape  # (T, C, Z, Y, X)
            label_dims = (dims[0], dims[2], dims[3], dims[4])  # (T, Z, Y, X)
            
            # --- Reload existing user labels and predictions from disk ---
            # 1. User Provided Labels
            if self.has_death:
                dead_labels_path = pixel_class_outdir / 'PixelClassifier_UserDeadLabels.zarr'
                if dead_labels_path.exists():
                    self.log("  Loading existing dead labels from disk...")
                    dead_labels = np.asarray(load_zarr(dead_labels_path))
                else:
                    dead_labels = np.zeros(label_dims, dtype=np.uint16)
                self.viewer.add_labels(dead_labels, name='User Provided Labels (Dead)', opacity=0.3)

            for cell_type in self.all_cell_types:
                user_label_path = pixel_class_outdir / f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr'
                if user_label_path.exists():
                    self.log(f"  Loading existing {cell_type} user labels from disk...")
                    user_labels = np.asarray(load_zarr(user_label_path))
                else:
                    user_labels = np.zeros(label_dims, dtype=np.uint16)
                self.viewer.add_labels(
                    user_labels,
                    name=f'User Provided Labels ({cell_type.capitalize()})',
                    opacity=0.3
                )

            # 2. Pixel Classification predictions
            if self.has_death:
                pred_death_path = pixel_class_outdir / 'PixelClassifier_Death_PredictedLabels.zarr'
                if pred_death_path.exists():
                    self.log("  Loading existing dead predictions from disk...")
                    pred_dead = np.asarray(load_zarr(pred_death_path))
                else:
                    pred_dead = np.zeros(label_dims, dtype=np.uint16)
                self.viewer.add_labels(pred_dead, name='Pixel Classification (Dead)', opacity=0.3, visible=False)

            for cell_type in self.all_cell_types:
                pred_label_path = pixel_class_outdir / f'PixelClassifier_{cell_type.capitalize()}_PredictedLabels.zarr'
                if pred_label_path.exists():
                    self.log(f"  Loading existing {cell_type} predictions from disk...")
                    pred_data = np.asarray(load_zarr(pred_label_path))
                else:
                    pred_data = np.zeros(label_dims, dtype=np.uint16)
                self.viewer.add_labels(
                    pred_data,
                    name=f'Pixel Classification ({cell_type.capitalize()})',
                    opacity=0.3,
                    visible=False
                )

            # 3. Segments
            for cell_type in self.all_cell_types:
                 self.viewer.add_labels(
                    np.zeros(label_dims, dtype=np.uint16),
                    name=f'{cell_type.capitalize()} Segments',
                    opacity=0.7,
                    visible=False
                )

            self.log("Training data loaded successfully!")
            
        except Exception as e:
            traceback.print_exc()
            self.log(f"Failed to load training data: {e}")

    def _on_train_clicked(self):
        """Train classifiers, predict pixels, and segment — mirrors original segment_and_update()."""
        self.log("Training classifier & segmenting...")
        self._persist_params()
        self._save_user_labels()  # Always persist labels before training
        
        # ── Load features ───────────────────────────────────────────
        if self.all_features is None:
            output_dir = Path(self.metadata_loader.output_dir)
            pixel_class_outdir = output_dir / "images" / "PixelClassification"
            features_outpath = pixel_class_outdir / 'PixelClassifier_Features.zarr'
            if features_outpath.exists():
                self.log("Loading features from disk...")
                self.all_features = np.asarray(load_zarr(features_outpath))
            else:
                self.log("Error: No features found. Please load training data first.")
                return

        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)
        n_workers = self.spin_workers.value()

        # all_features shape: (T, Z, Y, X, C_feat) — channel LAST
        all_features = np.asarray(self.all_features)
        self.log(f"Features shape: {all_features.shape}")

        # ── Collect targets ─────────────────────────────────────────
        targets = []
        if self.has_death:
            targets.append('Dead')
        targets.extend([t.capitalize() for t in self.all_cell_types])

        try:
            # ── 1. Save & postprocess user labels, then train ───────
            classifiers = {}
            cell_type_labels = {}

            for target in targets:
                label_layer_name = f'User Provided Labels ({target})'
                if label_layer_name not in self.viewer.layers:
                    self.log(f"  Layer {label_layer_name} not found. Skipping.")
                    continue
                
                labels = self.viewer.layers[label_layer_name].data.copy()

                if labels.max() == 0:
                    self.log(f"  No labels for {target}. Skipping.")
                    continue

                # Get postprocessing params for this target
                cell_type_key = target.lower()
                if cell_type_key in self.param_widgets:
                    w = self.param_widgets[cell_type_key]
                    opening = int(w['opening'].value())
                    fill = bool(w['fill_holes'].isChecked())
                elif cell_type_key == 'dead':
                    opening, fill = 0, True
                else:
                    opening, fill = 0, True

                # Postprocess user labels before training (match original lines 594-618)
                fg_mask = (labels == 2).astype(bool)
                if np.any(fg_mask):
                    processed_fg = postprocess_mask(fg_mask, fill_holes=fill, opening_nr_pixels=opening)
                    labels = np.where(processed_fg, 2,
                                      np.where(labels == 1, 1, 0)).astype(labels.dtype)

                # Save user labels to zarr
                labels_outpath = pixel_class_outdir / f'PixelClassifier_User{target}Labels.zarr'
                if labels_outpath.exists():
                    shutil.rmtree(labels_outpath)
                save_as_zarr(labels, labels_outpath)
                cell_type_labels[target] = labels
                self.log(f"  Saved {target} user labels (postprocessed)")

                # Train classifier (match original train_classifier)
                self.log(f"  Training {target}...")
                flat_labels = labels.ravel()
                flat_features = all_features.reshape(-1, all_features.shape[-1])

                label_indices = np.flatnonzero(flat_labels > 0)
                selected_features = flat_features[label_indices]
                selected_labels = flat_labels[label_indices]

                nr_bg = int(np.sum(selected_labels == 1))
                nr_fg = int(np.sum(selected_labels == 2))
                total = nr_bg + nr_fg
                self.log(f"    BG pixels: {nr_bg}, FG pixels: {nr_fg}")

                if total == 0:
                    self.log(f"  No labeled pixels for {target}. Skipping.")
                    continue

                class_weights = {1: nr_bg / total, 2: nr_fg / total}

                clf = RandomForestClassifier(
                    n_estimators=50,
                    n_jobs=-1,
                    max_depth=20,
                    class_weight=class_weights
                )
                from skimage import future
                clf = future.fit_segmenter(selected_labels, selected_features, clf)

                clf_path = pixel_class_outdir / f'PixelClassifier_{target}.joblib'
                joblib.dump(clf, clf_path)
                self.log(f"  Saved classifier: {clf_path.name}")
                classifiers[target] = clf

            # ── 2. Predict all pixels ───────────────────────────────
            pred_masks = {}
            n_timepoints = all_features.shape[0]

            for target, clf in classifiers.items():
                self.log(f"  Predicting {target} pixels...")
                prediction_stack = []

                for t in range(n_timepoints):
                    # all_features[t] shape: (Z, Y, X, C_feat) — already channel-last
                    feat_t = all_features[t]
                    spatial_shape = feat_t.shape[:-1]  # (Z, Y, X)
                    X_t = feat_t.reshape(-1, feat_t.shape[-1])

                    # Predict in batches
                    batch_size = 100_000
                    y_pred_t = np.zeros(X_t.shape[0], dtype=np.uint8)
                    for i in range(0, X_t.shape[0], batch_size):
                        y_pred_t[i:i+batch_size] = clf.predict(X_t[i:i+batch_size])

                    y_pred_t = y_pred_t.reshape(spatial_shape)
                    prediction_stack.append(y_pred_t)

                full_prediction = np.stack(prediction_stack)  # (T, Z, Y, X)

                # Postprocess predictions (match original lines 728-738)
                cell_type_key = target.lower()
                if cell_type_key in self.param_widgets:
                    w = self.param_widgets[cell_type_key]
                    opening = int(w['opening'].value())
                    fill = bool(w['fill_holes'].isChecked())
                elif cell_type_key == 'dead':
                    opening, fill = 0, True
                else:
                    opening, fill = 0, True

                fg_mask = (full_prediction == 2).astype(bool)
                processed_fg = postprocess_mask(fg_mask, fill_holes=fill, opening_nr_pixels=opening)
                full_prediction = np.where(processed_fg, 2,
                                           np.where(full_prediction == 1, 1, 0)).astype(np.uint8)

                pred_masks[target] = full_prediction

                # Save predictions to zarr
                pred_outpath = pixel_class_outdir / f'PixelClassifier_{target}_PredictedLabels.zarr'
                if pred_outpath.exists():
                    shutil.rmtree(pred_outpath)
                save_as_zarr(full_prediction, pred_outpath)

                # Update viewer layer
                pred_layer_name = f'Pixel Classification ({target})'
                if pred_layer_name in self.viewer.layers:
                    self.viewer.layers[pred_layer_name].data = full_prediction
                    # self.viewer.layers[pred_layer_name].visible = True
                    self.viewer.layers[pred_layer_name].refresh()
                else:
                    self.viewer.add_labels(full_prediction, name=pred_layer_name, opacity=0.3)
                self.log(f"  {target} prediction updated.")

            # ── 3. Segment instances (match original lines 757-819) ──
            self.log("Segmenting cell instances...")
            for cell_type in self.all_cell_types:
                target = cell_type.capitalize()
                if target not in pred_masks:
                    continue

                cell_type_key = cell_type.lower()
                if cell_type_key in self.param_widgets:
                    w = self.param_widgets[cell_type_key]
                    edt_threshold = float(w['edt'].value())
                    segment_size_min = int(w['min_size'].value())
                else:
                    if cell_type in self.organoid_types:
                        edt_threshold, segment_size_min = 12.0, 1000
                    elif cell_type in self.immune_types:
                        edt_threshold, segment_size_min = 2.5, 10
                    else:
                        edt_threshold, segment_size_min = 1.0, 10

                self.log(f"  Segmenting {target} (EDT={edt_threshold}, min_size={segment_size_min})...")
                pred_mask = pred_masks[target]

                segmented_timepoints = []
                for t_idx in range(pred_mask.shape[0]):
                    mask_t = (pred_mask[t_idx] == 2)
                    fg_pixels = int(mask_t.sum())

                    if fg_pixels == 0:
                        segmented = np.zeros_like(mask_t, dtype=np.uint16)
                    else:
                        segmented = segment_mask(
                            mask=mask_t,
                            edt_thr=edt_threshold,
                            edt_thr_refined=None,
                            segment_size_min=segment_size_min,
                            use_dims=2,
                            n_workers=1,
                        ).astype(np.uint16, copy=False)

                    segmented_timepoints.append(segmented)

                full_seg = np.stack(segmented_timepoints, axis=0)
                self.log(f"  {target}: max label = {int(full_seg.max())}")

                seg_layer_name = f'{target} Segments'
                if seg_layer_name in self.viewer.layers:
                    self.viewer.layers[seg_layer_name].data = full_seg
                    self.viewer.layers[seg_layer_name].visible = True
                    self.viewer.layers[seg_layer_name].refresh()
                else:
                    self.viewer.add_labels(full_seg, name=seg_layer_name, opacity=0.7)

            self.log("Training, prediction, and segmentation complete!")

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during training/segmentation: {e}")

    def _on_resegment_clicked(self):
        self.log("Testing segmentation parameters on current view...")
        try:
            # 1. Identify active timepoint
            current_step = self.viewer.dims.current_step
            # Assuming (T, Z, Y, X) or (T, C, Z, Y, X)
            # We need to find T index. 
            # If dims order is standard BEHAV3D: T is 0
            t_idx = current_step[0]
            self.log(f"Current timepoint index: {t_idx}")

            # 2. Iterate over cell types
            targets = []
            if self.has_death:
                targets.append('Dead')
            targets.extend([t.capitalize() for t in self.all_cell_types])

            for target in targets:
                # Get Pixel Classification layer (Prediction)
                pred_layer_name = f'Pixel Classification ({target})'
                if pred_layer_name not in self.viewer.layers:
                    self.log(f"  Warning: Layer '{pred_layer_name}' not found. skipping {target}.")
                    continue
                
                pred_layer = self.viewer.layers[pred_layer_name]
                pred_data = pred_layer.data # (T, Z, Y, X)
                
                if t_idx >= pred_data.shape[0]:
                    self.log(f"  Error: Current timepoint {t_idx} out of range for {target}.")
                    continue

                # Get prediction for current timepoint
                # Expected shape (Z, Y, X)
                mask_t = pred_data[t_idx]
                
                # Check for content
                # Background=1, Foreground=2 (usually)
                # But let's check values
                
                # Get parameters for this target
                # If 'Dead', use defaults or specific widget if added (currently Dead doesn't have params widget in UI logic above, need check)
                # The UI logic in _refresh_params_ui iterates over self.all_cell_types only.
                # So 'Dead' might not have adjustable params in UI. 
                # Let's skip 'Dead' if not in param_widgets or use defaults.
                
                cell_type_key = target.lower()
                if cell_type_key == 'dead':
                     # Use defaults for dead
                     edt, size, opening, fill = 2.5, 10, 0, True
                else:
                    if cell_type_key in self.param_widgets:
                        widgets = self.param_widgets[cell_type_key]
                        edt = widgets['edt'].value()
                        size = widgets['min_size'].value()
                        opening = widgets['opening'].value()
                        fill = widgets['fill_holes'].isChecked()
                    else:
                        self.log(f"  Warning: No parameters found for {target}. Using defaults.")
                        edt, size, opening, fill = 2.5, 10, 0, True

                self.log(f"  Segmenting {target} (EDT={edt}, Size={size}, Open={opening})...")
                
                # Preprocess mask (Postprocess prediction)
                # logic from napari_pixelclassifier.segment_and_update:
                # mask_t has labels: 0=background, 1=background_label, 2=foreground
                # extract fg (2)
                
                # Postprocess the prediction mask
                fg_mask = (mask_t == 2).astype(bool)
                
                if np.any(fg_mask):
                    fg_mask = postprocess_mask(fg_mask, fill_holes=fill, opening_nr_pixels=opening)
                    
                # Now Segment
                if np.sum(fg_mask) == 0:
                    segmented = np.zeros_like(fg_mask, dtype=np.uint16)
                    self.log(f"  {target}: Empty mask.")
                else:
                    segmented = segment_mask(
                        fg_mask, 
                        edt_thr=edt, 
                        segment_size_min=size
                    )
                    self.log(f"  {target}: {segmented.max()} objects found.")

                # Update 'Segments' layer
                seg_layer_name = f'{target} Segments' # e.g. "Organoid1 Segments"
                # "Dead" usually doesn't have a segments layer in the original code snippet 
                # (only 'User Provided' and 'Pixel Classification' for Dead).
                # But let's check create_layers logic. 
                # In _on_load_training_clicked:
                # self.viewer.add_labels(..., name=f'{cell_type.capitalize()} Segments', ...) for cell types
                # Not for Dead.
                
                if target == 'Dead':
                     continue 
                
                if seg_layer_name in self.viewer.layers:
                    layer = self.viewer.layers[seg_layer_name]
                    # Check shape
                    if layer.data.shape != pred_data.shape:
                        self.log(f"  Warning: Shape mismatch for {seg_layer_name}. Re-initializing layer data.")
                        layer.data = np.zeros_like(pred_data, dtype=np.uint16)
                    
                    # Update only current timepoint
                    # We need to write to the layer.data
                    # Accessing layer.data returns numpy-like (or tensorstore/dask)
                    # If it's in-memory numpy:
                    try:
                        layer.data[t_idx] = segmented
                        layer.refresh()
                    except Exception as e:
                         self.log(f"  Error updating layer {seg_layer_name}: {e}")
                else:
                    self.log(f"  Layer {seg_layer_name} not found. Creating...")
                    # Create full empty array
                    new_data = np.zeros_like(pred_data, dtype=np.uint16)
                    new_data[t_idx] = segmented
                    self.viewer.add_labels(new_data, name=seg_layer_name, opacity=0.7)

            self.log("Segmentation preview updated.")

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error resegmenting: {e}")

    def _save_user_labels(self):
        """Save all user-provided label layers to zarr on disk."""
        if not self.metadata_loader.output_dir:
            self.log("Cannot save labels: no output directory set.")
            return

        output_dir = Path(self.metadata_loader.output_dir)
        pixel_class_outdir = output_dir / "images" / "PixelClassification"
        pixel_class_outdir.mkdir(parents=True, exist_ok=True)

        saved_any = False

        if self.has_death:
            layer_name = 'User Provided Labels (Dead)'
            if layer_name in self.viewer.layers:
                outpath = pixel_class_outdir / 'PixelClassifier_UserDeadLabels.zarr'
                if outpath.exists():
                    shutil.rmtree(outpath)
                save_as_zarr(self.viewer.layers[layer_name].data, outpath)
                self.log(f"Saved Dead labels to {outpath.name}")
                saved_any = True

        for cell_type in self.all_cell_types:
            layer_name = f'User Provided Labels ({cell_type.capitalize()})'
            if layer_name in self.viewer.layers:
                outpath = pixel_class_outdir / f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr'
                if outpath.exists():
                    shutil.rmtree(outpath)
                save_as_zarr(self.viewer.layers[layer_name].data, outpath)
                self.log(f"Saved {cell_type} labels to {outpath.name}")
                saved_any = True

        if not saved_any:
            self.log("No user label layers found in viewer to save.")
        else:
            self.log("All user labels saved.")

    def _on_process_all_toggled(self, checked):
        self.spin_t_start.setEnabled(not checked)
        self.spin_t_end.setEnabled(not checked)

    def _on_run_segmentation_clicked(self):
        self.log("Starting batch segmentation...")
        self._persist_params()
        try:
            # Get parameters
            overwrite = self.check_overwrite.isChecked()
            
            # Timepoint range
            if self.check_process_all.isChecked():
                timepoint_range = None
                self.log("Processing ALL timepoints.")
            else:
                t_start = self.spin_t_start.value()
                t_end = self.spin_t_end.value()
                if t_start > t_end:
                    self.log("Error: Start timepoint must be <= End timepoint.")
                    return
                timepoint_range = (t_start, t_end)
                self.log(f"Processing timepoints: {t_start} - {t_end}")
            
            # Get cell types and parameters
            targets = []
            if self.has_death:
                targets.append('Dead')
            targets.extend([t.capitalize() for t in self.all_cell_types])
            
            # Extract segmentation parameters for each cell type
            segmentation_params = {}
            for target in targets:
                cell_type_key = target.lower()
                params = {}
                if cell_type_key == 'dead':
                    # Use defaults for dead channel if not in UI
                    params['edt_threshold'] = 2.5
                    params['segment_size_min'] = 10
                    params['opening_nr_pixels'] = 0
                    params['fill_holes'] = True
                else:
                    if cell_type_key in self.param_widgets:
                        widgets = self.param_widgets[cell_type_key]
                        params['edt_threshold'] = widgets['edt'].value()
                        params['segment_size_min'] = widgets['min_size'].value()
                        params['opening_nr_pixels'] = widgets['opening'].value()
                        params['fill_holes'] = widgets['fill_holes'].isChecked()
                    else:
                        self.log(f"Warning: No parameters found for {target}. Using defaults.")
                        params['edt_threshold'] = 2.5
                        params['segment_size_min'] = 10
                        params['opening_nr_pixels'] = 0
                        params['fill_holes'] = True
                segmentation_params[target] = params

            # Prepare kwargs for segment_and_update
            kwargs = {}
            for target, params in segmentation_params.items():
                for key, value in params.items():
                    kwargs[f'{target.lower()}_{key}'] = value
                    
            # Collect classifier paths
            output_dir = Path(self.metadata_loader.output_dir)
            pixel_class_outdir = output_dir / "images" / "PixelClassification"
            
            clf_organoid_paths = {}
            clf_immune_paths = {}
            clf_other_paths = {}
            clf_death_path = None
            
            # Helper to check clf exists
            def get_clf_path(target):
                 p = pixel_class_outdir / f'PixelClassifier_{target.capitalize()}.joblib'
                 return p if p.exists() else None

            for ct in self.organoid_types:
                clf_organoid_paths[ct] = get_clf_path(ct)
            for ct in self.immune_types:
                clf_immune_paths[ct] = get_clf_path(ct)
            for ct in self.other_types:
                clf_other_paths[ct] = get_clf_path(ct)
            if self.has_death:
                clf_death_path = get_clf_path('Dead')

            # Call the segmentation function
            updated_metadata = run_pixel_classifier_segmentation(
                metadata=self.metadata_loader.metadata,
                output_dir=self.metadata_loader.output_dir,
                only_segment=False,  # False means Predict + Segment
                n_workers=self.spin_workers.value(),
                log_callback=self.log, # Pass the logger callback
                overwrite_existing=overwrite,
                timepoint_range=timepoint_range,
                clf_organoid_paths=clf_organoid_paths,
                clf_immune_paths=clf_immune_paths,
                clf_other_paths=clf_other_paths,
                clf_death_path=clf_death_path,
                **kwargs
            )
            # Update metadata in loader
            self.metadata_loader.metadata = updated_metadata
            self.log("Batch segmentation finished successfully!")

        except Exception as e:
            traceback.print_exc()
            self.log(f"Error during batch segmentation: {e}")


class CellposeWidget(QWidget):
    def __init__(self, viewer, metadata_loader):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._init_ui()
    
    def _init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)
        layout.addWidget(QLabel("Cellpose Configuration (Coming Soon)"))
        layout.addStretch()

class ImportWidget(QWidget):
    def __init__(self, viewer, metadata_loader):
        super().__init__()
        self.viewer = viewer
        self.metadata_loader = metadata_loader
        self._init_ui()
    
    def _init_ui(self):
        layout = QVBoxLayout()
        self.setLayout(layout)
        layout.addWidget(QLabel("Import pre-existing segmentation labels (Coming Soon)"))
        layout.addStretch()