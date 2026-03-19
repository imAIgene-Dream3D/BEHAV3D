"""
Custom APOC (Accelerated Pixel and Object Classification) widget for BEHAV3D.

Replicates the original napari-accelerated-pixel-and-object-classification GUI
using only the `apoc` Python backend, avoiding the plugin dependency that
conflicts with NumPy 2.0.

This file is self-contained: it handles image loading, the Qt training widget,
and APOC train/predict calls. The original pixel classifier in
napari_pixelclassifier.py is left completely untouched.
"""

import os
import gc
import time
import shutil
from pathlib import Path

import numpy as np
import napari
from magicgui.widgets import PushButton
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QSpinBox, QDoubleSpinBox, QLineEdit, QPushButton as QtPushButton,
    QCheckBox, QGroupBox, QPlainTextEdit, QApplication, QScrollArea,
    QSizePolicy,
)
from qtpy.QtCore import Qt

import dask.array as da
import zarr

from behav3d.io.images import load_image, load_zarr, save_as_zarr, append_to_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape

# ---------------------------------------------------------------------------
# Feature set presets (matching the original APOC plugin)
# ---------------------------------------------------------------------------
FEATURE_PRESETS = {
    "small_quick": "original gaussian_blur=1 sobel_of_gaussian_blur=1",
    "medium_quick": "original gaussian_blur=1 gaussian_blur=5 sobel_of_gaussian_blur=1 sobel_of_gaussian_blur=5",
    "large_quick": (
        "original gaussian_blur=1 gaussian_blur=5 gaussian_blur=25 "
        "sobel_of_gaussian_blur=1 sobel_of_gaussian_blur=5 sobel_of_gaussian_blur=25"
    ),
    "custom": "",
}


def _load_training_images(
    metadata, output_dir, examples_per_sample, organoid_types, immune_types, other_types
):
    """
    Load training images from metadata, reusing the same logic as the original
    pixel classifier.  Returns a list of (C, Z, Y, X) arrays — one per sample
    timepoint — plus the pixel_class_outdir path and cell type lists.
    """
    from behav3d.core.metadata import has_dead_channel

    has_death = has_dead_channel(metadata)
    all_cell_types = organoid_types + immune_types + other_types

    pixel_class_outdir = Path(output_dir, "images", "PixelClassification")
    pixel_class_outdir.mkdir(exist_ok=True, parents=True)

    image_outpath = Path(pixel_class_outdir, "PixelClassifier_Images.zarr")

    # If images were already saved from a previous session, reload them
    if image_outpath.exists():
        print("Loading cached training images from previous session...")
        cached = np.asarray(load_image(image_outpath))  # (C, T, Z, Y, X)
        # Split T back into a list of (C, Z, Y, X)
        all_images = [cached[:, t, ...] for t in range(cached.shape[1])]
        return all_images, pixel_class_outdir, has_death, all_cell_types

    # Otherwise load fresh from raw files
    all_images = []
    max_shape = None

    for _, sample in metadata.iterrows():
        sample_name = sample.get("sample_name", "unknown")
        raw_image_path = sample.get("raw_image_path", "")

        if not raw_image_path or not Path(raw_image_path).exists():
            print(f"⚠️ Skipping {sample_name}: raw image not found")
            continue

        img_dir = Path(output_dir, "images", sample_name)
        img_dir.mkdir(parents=True, exist_ok=True)

        raw_zarr = Path(img_dir, "RawImage.zarr")
        if raw_zarr.exists():
            img = load_image(raw_zarr)
        else:
            img = load_image(raw_image_path)

        img = np.asarray(img)
        n_timepoints = img.shape[0]

        # Select random timepoints
        np.random.seed(42)
        t_indices = sorted(
            np.random.choice(n_timepoints, min(examples_per_sample, n_timepoints), replace=False)
        )
        print(f"  {sample_name}: selected timepoints {t_indices}")

        for t_idx in t_indices:
            frame = img[t_idx]  # (C, Z, Y, X)
            all_images.append(frame)
            frame_shape = frame.shape
            if max_shape is None:
                max_shape = list(frame_shape)
            else:
                for i in range(len(max_shape)):
                    max_shape[i] = max(max_shape[i], frame_shape[i])

    if max_shape is not None:
        all_images = [zeropad_image_to_match_shape(img, max_shape) for img in all_images]

    return all_images, pixel_class_outdir, has_death, all_cell_types


class APOCTrainingWidget(QWidget):
    """
    A Qt widget that replicates the original APOC Napari plugin interface.
    Docks directly into the Napari viewer.
    """

    def __init__(self, viewer, pixel_class_outdir, all_cell_types, has_death, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.pixel_class_outdir = pixel_class_outdir
        self.all_cell_types = all_cell_types
        self.has_death = has_death

        self._build_ui()
        self._connect_signals()
        # Populate layer dropdowns
        self._refresh_layers()
        # Listen for layer changes
        self.viewer.layers.events.inserted.connect(lambda _: self._refresh_layers())
        self.viewer.layers.events.removed.connect(lambda _: self._refresh_layers())

    def _build_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # Title
        title = QLabel("<h3>APOC Pixel Classification</h3>")
        layout.addWidget(title)

        # Labeling guidance
        info = QLabel("Draw labels: <b>1</b> = background, <b>2</b> = foreground")
        info.setStyleSheet("color: #888; font-style: italic; padding: 2px 0;")
        layout.addWidget(info)

        # --- Image layers section ---
        img_group = QGroupBox("Image Layers (select inputs)")
        img_layout = QVBoxLayout()
        img_layout.setSpacing(2)
        self.image_checkboxes = []
        # Will be populated dynamically when layers change
        self.img_checkbox_container = QWidget()
        self.img_checkbox_layout = QVBoxLayout()
        self.img_checkbox_layout.setContentsMargins(0, 0, 0, 0)
        self.img_checkbox_container.setLayout(self.img_checkbox_layout)
        img_layout.addWidget(self.img_checkbox_container)
        img_group.setLayout(img_layout)
        layout.addWidget(img_group)

        # --- Feature set ---
        feat_layout = QHBoxLayout()
        feat_layout.addWidget(QLabel("Feature set:"))
        self.feature_combo = QComboBox()
        self.feature_combo.addItems(list(FEATURE_PRESETS.keys()))
        self.feature_combo.setCurrentText("medium_quick")
        feat_layout.addWidget(self.feature_combo)
        layout.addLayout(feat_layout)

        # Custom features text
        self.custom_features_label = QLabel("Custom features:")
        self.custom_features_edit = QLineEdit()
        self.custom_features_edit.setPlaceholderText(
            "e.g. original gaussian_blur=1 sobel_of_gaussian_blur=1"
        )
        self.custom_features_edit.setText(FEATURE_PRESETS["medium_quick"])
        layout.addWidget(self.custom_features_label)
        layout.addWidget(self.custom_features_edit)
        # Initially hide custom field
        self._toggle_custom_features()

        # --- RF parameters ---
        rf_layout = QHBoxLayout()
        rf_layout.addWidget(QLabel("Max depth:"))
        self.max_depth_spin = QSpinBox()
        self.max_depth_spin.setRange(1, 20)
        self.max_depth_spin.setValue(2)
        rf_layout.addWidget(self.max_depth_spin)

        rf_layout.addWidget(QLabel("Trees:"))
        self.num_ensembles_spin = QSpinBox()
        self.num_ensembles_spin.setRange(10, 1000)
        self.num_ensembles_spin.setSingleStep(10)
        self.num_ensembles_spin.setValue(100)
        rf_layout.addWidget(self.num_ensembles_spin)
        layout.addLayout(rf_layout)

        # --- Buttons ---
        btn_layout = QHBoxLayout()
        self.run_btn = QtPushButton("Train & Apply APOC Classifiers")
        self.run_btn.setStyleSheet("background-color: #2196F3; color: white; font-weight: bold; padding: 8px 16px;")
        btn_layout.addWidget(self.run_btn)
        layout.addLayout(btn_layout)

        # --- Status ---
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        layout.addStretch()
        self.setLayout(layout)

    def _connect_signals(self):
        self.feature_combo.currentTextChanged.connect(self._toggle_custom_features)
        self.run_btn.clicked.connect(self._on_run_all)

    def _toggle_custom_features(self, *_):
        is_custom = self.feature_combo.currentText() == "custom"
        self.custom_features_label.setVisible(is_custom)
        self.custom_features_edit.setVisible(is_custom)

    def _refresh_layers(self):
        """Rebuild image checkboxes from viewer layers."""
        # Clear old checkboxes
        for cb in self.image_checkboxes:
            self.img_checkbox_layout.removeWidget(cb)
            cb.deleteLater()
        self.image_checkboxes = []

        # Add checkboxes for Image layers
        for layer in self.viewer.layers:
            if isinstance(layer, napari.layers.Image) and not layer.name.startswith("Pixel Classification"):
                cb = QCheckBox(layer.name)
                cb.setChecked(True)
                self.img_checkbox_layout.addWidget(cb)
                self.image_checkboxes.append(cb)

    def _get_feature_string(self):
        preset = self.feature_combo.currentText()
        if preset == "custom":
            return self.custom_features_edit.text().strip()
        return FEATURE_PRESETS.get(preset, FEATURE_PRESETS["medium_quick"])

    def _get_selected_images(self):
        """Return a list of numpy arrays for checked image layers."""
        images = []
        for cb in self.image_checkboxes:
            if cb.isChecked():
                try:
                    layer = self.viewer.layers[cb.text()]
                    images.append(np.asarray(layer.data))
                except KeyError:
                    pass
        return images

    def _on_run_all(self):
        import apoc

        images = self._get_selected_images()
        if not images:
            self.status_label.setText("⚠️ No image layers selected!")
            return

        feature_string = self._get_feature_string()
        max_depth = self.max_depth_spin.value()
        num_ensembles = self.num_ensembles_spin.value()
        
        self.status_label.setText("Starting training/prediction run...")
        QApplication.processEvents()

        # Build list of tasks
        tasks = []
        for ct in self.all_cell_types:
            layer_name = f"User Provided Labels ({ct.capitalize()})"
            clf_name = f"PixelClassifier_{ct.capitalize()}.cl"
            pred_layer_name = f"Pixel Classification ({ct.capitalize()})"
            tasks.append((layer_name, clf_name, pred_layer_name))
            
        if self.has_death:
            tasks.append(("User Provided Labels (Dead)", "PixelClassifier_Death.cl", "Pixel Classification (Dead)"))

        successes = []
        try:
            for layer_name, clf_name, pred_layer_name in tasks:
                self.status_label.setText(f"Processing {clf_name}...")
                QApplication.processEvents()
                
                # 1. Fetch annotation
                try:
                    annotation = np.asarray(self.viewer.layers[layer_name].data)
                except KeyError:
                    # Layer deleted or not found
                    continue
                    
                # 2. Check if there are any labels drawn
                if not np.any(annotation):
                    print(f"Skipping {clf_name} (no user labels found)")
                    continue
                    
                clf_path = str(Path(self.pixel_class_outdir, clf_name))
                
                # 3. Train
                n_timepoints = annotation.shape[0] if annotation.ndim == 4 else 1
                has_trained = False
                
                # Erase any preexisting classifier to guarantee we train from scratch as per APOC documentation
                if Path(clf_path).exists():
                    apoc.erase_classifier(clf_path)
                
                # Use ObjectSegmenter for native APOC instance segmentation
                clf = apoc.ObjectSegmenter(
                    opencl_filename=clf_path, 
                    max_depth=max_depth, 
                    num_ensembles=num_ensembles,
                    positive_class_identifier=2
                )
                
                if annotation.ndim == 4:
                    for t in range(n_timepoints):
                        ann_t = annotation[t]
                        if not np.any(ann_t): continue
                        imgs_t = [img[t] for img in images]
                        if len(imgs_t) == 1: imgs_t = imgs_t[0]
                        clf.train(feature_string, ann_t, imgs_t, continue_training=has_trained)
                        has_trained = True
                else:
                    imgs = images[0] if len(images) == 1 else images
                    clf.train(feature_string, annotation, imgs)
                    has_trained = True
                    
                if not has_trained:
                    continue
                    
                # 4. Predict / Apply
                ndim = images[0].ndim
                if ndim == 4:
                    results = []
                    for t in range(n_timepoints):
                        imgs_t = [img[t] for img in images]
                        if len(imgs_t) == 1: imgs_t = imgs_t[0]
                        res_t = clf.predict(image=imgs_t)
                        results.append(np.asarray(res_t).astype(np.int16))
                    result = np.stack(results, axis=0)
                else:
                    imgs = images[0] if len(images) == 1 else images
                    result = clf.predict(image=imgs)
                    result = np.asarray(result).astype(np.int16)
                    
                # 5. Add or update result layer (Predicted instance labels)
                # We name the layer Segments directly since ObjectSegmenter runs connected components
                layer_name_final = pred_layer_name if "Death" in clf_name else f"{clf_name.replace('PixelClassifier_', '').replace('.cl', '')} Segments"
                
                if layer_name_final in [l.name for l in self.viewer.layers]:
                    self.viewer.layers[layer_name_final].data = result
                else:
                    self.viewer.add_labels(result, name=layer_name_final, opacity=0.8)
                    
                # Save to disk to match original behavior
                from behav3d.io.images import save_as_zarr
                
                # Save as PredictedLabels.zarr so the rest of the BEHAV3D pipeline can find it
                arr_path = Path(self.pixel_class_outdir, f"{Path(clf_path).stem}_PredictedLabels.zarr")
                if arr_path.exists():
                    import shutil
                    shutil.rmtree(arr_path)
                save_as_zarr(result, arr_path)
                
                # Also save the raw segments if not death
                if "Death" not in clf_name:
                    seg_path = Path(self.pixel_class_outdir, f"segmentation_{Path(clf_path).stem.replace('PixelClassifier_','')}.zarr")
                    if seg_path.exists():
                        import shutil
                        shutil.rmtree(seg_path)
                    save_as_zarr(result, seg_path)
                
                successes.append(clf_name)

            if successes:
                self.status_label.setText(f"✅ Successfully trained & applied:\n{', '.join(successes)}")
            else:
                self.status_label.setText("⚠️ Finished, but found NO drawn labels in any layer!")

        except Exception as e:
            self.status_label.setText(f"❌ Failed: {e}")
            import traceback
            traceback.print_exc()


def train_pixel_classifier_apoc(
    output_dir,
    metadata=None,
    examples_per_sample=3,
    n_workers=None,
    organoid_types=None,
    immune_types=None,
    other_types=None,
    initial_params=None,
    on_params_changed=None,
):
    """
    APOC version of the pixel classifier training.
    Opens a Napari viewer with the custom APOC widget docked on the right.
    The original pixel classifier code is NOT called.
    """
    from behav3d.core.metadata import has_dead_channel
    from functools import partial

    organoid_types = organoid_types or []
    immune_types = immune_types or []
    other_types = other_types or []
    all_cell_types = organoid_types + immune_types + other_types
    has_death = has_dead_channel(metadata) if metadata is not None else False

    # --- Load images ---
    print("Loading training images for APOC...")
    image_list, pixel_class_outdir, has_death, all_cell_types = _load_training_images(
        metadata=metadata,
        output_dir=output_dir,
        examples_per_sample=examples_per_sample,
        organoid_types=organoid_types,
        immune_types=immune_types,
        other_types=other_types,
    )

    if not image_list:
        print("⚠️ No training images found!")
        return None

    # --- Stack along Time ---
    max_shape = list(image_list[0].shape)
    for img in image_list[1:]:
        for i in range(len(max_shape)):
            max_shape[i] = max(max_shape[i], img.shape[i])
    image_list = [zeropad_image_to_match_shape(img, max_shape) for img in image_list]

    stacked = np.stack(image_list, axis=0)  # Shape depends on disk cache order
    
    # Heuristic constraint: Ensure shape is (Channels, Timepoints, Z, Y, X)
    # Usually Timepoints > Channels. If it's (Timepoints, Channels, Z, Y, X), flip the first two.
    if stacked.shape[0] > stacked.shape[1]:
        print(f"Flipping dimensions {stacked.shape[:2]} to match (C, T, Z, Y, X)")
        stacked = stacked.transpose(1, 0, 2, 3, 4)
        
    print(f"Final normalized shape {stacked.shape} (C, T, Z, Y, X)")

    # --- Create Napari viewer ---
    viewer = napari.Viewer()

    # Add each channel as a separate image layer
    n_channels = stacked.shape[0]
    channel_colors = [
        'cyan', 'yellow', 'red', 'green', 'magenta', 'blue',
        'gray', 'turbo', 'viridis', 'plasma', 'inferno', 'twilight'
    ]
    for ch in range(n_channels):
        channel_data = stacked[ch, ...] # (T, Z, Y, X)
        nonzero = channel_data[channel_data > 0]
        if nonzero.size > 0:
            clim = (0, float(np.percentile(nonzero, 99.8)))
        else:
            clim = (0, 1e-3)

        viewer.add_image(
            channel_data,
            name=f"Channel {ch}",
            contrast_limits=clim,
            colormap=channel_colors[ch % len(channel_colors)],
            blending='additive',
            opacity=0.8,
        )

    # Add annotation labels layers per cell type
    spatial_shape = stacked.shape[2:]  # (Z, Y, X)
    time_shape = stacked.shape[1]      # T
    full_shape = (time_shape,) + spatial_shape # (T, Z, Y, X)
    
    for cell_type in all_cell_types:
        labels_outpath = Path(pixel_class_outdir, f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr')
        if labels_outpath.exists():
            print(f"Loading existing labels for {cell_type}")
            existing = np.asarray(load_zarr(labels_outpath))
            
            # If stored as mosaiced, we can't easily map it perfectly if we don't know N_samples
            # So let's just initialize an empty array of the correct shape if it mismatches
            if existing.shape != full_shape:
                print(f"⚠️ Existing labels shape {existing.shape} doesn't match {full_shape}. Initializing empty.")
                user_labels = np.zeros(full_shape, dtype=np.int16)
            else:
                user_labels = existing
        else:
            user_labels = np.zeros(full_shape, dtype=np.int16)

        viewer.add_labels(
            user_labels,
            name=f"User Provided Labels ({cell_type.capitalize()})",
            opacity=0.5,
        )

        # Create segment layers to mimic original widget
        seg_layer_name = f"{cell_type.capitalize()} Segments"
        seg_labels = np.zeros(full_shape, dtype=np.int16)
        viewer.add_labels(
            seg_labels,
            name=seg_layer_name,
            opacity=0.8,
            visible=False
        )

    if has_death:
        dead_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
        if dead_outpath.exists():
            dead_labels = np.asarray(load_zarr(dead_outpath))
            if dead_labels.shape != full_shape:
                dead_labels = np.zeros(full_shape, dtype=np.int16)
        else:
            dead_labels = np.zeros(full_shape, dtype=np.int16)
        viewer.add_labels(dead_labels, name="User Provided Labels (Dead)", opacity=0.5)

    # --- Dock the APOC widget ---
    apoc_widget = APOCTrainingWidget(
        viewer=viewer,
        pixel_class_outdir=pixel_class_outdir,
        all_cell_types=all_cell_types,
        has_death=has_death,
    )
    viewer.window.add_dock_widget(apoc_widget, area="right", name="APOC Pixel Classification")

    # --- Save labels button ---
    def save_labels(log=print):
        """Save user labels to Zarr."""
        if has_death:
            dead_layer = viewer.layers['User Provided Labels (Dead)']
            dead_outpath = Path(pixel_class_outdir, 'PixelClassifier_UserDeadLabels.zarr')
            if dead_outpath.exists():
                shutil.rmtree(dead_outpath)
            save_as_zarr(dead_layer.data, dead_outpath)
            log(f"Saved death labels to: {dead_outpath}")

        for cell_type in all_cell_types:
            layer_name = f'User Provided Labels ({cell_type.capitalize()})'
            label_layer = viewer.layers[layer_name]
            labels_outpath = Path(pixel_class_outdir, f'PixelClassifier_User{cell_type.capitalize()}Labels.zarr')
            if labels_outpath.exists():
                shutil.rmtree(labels_outpath)
            save_as_zarr(label_layer.data, labels_outpath)
            log(f"Saved {cell_type} labels to: {labels_outpath}")

        log("All user labels saved!")

    # Add log widget + save button
    log_output = QPlainTextEdit()
    log_output.setReadOnly(True)
    log_widget = QWidget()
    log_layout = QVBoxLayout()
    log_layout.addWidget(log_output)
    log_widget.setLayout(log_layout)
    viewer.window.add_dock_widget(log_widget, area="right", name="Log Output")

    save_button = QtPushButton("Save User Labels")
    save_button.setStyleSheet("background-color: #FF9800; color: white; font-weight: bold; padding: 6px;")
    save_button.clicked.connect(lambda: save_labels(log=log_output.appendPlainText))
    apoc_widget.layout().addWidget(save_button)

    return viewer
