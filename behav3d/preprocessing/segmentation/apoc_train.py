"""
Custom APOC (Accelerated Pixel and Object Classification) widget for BEHAV3D.

This file is self-contained: it handles image loading, the Qt training widget
(with per-cell-type tabs), and APOC train/predict calls.
The original pixel classifier in napari_pixelclassifier.py is left completely untouched.
"""

import os
import gc
import json
import time
import shutil
from pathlib import Path

import numpy as np
import napari
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QSpinBox, QLineEdit, QPushButton as QtPushButton,
    QCheckBox, QGroupBox, QPlainTextEdit, QApplication, QScrollArea,
    QSizePolicy, QTabWidget, QFrame, QMessageBox,
)
from qtpy.QtCore import Qt

import dask.array as da
import zarr
import pyclesperanto_prototype as cle

from behav3d.io.images import load_image, load_zarr, save_as_zarr, append_to_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape

# ---------------------------------------------------------------------------
# Feature set presets (matching APOC convention)
# ---------------------------------------------------------------------------

FEATURE_PRESETS = {
    "small_quick": {
        "features": ["gaussian_blur", "sobel_of_gaussian_blur"],
        "default_sigmas": "1",
    },
    "medium_quick": {
        "features": ["original", "gaussian_blur", "sobel_of_gaussian_blur"],
        "default_sigmas": "1, 5",
    },
    "large_quick": {
        "features": ["original", "gaussian_blur", "sobel_of_gaussian_blur",
                      "laplace_box_of_gaussian_blur"],
        "default_sigmas": "1, 5, 25",
    },
    "custom": {
        "features": [],
        "default_sigmas": "",
    },
}


def _build_feature_string(preset_name, sigmas_str):
    """Build a full APOC feature string from a preset name + comma-separated sigmas.
    For 'custom' preset, sigmas_str is used as the raw feature string directly.
    """
    if preset_name == "custom":
        return sigmas_str.strip()

    preset = FEATURE_PRESETS.get(preset_name, FEATURE_PRESETS["medium_quick"])
    features = preset["features"]

    try:
        sigmas = [s.strip() for s in sigmas_str.split(",") if s.strip()]
    except Exception:
        sigmas = ["1"]

    if not sigmas:
        sigmas = ["1"]

    parts = []
    for feat in features:
        if feat == "original":
            parts.append("original")
        else:
            for s in sigmas:
                parts.append(f"{feat}={s}")
    return " ".join(parts)


# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Image loading (same logic as CPU pipeline for consistency)
# ---------------------------------------------------------------------------

def _load_training_images(
    metadata, output_dir, examples_per_sample, organoid_types, immune_types, other_types, overwrite_images=False
):
    """
    Load training images from metadata.
    Returns a list of (C, Z, Y, X) arrays — one per selected timepoint — plus
    the pixel_class_outdir path, has_death flag, and all cell types list.
    """
    from behav3d.core.metadata import has_dead_channel

    has_death = has_dead_channel(metadata)
    all_cell_types = organoid_types + immune_types + other_types

    pixel_class_outdir = Path(output_dir, "images", "PixelClassification")
    pixel_class_outdir.mkdir(exist_ok=True, parents=True)

    image_outpath = Path(pixel_class_outdir, "PixelClassifier_Images.zarr")

    # --- Check cached dataset ---
    if image_outpath.exists():
        if overwrite_images:
            print("Overwrite sample images requested — deleting cached data...")
            shutil.rmtree(image_outpath)
        else:
            print("Loading cached training images from previous session...")
            # We assume the cache matches the metadata order of the experimental setup
            # Find the channel axis based on metadata of the first sample
            first_sample = metadata.iloc[0]
            axis_order = first_sample.get("dimension_order", "TCZYX")
            if not isinstance(axis_order, str) or not axis_order:
                axis_order = "TCZYX"
            
            remaining_axes = axis_order.replace('T', '')
            c_axis_in_stacked = 1 + remaining_axes.find('C')
            
            cached = load_image(image_outpath) 
            n_channels = cached.shape[c_axis_in_stacked]
            all_images = [np.take(cached, ch, axis=c_axis_in_stacked) for ch in range(n_channels)]
            all_images = [np.asarray(img) for img in all_images] # Load into mem for viewer
            return all_images, pixel_class_outdir, has_death, all_cell_types

    # --- Load fresh from raw files ---
    all_images = []
    max_shape = None
    n_samples = 0

    for _, sample in metadata.iterrows():
        sample_name = sample.get("sample_name", "unknown")
        raw_image_path = sample.get("raw_image_path", "")

        if not raw_image_path or not Path(raw_image_path).exists():
            print(f"⚠️ Skipping {sample_name}: raw image not found")
            continue

        n_samples += 1
        img_dir = Path(output_dir, "images", sample_name)
        img_dir.mkdir(parents=True, exist_ok=True)

        # Get dimension order from metadata if available
        axis_order = sample.get("dimension_order", "TCZYX")
        if not isinstance(axis_order, str) or not axis_order:
            axis_order = "TCZYX"

        # Use raw_image_path directly as it is the source of truth in metadata
        # (and is updated to the .zarr path if conversion was run)
        img = load_image(raw_image_path, axis_order=axis_order)

        # Find which axis is Time in this specific image
        t_axis = axis_order.find('T')
        if t_axis == -1:
            t_axis = 0 # Fallback
            
        n_timepoints = img.shape[t_axis]
        n_to_select = min(examples_per_sample, n_timepoints)

        # Select equally spaced timepoints (first, last, and middle)
        if n_to_select <= 1:
            t_indices = [0]
        else:
            t_indices = np.round(np.linspace(0, n_timepoints - 1, n_to_select)).astype(int)
            t_indices = sorted(list(set(t_indices)))
            
        print(f"  {sample_name}: selected {len(t_indices)} equidistant timepoints {t_indices}")

        for t_idx in t_indices:
            # Fetch only the specific timepoint slice into memory
            frame = np.asarray(np.take(img, t_idx, axis=t_axis)) 
            all_images.append(frame)
            frame_shape = frame.shape
            if max_shape is None:
                max_shape = list(frame_shape)
            else:
                for i in range(len(max_shape)):
                    max_shape[i] = max(max_shape[i], frame_shape[i])

    if max_shape is not None:
        all_images = [zeropad_image_to_match_shape(img, max_shape) for img in all_images]

    # Stash on disk for next time
    import dask.array as da
    import gc
    all_images_stack = da.stack(all_images)
    save_as_zarr(all_images_stack, image_outpath)
    del all_images_stack
    gc.collect()

    # Get dimension info from the last processed sample's metadata
    # (assuming all samples in the same training set share the same dimension order)
    if 'axis_order' not in locals():
        # Fallback if the loop didn't run for some reason
        axis_order = "TCZYX"
        
    remaining_axes = axis_order.replace('T', '')
    c_axis_in_stacked = 1 + remaining_axes.find('C')

    all_images_cached = load_image(image_outpath)
    n_channels = all_images_cached.shape[c_axis_in_stacked]
    all_images = [np.take(all_images_cached, ch, axis=c_axis_in_stacked) for ch in range(n_channels)]
    all_images = [np.asarray(img) for img in all_images]

    return all_images, pixel_class_outdir, has_death, all_cell_types


# ---------------------------------------------------------------------------
# Per-cell-type tab panel
# ---------------------------------------------------------------------------

class CellTypeTab(QWidget):
    """
    Widget for a single cell type tab.
    Contains: channel selection, feature preset, sigma values, read-only
    feature preview, max_depth, num_ensembles.
    """

    def __init__(self, cell_type, viewer, initial_params=None, parent=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.viewer = viewer
        ip = initial_params or {}

        layout = QVBoxLayout()
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(5)

        # --- Channel selection ---
        chan_group = QGroupBox("Image Channel Inputs")
        chan_layout = QVBoxLayout()
        chan_layout.setSpacing(2)
        self.channel_checkboxes = []
        self.chan_checkbox_container = QWidget()
        self.chan_checkbox_layout = QVBoxLayout()
        self.chan_checkbox_layout.setContentsMargins(0, 0, 0, 0)
        self.chan_checkbox_container.setLayout(self.chan_checkbox_layout)
        chan_layout.addWidget(self.chan_checkbox_container)
        chan_group.setLayout(chan_layout)
        layout.addWidget(chan_group)

        # Labeling hint
        hint = QLabel("Labels: <b>1</b> = background&nbsp;&nbsp; <b>2</b> = foreground")
        hint.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(hint)

        # --- Feature preset ---
        feat_row = QHBoxLayout()
        feat_row.addWidget(QLabel("Preset:"))
        self.feature_combo = QComboBox()
        self.feature_combo.addItems(list(FEATURE_PRESETS.keys()))
        saved_preset = ip.get(f"apoc_{cell_type}_feature_preset", "medium_quick")
        self.feature_combo.setCurrentText(saved_preset if saved_preset in FEATURE_PRESETS else "medium_quick")
        feat_row.addWidget(self.feature_combo)
        layout.addLayout(feat_row)

        # --- Sigma values (comma-separated) ---
        sigma_row = QHBoxLayout()
        self.sigma_label = QLabel("Sigma values:")
        sigma_row.addWidget(self.sigma_label)
        self.sigma_edit = QLineEdit()
        default_sigmas = FEATURE_PRESETS.get(saved_preset, FEATURE_PRESETS["medium_quick"])["default_sigmas"]
        saved_sigmas = ip.get(f"apoc_{cell_type}_sigmas", default_sigmas)
        self.sigma_edit.setText(str(saved_sigmas))
        self.sigma_edit.setPlaceholderText("e.g. 1, 5, 25")
        self.sigma_edit.setToolTip("Comma-separated sigma (scale) values for Gaussian filters.\nSmall σ (1-2) = fine detail; Large σ (5-25) = coarse/shape.")
        sigma_row.addWidget(self.sigma_edit)
        layout.addLayout(sigma_row)

        # --- Custom feature string (only visible when preset = custom) ---
        self.custom_feat_label = QLabel("Custom feature string:")
        self.custom_feat_edit = QLineEdit()
        saved_custom = ip.get(f"apoc_{cell_type}_custom_feature_string", "")
        self.custom_feat_edit.setText(saved_custom)
        self.custom_feat_edit.setPlaceholderText("e.g. original gaussian_blur=1 sobel_of_gaussian_blur=5")
        layout.addWidget(self.custom_feat_label)
        layout.addWidget(self.custom_feat_edit)

        # --- Read-only feature preview ---
        self.preview_label = QLabel("")
        self.preview_label.setWordWrap(True)
        self.preview_label.setStyleSheet("color: #888; font-size: 11px; padding: 2px 4px; background: #f5f5f5; border-radius: 3px;")
        layout.addWidget(self.preview_label)

        # Connect signals to update preview
        self.feature_combo.currentTextChanged.connect(self._on_preset_changed)
        self.sigma_edit.textChanged.connect(self._update_preview)
        self.custom_feat_edit.textChanged.connect(self._update_preview)

        # Initial visibility / preview
        self._on_preset_changed(self.feature_combo.currentText())

        # --- RF parameters ---
        rf_row = QHBoxLayout()
        rf_row.addWidget(QLabel("Max depth:"))
        self.max_depth_spin = QSpinBox()
        self.max_depth_spin.setRange(1, 20)
        self.max_depth_spin.setValue(int(ip.get(f"apoc_{cell_type}_max_depth", 2)))
        rf_row.addWidget(self.max_depth_spin)

        rf_row.addWidget(QLabel("Trees:"))
        self.num_ensembles_spin = QSpinBox()
        self.num_ensembles_spin.setRange(10, 1000)
        self.num_ensembles_spin.setSingleStep(10)
        self.num_ensembles_spin.setValue(int(ip.get(f"apoc_{cell_type}_num_ensembles", 100)))
        rf_row.addWidget(self.num_ensembles_spin)
        layout.addLayout(rf_row)

        layout.addStretch()
        self.setLayout(layout)

    def _on_preset_changed(self, preset_name):
        """Toggle visibility of sigma vs custom field, update sigma defaults."""
        is_custom = preset_name == "custom"
        self.sigma_edit.setVisible(not is_custom)
        self.sigma_label.setVisible(not is_custom)
        self.custom_feat_label.setVisible(is_custom)
        self.custom_feat_edit.setVisible(is_custom)

        if not is_custom:
            # Auto-fill default sigmas for the new preset
            preset = FEATURE_PRESETS.get(preset_name, FEATURE_PRESETS["medium_quick"])
            self.sigma_edit.setText(preset["default_sigmas"])

        self._update_preview()

    def _update_preview(self):
        """Rebuild the feature string preview from current preset + sigmas."""
        feat_str = self.get_feature_string()
        self.preview_label.setText(f"<b>Features:</b> {feat_str}")

    def get_feature_string(self):
        """Return the full APOC feature string based on current widget values."""
        preset = self.feature_combo.currentText()
        if preset == "custom":
            return self.custom_feat_edit.text().strip()
        return _build_feature_string(preset, self.sigma_edit.text())

    def refresh_channel_checkboxes(self):
        """Rebuild channel checkboxes from current Napari image layers."""
        for cb in self.channel_checkboxes:
            self.chan_checkbox_layout.removeWidget(cb)
            cb.deleteLater()
        self.channel_checkboxes = []

        for layer in self.viewer.layers:
            if isinstance(layer, napari.layers.Image) and not layer.name.startswith("Pixel Classification"):
                cb = QCheckBox(layer.name)
                cb.setChecked(True)
                self.chan_checkbox_layout.addWidget(cb)
                self.channel_checkboxes.append(cb)

    def get_config(self):
        """Return a dict with all current widget values."""
        return {
            "feature_preset":          self.feature_combo.currentText(),
            "sigmas":                  self.sigma_edit.text().strip(),
            "custom_feature_string":   self.custom_feat_edit.text().strip(),
            "feature_string":          self.get_feature_string(),  # computed, for training
            "max_depth":               self.max_depth_spin.value(),
            "num_ensembles":           self.num_ensembles_spin.value(),
            "channels":                [cb.text() for cb in self.channel_checkboxes if cb.isChecked()],
        }

    def apply_config(self, cfg):
        """Restore widget values from a config dict."""
        if "feature_preset" in cfg:
            self.feature_combo.setCurrentText(cfg["feature_preset"])
        if "sigmas" in cfg:
            self.sigma_edit.setText(cfg["sigmas"])
        if "custom_feature_string" in cfg:
            self.custom_feat_edit.setText(cfg["custom_feature_string"])
        if "max_depth" in cfg:
            self.max_depth_spin.setValue(int(cfg["max_depth"]))
        if "num_ensembles" in cfg:
            self.num_ensembles_spin.setValue(int(cfg["num_ensembles"]))
        if "channels" in cfg:
            for cb in self.channel_checkboxes:
                cb.setChecked(cb.text() in cfg["channels"])


# ---------------------------------------------------------------------------
# Main APOC training widget
# ---------------------------------------------------------------------------

class APOCTrainingWidget(QWidget):
    """
    Tabbed APOC training widget docked inside the Napari viewer.
    One tab per cell type (+ Death if applicable), each with:
      - channel selection
      - feature preset + sigma values + auto-preview
      - max_depth / num_ensembles
    Global controls: Apply-to-All, Train Current, Train-All, Save Labels.
    """

    def __init__(
        self,
        viewer,
        pixel_class_outdir,
        all_cell_types,
        has_death,
        initial_params=None,
        on_params_changed=None,
        parent=None,
    ):
        super().__init__(parent)
        self.viewer = viewer
        self.pixel_class_outdir = pixel_class_outdir
        self.all_cell_types = all_cell_types
        self.has_death = has_death
        self._initial_params = initial_params or {}
        self._on_params_changed = on_params_changed

        # All tab labels = cell types + optional Death
        self._tab_cell_types = list(all_cell_types)
        if has_death:
            self._tab_cell_types.append("dead")

        self._build_ui()
        self._connect_signals()
        self._refresh_all_channels()

        # Listen for layer changes
        self.viewer.layers.events.inserted.connect(lambda _: self._refresh_all_channels())
        self.viewer.layers.events.removed.connect(lambda _: self._refresh_all_channels())

    def _build_ui(self):
        layout = QVBoxLayout()
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)

        title = QLabel("<h3>APOC Pixel Classification</h3>")
        layout.addWidget(title)

        # --- Tab widget ---
        self.tab_widget = QTabWidget()
        self.tabs = {}  # cell_type -> CellTypeTab

        for ct in self._tab_cell_types:
            tab = CellTypeTab(ct, self.viewer, initial_params=self._initial_params)
            self.tabs[ct] = tab
            self.tab_widget.addTab(tab, ct.capitalize())

        layout.addWidget(self.tab_widget)

        # --- Separator ---
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        layout.addWidget(line)

        # --- Global buttons ---
        global_row1 = QHBoxLayout()

        self.apply_all_btn = QtPushButton("⬇ Apply config to all tabs")
        self.apply_all_btn.setStyleSheet("padding: 5px 10px;")
        self.apply_all_btn.setToolTip("Copy the current tab's preset, sigmas, depth and trees to ALL other tabs")
        global_row1.addWidget(self.apply_all_btn)
        layout.addLayout(global_row1)

        global_row2 = QHBoxLayout()

        self.train_current_btn = QtPushButton("▶ Train current tab")
        self.train_current_btn.setStyleSheet(
            "background-color: #1976D2; color: white; font-weight: bold; padding: 6px 12px;"
        )
        global_row2.addWidget(self.train_current_btn)

        self.train_all_btn = QtPushButton("▶▶ Train ALL classifiers")
        self.train_all_btn.setStyleSheet(
            "background-color: #2e7d32; color: white; font-weight: bold; padding: 6px 12px;"
        )
        global_row2.addWidget(self.train_all_btn)
        layout.addLayout(global_row2)

        # --- Status label ---
        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)

        layout.addStretch()
        self.setLayout(layout)

    def _connect_signals(self):
        self.apply_all_btn.clicked.connect(self._on_apply_to_all)
        self.train_current_btn.clicked.connect(self._on_train_current)
        self.train_all_btn.clicked.connect(self._on_train_all)

    def _refresh_all_channels(self):
        """Refresh channel checkboxes in all tabs."""
        for tab in self.tabs.values():
            tab.refresh_channel_checkboxes()

    # ------------------------------------------------------------------
    # Button handlers
    # ------------------------------------------------------------------

    def _on_apply_to_all(self):
        """Copy current tab's config to all other tabs."""
        current_ct = self._tab_cell_types[self.tab_widget.currentIndex()]
        cfg = self.tabs[current_ct].get_config()
        for ct, tab in self.tabs.items():
            if ct != current_ct:
                tab.apply_config(cfg)
        self.status_label.setText(f"↪ Config applied from '{current_ct}' to all other tabs.")

    def _on_train_current(self):
        """Train only the classifier for the currently visible tab."""
        current_ct = self._tab_cell_types[self.tab_widget.currentIndex()]
        self._run_training([current_ct])

    def _on_train_all(self):
        """Train classifiers for ALL cell types using each tab's individual config."""
        self._run_training(self._tab_cell_types)

    # ------------------------------------------------------------------
    # Core training logic
    # ------------------------------------------------------------------

    def _get_images_for_tab(self, ct):
        """Return numpy arrays for the image layers checked in a tab."""
        tab = self.tabs[ct]
        images = []
        for cb in tab.channel_checkboxes:
            if cb.isChecked():
                try:
                    layer = self.viewer.layers[cb.text()]
                    # Keep data lazy (dask/zarr) to save memory
                    images.append(layer.data)
                except KeyError:
                    pass
        return images

    def save_user_labels(self, log=print):
        """Save all user-provided labels for all cell types and Dead (if present)."""
        import shutil
        for cell_type in self.all_cell_types:
            lname = f"User Provided Labels ({cell_type.capitalize()})"
            if lname not in [l.name for l in self.viewer.layers]:
                continue
            label_layer = self.viewer.layers[lname]
            outpath = Path(self.pixel_class_outdir, f"PixelClassifier_User{cell_type.capitalize()}Labels.zarr")
            if outpath.exists():
                shutil.rmtree(outpath)
            save_as_zarr(label_layer.data, outpath)
            log(f"Saved {cell_type} labels → {outpath}")

        if self.has_death and "User Provided Labels (Dead)" in [l.name for l in self.viewer.layers]:
            dead_layer = self.viewer.layers["User Provided Labels (Dead)"]
            dead_outpath = Path(self.pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
            if dead_outpath.exists():
                shutil.rmtree(dead_outpath)
            save_as_zarr(dead_layer.data, dead_outpath)
            log(f"Saved Death labels → {dead_outpath}")

        log("✅ All user labels saved!")

    def _run_training(self, cell_types_to_train):
        """Train (and apply) APOC classifiers for the given list of cell types."""
        import apoc

        # Auto-save labels before training
        print("Auto-saving user labels before training...")
        self.save_user_labels(log=print)

        successes = []
        try:
            for ct in cell_types_to_train:
                tab = self.tabs[ct]
                cfg = tab.get_config()
                feature_string = cfg["feature_string"]
                max_depth = cfg["max_depth"]
                num_ensembles = cfg["num_ensembles"]

                self.status_label.setText(f"Processing {ct}...")
                QApplication.processEvents()

                images = self._get_images_for_tab(ct)
                if not images:
                    self.status_label.setText(f"⚠️ No image layers selected for '{ct}'!")
                    continue

                # Get annotation layer (keep lazy if possible)
                if ct == "dead":
                    layer_name = "User Provided Labels (Dead)"
                else:
                    layer_name = f"User Provided Labels ({ct.capitalize()})"

                try:
                    annotation = self.viewer.layers[layer_name].data
                except KeyError:
                    print(f"Skipping '{ct}': annotation layer not found")
                    continue

                if not np.any(annotation):
                    print(f"Skipping '{ct}': no labels drawn")
                    continue

                if ct == "dead":
                    clf_name = "PixelClassifier_Death.cl"
                else:
                    clf_name = f"PixelClassifier_{ct.capitalize()}.cl"

                clf_path = str(Path(self.pixel_class_outdir, clf_name))

                # Erase existing classifier
                if Path(clf_path).exists():
                    apoc.erase_classifier(clf_path)

                # Train with ObjectSegmenter
                clf = apoc.ObjectSegmenter(
                    opencl_filename=clf_path,
                    max_depth=max_depth,
                    num_ensembles=num_ensembles,
                    positive_class_identifier=2,
                )

                has_trained = False
                n_timepoints = annotation.shape[0] if annotation.ndim == 4 else 1

                if annotation.ndim == 4:
                    for t in range(n_timepoints):
                        # Load only the current timepoint slice into memory
                        ann_t = np.asarray(annotation[t])
                        if not np.any(ann_t):
                            continue
                        imgs_t = [np.asarray(img[t]) for img in images]
                        imgs_to_pass = imgs_t[0] if len(imgs_t) == 1 else imgs_t
                        clf.train(feature_string, ann_t, imgs_to_pass, continue_training=has_trained)
                        has_trained = True
                else:
                    # Single timepoint
                    ann_np = np.asarray(annotation)
                    imgs_np = [np.asarray(img) for img in images]
                    imgs_to_pass = imgs_np[0] if len(imgs_np) == 1 else imgs_np
                    clf.train(feature_string, ann_np, imgs_to_pass)
                    has_trained = True

                if not has_trained:
                    continue

                # Predict (visual confirmation)
                if images[0].ndim == 4:
                    results = []
                    for t in range(n_timepoints):
                        # Load only current timepoint for prediction
                        imgs_t = [np.asarray(img[t]) for img in images]
                        imgs_to_pass = imgs_t[0] if len(imgs_t) == 1 else imgs_t
                        res_t = clf.predict(image=imgs_to_pass)
                        results.append(np.asarray(res_t).astype(np.int16))
                    result = np.stack(results, axis=0)
                else:
                    imgs = images[0] if len(images) == 1 else images
                    result = np.asarray(clf.predict(image=imgs)).astype(np.int16)

                # Show in viewer
                seg_layer_name = (
                    "Pixel Classification (Dead)" if ct == "dead"
                    else f"{ct.capitalize()} Segments"
                )
                if seg_layer_name in [l.name for l in self.viewer.layers]:
                    self.viewer.layers[seg_layer_name].data = result
                    self.viewer.layers[seg_layer_name].visible = True
                else:
                    self.viewer.add_labels(result, name=seg_layer_name, opacity=0.8, visible=True)

                # Save to disk
                arr_path = Path(self.pixel_class_outdir, f"{Path(clf_path).stem}_PredictedLabels.zarr")
                if arr_path.exists():
                    shutil.rmtree(arr_path)
                save_as_zarr(result, arr_path)

                successes.append(ct)

            if successes:
                self.status_label.setText(f"✅ Trained: {', '.join(successes)}")
            else:
                self.status_label.setText("⚠️ No cell types were trained (check labels).")

            # Persist all params back to caller
            self._persist_params()

        except Exception as e:
            self.status_label.setText(f"❌ Error: {e}")
            import traceback
            traceback.print_exc()

    # ------------------------------------------------------------------
    # Config persistence
    # ------------------------------------------------------------------

    def _collect_all_params(self):
        """Return a flat param dict for all tabs, suitable for storing in BEHAV3D config."""
        params = {}
        for ct, tab in self.tabs.items():
            cfg = tab.get_config()
            params[f"apoc_{ct}_feature_preset"]        = cfg["feature_preset"]
            params[f"apoc_{ct}_sigmas"]                = cfg["sigmas"]
            params[f"apoc_{ct}_custom_feature_string"] = cfg["custom_feature_string"]
            params[f"apoc_{ct}_feature_string"]        = cfg["feature_string"]
            params[f"apoc_{ct}_max_depth"]             = cfg["max_depth"]
            params[f"apoc_{ct}_num_ensembles"]         = cfg["num_ensembles"]
            params[f"apoc_{ct}_channels"]              = cfg["channels"]
        return params

    def _persist_params(self):
        """Push current params back to the notebook via the on_params_changed callback."""
        if callable(self._on_params_changed):
            try:
                self._on_params_changed(self._collect_all_params())
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Entry-point called from segmentation.py
# ---------------------------------------------------------------------------

def train_pixel_classifier_apoc(
    output_dir,
    metadata,
    examples_per_sample=3,
    overwrite_images=False,
    organoid_types=None,
    immune_types=None,
    other_types=None,
    initial_params=None,
    on_params_changed=None,
):
    """
    Open a Napari viewer with the tabbed APOC classification widget.
    The original napari_pixelclassifier.py is NOT called.
    """
    ip = initial_params or {}
    gpu_device = ip.get("gpu_device_name")
    if gpu_device:
        print(f"APOC Training: Selecting device {gpu_device}")
        cle.select_device(gpu_device)

    organoid_types = organoid_types or []
    immune_types   = immune_types   or []
    other_types    = other_types    or []
    all_cell_types = organoid_types + immune_types + other_types

    # --- Load images ---
    print("Loading training images for APOC...")
    image_list, pixel_class_outdir, has_death, all_cell_types = _load_training_images(
        metadata=metadata,
        output_dir=output_dir,
        examples_per_sample=examples_per_sample,
        organoid_types=organoid_types,
        immune_types=immune_types,
        other_types=other_types,
        overwrite_images=overwrite_images,
    )

    if not image_list:
        print("⚠️ No training images found!")
        return None

    # --- Pad images to the same shape and stack along T ---
    max_shape = list(image_list[0].shape)
    for img in image_list[1:]:
        for i in range(len(max_shape)):
            max_shape[i] = max(max_shape[i], img.shape[i])
    image_list = [zeropad_image_to_match_shape(img, max_shape) for img in image_list]

    stacked = np.stack(image_list, axis=0)  # (T_selected, C, Z, Y, X)

    print(f"Training data shape (TCZYX): {stacked.shape}")

    # --- Create Napari viewer ---
    viewer = napari.Viewer()

    n_channels = stacked.shape[1]
    channel_colors = [
        "cyan", "yellow", "red", "green", "magenta", "blue",
        "gray", "turbo", "viridis", "plasma", "inferno", "twilight",
    ]

    for ch in range(n_channels):
        channel_data = stacked[:, ch, ...]  # (T, Z, Y, X)
        nonzero = channel_data[channel_data > 0]
        clim = (0, float(np.percentile(nonzero, 99.8))) if nonzero.size > 0 else (0, 1e-3)
        viewer.add_image(
            channel_data,
            name=f"Channel {ch}",
            contrast_limits=clim,
            colormap=channel_colors[ch % len(channel_colors)],
            blending="additive",
            opacity=0.8,
        )

    # --- Annotation label layers per cell type ---
    time_shape    = stacked.shape[0]
    spatial_shape = stacked.shape[2:]
    full_shape    = (time_shape,) + spatial_shape  # (T, Z, Y, X)

    ip = initial_params or {}

    # 1. User Provided Labels (all cell types)
    for cell_type in all_cell_types:
        saved_path = Path(pixel_class_outdir, f"PixelClassifier_User{cell_type.capitalize()}Labels.zarr")
        if saved_path.exists():
            existing = np.asarray(load_zarr(saved_path))
            user_labels = existing if existing.shape == full_shape else np.zeros(full_shape, dtype=np.int16)
        else:
            user_labels = np.zeros(full_shape, dtype=np.int16)

        viewer.add_labels(
            user_labels,
            name=f"User Provided Labels ({cell_type.capitalize()})",
            opacity=0.5,
        )

    if has_death:
        dead_path = Path(pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
        if dead_path.exists():
            dead_labels = np.asarray(load_zarr(dead_path))
            if dead_labels.shape != full_shape:
                dead_labels = np.zeros(full_shape, dtype=np.int16)
        else:
            dead_labels = np.zeros(full_shape, dtype=np.int16)
        viewer.add_labels(dead_labels, name="User Provided Labels (Dead)", opacity=0.5)

    # 2. Results (Pixel Classification / Segments) on top
    # Dead result
    if has_death:
        pred_death_path = Path(pixel_class_outdir, "PixelClassifier_Death_PredictedLabels.zarr")
        if pred_death_path.exists():
            pred_death = np.asarray(load_zarr(pred_death_path))
            if pred_death.shape != full_shape:
                pred_death = np.zeros(full_shape, dtype=np.int16)
        else:
            pred_death = np.zeros(full_shape, dtype=np.int16)
        viewer.add_labels(pred_death, name="Pixel Classification (Dead)", opacity=0.8, visible=False)

    # Cell type segments
    for cell_type in all_cell_types:
        pred_path = Path(pixel_class_outdir, f"PixelClassifier_{cell_type.capitalize()}_PredictedLabels.zarr")
        if pred_path.exists():
            pred_data = np.asarray(load_zarr(pred_path))
            if pred_data.shape != full_shape:
                pred_data = np.zeros(full_shape, dtype=np.int16)
        else:
            pred_data = np.zeros(full_shape, dtype=np.int16)

        viewer.add_labels(
            pred_data,
            name=f"{cell_type.capitalize()} Segments",
            opacity=0.8,
            visible=False,
        )

    # --- Dock the APOC widget ---
    apoc_widget = APOCTrainingWidget(
        viewer=viewer,
        pixel_class_outdir=pixel_class_outdir,
        all_cell_types=all_cell_types,
        has_death=has_death,
        initial_params=ip,
        on_params_changed=on_params_changed,
    )
    viewer.window.add_dock_widget(apoc_widget, area="right", name="APOC Pixel Classification")

    log_output = QPlainTextEdit()
    log_output.setReadOnly(True)
    log_output.setMaximumHeight(120)
    log_widget = QWidget()
    log_layout = QVBoxLayout()
    log_layout.addWidget(log_output)
    log_widget.setLayout(log_layout)
    viewer.window.add_dock_widget(log_widget, area="right", name="Log Output")

    save_button = QtPushButton("💾 Save User Labels")
    save_button.setStyleSheet("background-color: #FF9800; color: white; font-weight: bold; padding: 6px;")
    save_button.clicked.connect(lambda: apoc_widget.save_user_labels(log=log_output.appendPlainText))
    apoc_widget.layout().addWidget(save_button)

    return viewer
