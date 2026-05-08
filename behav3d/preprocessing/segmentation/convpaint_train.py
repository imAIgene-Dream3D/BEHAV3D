"""
ConvPaint training widget for BEHAV3D.

Single global config panel (shared across all cell types).
All channels are always used. Excludes *_merge cell types from training.
"""

import gc, shutil, time
from pathlib import Path
import numpy as np
import napari
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QSpinBox, QDoubleSpinBox, QPushButton as QtPushButton,
    QGroupBox, QPlainTextEdit, QApplication, QScrollArea,
)
from behav3d.io.images import load_image, load_zarr, save_as_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape

CONVPAINT_FE_OPTIONS = [
    ("VGG16 (default)", "vgg"), ("VGG16 Medium", "vgg-m"), ("VGG16 Large", "vgg-l"),
    ("DINOv2", "dino"), ("JAFAR DINOv2", "dino-jafar"),
    ("Gaussian", "gaussian"), ("Cellpose", "cellpose"), ("Ilastik", "ilastik"),
]
CHANNEL_MODE_OPTIONS = [("Single channel", "single"), ("Multichannel", "multi")]
NORMALIZE_OPTIONS = [("No normalization", 1), ("Normalize stack", 2), ("Normalize per plane", 3)]
CHANNEL_COLORS = [
    "cyan", "yellow", "red", "green", "magenta", "blue",
    "gray", "turbo", "viridis", "plasma", "inferno", "twilight",
]


def _detect_torch_devices():
    devices = [("CPU", "cpu"), ("Auto (best available)", "auto")]
    try:
        import torch
        if torch.cuda.is_available():
            for i in range(torch.cuda.device_count()):
                devices.append((f"GPU {i}: {torch.cuda.get_device_name(i)}", f"cuda:{i}"))
    except Exception:
        pass
    return devices


def _set_torch_device(device_str):
    if not device_str or device_str in ("auto", "cpu"):
        return
    try:
        import torch
        if device_str.startswith("cuda:"):
            idx = int(device_str.split(":")[1])
            torch.cuda.set_device(idx)
            print(f"ConvPaint: Selected GPU {idx} ({torch.cuda.get_device_name(idx)})")
    except Exception as e:
        print(f"\u26a0\ufe0f Could not set CUDA device '{device_str}': {e}")


def _filter_merge_types(cell_types):
    """Remove derived output types (*_merged, *_grouped). Keep *_N_multicolor."""
    from behav3d.core.metadata import is_combined_multicolor_celltype
    return [ct for ct in cell_types
            if not is_combined_multicolor_celltype(ct)]


# ---------------------------------------------------------------------------
# Image loading
# ---------------------------------------------------------------------------

def _load_training_images(metadata, output_dir, examples_per_sample,
                          organoid_types, immune_types, other_types,
                          overwrite_images=False):
    from behav3d.core.metadata import has_dead_channel
    has_death = has_dead_channel(metadata)
    all_cell_types = organoid_types + immune_types + other_types
    pixel_class_outdir = Path(output_dir, "images", "PixelClassification")
    pixel_class_outdir.mkdir(exist_ok=True, parents=True)
    image_outpath = Path(pixel_class_outdir, "PixelClassifier_Images.zarr")

    if image_outpath.exists():
        if overwrite_images:
            shutil.rmtree(image_outpath)
        else:
            print("Loading cached training images...")
            cached = load_zarr(image_outpath)
            all_images = [np.asarray(cached[:, t]) for t in range(cached.shape[1])]
            return all_images, pixel_class_outdir, has_death, all_cell_types

    all_images, max_shape = [], None
    for _, sample in metadata.iterrows():
        sample_name = sample.get("sample_name", "unknown")
        raw_path = sample.get("raw_image_path", "")
        if not raw_path or not Path(raw_path).exists():
            continue
        axis_order = sample.get("dimension_order", "TCZYX")
        if not isinstance(axis_order, str) or not axis_order:
            axis_order = "TCZYX"
        img = load_image(raw_path, axis_order=axis_order)
        n_tp = img.shape[0]
        n_sel = min(examples_per_sample, n_tp)
        t_indices = [0] if n_sel <= 1 else sorted(set(np.round(np.linspace(0, n_tp-1, n_sel)).astype(int)))
        print(f"  {sample_name}: selected {len(t_indices)} timepoints {list(t_indices)}")
        for t_idx in t_indices:
            frame = np.asarray(np.take(img, t_idx, axis=0))
            all_images.append(frame)
            if max_shape is None:
                max_shape = list(frame.shape)
            else:
                for d in range(len(max_shape)):
                    max_shape[d] = max(max_shape[d], frame.shape[d])

    if max_shape is not None:
        all_images = [zeropad_image_to_match_shape(im, max_shape) for im in all_images]

    import dask.array as da
    stack = da.stack(all_images).transpose(1, 0, 2, 3, 4)
    save_as_zarr(stack, image_outpath)
    del stack; gc.collect()
    cached = load_zarr(image_outpath)
    all_images = [np.asarray(cached[:, t]) for t in range(cached.shape[1])]
    return all_images, pixel_class_outdir, has_death, all_cell_types


# ---------------------------------------------------------------------------
# Widget
# ---------------------------------------------------------------------------

class ConvPaintTrainingWidget(QWidget):
    def __init__(self, viewer, all_images, pixel_class_outdir,
                 all_cell_types, has_death, initial_params=None,
                 on_params_changed=None, convpaint_strategy=None, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.all_images = all_images
        self.pixel_class_outdir = Path(pixel_class_outdir)
        self.all_cell_types = all_cell_types
        self.has_death = has_death
        self._on_params_changed = on_params_changed
        self.ip = initial_params or {}
        self.convpaint_strategy = convpaint_strategy or self.ip.get(
            "convpaint_strategy", "ConvPaint (Direct Segmentation)")
        self._build_ui()

    def _build_ui(self):
        root = QVBoxLayout()
        root.setContentsMargins(4, 4, 4, 4)
        root.addWidget(QLabel("<h3>ConvPaint Training</h3>"))

        # ── Active cell type ─────────────────────────────────────────
        ct_row = QHBoxLayout()
        ct_row.addWidget(QLabel("Active cell type:"))
        self.cell_type_combo = QComboBox()
        tab_types = list(self.all_cell_types)
        if self.has_death:
            tab_types.append("dead")
        for ct in tab_types:
            self.cell_type_combo.addItem(ct.capitalize(), ct)
        self.cell_type_combo.currentIndexChanged.connect(self._on_cell_type_changed)
        ct_row.addWidget(self.cell_type_combo)
        root.addLayout(ct_row)

        # ── Train buttons ────────────────────────────────────────────
        btn_row = QHBoxLayout()
        self.btn_train_current = QtPushButton("▶ Train Current")
        self.btn_train_current.setStyleSheet(
            "background-color: #1976D2; color: white; font-weight: bold; padding: 6px 12px;"
        )
        self.btn_train_current.clicked.connect(self._on_train_current)
        btn_row.addWidget(self.btn_train_current)
        self.btn_train_all = QtPushButton("▶▶ Train All")
        self.btn_train_all.setStyleSheet(
            "background-color: #2e7d32; color: white; font-weight: bold; padding: 6px 12px;"
        )
        self.btn_train_all.clicked.connect(self._on_train_all)
        btn_row.addWidget(self.btn_train_all)
        root.addLayout(btn_row)

        # ── Segmentation params (for EDT/Probability strategies) ─────
        is_edt = "EDT" in self.convpaint_strategy
        is_prob = "Probability" in self.convpaint_strategy
        show_seg_params = is_edt or is_prob

        self.seg_group = QGroupBox("Segmentation Parameters")
        seg_layout = QVBoxLayout()

        # EDT threshold
        r1 = QHBoxLayout()
        r1.addWidget(QLabel("EDT threshold:"))
        self.edt_thr_spin = QDoubleSpinBox()
        self.edt_thr_spin.setRange(0.0, 50.0)
        self.edt_thr_spin.setSingleStep(0.5)
        self.edt_thr_spin.setDecimals(1)
        r1.addWidget(self.edt_thr_spin)
        seg_layout.addLayout(r1)

        # Min segment size
        r2 = QHBoxLayout()
        r2.addWidget(QLabel("Min segment size:"))
        self.seg_min_spin = QSpinBox()
        self.seg_min_spin.setRange(0, 100000)
        r2.addWidget(self.seg_min_spin)
        seg_layout.addLayout(r2)

        # Opening
        r3 = QHBoxLayout()
        r3.addWidget(QLabel("Opening pixels:"))
        self.opening_spin = QSpinBox()
        self.opening_spin.setRange(0, 50)
        r3.addWidget(self.opening_spin)
        seg_layout.addLayout(r3)

        # Fill holes (checkbox-like spin: 0=no, 1=yes)
        r4 = QHBoxLayout()
        r4.addWidget(QLabel("Fill holes:"))
        self.fill_holes_combo = QComboBox()
        self.fill_holes_combo.addItems(["No", "Yes"])
        r4.addWidget(self.fill_holes_combo)
        seg_layout.addLayout(r4)

        # Probability-specific
        self.prob_label_mask = QLabel("Prob mask threshold:")
        self.prob_mask_spin = QDoubleSpinBox()
        self.prob_mask_spin.setRange(0.0, 1.0)
        self.prob_mask_spin.setSingleStep(0.05)
        self.prob_mask_spin.setDecimals(2)
        r5 = QHBoxLayout()
        r5.addWidget(self.prob_label_mask)
        r5.addWidget(self.prob_mask_spin)
        seg_layout.addLayout(r5)

        self.prob_label_seed = QLabel("Prob seed threshold:")
        self.prob_seed_spin = QDoubleSpinBox()
        self.prob_seed_spin.setRange(0.0, 1.0)
        self.prob_seed_spin.setSingleStep(0.05)
        self.prob_seed_spin.setDecimals(2)
        r6 = QHBoxLayout()
        r6.addWidget(self.prob_label_seed)
        r6.addWidget(self.prob_seed_spin)
        seg_layout.addLayout(r6)

        self.seg_group.setLayout(seg_layout)
        root.addWidget(self.seg_group)

        # Show/hide prob-specific widgets
        self.prob_label_mask.setVisible(is_prob)
        self.prob_mask_spin.setVisible(is_prob)
        self.prob_label_seed.setVisible(is_prob)
        self.prob_seed_spin.setVisible(is_prob)
        self.seg_group.setVisible(show_seg_params)

        # Load initial values for first cell type
        self._load_seg_params_for_cell_type()

        # Connect changes to persist
        for w in [self.edt_thr_spin, self.seg_min_spin, self.opening_spin,
                  self.prob_mask_spin, self.prob_seed_spin]:
            if isinstance(w, QDoubleSpinBox):
                w.valueChanged.connect(self._on_seg_param_changed)
            else:
                w.valueChanged.connect(self._on_seg_param_changed)
        self.fill_holes_combo.currentIndexChanged.connect(self._on_seg_param_changed)

        # ── Predict button ───────────────────────────────────────────
        self.btn_predict = QtPushButton("🔮 Predict All Cell Types")
        self.btn_predict.setStyleSheet(
            "background-color: #7B1FA2; color: white; font-weight: bold; padding: 6px 12px;"
        )
        self.btn_predict.clicked.connect(self._on_predict_all)
        root.addWidget(self.btn_predict)

        # ── Save Labels button ───────────────────────────────────────
        self.btn_save_labels = QtPushButton("💾 Save User Labels")
        self.btn_save_labels.setStyleSheet(
            "background-color: #FF9800; color: white; font-weight: bold; padding: 6px;"
        )
        self.btn_save_labels.clicked.connect(self._on_save_labels)
        root.addWidget(self.btn_save_labels)

        # ── Feature Extractor config ─────────────────────────────────
        fe_group = QGroupBox("Feature Extractor")
        fe_layout = QVBoxLayout()
        r = QHBoxLayout()
        r.addWidget(QLabel("Model:"))
        self.fe_combo = QComboBox()
        for label, alias in CONVPAINT_FE_OPTIONS:
            self.fe_combo.addItem(label, alias)
        idx = self.fe_combo.findData(self.ip.get("convpaint_fe_alias", "vgg"))
        if idx >= 0: self.fe_combo.setCurrentIndex(idx)
        r.addWidget(self.fe_combo)
        fe_layout.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Channel mode:"))
        self.channel_mode_combo = QComboBox()
        for label, mode in CHANNEL_MODE_OPTIONS:
            self.channel_mode_combo.addItem(label, mode)
        idx = self.channel_mode_combo.findData(self.ip.get("convpaint_channel_mode", "multi"))
        if idx >= 0: self.channel_mode_combo.setCurrentIndex(idx)
        r.addWidget(self.channel_mode_combo)
        fe_layout.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Normalization:"))
        self.normalize_combo = QComboBox()
        for label, val in NORMALIZE_OPTIONS:
            self.normalize_combo.addItem(label, val)
        idx = self.normalize_combo.findData(int(self.ip.get("convpaint_normalize", 2)))
        if idx >= 0: self.normalize_combo.setCurrentIndex(idx)
        r.addWidget(self.normalize_combo)
        fe_layout.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Downsample:"))
        self.downsample_spin = QSpinBox(); self.downsample_spin.setRange(-4, 16)
        self.downsample_spin.setValue(int(self.ip.get("convpaint_image_downsample", 1)))
        r.addWidget(self.downsample_spin)
        r.addWidget(QLabel("Smoothen:"))
        self.smooth_spin = QSpinBox(); self.smooth_spin.setRange(0, 20)
        self.smooth_spin.setValue(int(self.ip.get("convpaint_seg_smoothening", 1)))
        r.addWidget(self.smooth_spin)
        fe_layout.addLayout(r)
        fe_group.setLayout(fe_layout)
        root.addWidget(fe_group)

        # ── Classifier params ────────────────────────────────────────
        clf_group = QGroupBox("Classifier (CatBoost)")
        clf_layout = QVBoxLayout()
        r = QHBoxLayout()
        r.addWidget(QLabel("Iterations:"))
        self.iterations_spin = QSpinBox(); self.iterations_spin.setRange(10, 2000)
        self.iterations_spin.setValue(int(self.ip.get("convpaint_clf_iterations", 100)))
        r.addWidget(self.iterations_spin)
        r.addWidget(QLabel("Depth:"))
        self.depth_spin = QSpinBox(); self.depth_spin.setRange(1, 16)
        self.depth_spin.setValue(int(self.ip.get("convpaint_clf_depth", 5)))
        r.addWidget(self.depth_spin)
        clf_layout.addLayout(r)
        r = QHBoxLayout()
        r.addWidget(QLabel("Learning rate:"))
        self.lr_spin = QDoubleSpinBox(); self.lr_spin.setRange(0.001, 1.0)
        self.lr_spin.setSingleStep(0.01); self.lr_spin.setDecimals(3)
        self.lr_spin.setValue(float(self.ip.get("convpaint_clf_learning_rate", 0.1)))
        r.addWidget(self.lr_spin)
        clf_layout.addLayout(r)
        clf_group.setLayout(clf_layout)
        root.addWidget(clf_group)

        # ── Device ───────────────────────────────────────────────────
        r = QHBoxLayout()
        r.addWidget(QLabel("Device:"))
        self.device_combo = QComboBox()
        for label, dev_str in _detect_torch_devices():
            self.device_combo.addItem(label, dev_str)
        idx = self.device_combo.findData(self.ip.get("convpaint_device", "auto"))
        if idx >= 0: self.device_combo.setCurrentIndex(idx)
        r.addWidget(self.device_combo)
        root.addLayout(r)

        hint = QLabel("<i>Labels: <b>1</b>=background <b>2</b>=foreground. All channels used.</i>")
        hint.setStyleSheet("color: #666; margin-top: 4px;")
        root.addWidget(hint)

        # ── Log ──────────────────────────────────────────────────────
        self.log_box = QPlainTextEdit(); self.log_box.setReadOnly(True); self.log_box.setMaximumHeight(150)
        root.addWidget(self.log_box)
        root.addStretch()

        scroll_content = QWidget(); scroll_content.setLayout(root)
        scroll = QScrollArea(); scroll.setWidget(scroll_content); scroll.setWidgetResizable(True)
        outer = QVBoxLayout(); outer.setContentsMargins(0, 0, 0, 0); outer.addWidget(scroll)
        self.setLayout(outer)

    # ── Seg params per cell type ─────────────────────────────────────

    def _load_seg_params_for_cell_type(self):
        ct = self.cell_type_combo.currentData() or ""
        ip = self.ip
        self.edt_thr_spin.setValue(float(ip.get(f"{ct}_edt_threshold", 1.0)))
        self.seg_min_spin.setValue(int(ip.get(f"{ct}_segment_size_min", 10)))
        self.opening_spin.setValue(int(ip.get(f"{ct}_opening_nr_pixels", 0)))
        self.fill_holes_combo.setCurrentIndex(1 if ip.get(f"{ct}_fill_holes", True) else 0)
        self.prob_mask_spin.setValue(float(ip.get(f"{ct}_prob_mask_threshold", 0.5)))
        self.prob_seed_spin.setValue(float(ip.get(f"{ct}_prob_seed_threshold", 0.8)))

    def _on_cell_type_changed(self, _idx):
        self._load_seg_params_for_cell_type()

    def _on_seg_param_changed(self, _val=None):
        """Persist seg params back to notebook widget."""
        ct = self.cell_type_combo.currentData() or ""
        params = {
            f"{ct}_edt_threshold": self.edt_thr_spin.value(),
            f"{ct}_segment_size_min": self.seg_min_spin.value(),
            f"{ct}_opening_nr_pixels": self.opening_spin.value(),
            f"{ct}_fill_holes": (self.fill_holes_combo.currentIndex() == 1),
            f"{ct}_prob_mask_threshold": self.prob_mask_spin.value(),
            f"{ct}_prob_seed_threshold": self.prob_seed_spin.value(),
        }
        # Also update our ip cache
        self.ip.update(params)
        if self._on_params_changed:
            self._on_params_changed(params)

    # ── Model building ───────────────────────────────────────────────

    def _log(self, msg):
        self.log_box.appendPlainText(msg); QApplication.processEvents()

    def _build_model(self):
        from napari_convpaint import ConvpaintModel
        model = ConvpaintModel(self.fe_combo.currentData())
        model.set_params(
            channel_mode=self.channel_mode_combo.currentData(),
            normalize=self.normalize_combo.currentData(),
            image_downsample=self.downsample_spin.value(),
            seg_smoothening=self.smooth_spin.value(),
            clf_iterations=self.iterations_spin.value(),
            clf_learning_rate=self.lr_spin.value(),
            clf_depth=self.depth_spin.value(),
        )
        device = self.device_combo.currentData()
        if device and device not in ("auto", "cpu"):
            _set_torch_device(device); fe_device = "gpu"
        elif device == "cpu":
            fe_device = "cpu"
        else:
            fe_device = None
        return model, fe_device

    def _get_config(self):
        return {
            "convpaint_fe_alias": self.fe_combo.currentData(),
            "convpaint_channel_mode": self.channel_mode_combo.currentData(),
            "convpaint_normalize": self.normalize_combo.currentData(),
            "convpaint_image_downsample": self.downsample_spin.value(),
            "convpaint_seg_smoothening": self.smooth_spin.value(),
            "convpaint_clf_iterations": self.iterations_spin.value(),
            "convpaint_clf_depth": self.depth_spin.value(),
            "convpaint_clf_learning_rate": self.lr_spin.value(),
            "convpaint_device": self.device_combo.currentData(),
        }

    # ── Training ─────────────────────────────────────────────────────

    def _train_cell_type(self, cell_type):
        layer_name = (f"User Provided Labels ({cell_type.capitalize()})"
                      if cell_type != "dead" else "User Provided Labels (Dead)")
        if layer_name not in self.viewer.layers:
            self._log(f"\u274c Label layer '{layer_name}' not found!"); return
        annotations = np.asarray(self.viewer.layers[layer_name].data)
        train_images, train_annots = [], []
        for i, img in enumerate(self.all_images):
            if np.any(annotations[i] > 0):
                train_images.append(img); train_annots.append(annotations[i])
        if not train_images:
            self._log(f"\u26a0\ufe0f No annotations for {cell_type}."); return
        self._log(f"Training {cell_type} on {len(train_images)} frames...")
        try:
            model, fe_device = self._build_model()
            model.train(train_images, train_annots, fe_use_device=fe_device)
            model_path = self.pixel_class_outdir / f"ConvPaintModel_{cell_type.capitalize()}.pkl"
            model.save(str(model_path))
            self._log(f"\u2705 Saved: {model_path.name}")
            if self._on_params_changed:
                self._on_params_changed(self._get_config())
        except Exception as e:
            self._log(f"\u274c Training failed: {e}")
            import traceback; traceback.print_exc()

    def _on_train_current(self):
        # Auto-save labels before training
        self._log("Auto-saving user labels before training...")
        self.save_user_labels()
        ct = self.cell_type_combo.currentData()
        if ct: self._train_cell_type(ct)

    def _on_train_all(self):
        # Auto-save labels before training
        self._log("Auto-saving user labels before training...")
        self.save_user_labels()
        all_types = list(self.all_cell_types)
        if self.has_death: all_types.append("dead")
        for ct in all_types:
            self._train_cell_type(ct)

    # ── Prediction ───────────────────────────────────────────────────

    def _on_predict_all(self):
        from napari_convpaint import ConvpaintModel
        device = self.device_combo.currentData()
        if device and device not in ("auto", "cpu"):
            _set_torch_device(device); fe_device = "gpu"
        elif device == "cpu":
            fe_device = "cpu"
        else:
            fe_device = None

        all_types = list(self.all_cell_types)
        if self.has_death: all_types.append("dead")

        for ct in all_types:
            model_path = self.pixel_class_outdir / f"ConvPaintModel_{ct.capitalize()}.pkl"
            if not model_path.exists():
                self._log(f"\u26a0\ufe0f No model for {ct}, skipping."); continue
            try:
                model = ConvpaintModel(model_path=str(model_path))
                preds = [np.asarray(model.segment(img, fe_use_device=fe_device))
                         for img in self.all_images]
                pred_stack = np.stack(preds, axis=0)
                seg_name = f"{ct.capitalize()} Segments"
                if seg_name in self.viewer.layers:
                    self.viewer.layers[seg_name].data = pred_stack
                else:
                    self.viewer.add_labels(pred_stack, name=seg_name)
                self._log(f"\u2705 Predicted {ct}")
            except Exception as e:
                self._log(f"\u274c Prediction failed for {ct}: {e}")

    # ── Save / load labels ───────────────────────────────────────

    def save_user_labels(self):
        """Save all user-provided labels for all cell types and Dead (if present).
        Uses the same path convention as APOC so annotations are reusable."""
        for cell_type in self.all_cell_types:
            lname = f"User Provided Labels ({cell_type.capitalize()})"
            if lname not in [l.name for l in self.viewer.layers]:
                continue
            label_layer = self.viewer.layers[lname]
            label_data = np.asarray(label_layer.data)
            outpath = Path(self.pixel_class_outdir, f"PixelClassifier_User{cell_type.capitalize()}Labels.zarr")
            if outpath.exists():
                import shutil; shutil.rmtree(outpath)
            save_as_zarr(label_data, outpath)
            self._log(f"Saved {cell_type} labels \u2192 {outpath.name}")

        if self.has_death and "User Provided Labels (Dead)" in [l.name for l in self.viewer.layers]:
            dead_layer = self.viewer.layers["User Provided Labels (Dead)"]
            dead_data = np.asarray(dead_layer.data)
            dead_outpath = Path(self.pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
            if dead_outpath.exists():
                import shutil; shutil.rmtree(dead_outpath)
            save_as_zarr(dead_data, dead_outpath)
            self._log(f"Saved Dead labels \u2192 {dead_outpath.name}")

        self._log("\u2705 All user labels saved!")

    def _on_save_labels(self):
        self.save_user_labels()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def train_pixel_classifier_convpaint(
    output_dir, metadata, examples_per_sample=3, overwrite_images=False,
    organoid_types=None, immune_types=None, other_types=None,
    initial_params=None, on_params_changed=None,
):
    organoid_types = _filter_merge_types(organoid_types or [])
    immune_types = _filter_merge_types(immune_types or [])
    other_types = _filter_merge_types(other_types or [])

    all_images, pixel_class_outdir, has_death, all_cell_types = _load_training_images(
        metadata, output_dir, examples_per_sample,
        organoid_types, immune_types, other_types,
        overwrite_images=overwrite_images,
    )
    # Filter merge types from detected cell types too
    all_cell_types = _filter_merge_types(all_cell_types)

    if not all_images:
        print("\u274c No training images loaded."); return None

    n_channels = all_images[0].shape[0]
    print(f"Loaded {len(all_images)} training images with {n_channels} channels")

    stacked = np.stack(all_images, axis=0)
    T_total = stacked.shape[0]
    viewer = napari.Viewer(title="BEHAV3D - ConvPaint Training")

    # Add channels (matching APOC: additive blending, per-channel colors)
    for ch in range(n_channels):
        ch_data = stacked[:, ch, :, :, :]
        nonzero = ch_data[ch_data > 0]
        clim = (0, float(np.percentile(nonzero, 99.8))) if nonzero.size > 0 else (0, 1e-3)
        layer = viewer.add_image(ch_data, name=f"Channel {ch}",
            contrast_limits=clim, colormap=CHANNEL_COLORS[ch % len(CHANNEL_COLORS)],
            blending="additive", opacity=0.8)
        layer.contrast_limits_range = (0, float(ch_data.max()))

    # Annotation layers
    label_shape = (T_total,) + stacked.shape[2:]
    ip = initial_params or {}

    def _restore(path, fallback_path, shape, name):
        """Try to restore labels from primary path, then fallback."""
        for p in [path, fallback_path]:
            if p is not None and p.exists():
                try:
                    existing = np.asarray(load_zarr(p))
                    if existing.shape == shape:
                        print(f"  \u2705 Restored saved labels for '{name}' from {p.name}")
                        return existing
                    print(f"  \u26a0\ufe0f Shape mismatch for '{name}': {existing.shape} vs {shape}")
                except Exception as exc:
                    print(f"  \u26a0\ufe0f Could not restore '{name}' from {p.name}: {exc}")
        return np.zeros(shape, dtype=np.int16)

    for ct in all_cell_types:
        # Primary: APOC-compatible path; Fallback: old ConvPaint path
        primary = Path(pixel_class_outdir, f"PixelClassifier_User{ct.capitalize()}Labels.zarr")
        fallback = Path(pixel_class_outdir, f"ConvPaintLabels_{ct.capitalize()}.zarr")
        viewer.add_labels(_restore(primary, fallback, label_shape, ct),
            name=f"User Provided Labels ({ct.capitalize()})", opacity=0.5)

    if has_death:
        primary = Path(pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
        fallback = Path(pixel_class_outdir, "ConvPaintLabels_Dead.zarr")
        viewer.add_labels(_restore(primary, fallback, label_shape, "dead"),
            name="User Provided Labels (Dead)", opacity=0.5)

    # Resolve strategy from initial_params
    cp_strategy = ip.get("convpaint_strategy", "ConvPaint (Direct Segmentation)")

    widget = ConvPaintTrainingWidget(
        viewer=viewer, all_images=all_images, pixel_class_outdir=pixel_class_outdir,
        all_cell_types=all_cell_types, has_death=has_death,
        initial_params=ip, on_params_changed=on_params_changed,
        convpaint_strategy=cp_strategy,
    )
    viewer.window.add_dock_widget(widget, name="ConvPaint", area="right")
    return viewer
