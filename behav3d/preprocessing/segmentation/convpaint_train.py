"""
ConvPaint training widget for BEHAV3D — UNIFIED MULTI-CLASS MODE.

A single ConvPaint classifier predicts every voxel into one of N+1 classes:

* class 1 = background
* classes 2..N+1 = one per cell type (organoid → immune → other)

The mapping between class indices and cell-type names lives in
``ConvPaintModel_Unified_LabelMap.json``, written next to the trained
``ConvPaintModel_Unified.pkl`` file. See
:mod:`behav3d.preprocessing.segmentation.convpaint_label_map`.

Architecture
------------

* A global panel for the Feature Extractor, CatBoost classifier and Device.
* A ``QTabWidget`` with:
    1. **Annotation Legend** — first tab; shows index ↔ color ↔ cell-type with
       voxel counts. Click a row to switch the active paint label on the
       unified labels layer.
    2. One tab per cell type with the strategy-specific instance-segmentation
       parameters (EDT / Probability) and a *Run instance segmentation*
       button. Pressing it does **post-hoc only** processing on the cached
       unified-classifier output (no feature extraction).
    3. *Dead* tab (when applicable) for a separate binary ConvPaint death
       model — death is a state, not a cell-type identity, so it stays as
       its own classifier writing ``mask_dead.zarr``.

Training the unified classifier runs feature extraction exactly once per
training frame, regardless of how many cell types are configured. Multi-class
prediction also runs once per frame to produce the unified labels volume and
per-cell-type probability slices, which feed the per-tab instance previews.

Death model (binary): retains the legacy per-classifier behaviour.
"""

import gc
import shutil
import traceback
from pathlib import Path

import numpy as np
import napari
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QComboBox,
    QSpinBox, QDoubleSpinBox, QPushButton as QtPushButton,
    QGroupBox, QPlainTextEdit, QApplication, QScrollArea,
    QTabWidget, QFrame, QCheckBox, QSizePolicy, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView,
)
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

from behav3d.io.images import load_image, load_zarr, save_as_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape
from behav3d.preprocessing.segmentation.convpaint_postprocess import (
    segment_per_z,
    predict_probas_per_z,
    mask_to_instances,
    probability_to_instances,
)
from behav3d.preprocessing.segmentation.convpaint_label_map import (
    BACKGROUND_LABEL,
    build_label_map,
    save_label_map,
    load_label_map,
    cell_type_order,
    n_classes,
    label_map_matches,
    unified_model_path,
    unified_label_map_path,
    unified_input_channels_path,
    death_input_channels_path,
    save_input_channels,
    load_input_channels,
    unified_user_labels_path,
    unified_predicted_labels_path,
)


# ---------------------------------------------------------------------------
# Strategy + UI option constants.
# ---------------------------------------------------------------------------

STRATEGY_EDT = "ConvPaint Mask + EDT/Watershed"
STRATEGY_PEAK_EDT = "ConvPaint Mask + Peak EDT/Watershed"
STRATEGY_PROB = "ConvPaint Probability + Watershed"
DEFAULT_STRATEGY = STRATEGY_EDT


CONVPAINT_FE_OPTIONS = [
    ("VGG16 (default)", "vgg"), ("VGG16 Medium", "vgg-m"), ("VGG16 Large", "vgg-l"),
    ("DINOv2", "dino"), ("JAFAR DINOv2", "dino-jafar"),
    ("Gaussian", "gaussian"), ("Cellpose", "cellpose"), ("Ilastik", "ilastik"),
]
CHANNEL_MODE_OPTIONS = [
    ("Single channel", "single"),
    ("Multichannel", "multi"),
    ("RGB (3-ch color images)", "rgb"),
]
NORMALIZE_OPTIONS = [
    ("No normalization", 1), ("Normalize stack", 2), ("Normalize per plane", 3),
]
CLASS_WEIGHTS_OPTIONS = [
    ("Balanced (default)", "Balanced"),
    ("SqrtBalanced", "SqrtBalanced"),
    ("None", "None"),
]
DEFAULT_CLASS_WEIGHTS = "Balanced"

CHANNEL_COLORS = [
    "cyan", "yellow", "red", "green", "magenta", "blue",
    "gray", "turbo", "viridis", "plasma", "inferno", "twilight",
]


# ---------------------------------------------------------------------------
# Napari layer names (unified mode).
# ---------------------------------------------------------------------------

UNIFIED_LABELS_LAYER_NAME = "User Provided Labels"
UNIFIED_PREDICTED_LAYER_NAME = "Pixel Classification (Unified)"
DEAD_LABELS_LAYER_NAME = "User Provided Labels (Dead)"
DEAD_PREDICTED_LAYER_NAME = "Pixel Classification (Dead)"


def _segments_layer_name(cell_type):
    return f"{cell_type.capitalize()} Segments"


def _probability_layer_name(cell_type):
    return f"Probability Map ({cell_type.capitalize()})"


# ---------------------------------------------------------------------------
# Small helpers.
# ---------------------------------------------------------------------------

def _compact_combo(combo, min_chars=8):
    """Prevent a QComboBox from expanding to fit its longest item text."""
    combo.setSizeAdjustPolicy(QComboBox.AdjustToMinimumContentsLengthWithIcon)
    combo.setMinimumContentsLength(min_chars)
    combo.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)


def _normalize_strategy(value):
    """Return the matching strategy constant, defaulting to EDT."""
    val = str(value or "").strip()
    if val == STRATEGY_PROB:
        return STRATEGY_PROB
    if val == STRATEGY_PEAK_EDT:
        return STRATEGY_PEAK_EDT
    return STRATEGY_EDT


def _normalize_class_weights(value):
    """Coerce config-supplied class-weights value to a CatBoost-friendly form."""
    val = str(value or "").strip()
    if val.lower() in ("", "none", "null"):
        return None
    if val.lower() == "balanced":
        return "Balanced"
    if val.lower() == "sqrtbalanced":
        return "SqrtBalanced"
    return val


def _detect_torch_devices():
    devices = [("CPU", "cpu"), ("Auto (best available)", "auto")]
    try:
        import torch
        if torch.cuda.is_available():
            for i in range(torch.cuda.device_count()):
                devices.append(
                    (f"GPU {i}: {torch.cuda.get_device_name(i)}", f"cuda:{i}")
                )
    except Exception:
        pass
    return devices


def _default_device():
    """Return ``'cuda:0'`` when an NVIDIA GPU is present, otherwise ``'auto'``."""
    try:
        import torch
        if torch.cuda.is_available():
            return "cuda:0"
    except Exception:
        pass
    return "auto"


def _set_torch_device(device_str):
    if not device_str or device_str in ("auto", "cpu"):
        return
    try:
        import torch
        if device_str.startswith("cuda:"):
            idx = int(device_str.split(":")[1])
            torch.cuda.set_device(idx)
            print(
                f"ConvPaint: Selected GPU {idx} ({torch.cuda.get_device_name(idx)})"
            )
    except Exception as e:
        print(f"\u26a0\ufe0f Could not set CUDA device '{device_str}': {e}")


def _filter_merge_types(cell_types):
    """Remove derived output types (*_merged, *_grouped). Keep *_N_multicolor."""
    from behav3d.core.metadata import is_combined_multicolor_celltype
    return [ct for ct in cell_types
            if not is_combined_multicolor_celltype(ct)]


# ---------------------------------------------------------------------------
# Disk paths for cached widget previews.
# Per cell type we still keep one probability map and one predicted-labels
# zarr (the latter caches the *post-processed* instance segmentation), so the
# per-tab live preview can repaint without re-extracting features.
# ---------------------------------------------------------------------------

def _predicted_labels_path(pixel_class_outdir, cell_type):
    if not pixel_class_outdir:
        return None
    fname = f"ConvPaintModel_{cell_type.capitalize()}_PredictedLabels.zarr"
    return Path(pixel_class_outdir) / fname


def _probability_map_path(pixel_class_outdir, cell_type):
    if not pixel_class_outdir:
        return None
    fname = f"ConvPaintModel_{cell_type.capitalize()}_ProbabilityMap.zarr"
    return Path(pixel_class_outdir) / fname


def _death_model_path(pixel_class_outdir):
    """Death model retains its own binary classifier (separate from the unified one)."""
    return Path(pixel_class_outdir) / "ConvPaintModel_Dead.pkl"


def _unique_valid_dead_channels(metadata, n_channels):
    """Return unique in-range dead-channel indices from metadata."""
    if metadata is None or "dead_channel" not in metadata.columns:
        return []
    channels = []
    for value in metadata["dead_channel"].dropna().tolist():
        try:
            idx = int(value)
        except (TypeError, ValueError):
            continue
        if 0 <= idx < int(n_channels) and idx not in channels:
            channels.append(idx)
    return channels


def _derive_convpaint_input_channels(metadata, n_channels):
    """Derive model-specific channel subsets from metadata.

    We only create channel-specific models when the metadata has one consistent
    dead-channel index. If channel configuration is ambiguous, return all
    channels so behavior remains identical to the legacy all-channel training.
    """
    n_channels = int(n_channels)
    all_channels = list(range(n_channels))
    dead_channels = _unique_valid_dead_channels(metadata, n_channels)
    if len(dead_channels) != 1:
        return all_channels, all_channels

    dead_channel = dead_channels[0]
    unified_channels = [idx for idx in all_channels if idx != dead_channel]
    if not unified_channels:
        unified_channels = all_channels
    return unified_channels, [dead_channel]


def _slice_image_channels(image, channel_indices):
    """Slice a ``(C, Z, Y, X)`` image by saved input channels with safe fallback."""
    if not channel_indices:
        return image
    valid = []
    n_channels = int(image.shape[0])
    for idx in channel_indices:
        try:
            idx = int(idx)
        except (TypeError, ValueError):
            continue
        if 0 <= idx < n_channels and idx not in valid:
            valid.append(idx)
    if not valid:
        return image
    return image[valid]


def _save_preview_array(path, data):
    """Persist a preview array to zarr, replacing the old preview if needed."""
    path = Path(path)
    if path.exists():
        shutil.rmtree(path)
    save_as_zarr(np.asarray(data), path)


def _reorder_convpaint_layers(viewer, all_cell_types, has_death=False):
    """Stable layer ordering for the unified-mode viewer."""
    layer_names = [layer.name for layer in viewer.layers]
    channel_names = [name for name in layer_names if name.startswith("Channel ")]
    ordered_names = list(channel_names)

    ordered_names.append(UNIFIED_LABELS_LAYER_NAME)
    ordered_names.append(UNIFIED_PREDICTED_LAYER_NAME)

    for cell_type in all_cell_types:
        ordered_names.append(_probability_layer_name(cell_type))
        ordered_names.append(_segments_layer_name(cell_type))

    if has_death:
        ordered_names.append(DEAD_LABELS_LAYER_NAME)
        ordered_names.append(DEAD_PREDICTED_LAYER_NAME)

    target_names = [name for name in ordered_names if name in layer_names]
    trailing_names = [name for name in layer_names if name not in target_names]
    final_order = target_names + trailing_names

    for target_idx, name in enumerate(final_order):
        current_idx = viewer.layers.index(name)
        if current_idx != target_idx:
            viewer.layers.move(current_idx, target_idx)


def _napari_label_color(layer, label_idx):
    """Return the napari labels-layer color for ``label_idx`` as a QColor.

    Falls back to a neutral gray if napari's API does not expose the color.
    """
    try:
        rgba = None
        # Modern napari (>=0.4.16) exposes get_color
        if hasattr(layer, "get_color"):
            rgba = layer.get_color(int(label_idx))
        # Older napari: use colormap.map directly
        if rgba is None and hasattr(layer, "colormap"):
            cm = layer.colormap
            if hasattr(cm, "map"):
                rgba = np.asarray(cm.map(np.array([int(label_idx)]))).reshape(-1)
        if rgba is None:
            return QColor(160, 160, 160)
        rgba = np.asarray(rgba).reshape(-1)
        if rgba.size < 3:
            return QColor(160, 160, 160)
        r, g, b = (int(round(float(rgba[i]) * 255)) for i in range(3))
        a = int(round(float(rgba[3]) * 255)) if rgba.size >= 4 else 255
        return QColor(r, g, b, a)
    except Exception:
        return QColor(160, 160, 160)


# ---------------------------------------------------------------------------
# Image loading (unchanged from per-cell-type version).
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
        t_indices = (
            [0] if n_sel <= 1
            else sorted(set(np.round(np.linspace(0, n_tp - 1, n_sel)).astype(int)))
        )
        print(
            f"  {sample_name}: selected {len(t_indices)} timepoints {list(t_indices)}"
        )
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
    del stack
    gc.collect()
    cached = load_zarr(image_outpath)
    all_images = [np.asarray(cached[:, t]) for t in range(cached.shape[1])]
    return all_images, pixel_class_outdir, has_death, all_cell_types


# ---------------------------------------------------------------------------
# Annotation Legend tab — first tab in the QTabWidget.
# ---------------------------------------------------------------------------

class AnnotationLegendTab(QWidget):
    """Legend table for the unified User Provided Labels layer.

    Each row shows ``[Index | Color swatch | Cell type | Voxel count]``.
    Clicking a row sets the labels layer's ``selected_label`` to that index
    and switches the layer to paint mode, so the user can immediately start
    annotating that class.
    """

    def __init__(self, viewer, label_map, has_death=False, parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.label_map = label_map
        self.has_death = has_death

        layout = QVBoxLayout()
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        intro = QLabel(
            "<b>Annotation legend.</b> Click a row to switch the active paint "
            "label on the unified <i>User Provided Labels</i> layer. "
            "Background must always be label <b>1</b>."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        self.table = QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(
            ["Idx", "Color", "Cell type", "Annotated voxels"]
        )
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.Stretch)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.table.cellClicked.connect(self._on_cell_clicked)
        layout.addWidget(self.table)

        button_row = QHBoxLayout()
        self.btn_refresh = QtPushButton("Refresh counts & colors")
        self.btn_refresh.clicked.connect(self.refresh)
        button_row.addWidget(self.btn_refresh)
        self.note = QLabel(
            "<i>Death is a separate binary classifier — annotate it on the "
            "<b>User Provided Labels (Dead)</b> layer.</i>"
        )
        self.note.setWordWrap(True)
        if not has_death:
            self.note.hide()
        layout.addWidget(self.note)
        layout.addLayout(button_row)
        layout.addStretch()
        self.setLayout(layout)

        self.refresh()

    # ── Refresh ─────────────────────────────────────────────────────

    def refresh(self):
        """Rebuild the legend table from the current label map + label layer."""
        layer = self.viewer.layers.get(UNIFIED_LABELS_LAYER_NAME) \
            if hasattr(self.viewer.layers, "get") else None
        if layer is None:
            for lyr in self.viewer.layers:
                if lyr.name == UNIFIED_LABELS_LAYER_NAME:
                    layer = lyr
                    break

        # Voxel-count tally over the whole annotation volume.
        counts = {}
        if layer is not None:
            data = np.asarray(layer.data)
            if data.size:
                values, freq = np.unique(data, return_counts=True)
                counts = {int(v): int(c) for v, c in zip(values, freq) if int(v) > 0}

        ordered_cts = cell_type_order(self.label_map)
        bg_label = int(self.label_map.get("background_label", BACKGROUND_LABEL))
        rows = [(bg_label, "background")]
        for ct in ordered_cts:
            rows.append((self.label_map["celltype_to_label"][ct], ct))

        self.table.setRowCount(len(rows))
        for row, (idx, name) in enumerate(rows):
            # Index
            idx_item = QTableWidgetItem(str(idx))
            idx_item.setData(Qt.UserRole, int(idx))
            idx_item.setTextAlignment(Qt.AlignCenter)
            self.table.setItem(row, 0, idx_item)

            # Color swatch
            swatch = QWidget()
            swatch.setFixedSize(24, 24)
            color = _napari_label_color(layer, idx) if layer is not None else QColor(160, 160, 160)
            swatch.setStyleSheet(
                f"background-color: rgba({color.red()}, {color.green()}, "
                f"{color.blue()}, {color.alpha()}); border: 1px solid #444;"
            )
            self.table.setCellWidget(row, 1, swatch)

            # Cell type
            name_item = QTableWidgetItem(name)
            name_item.setData(Qt.UserRole, int(idx))
            self.table.setItem(row, 2, name_item)

            # Voxel count
            count = counts.get(int(idx), 0)
            count_item = QTableWidgetItem(f"{count:,}")
            count_item.setData(Qt.UserRole, int(idx))
            count_item.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
            self.table.setItem(row, 3, count_item)

        self.table.resizeColumnsToContents()
        # Re-stretch the cell-type column after resizeColumnsToContents reset it.
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)

    # ── Click handler ───────────────────────────────────────────────

    def _on_cell_clicked(self, row, _col):
        idx_item = self.table.item(row, 0)
        if idx_item is None:
            return
        try:
            label_idx = int(idx_item.data(Qt.UserRole))
        except (TypeError, ValueError):
            return

        layer = None
        for lyr in self.viewer.layers:
            if lyr.name == UNIFIED_LABELS_LAYER_NAME:
                layer = lyr
                break
        if layer is None:
            return

        layer.selected_label = label_idx
        try:
            layer.mode = "paint"
        except Exception:
            pass
        try:
            self.viewer.layers.selection.active = layer
        except Exception:
            try:
                self.viewer.layers.selection = {layer}
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Per-cell-type segmentation tab (instance-segmentation parameters only).
# Identical in shape to the previous implementation; the *callback* now runs
# post-hoc work only.
# ---------------------------------------------------------------------------

class CellTypeConvPaintTab(QWidget):
    """One tab per cell type with strategy-specific spinners + preview button."""

    def __init__(self, cell_type, strategy, initial_params=None,
                 on_params_changed=None, run_preview_callback=None,
                 parent=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.strategy = _normalize_strategy(strategy)
        self._on_params_changed = on_params_changed
        self._run_preview_callback = run_preview_callback
        self._is_dead = (cell_type == "dead")

        ip = initial_params or {}

        layout = QVBoxLayout()
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        if self._is_dead:
            layout.addWidget(QLabel(
                "<i>Death classifier is a separate binary ConvPaint model. "
                "Use the <b>User Provided Labels (Dead)</b> layer to annotate it. "
                "No instance-segmentation parameters required.</i>"
            ))
            self.btn_preview = QtPushButton("Run preview (death mask)")
            self.btn_preview.clicked.connect(self._on_run_preview)
            layout.addWidget(self.btn_preview)
            layout.addStretch()
            self.setLayout(layout)
            self.edt_threshold_spin = None
            self.opening_nr_pixels_spin = None
            self.segment_size_min_spin = None
            self.fill_holes_cb = None
            self.prob_mask_threshold_spin = None
            self.prob_seed_threshold_spin = None
            self.peak_min_distance_spin = None
            self.peak_min_ratio_spin = None
            return

        strategy_lbl = QLabel(f"<i>Strategy:</i> <b>{self.strategy}</b>")
        strategy_lbl.setWordWrap(True)
        layout.addWidget(strategy_lbl)

        self.edt_threshold_spin = None
        self.opening_nr_pixels_spin = None
        self.segment_size_min_spin = None
        self.fill_holes_cb = None
        self.prob_mask_threshold_spin = None
        self.prob_seed_threshold_spin = None
        self.peak_min_distance_spin = None
        self.peak_min_ratio_spin = None

        if self.strategy == STRATEGY_PROB:
            row1 = QHBoxLayout()
            row1.addWidget(QLabel("Mask threshold:"))
            self.prob_mask_threshold_spin = QDoubleSpinBox()
            self.prob_mask_threshold_spin.setRange(0.0, 1.0)
            self.prob_mask_threshold_spin.setSingleStep(0.05)
            self.prob_mask_threshold_spin.setDecimals(2)
            self.prob_mask_threshold_spin.setValue(
                float(ip.get(f"{cell_type}_prob_mask_threshold", 0.5))
            )
            row1.addWidget(self.prob_mask_threshold_spin)

            row1.addWidget(QLabel("Seed threshold:"))
            self.prob_seed_threshold_spin = QDoubleSpinBox()
            self.prob_seed_threshold_spin.setRange(0.0, 1.0)
            self.prob_seed_threshold_spin.setSingleStep(0.05)
            self.prob_seed_threshold_spin.setDecimals(2)
            self.prob_seed_threshold_spin.setValue(
                float(ip.get(f"{cell_type}_prob_seed_threshold", 0.8))
            )
            row1.addWidget(self.prob_seed_threshold_spin)
            layout.addLayout(row1)

            row2 = QHBoxLayout()
            row2.addWidget(QLabel("Opening px:"))
            self.opening_nr_pixels_spin = QSpinBox()
            self.opening_nr_pixels_spin.setRange(0, 50)
            self.opening_nr_pixels_spin.setValue(
                int(ip.get(f"{cell_type}_opening_nr_pixels", 0))
            )
            row2.addWidget(self.opening_nr_pixels_spin)

            row2.addWidget(QLabel("Min size:"))
            self.segment_size_min_spin = QSpinBox()
            self.segment_size_min_spin.setRange(0, 100000)
            self.segment_size_min_spin.setSingleStep(10)
            self.segment_size_min_spin.setValue(
                int(ip.get(f"{cell_type}_segment_size_min", 10))
            )
            row2.addWidget(self.segment_size_min_spin)
            layout.addLayout(row2)
        else:
            row1 = QHBoxLayout()
            row1.addWidget(QLabel("EDT threshold:"))
            self.edt_threshold_spin = QDoubleSpinBox()
            self.edt_threshold_spin.setRange(0.0, 50.0)
            self.edt_threshold_spin.setSingleStep(0.5)
            self.edt_threshold_spin.setDecimals(2)
            self.edt_threshold_spin.setValue(
                float(ip.get(f"{cell_type}_edt_threshold", 1.0))
            )
            row1.addWidget(self.edt_threshold_spin)

            row1.addWidget(QLabel("Opening px:"))
            self.opening_nr_pixels_spin = QSpinBox()
            self.opening_nr_pixels_spin.setRange(0, 50)
            self.opening_nr_pixels_spin.setValue(
                int(ip.get(f"{cell_type}_opening_nr_pixels", 0))
            )
            row1.addWidget(self.opening_nr_pixels_spin)
            layout.addLayout(row1)

            row2 = QHBoxLayout()
            row2.addWidget(QLabel("Min size:"))
            self.segment_size_min_spin = QSpinBox()
            self.segment_size_min_spin.setRange(0, 100000)
            self.segment_size_min_spin.setSingleStep(10)
            self.segment_size_min_spin.setValue(
                int(ip.get(f"{cell_type}_segment_size_min", 10))
            )
            row2.addWidget(self.segment_size_min_spin)

            self.fill_holes_cb = QCheckBox("Fill holes")
            self.fill_holes_cb.setChecked(
                bool(ip.get(f"{cell_type}_fill_holes", True))
            )
            row2.addWidget(self.fill_holes_cb)
            layout.addLayout(row2)

            if self.strategy == STRATEGY_PEAK_EDT:
                row3 = QHBoxLayout()
                row3.addWidget(QLabel("Peak min dist:"))
                self.peak_min_distance_spin = QDoubleSpinBox()
                self.peak_min_distance_spin.setRange(0.0, 50.0)
                self.peak_min_distance_spin.setSingleStep(0.5)
                self.peak_min_distance_spin.setDecimals(2)
                self.peak_min_distance_spin.setValue(
                    float(ip.get(f"{cell_type}_peak_min_distance", 0.0))
                )
                self.peak_min_distance_spin.setToolTip(
                    "Minimum distance (µm) between local EDT peaks used as watershed seeds"
                )
                row3.addWidget(self.peak_min_distance_spin)

                row3.addWidget(QLabel("Peak min ratio:"))
                self.peak_min_ratio_spin = QDoubleSpinBox()
                self.peak_min_ratio_spin.setRange(0.0, 1.0)
                self.peak_min_ratio_spin.setSingleStep(0.05)
                self.peak_min_ratio_spin.setDecimals(2)
                self.peak_min_ratio_spin.setValue(
                    float(ip.get(f"{cell_type}_peak_min_ratio", 0.35))
                )
                self.peak_min_ratio_spin.setToolTip(
                    "Minimum EDT peak height as a fraction of the local maximum (0–1)"
                )
                row3.addWidget(self.peak_min_ratio_spin)
                layout.addLayout(row3)

        hint = QLabel(
            "<i>Re-runs only post-processing on the cached classifier output "
            "(no feature extraction).</i>"
        )
        hint.setStyleSheet("color: #666;")
        hint.setWordWrap(True)
        layout.addWidget(hint)

        self.btn_preview = QtPushButton("Run instance segmentation")
        self.btn_preview.setStyleSheet(
            "background-color: #6A1B9A; color: white; font-weight: bold; padding: 5px 10px;"
        )
        self.btn_preview.clicked.connect(self._on_run_preview)
        layout.addWidget(self.btn_preview)
        layout.addStretch()
        self.setLayout(layout)

        for spin in [
            self.edt_threshold_spin, self.opening_nr_pixels_spin,
            self.segment_size_min_spin,
            self.prob_mask_threshold_spin, self.prob_seed_threshold_spin,
            self.peak_min_distance_spin, self.peak_min_ratio_spin,
        ]:
            if spin is not None:
                spin.valueChanged.connect(self._emit_params_changed)
        if self.fill_holes_cb is not None:
            self.fill_holes_cb.stateChanged.connect(self._emit_params_changed)

    def _on_run_preview(self):
        if callable(self._run_preview_callback):
            self._run_preview_callback(self.cell_type)

    def _emit_params_changed(self, *_args):
        if callable(self._on_params_changed):
            self._on_params_changed(self.collect_params())

    def collect_params(self):
        """Return a flat dict of this tab's strategy-relevant params."""
        ct = self.cell_type
        params = {}
        if self.edt_threshold_spin is not None:
            params[f"{ct}_edt_threshold"] = float(self.edt_threshold_spin.value())
        if self.opening_nr_pixels_spin is not None:
            params[f"{ct}_opening_nr_pixels"] = int(self.opening_nr_pixels_spin.value())
        if self.segment_size_min_spin is not None:
            params[f"{ct}_segment_size_min"] = int(self.segment_size_min_spin.value())
        if self.fill_holes_cb is not None:
            params[f"{ct}_fill_holes"] = bool(self.fill_holes_cb.isChecked())
        if self.prob_mask_threshold_spin is not None:
            params[f"{ct}_prob_mask_threshold"] = float(
                self.prob_mask_threshold_spin.value()
            )
        if self.prob_seed_threshold_spin is not None:
            params[f"{ct}_prob_seed_threshold"] = float(
                self.prob_seed_threshold_spin.value()
            )
        if self.peak_min_distance_spin is not None:
            params[f"{ct}_peak_min_distance"] = float(
                self.peak_min_distance_spin.value()
            )
        if self.peak_min_ratio_spin is not None:
            params[f"{ct}_peak_min_ratio"] = float(
                self.peak_min_ratio_spin.value()
            )
        return params


# ---------------------------------------------------------------------------
# Main widget
# ---------------------------------------------------------------------------

class ConvPaintTrainingWidget(QWidget):
    def __init__(self, viewer, all_images, pixel_class_outdir,
                 all_cell_types, has_death, initial_params=None,
                 on_params_changed=None, convpaint_strategy=None,
                 unified_input_channels=None, death_input_channels=None,
                 parent=None):
        super().__init__(parent)
        self.viewer = viewer
        self.all_images = all_images
        self.pixel_class_outdir = Path(pixel_class_outdir)
        self.all_cell_types = list(all_cell_types)
        self.has_death = has_death
        self.unified_input_channels = list(unified_input_channels or [])
        self.death_input_channels = list(death_input_channels or [])
        self._on_params_changed = on_params_changed
        self.ip = initial_params or {}
        self.convpaint_strategy = _normalize_strategy(
            convpaint_strategy or self.ip.get("convpaint_strategy", DEFAULT_STRATEGY)
        )

        # Active label map (always built fresh from the current cell-type set).
        self.label_map = build_label_map(self.all_cell_types)

        # Per-cell-type segmentation tabs + (optional) Dead tab.
        self._tab_cell_types = list(self.all_cell_types)
        if self.has_death:
            self._tab_cell_types.append("dead")

        self.tabs = {}
        self.legend_tab = None
        self._build_ui()

    # ── UI construction ─────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout()
        root.setContentsMargins(4, 4, 4, 4)
        root.addWidget(QLabel("<h3>ConvPaint Training (unified)</h3>"))

        legend_intro = QLabel(
            "<i>One classifier predicts every cell type at once. "
            "Annotate every class on the single <b>User Provided Labels</b> layer "
            "(see legend tab below for index ↔ cell type).</i>"
        )
        legend_intro.setWordWrap(True)
        legend_intro.setStyleSheet("color: #444; margin-bottom: 4px;")
        root.addWidget(legend_intro)

        # Feature Extractor (global)
        fe_group = QGroupBox("Feature Extractor")
        fe_layout = QVBoxLayout()
        r = QHBoxLayout()
        r.addWidget(QLabel("Model:"))
        self.fe_combo = QComboBox()
        for label, alias in CONVPAINT_FE_OPTIONS:
            self.fe_combo.addItem(label, alias)
        idx = self.fe_combo.findData(self.ip.get("convpaint_fe_alias", "vgg"))
        if idx >= 0:
            self.fe_combo.setCurrentIndex(idx)
        _compact_combo(self.fe_combo, min_chars=10)
        r.addWidget(self.fe_combo)
        fe_layout.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Channels:"))
        self.channel_mode_combo = QComboBox()
        for label, mode in CHANNEL_MODE_OPTIONS:
            self.channel_mode_combo.addItem(label, mode)
        idx = self.channel_mode_combo.findData(
            self.ip.get("convpaint_channel_mode", "multi")
        )
        if idx >= 0:
            self.channel_mode_combo.setCurrentIndex(idx)
        _compact_combo(self.channel_mode_combo, min_chars=8)
        r.addWidget(self.channel_mode_combo)
        fe_layout.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Normalize:"))
        self.normalize_combo = QComboBox()
        for label, val in NORMALIZE_OPTIONS:
            self.normalize_combo.addItem(label, val)
        idx = self.normalize_combo.findData(int(self.ip.get("convpaint_normalize", 2)))
        if idx >= 0:
            self.normalize_combo.setCurrentIndex(idx)
        _compact_combo(self.normalize_combo, min_chars=8)
        r.addWidget(self.normalize_combo)
        fe_layout.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Downsample:"))
        self.downsample_spin = QSpinBox()
        self.downsample_spin.setRange(-4, 16)
        self.downsample_spin.setValue(int(self.ip.get("convpaint_image_downsample", 1)))
        r.addWidget(self.downsample_spin)
        r.addWidget(QLabel("Smoothen:"))
        self.smooth_spin = QSpinBox()
        self.smooth_spin.setRange(0, 20)
        self.smooth_spin.setValue(int(self.ip.get("convpaint_seg_smoothening", 1)))
        r.addWidget(self.smooth_spin)
        fe_layout.addLayout(r)
        fe_group.setLayout(fe_layout)
        root.addWidget(fe_group)

        # Classifier (global)
        clf_group = QGroupBox("Classifier (CatBoost)")
        clf_layout = QVBoxLayout()
        r = QHBoxLayout()
        r.addWidget(QLabel("Iterations:"))
        self.iterations_spin = QSpinBox()
        self.iterations_spin.setRange(10, 2000)
        self.iterations_spin.setValue(int(self.ip.get("convpaint_clf_iterations", 100)))
        r.addWidget(self.iterations_spin)
        r.addWidget(QLabel("Depth:"))
        self.depth_spin = QSpinBox()
        self.depth_spin.setRange(1, 16)
        self.depth_spin.setValue(int(self.ip.get("convpaint_clf_depth", 5)))
        r.addWidget(self.depth_spin)
        clf_layout.addLayout(r)
        r = QHBoxLayout()
        r.addWidget(QLabel("Learning rate:"))
        self.lr_spin = QDoubleSpinBox()
        self.lr_spin.setRange(0.001, 1.0)
        self.lr_spin.setSingleStep(0.01)
        self.lr_spin.setDecimals(3)
        self.lr_spin.setValue(float(self.ip.get("convpaint_clf_learning_rate", 0.1)))
        r.addWidget(self.lr_spin)
        clf_layout.addLayout(r)
        r = QHBoxLayout()
        r.addWidget(QLabel("Class weights:"))
        self.class_weights_combo = QComboBox()
        for label, val in CLASS_WEIGHTS_OPTIONS:
            self.class_weights_combo.addItem(label, val)
        saved_cw = str(self.ip.get("convpaint_class_weights", DEFAULT_CLASS_WEIGHTS))
        idx = self.class_weights_combo.findData(saved_cw)
        if idx < 0:
            idx = self.class_weights_combo.findData(DEFAULT_CLASS_WEIGHTS)
        if idx >= 0:
            self.class_weights_combo.setCurrentIndex(idx)
        _compact_combo(self.class_weights_combo, min_chars=10)
        r.addWidget(self.class_weights_combo)
        clf_layout.addLayout(r)
        clf_group.setLayout(clf_layout)
        root.addWidget(clf_group)

        # Device (global)
        r = QHBoxLayout()
        r.addWidget(QLabel("Device:"))
        self.device_combo = QComboBox()
        for label, dev_str in _detect_torch_devices():
            self.device_combo.addItem(label, dev_str)
        idx = self.device_combo.findData(
            self.ip.get("convpaint_device", _default_device())
        )
        if idx >= 0:
            self.device_combo.setCurrentIndex(idx)
        _compact_combo(self.device_combo, min_chars=10)
        r.addWidget(self.device_combo)
        root.addLayout(r)

        # Tiling options (global)
        tile_row = QHBoxLayout()
        self.tile_annotations_cb = QCheckBox("Tile annotations")
        self.tile_annotations_cb.setToolTip(
            "Only extract features inside annotated regions during training "
            "(faster on sparse annotations)."
        )
        self.tile_annotations_cb.setChecked(
            bool(self.ip.get("convpaint_tile_annotations", False))
        )
        tile_row.addWidget(self.tile_annotations_cb)
        self.tile_image_cb = QCheckBox("Tile image")
        self.tile_image_cb.setToolTip(
            "Extract features in tiles during prediction "
            "(reduces peak memory on large images)."
        )
        self.tile_image_cb.setChecked(
            bool(self.ip.get("convpaint_tile_image", False))
        )
        tile_row.addWidget(self.tile_image_cb)
        self.use_dask_cb = QCheckBox("Use Dask")
        self.use_dask_cb.setToolTip(
            "Parallelize prediction across tiles with Dask. Only takes effect "
            "when 'Tile image' is also enabled (per the napari-convpaint API). "
            "Beta feature \u2014 may not be fully optimized yet."
        )
        self.use_dask_cb.setChecked(
            bool(self.ip.get("convpaint_use_dask", False))
        )
        tile_row.addWidget(self.use_dask_cb)
        root.addLayout(tile_row)

        # Tabs: legend first, then per-cell-type segmentation params, then Dead.
        seg_group = QGroupBox("Annotation & Segmentation")
        seg_layout = QVBoxLayout()
        self.tab_widget = QTabWidget()

        self.legend_tab = AnnotationLegendTab(
            self.viewer, self.label_map, has_death=self.has_death,
        )
        self.tab_widget.addTab(self.legend_tab, "Legend")

        for ct in self._tab_cell_types:
            tab = CellTypeConvPaintTab(
                cell_type=ct,
                strategy=self.convpaint_strategy,
                initial_params=self.ip,
                on_params_changed=self._on_tab_params_changed,
                run_preview_callback=self._run_instance_preview,
            )
            self.tabs[ct] = tab
            self.tab_widget.addTab(tab, ct.capitalize())
        seg_layout.addWidget(self.tab_widget)
        seg_group.setLayout(seg_layout)
        root.addWidget(seg_group)

        # Train button (one button, trains the unified model + dead model if present)
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        root.addWidget(line)

        self.btn_train = QtPushButton("\u25B6 Train Classifier")
        self.btn_train.setStyleSheet(
            "background-color: #2e7d32; color: white; font-weight: bold; padding: 6px 12px;"
        )
        self.btn_train.setToolTip(
            "Train the unified multi-class classifier on User Provided Labels "
            "(and the binary death classifier if applicable)."
        )
        self.btn_train.clicked.connect(self._on_train_clicked)
        root.addWidget(self.btn_train)

        # Re-run preview (uses cached unified labels / probability map, no FE)
        self.btn_rerun_preview = QtPushButton("\U0001F501 Re-run preview")
        self.btn_rerun_preview.setStyleSheet(
            "background-color: #7B1FA2; color: white; font-weight: bold; padding: 6px 12px;"
        )
        self.btn_rerun_preview.setToolTip(
            "Re-run the post-processing for the current tab using the cached "
            "classifier output (no feature extraction)."
        )
        self.btn_rerun_preview.clicked.connect(self._on_rerun_preview)
        root.addWidget(self.btn_rerun_preview)

        # Save labels
        self.btn_save_labels = QtPushButton("\U0001F4BE Save User Labels")
        self.btn_save_labels.setStyleSheet(
            "background-color: #FF9800; color: white; font-weight: bold; padding: 6px;"
        )
        self.btn_save_labels.clicked.connect(self._on_save_labels)
        root.addWidget(self.btn_save_labels)

        hint = QLabel(
            "<i>Labels: 1=background, 2..N+1=one per cell type "
            "(see Legend tab). Dead labels live on a separate layer.</i>"
        )
        hint.setStyleSheet("color: #666; margin-top: 4px;")
        hint.setWordWrap(True)
        root.addWidget(hint)

        self.log_box = QPlainTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(150)
        root.addWidget(self.log_box)
        root.addStretch()

        scroll_content = QWidget()
        scroll_content.setLayout(root)
        scroll_content.setMaximumWidth(420)

        scroll = QScrollArea()
        scroll.setWidget(scroll_content)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setMaximumWidth(420)

        outer = QVBoxLayout()
        outer.setContentsMargins(0, 0, 0, 0)
        outer.addWidget(scroll)
        self.setLayout(outer)
        self.setMaximumWidth(420)

    # ── Callbacks ───────────────────────────────────────────────────

    def _on_tab_params_changed(self, params):
        self.ip.update(params)
        if self._on_params_changed:
            self._on_params_changed(params)

    def _on_save_labels(self):
        self.save_user_labels()

    def _log(self, msg):
        self.log_box.appendPlainText(str(msg))
        QApplication.processEvents()

    # ── Model building ──────────────────────────────────────────────

    def _build_model(self):
        from napari_convpaint import ConvpaintModel
        model = ConvpaintModel(self.fe_combo.currentData())
        params = dict(
            channel_mode=self.channel_mode_combo.currentData(),
            normalize=self.normalize_combo.currentData(),
            image_downsample=self.downsample_spin.value(),
            seg_smoothening=self.smooth_spin.value(),
            clf_iterations=self.iterations_spin.value(),
            clf_learning_rate=self.lr_spin.value(),
            clf_depth=self.depth_spin.value(),
            tile_annotations=self.tile_annotations_cb.isChecked(),
            tile_image=self.tile_image_cb.isChecked(),
            use_dask=self.use_dask_cb.isChecked(),
        )
        model.set_params(**params)

        # use_dask only takes effect when tile_image is also enabled.
        if self.use_dask_cb.isChecked() and not self.tile_image_cb.isChecked():
            self._log(
                "\u2139\ufe0f 'Use Dask' is enabled but 'Tile image' is not \u2014 "
                "Dask only parallelizes tiled predictions, so it will have no "
                "effect until 'Tile image' is enabled."
            )

        # Class weights are forwarded only when the underlying ConvPaint
        # version supports the parameter; otherwise we log and continue.
        cw = _normalize_class_weights(self.class_weights_combo.currentData())
        if cw is not None:
            try:
                model.set_params(clf_class_weights=cw)
            except Exception as exc:
                self._log(
                    "\u26a0\ufe0f Could not apply class weights "
                    f"({cw!r}); ConvPaint may be too old. ({exc})"
                )
        return model, self._fe_device()

    def _fe_device(self):
        device = self.device_combo.currentData()
        if device and device not in ("auto", "cpu"):
            _set_torch_device(device)
            return "gpu"
        if device == "cpu":
            return "cpu"
        return "auto"

    def _get_global_config(self):
        return {
            "convpaint_fe_alias": self.fe_combo.currentData(),
            "convpaint_channel_mode": self.channel_mode_combo.currentData(),
            "convpaint_normalize": self.normalize_combo.currentData(),
            "convpaint_image_downsample": self.downsample_spin.value(),
            "convpaint_seg_smoothening": self.smooth_spin.value(),
            "convpaint_clf_iterations": self.iterations_spin.value(),
            "convpaint_clf_depth": self.depth_spin.value(),
            "convpaint_clf_learning_rate": self.lr_spin.value(),
            "convpaint_class_weights": _normalize_class_weights(
                self.class_weights_combo.currentData()
            ) or "None",
            "convpaint_device": self.device_combo.currentData(),
            "convpaint_strategy": self.convpaint_strategy,
            "convpaint_tile_annotations": self.tile_annotations_cb.isChecked(),
            "convpaint_tile_image": self.tile_image_cb.isChecked(),
            "convpaint_use_dask": self.use_dask_cb.isChecked(),
        }

    def _persist_all_params(self):
        if not callable(self._on_params_changed):
            return
        params = dict(self._get_global_config())
        for tab in self.tabs.values():
            params.update(tab.collect_params())
        self._on_params_changed(params)

    # ── Training ────────────────────────────────────────────────────

    def _on_train_clicked(self):
        self._log("Auto-saving user labels before training...")
        self.save_user_labels()
        self._train_unified_classifier()
        if self.has_death:
            self._train_death_classifier()
        self._persist_all_params()
        if self.legend_tab is not None:
            self.legend_tab.refresh()

    def _train_unified_classifier(self):
        """Train the single multi-class ConvPaint model."""
        layer = None
        for lyr in self.viewer.layers:
            if lyr.name == UNIFIED_LABELS_LAYER_NAME:
                layer = lyr
                break
        if layer is None:
            self._log(f"\u274c Label layer '{UNIFIED_LABELS_LAYER_NAME}' not found!")
            return

        annotations = np.asarray(layer.data)

        # Determine which cell types actually have annotations (≥ 1 voxel painted).
        unique_vals = set(int(v) for v in np.unique(annotations) if int(v) > 0)
        if BACKGROUND_LABEL not in unique_vals:
            self._log(
                "\u26a0\ufe0f No background labels (1) painted. "
                "Training will likely produce poor results — paint some background."
            )

        active_cell_types = []
        skipped_cell_types = []
        for ct in self.all_cell_types:
            idx = self.label_map["celltype_to_label"][ct]
            if int(idx) in unique_vals:
                active_cell_types.append(ct)
            else:
                skipped_cell_types.append(ct)
        if skipped_cell_types:
            self._log(
                "\u26a0\ufe0f No annotations for: "
                f"{', '.join(skipped_cell_types)}. "
                "Those cell types will not be predicted by this model."
            )
        if not active_cell_types:
            self._log("\u274c No annotated cell types — aborting training.")
            return

        # Build a label map restricted to *trained* cell types so that downstream
        # decoders never request a class index that the model didn't learn.
        active_label_map = build_label_map(active_cell_types)

        # Remap annotations to dense indices [1..len(active_cell_types)+1].
        # All non-active foreground labels get demoted to 0 (unpainted).
        remap = {0: 0, BACKGROUND_LABEL: BACKGROUND_LABEL}
        for ct in active_cell_types:
            old_idx = int(self.label_map["celltype_to_label"][ct])
            new_idx = int(active_label_map["celltype_to_label"][ct])
            remap[old_idx] = new_idx

        remapped = np.zeros_like(annotations)
        for old_v, new_v in remap.items():
            if old_v == 0:
                continue
            remapped[annotations == old_v] = new_v

        # Build per-frame training pairs.
        train_images, train_annots = [], []
        for i, img in enumerate(self.all_images):
            frame_annots = remapped[i]
            if np.any(frame_annots > 0):
                train_images.append(
                    _slice_image_channels(img, self.unified_input_channels)
                )
                train_annots.append(frame_annots)
        if not train_images:
            self._log("\u26a0\ufe0f No annotated frames — aborting training.")
            return

        self._log(
            f"Training unified classifier on {len(train_images)} frames "
            f"({len(active_cell_types)} cell types: {active_cell_types})..."
        )
        self._log(f"Unified model input channels: {self.unified_input_channels}")
        try:
            model, fe_device = self._build_model()
            model.train(
                train_images, train_annots,
                fe_use_device=fe_device, clf_use_device=fe_device,
            )

            # Save model + label-map sidecar.
            model_path = unified_model_path(self.pixel_class_outdir)
            map_path = unified_label_map_path(self.pixel_class_outdir)
            model.save(str(model_path))
            save_label_map(map_path, active_label_map)
            save_input_channels(
                unified_input_channels_path(self.pixel_class_outdir),
                self.unified_input_channels,
                model_name="unified",
            )
            self._log(f"\u2705 Saved model: {model_path.name}")
            self._log(f"\u2705 Saved label map: {map_path.name}")

            # Update the live label map so downstream tabs decode correctly.
            self.label_map = active_label_map

            # Single-pass inference on the training stack: produce unified
            # labels + per-cell-type probability maps, cache to disk, refresh
            # the per-cell-type instance previews.
            self._run_unified_inference(model, fe_device, active_cell_types)
        except Exception as e:
            self._log(f"\u274c Unified training failed: {e}")
            traceback.print_exc()

    def _run_unified_inference(self, model, fe_device, active_cell_types):
        """One inference pass on the training stack → caches + previews."""
        self._log("Running unified inference on the training stack...")

        # Multi-class label volume: stack `segment` outputs over T.
        seg_stack = np.stack(
            [
                segment_per_z(
                    model,
                    _slice_image_channels(img, self.unified_input_channels),
                    fe_use_device=fe_device,
                )
                for img in self.all_images
            ],
            axis=0,
        ).astype(np.int16)
        _save_preview_array(
            unified_predicted_labels_path(self.pixel_class_outdir), seg_stack,
        )
        self._set_labels_layer(
            UNIFIED_PREDICTED_LAYER_NAME, seg_stack, visible=False,
        )

        # Per-cell-type probability slice: one inference call per frame, stack
        # along T. Keep only the foreground-class slice for each active cell
        # type (probas axis 0 indexes class labels in ascending order, mapped
        # to the *active* label map).
        if self.convpaint_strategy == STRATEGY_PROB or any(
            self.tabs[ct].prob_mask_threshold_spin is not None
            for ct in active_cell_types if ct in self.tabs
        ):
            self._cache_per_celltype_probabilities(model, fe_device, active_cell_types)
        else:
            # EDT mode does not strictly need probabilities; skip them to save
            # disk + time.
            for ct in active_cell_types:
                p = _probability_map_path(self.pixel_class_outdir, ct)
                if p is not None and p.exists():
                    shutil.rmtree(p)
                lname = _probability_layer_name(ct)
                if lname in self.viewer.layers:
                    self.viewer.layers[lname].visible = False

        # Per-cell-type post-hoc instance segmentation.
        for ct in active_cell_types:
            self._run_instance_preview(ct)

        _reorder_convpaint_layers(
            self.viewer, self.all_cell_types, has_death=self.has_death,
        )
        self._log("\u2705 Unified inference + previews complete.")

    def _cache_per_celltype_probabilities(self, model, fe_device, active_cell_types):
        """Run predict_probas once per frame; slice & cache per cell type."""
        per_ct_stacks = {ct: [] for ct in active_cell_types}
        for img in self.all_images:
            probas = predict_probas_per_z(
                model,
                _slice_image_channels(img, self.unified_input_channels),
                fe_use_device=fe_device,
            )
            probas = np.asarray(probas)
            # probas shape: (n_classes, Z, Y, X). Class indices in the active
            # label map run [1, 2, ..., n_classes]; axis 0 is in the same order.
            for ct in active_cell_types:
                k = int(self.label_map["celltype_to_label"][ct])
                # k is 2..N+1; channel index along axis 0 is k - 1.
                if 0 <= (k - 1) < probas.shape[0]:
                    per_ct_stacks[ct].append(
                        probas[k - 1].astype(np.float32)
                    )
                else:
                    per_ct_stacks[ct].append(
                        np.zeros(probas.shape[1:], dtype=np.float32)
                    )

        for ct, slices in per_ct_stacks.items():
            stack = np.stack(slices, axis=0)
            _save_preview_array(
                _probability_map_path(self.pixel_class_outdir, ct), stack,
            )
            self._set_image_layer(
                _probability_layer_name(ct), stack, visible=False,
            )

    def _train_death_classifier(self):
        """Train the binary death ConvPaint model from its own label layer."""
        layer = None
        for lyr in self.viewer.layers:
            if lyr.name == DEAD_LABELS_LAYER_NAME:
                layer = lyr
                break
        if layer is None:
            self._log(
                f"\u26a0\ufe0f Death layer '{DEAD_LABELS_LAYER_NAME}' not found, "
                "skipping death training."
            )
            return
        annotations = np.asarray(layer.data)
        train_images, train_annots = [], []
        for i, img in enumerate(self.all_images):
            if np.any(annotations[i] > 0):
                train_images.append(
                    _slice_image_channels(img, self.death_input_channels)
                )
                train_annots.append(annotations[i])
        if not train_images:
            self._log("\u26a0\ufe0f No annotations for death — skipping.")
            return
        self._log(f"Training death classifier on {len(train_images)} frames...")
        self._log(f"Death model input channels: {self.death_input_channels}")
        try:
            model, fe_device = self._build_model()
            model.train(
                train_images, train_annots,
                fe_use_device=fe_device, clf_use_device=fe_device,
            )
            mp = _death_model_path(self.pixel_class_outdir)
            model.save(str(mp))
            save_input_channels(
                death_input_channels_path(self.pixel_class_outdir),
                self.death_input_channels,
                model_name="death",
            )
            self._log(f"\u2705 Saved death model: {mp.name}")
            self._preview_death(model, fe_device, self.death_input_channels)
        except Exception as e:
            self._log(f"\u274c Death training failed: {e}")
            traceback.print_exc()

    # ── Instance-segmentation preview (post-hoc, no feature extraction) ──

    def _on_rerun_preview(self):
        idx = self.tab_widget.currentIndex()
        # Tabs are: [0]=Legend, [1..]=cell types (+ optional dead).
        if idx <= 0:
            self._log(
                "\u2139\ufe0f Switch to a cell-type tab to re-run its preview."
            )
            return
        ct_idx = idx - 1
        if 0 <= ct_idx < len(self._tab_cell_types):
            ct = self._tab_cell_types[ct_idx]
            self._run_instance_preview(ct)

    def _run_instance_preview(self, cell_type):
        """Re-run post-processing for ``cell_type`` from cached classifier output."""
        try:
            if cell_type == "dead":
                self._rerun_death_preview()
                return

            tab = self.tabs.get(cell_type)
            if tab is None:
                self._log(f"\u274c No tab for {cell_type}.")
                return

            if cell_type not in self.label_map.get("celltype_to_label", {}):
                self._log(
                    f"\u26a0\ufe0f '{cell_type}' is not in the trained label map; "
                    "annotate + train first."
                )
                return

            self._log(f"Running instance segmentation for {cell_type}...")

            if self.convpaint_strategy == STRATEGY_PROB:
                self._preview_probability(cell_type, tab)
            else:
                self._preview_edt(cell_type, tab)

            _reorder_convpaint_layers(
                self.viewer, self.all_cell_types, has_death=self.has_death,
            )
            self._log(f"\u2705 Preview updated for {cell_type}")
        except Exception as e:
            self._log(f"\u274c Preview failed for {cell_type}: {e}")
            traceback.print_exc()

    def _load_unified_predicted_labels(self):
        """Load the cached multi-class label volume for the training stack."""
        path = unified_predicted_labels_path(self.pixel_class_outdir)
        if not path.exists():
            return None
        try:
            return np.asarray(load_zarr(path))
        except Exception as e:
            self._log(f"\u26a0\ufe0f Could not load unified predicted labels: {e}")
            return None

    def _load_celltype_probabilities(self, cell_type):
        """Load the cached probability volume for ``cell_type``."""
        path = _probability_map_path(self.pixel_class_outdir, cell_type)
        if path is None or not path.exists():
            return None
        try:
            return np.asarray(load_zarr(path))
        except Exception as e:
            self._log(
                f"\u26a0\ufe0f Could not load probability map for {cell_type}: {e}"
            )
            return None

    def _preview_edt(self, cell_type, tab):
        unified = self._load_unified_predicted_labels()
        if unified is None:
            self._log(
                "\u26a0\ufe0f No cached unified labels — train the classifier first."
            )
            return
        k = int(self.label_map["celltype_to_label"][cell_type])
        mask_stack = (unified == k).astype(np.uint8)

        instances = mask_to_instances(
            mask_stack,
            edt_thr=float(tab.edt_threshold_spin.value()),
            opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
            fill_holes=bool(tab.fill_holes_cb.isChecked())
            if tab.fill_holes_cb is not None else True,
            segment_size_min=int(tab.segment_size_min_spin.value()),
            marker_strategy=(
                "peak" if self.convpaint_strategy == STRATEGY_PEAK_EDT
                else "threshold"
            ),
            peak_min_distance=float(tab.peak_min_distance_spin.value()) if tab.peak_min_distance_spin is not None else None,
            peak_min_ratio=float(tab.peak_min_ratio_spin.value()) if tab.peak_min_ratio_spin is not None else 0.35,
        )
        self._set_labels_layer(_segments_layer_name(cell_type), instances)

        # EDT strategy doesn't use the probability map; hide it if cached.
        prob_lname = _probability_layer_name(cell_type)
        if prob_lname in self.viewer.layers:
            self.viewer.layers[prob_lname].visible = False

        pred_path = _predicted_labels_path(self.pixel_class_outdir, cell_type)
        if pred_path is not None:
            _save_preview_array(pred_path, instances)

    def _preview_probability(self, cell_type, tab):
        prob_stack = self._load_celltype_probabilities(cell_type)
        if prob_stack is None:
            self._log(
                f"\u26a0\ufe0f No cached probability map for {cell_type} — "
                "train the classifier first."
            )
            return
        instances = probability_to_instances(
            prob_stack,
            mask_thr=float(tab.prob_mask_threshold_spin.value()),
            seed_thr=float(tab.prob_seed_threshold_spin.value()),
            opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
            segment_size_min=int(tab.segment_size_min_spin.value()),
        )
        self._set_labels_layer(_segments_layer_name(cell_type), instances)
        self._set_image_layer(
            _probability_layer_name(cell_type), prob_stack, visible=False,
        )
        pred_path = _predicted_labels_path(self.pixel_class_outdir, cell_type)
        if pred_path is not None:
            _save_preview_array(pred_path, instances)

    def _rerun_death_preview(self):
        """Reload the death model + re-run a binary preview."""
        from napari_convpaint import ConvpaintModel
        mp = _death_model_path(self.pixel_class_outdir)
        if not mp.exists():
            self._log("\u26a0\ufe0f No trained death model on disk.")
            return
        try:
            model = ConvpaintModel(model_path=str(mp))
        except Exception as e:
            self._log(f"\u274c Could not load death model: {e}")
            return
        input_channels = load_input_channels(
            death_input_channels_path(self.pixel_class_outdir)
        )
        self._preview_death(model, self._fe_device(), input_channels)

    def _preview_death(self, model, fe_device, input_channels=None):
        seg_stack = np.stack(
            [
                segment_per_z(
                    model,
                    _slice_image_channels(img, input_channels),
                    fe_use_device=fe_device,
                )
                for img in self.all_images
            ],
            axis=0,
        ).astype(np.int16)
        mask_stack = (seg_stack >= 2).astype(np.uint16)
        self._set_labels_layer(DEAD_PREDICTED_LAYER_NAME, mask_stack)
        pred_path = _predicted_labels_path(self.pixel_class_outdir, "dead")
        if pred_path is not None:
            _save_preview_array(pred_path, mask_stack)
        _reorder_convpaint_layers(
            self.viewer, self.all_cell_types, has_death=self.has_death,
        )
        self._log("\u2705 Death mask preview updated")

    # ── Napari layer helpers ────────────────────────────────────────

    def _set_labels_layer(self, name, data, visible=True, opacity=0.8):
        if name in self.viewer.layers:
            self.viewer.layers[name].data = np.asarray(data)
            self.viewer.layers[name].visible = visible
        else:
            self.viewer.add_labels(np.asarray(data), name=name, opacity=opacity,
                                   visible=visible)

    def _set_image_layer(self, name, data, visible=False, opacity=0.6):
        if name in self.viewer.layers:
            self.viewer.layers[name].data = np.asarray(data)
            self.viewer.layers[name].visible = visible
        else:
            self.viewer.add_image(
                np.asarray(data),
                name=name,
                opacity=opacity,
                blending="additive",
                colormap="magma",
                contrast_limits=(0.0, 1.0),
                visible=visible,
            )

    # ── Save / load labels ──────────────────────────────────────────

    def save_user_labels(self):
        # Unified labels.
        layer = None
        for lyr in self.viewer.layers:
            if lyr.name == UNIFIED_LABELS_LAYER_NAME:
                layer = lyr
                break
        if layer is not None:
            label_data = np.asarray(layer.data)
            outpath = unified_user_labels_path(self.pixel_class_outdir)
            if outpath.exists():
                shutil.rmtree(outpath)
            save_as_zarr(label_data, outpath)
            self._log(f"Saved unified labels \u2192 {outpath.name}")

            # Persist the canonical label map alongside the labels (without
            # requiring training) — useful for users iterating on annotations.
            try:
                save_label_map(
                    unified_label_map_path(self.pixel_class_outdir),
                    self.label_map,
                )
            except Exception:
                pass

        # Death labels (separate layer).
        if self.has_death:
            for lyr in self.viewer.layers:
                if lyr.name == DEAD_LABELS_LAYER_NAME:
                    dead_data = np.asarray(lyr.data)
                    dead_outpath = Path(
                        self.pixel_class_outdir,
                        "PixelClassifier_UserDeadLabels.zarr",
                    )
                    if dead_outpath.exists():
                        shutil.rmtree(dead_outpath)
                    save_as_zarr(dead_data, dead_outpath)
                    self._log(f"Saved Dead labels \u2192 {dead_outpath.name}")
                    break

        self._log("\u2705 All user labels saved.")
        if self.legend_tab is not None:
            self.legend_tab.refresh()


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
    all_cell_types = _filter_merge_types(all_cell_types)

    if not all_images:
        print("\u274c No training images loaded.")
        return None
    if not all_cell_types:
        print("\u274c No cell types detected — nothing to train.")
        return None

    n_channels = all_images[0].shape[0]
    print(f"Loaded {len(all_images)} training images with {n_channels} channels")
    unified_input_channels, death_input_channels = _derive_convpaint_input_channels(
        metadata, n_channels
    )
    print(f"Unified ConvPaint training input channels: {unified_input_channels}")
    if has_death:
        print(f"Death ConvPaint training input channels: {death_input_channels}")

    stacked = np.stack(all_images, axis=0)
    T_total = stacked.shape[0]
    viewer = napari.Viewer(title="BEHAV3D - ConvPaint Training (unified)")

    for ch in range(n_channels):
        ch_data = stacked[:, ch, :, :, :]
        nonzero = ch_data[ch_data > 0]
        clim = (
            (0, float(np.percentile(nonzero, 99.8)))
            if nonzero.size > 0 else (0, 1e-3)
        )
        layer = viewer.add_image(
            ch_data,
            name=f"Channel {ch}",
            contrast_limits=clim,
            colormap=CHANNEL_COLORS[ch % len(CHANNEL_COLORS)],
            blending="additive", opacity=0.8,
        )
        layer.contrast_limits_range = (0, float(ch_data.max()))

    label_shape = (T_total,) + stacked.shape[2:]
    ip = initial_params or {}

    # Build (or rebuild) the label map for the active cell-type set.
    label_map = build_label_map(all_cell_types)

    # Restore the unified User Provided Labels layer, validating that the
    # cached cell-type ordering matches the current one. If not, start fresh.
    unified_path = unified_user_labels_path(pixel_class_outdir)
    saved_map_path = unified_label_map_path(pixel_class_outdir)
    unified_data = np.zeros(label_shape, dtype=np.int16)
    if unified_path.exists():
        try:
            existing = np.asarray(load_zarr(unified_path))
            if existing.shape == label_shape:
                if saved_map_path.exists():
                    saved_map = load_label_map(saved_map_path)
                    if label_map_matches(saved_map, all_cell_types):
                        unified_data = existing
                        print(
                            "  \u2705 Restored unified User Provided Labels "
                            f"({unified_path.name})"
                        )
                    else:
                        print(
                            "  \u26a0\ufe0f Cached unified labels were trained "
                            "with a different cell-type set — starting fresh. "
                            "Old file kept on disk."
                        )
                else:
                    # No saved map; assume current ordering is the same.
                    unified_data = existing
                    print(
                        "  \u2705 Restored unified User Provided Labels (no map sidecar)"
                    )
            else:
                print(
                    "  \u26a0\ufe0f Unified labels shape mismatch "
                    f"({existing.shape} vs {label_shape}); starting fresh."
                )
        except Exception as exc:
            print(f"  \u26a0\ufe0f Could not restore unified labels: {exc}")

    viewer.add_labels(
        unified_data, name=UNIFIED_LABELS_LAYER_NAME, opacity=0.5,
    )

    # Death labels (separate binary layer).
    if has_death:
        dead_path = Path(pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
        dead_data = np.zeros(label_shape, dtype=np.int16)
        if dead_path.exists():
            try:
                existing = np.asarray(load_zarr(dead_path))
                if existing.shape == label_shape:
                    dead_data = existing
                    print(f"  \u2705 Restored death labels ({dead_path.name})")
            except Exception as exc:
                print(f"  \u26a0\ufe0f Could not restore death labels: {exc}")
        viewer.add_labels(dead_data, name=DEAD_LABELS_LAYER_NAME, opacity=0.5)

    # Restore cached unified predicted labels (for fast preview reruns).
    pred_path = unified_predicted_labels_path(pixel_class_outdir)
    if pred_path.exists():
        try:
            pred = np.asarray(load_zarr(pred_path))
            if pred.shape == label_shape:
                viewer.add_labels(
                    pred, name=UNIFIED_PREDICTED_LAYER_NAME, opacity=0.6,
                    visible=False,
                )
        except Exception as exc:
            print(f"  \u26a0\ufe0f Could not restore unified predicted labels: {exc}")

    # Restore per-cell-type cached previews (segments + probability maps).
    for ct in all_cell_types:
        seg_path = _predicted_labels_path(pixel_class_outdir, ct)
        if seg_path is not None and seg_path.exists():
            try:
                pred = np.asarray(load_zarr(seg_path))
                if pred.shape == label_shape:
                    viewer.add_labels(
                        pred, name=_segments_layer_name(ct),
                        opacity=0.8, visible=False,
                    )
            except Exception as exc:
                print(f"  \u26a0\ufe0f Could not restore segments for '{ct}': {exc}")

        prob_path = _probability_map_path(pixel_class_outdir, ct)
        if prob_path is not None and prob_path.exists():
            try:
                prob = np.asarray(load_zarr(prob_path))
                if prob.shape == label_shape:
                    viewer.add_image(
                        prob, name=_probability_layer_name(ct), opacity=0.6,
                        blending="additive", colormap="magma",
                        contrast_limits=(0.0, 1.0), visible=False,
                    )
            except Exception as exc:
                print(
                    f"  \u26a0\ufe0f Could not restore probability map for '{ct}': {exc}"
                )

    if has_death:
        seg_path = _predicted_labels_path(pixel_class_outdir, "dead")
        if seg_path is not None and seg_path.exists():
            try:
                pred = np.asarray(load_zarr(seg_path))
                if pred.shape == label_shape:
                    viewer.add_labels(
                        pred, name=DEAD_PREDICTED_LAYER_NAME, opacity=0.8,
                        visible=False,
                    )
            except Exception as exc:
                print(f"  \u26a0\ufe0f Could not restore death mask: {exc}")

    cp_strategy = _normalize_strategy(
        ip.get("convpaint_strategy", DEFAULT_STRATEGY)
    )

    widget = ConvPaintTrainingWidget(
        viewer=viewer,
        all_images=all_images,
        pixel_class_outdir=pixel_class_outdir,
        all_cell_types=all_cell_types,
        has_death=has_death,
        initial_params=ip,
        on_params_changed=on_params_changed,
        convpaint_strategy=cp_strategy,
        unified_input_channels=unified_input_channels,
        death_input_channels=death_input_channels,
    )
    viewer.window.add_dock_widget(widget, name="ConvPaint", area="right")
    _reorder_convpaint_layers(viewer, all_cell_types, has_death=has_death)
    return viewer
