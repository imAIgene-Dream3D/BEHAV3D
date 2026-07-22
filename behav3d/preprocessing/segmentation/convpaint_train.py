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
    QGroupBox, QGridLayout, QPlainTextEdit, QApplication, QScrollArea,
    QTabWidget, QFrame, QCheckBox, QSizePolicy, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView,
)
from qtpy.QtCore import Qt, Signal

from behav3d.napari._background_runner import BackgroundOperation, ThreadSafeLogger
from qtpy.QtGui import QColor

from behav3d.io.images import load_image, load_zarr, save_as_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape
from behav3d.core.qt_help import HelpButton, make_help_row
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


def _set_dead_mask_layer_color(layer, color=(0.9, 0.0, 0.0, 1.0)):
    """Render the death pixel-classification layer as a solid color (red by
    default) instead of napari's default per-label random colors — it is a
    binary mask (dead vs. not), not a multi-class label map. Any nonzero
    label (whichever convention is used for "dead") is painted ``color``;
    background (0) stays transparent."""
    try:
        from napari.utils.colormaps import DirectLabelColormap
        layer.colormap = DirectLabelColormap(color_dict={None: list(color), 0: [0, 0, 0, 0]})
    except Exception:
        try:
            layer.color = {1: color, 2: color}
        except Exception:
            pass


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
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.table.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self.table.setMinimumWidth(0)
        header = self.table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        header.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        header.setStretchLastSection(True)
        header.setMinimumSectionSize(20)
        self.table.cellClicked.connect(self._on_cell_clicked)
        layout.addWidget(self.table)

        button_row = QHBoxLayout()
        self.btn_refresh = QtPushButton("Refresh counts & colors")
        self.btn_refresh.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Fixed)
        self.btn_refresh.clicked.connect(self.refresh)
        button_row.addWidget(self.btn_refresh)
        button_row.addStretch()
        self.note = QLabel(
            "<i>Death is a separate binary classifier — annotate it on the "
            "<b>User Provided Labels (Dead)</b> layer.</i>"
        )
        self.note.setWordWrap(True)
        if not has_death:
            self.note.hide()
        layout.addWidget(self.note)

        if has_death:
            def create_button(text, color_hex, layer_name, label_idx):
                btn = QtPushButton(text)
                text_color = "white" if color_hex not in ("none", "white", "#ffffff", "transparent") else "black"
                bg_color = "transparent" if color_hex == "none" else color_hex
                btn.setStyleSheet(f"background-color: {bg_color}; color: {text_color}; font-size: 11px; padding: 4px; border: 1px solid #555; border-radius: 2px;")
                btn.clicked.connect(lambda _, ln=layer_name, idx=label_idx: self._select_brush(ln, idx))
                return btn

            dead_row = QHBoxLayout()
            dead_row.addWidget(QLabel("<b>Dead</b>"), stretch=1)
            layer_name = 'User Provided Labels (Dead)'
            dead_row.addWidget(create_button("0: Eraser", "none", layer_name, 0))
            dead_row.addWidget(create_button("1: Background", "#8b3a26", layer_name, 1))
            dead_row.addWidget(create_button("2: Foreground", "#00ffff", layer_name, 2))
            layout.addLayout(dead_row)

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

    def _select_brush(self, layer_name, label_idx):
        layer = None
        for lyr in self.viewer.layers:
            if lyr.name == layer_name:
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
    """One tab per cell type with strategy-specific spinners + preview button.

    Constructor flags ``show_strategy_combo`` and ``per_tab_strategies``
    (used by the napari GUI plugin) attach a ``\u21b3 Cell type strategy:``
    combo to the top of the tab so each cell type can override the global
    ConvPaint strategy. When changed, the post-processing controls are
    rebuilt to match the new strategy. The notebook entry point never sets
    these flags, so its UI is unchanged.
    """

    def __init__(self, cell_type, strategy, initial_params=None,
                 on_params_changed=None, run_preview_callback=None,
                 parent=None,
                 show_strategy_combo=False,
                 per_tab_strategies=None,
                 on_per_tab_strategy_changed=None,
                 pixel_sizes=None):
        super().__init__(parent)
        self.cell_type = cell_type
        self.strategy = _normalize_strategy(strategy)
        self._on_params_changed = on_params_changed
        self._run_preview_callback = run_preview_callback
        self._is_dead = (cell_type == "dead")
        self._show_strategy_combo = bool(show_strategy_combo) and not self._is_dead
        self._per_tab_strategies = list(per_tab_strategies or [])
        self._on_per_tab_strategy_changed = on_per_tab_strategy_changed
        self._per_tab_strategy_combo = None
        self._per_tab_strategy_widget = None
        # Optional resolved pixel sizes ({"xy_um":..., "z_um":...}) for the
        # physical(µm)/pixel unit toggle on EDT/min-size/peak-distance
        # controls (see behav3d.napari._units.UnitGroupManager).
        self._pixel_sizes = dict(pixel_sizes) if pixel_sizes else {}
        self._unit_mgr = None  # created lazily in _build_controls
        # Persists the user's px/µm display choice across strategy-change
        # rebuilds (a fresh UnitGroupManager is created each rebuild since
        # the spinbox instances themselves are discarded and recreated).
        self._units_physical = True

        # Cached so future strategy switches reuse persisted values.
        self._initial_params_cache = dict(initial_params or {})
        ip = self._initial_params_cache

        self._root_layout = QVBoxLayout()
        self._root_layout.setContentsMargins(8, 8, 8, 8)
        self._root_layout.setSpacing(6)

        if self._is_dead:
            _dead_lbl = QLabel(
                "<i>Death classifier is a separate binary ConvPaint model. "
                "Use the <b>User Provided Labels (Dead)</b> layer to annotate it. "
                "No instance-segmentation parameters required.</i>"
            )
            _dead_lbl.setWordWrap(True)
            self._root_layout.addWidget(_dead_lbl)
            self.btn_preview = QtPushButton("Run preview (death mask)")
            self.btn_preview.clicked.connect(self._on_run_preview)
            self._root_layout.addWidget(self.btn_preview)
            self._root_layout.addStretch()
            self.setLayout(self._root_layout)
            self.edt_threshold_spin = None
            self.opening_nr_pixels_spin = None
            self.segment_size_min_spin = None
            self.fill_holes_cb = None
            self.prob_mask_threshold_spin = None
            self.prob_seed_threshold_spin = None
            self.peak_min_distance_spin = None
            self.peak_min_ratio_spin = None
            return

        # Optional per-tab strategy combo at the very top of the tab.
        if self._show_strategy_combo and self._per_tab_strategies:
            saved_per_tab = ip.get(
                f"per_ct_convpaint_strategy_{cell_type}", self.strategy
            )
            saved_per_tab = _normalize_strategy(saved_per_tab)
            if saved_per_tab not in self._per_tab_strategies:
                saved_per_tab = self._per_tab_strategies[0]
            self.strategy = saved_per_tab
            self._per_tab_strategy_widget = QWidget()
            per_tab_lay = QHBoxLayout(self._per_tab_strategy_widget)
            per_tab_lay.setContentsMargins(0, 2, 0, 4)
            per_tab_lay.addWidget(QLabel("\u21b3 Cell type strategy:"))
            self._per_tab_strategy_combo = QComboBox()
            self._per_tab_strategy_combo.addItems(self._per_tab_strategies)
            self._per_tab_strategy_combo.setCurrentText(saved_per_tab)
            _compact_combo(self._per_tab_strategy_combo, min_chars=12)
            self._per_tab_strategy_combo.currentTextChanged.connect(
                self._on_per_tab_strategy_combo_changed
            )
            per_tab_lay.addWidget(self._per_tab_strategy_combo, stretch=1)
            self._root_layout.addWidget(self._per_tab_strategy_widget)

        # Placeholder layout where the strategy-dependent controls are
        # inserted. Lets us swap the controls on strategy change without
        # rebuilding the whole tab.
        self._controls_widget = QWidget()
        self._controls_layout = QVBoxLayout(self._controls_widget)
        self._controls_layout.setContentsMargins(0, 0, 0, 0)
        self._controls_layout.setSpacing(6)
        self._root_layout.addWidget(self._controls_widget)
        self._root_layout.addStretch()
        self.setLayout(self._root_layout)

        self.edt_threshold_spin = None
        self.opening_nr_pixels_spin = None
        self.segment_size_min_spin = None
        self.fill_holes_cb = None
        self.prob_mask_threshold_spin = None
        self.prob_seed_threshold_spin = None
        self.peak_min_distance_spin = None
        self.peak_min_ratio_spin = None
        self._build_controls(self.strategy)

    @staticmethod
    def _purge_layout(layout):
        """Recursively remove every item (widget or sub-layout) from *layout*.

        ``QLayout.takeAt`` removes the item from the layout but leaves
        sub-layout widgets orphaned on the parent widget, causing them to
        paint as ghost overlays.  This helper walks into nested layouts and
        calls ``setParent(None)`` + ``deleteLater()`` on every widget it
        finds so they are fully destroyed.
        """
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()
            else:
                child = item.layout()
                if child is not None:
                    CellTypeConvPaintTab._purge_layout(child)

    def _clear_controls_layout(self):
        self._purge_layout(self._controls_layout)

    @staticmethod
    def _pair_row(lbl1, w1, h1, lbl2, w2, h2):
        """Two label+widget+help groups stacked vertically (one per row).

        Previously placed both pairs on a single horizontal row, but the
        combined width forced a wide minimum on the cell-type tabs (and
        therefore the whole widget). Stacking them keeps each row narrow
        so the panel can shrink to the width of a single label+spinner+help.
        """
        for w in (w1, w2):
            if isinstance(w, (QSpinBox, QDoubleSpinBox)):
                w.setMaximumWidth(70)
                w.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Fixed)

        outer = QVBoxLayout()
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(2)

        for lbl, w, h in ((lbl1, w1, h1), (lbl2, w2, h2)):
            row = QHBoxLayout()
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(2)
            row.addWidget(lbl)
            row.addWidget(w)
            row.addWidget(h)
            row.addStretch()
            outer.addLayout(row)

        return outer

    @staticmethod
    def _solo_row(*widgets):
        """Single group of widgets on one row with a trailing stretch."""
        for w in widgets:
            if isinstance(w, (QSpinBox, QDoubleSpinBox)):
                w.setMaximumWidth(70)
                w.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Fixed)
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(2)
        for w in widgets:
            row.addWidget(w)
        row.addStretch()
        return row

    def get_native(self, widget):
        """Return the canonical native (px/voxel) value for a spinbox
        registered with this tab's unit manager, or its raw ``.value()``
        when no manager applies (dimensionless params, or no resolution
        available). Safe to call with ``widget=None``."""
        if widget is None:
            return None
        if self._unit_mgr is not None and widget in getattr(self._unit_mgr, "_entries", ()):
            return self._unit_mgr.get_native(widget)
        return widget.value()

    def set_native(self, widget, native_value):
        """Set a spinbox to a native value, converting to the currently
        displayed unit via this tab's unit manager when registered."""
        if widget is None or native_value is None:
            return
        if self._unit_mgr is not None and widget in getattr(self._unit_mgr, "_entries", ()):
            self._unit_mgr.set_native(widget, native_value)
        else:
            widget.setValue(native_value)

    def _on_units_toggled(self, checked):
        """Remember the px/µm display choice across strategy rebuilds.

        Native values are unaffected by this toggle — it only changes what
        the spinboxes *display* — so this does not need to trigger a param
        re-persist."""
        self._units_physical = bool(checked)

    def _build_controls(self, strategy):
        """(Re)create the strategy-dependent spinners + preview button."""
        self._clear_controls_layout()

        ip = self._initial_params_cache
        cell_type = self.cell_type

        # Reset references before rebuilding.
        self.edt_threshold_spin = None
        self.opening_nr_pixels_spin = None
        self.segment_size_min_spin = None
        self.fill_holes_cb = None
        self.prob_mask_threshold_spin = None
        self.prob_seed_threshold_spin = None
        self.peak_min_distance_spin = None
        self.peak_min_ratio_spin = None

        strategy = _normalize_strategy(strategy)
        self.strategy = strategy
        # A fresh manager is created below (rather than reused) since the
        # spinbox instances themselves are rebuilt from scratch here.
        self._unit_mgr = None

        layout = self._controls_layout

        strategy_lbl = QLabel(f"<i>Strategy:</i> <b>{strategy}</b>")
        strategy_lbl.setWordWrap(True)
        layout.addWidget(strategy_lbl)

        # Per-tab physical(µm)/pixel unit toggle for the distance/volume
        # controls below (EDT threshold, min size, peak distance). Mask/seed
        # thresholds and peak ratio are dimensionless (0-1) and excluded, as
        # is opening (fixed in pixels by definition).
        from behav3d.napari._units import UnitGroupManager
        self._unit_mgr = UnitGroupManager(
            xy_um=self._pixel_sizes.get("xy_um"),
            z_um=self._pixel_sizes.get("z_um"),
            default_physical=self._units_physical,
        )
        self._unit_mgr.switch.toggled.connect(self._on_units_toggled)
        layout.addWidget(self._unit_mgr.header_row("Units:"))

        if strategy == STRATEGY_PROB:
            self.prob_mask_threshold_spin = QDoubleSpinBox()
            self.prob_mask_threshold_spin.setRange(0.0, 1.0)
            self.prob_mask_threshold_spin.setSingleStep(0.05)
            self.prob_mask_threshold_spin.setDecimals(2)
            self.prob_mask_threshold_spin.setValue(
                float(ip.get(f"{cell_type}_prob_mask_threshold", 0.5))
            )

            self.prob_seed_threshold_spin = QDoubleSpinBox()
            self.prob_seed_threshold_spin.setRange(0.0, 1.0)
            self.prob_seed_threshold_spin.setSingleStep(0.05)
            self.prob_seed_threshold_spin.setDecimals(2)
            self.prob_seed_threshold_spin.setValue(
                float(ip.get(f"{cell_type}_prob_seed_threshold", 0.8))
            )
            layout.addLayout(self._pair_row(
                QLabel("Mask threshold:"), self.prob_mask_threshold_spin,
                HelpButton("Mask threshold",
                    "Foreground cutoff applied to the per-cell-type probability map.\n"
                    "Pixels above this value are kept as foreground (typical: 0.5)."),
                QLabel("Seed threshold:"), self.prob_seed_threshold_spin,
                HelpButton("Seed threshold",
                    "Higher cutoff used to place watershed seeds (≥ Mask threshold).\n"
                    "Higher values keep only each object's confident core as a "
                    "separate seed, splitting more touching objects. Lower values "
                    "(closer to Mask threshold) merge neighboring cores together, "
                    "splitting fewer objects.\n\n"
                    "Note: this is the SAME direction as the 'EDT threshold' on the "
                    "plain EDT/Watershed strategy (there, higher values also split "
                    "more), but the OPPOSITE direction from that same field on the "
                    "Peak EDT/Watershed strategy, where it acts as a peak-height "
                    "noise filter instead — same identical behaviour in APOC and "
                    "ConvPaint."),
            ))

            self.opening_nr_pixels_spin = QSpinBox()
            self.opening_nr_pixels_spin.setRange(0, 50)
            self.opening_nr_pixels_spin.setValue(
                int(ip.get(f"{cell_type}_opening_nr_pixels", 0))
            )

            self.segment_size_min_spin = QSpinBox()
            self.segment_size_min_spin.setRange(0, 1000000000)
            self.segment_size_min_spin.setSingleStep(10)
            self.segment_size_min_spin.setValue(
                int(ip.get(f"{cell_type}_segment_size_min", 10))
            )
            self._unit_mgr.register(
                self.segment_size_min_spin, "volume",
                ip.get(f"{cell_type}_segment_size_min", 10),
            )
            layout.addLayout(self._pair_row(
                QLabel("Opening px:"), self.opening_nr_pixels_spin,
                HelpButton("Morphological opening",
                    "Number of erosion-then-dilation iterations applied to the mask.\n"
                    "Smooths boundaries and removes small speckles. Set to 0 to disable."),
                QLabel("Min size:"), self.segment_size_min_spin,
                HelpButton("Minimum segment size",
                    "Segments with fewer voxels than this are discarded after watershed.\n"
                    "Use to filter out noise / debris."),
            ))
        else:
            self.edt_threshold_spin = QDoubleSpinBox()
            self.edt_threshold_spin.setRange(0.0, 1000.0)
            self.edt_threshold_spin.setSingleStep(0.5)
            self.edt_threshold_spin.setDecimals(2)
            self.edt_threshold_spin.setValue(
                float(ip.get(f"{cell_type}_edt_threshold", 1.0))
            )
            self._unit_mgr.register(
                self.edt_threshold_spin, "distance",
                ip.get(f"{cell_type}_edt_threshold", 1.0),
            )

            self.opening_nr_pixels_spin = QSpinBox()
            self.opening_nr_pixels_spin.setRange(0, 50)
            self.opening_nr_pixels_spin.setValue(
                int(ip.get(f"{cell_type}_opening_nr_pixels", 0))
            )
            _edt_help = (
                HelpButton(
                    "EDT threshold (peak height)",
                    "Minimum EDT peak height required for a local maximum to count as "
                    "a watershed seed — a noise filter, not an erosion cutoff.\n"
                    "Higher values suppress weaker peaks (e.g. a boundary bump on one "
                    "cell) → fewer seeds, less splitting. Lower values keep weak peaks "
                    "→ more seeds, more splitting.\n\n"
                    "Note: this is the OPPOSITE direction from this same field on the "
                    "plain 'Mask + EDT/Watershed' strategy, where higher values "
                    "increase splitting instead (there it's an erosion cutoff, not a "
                    "peak-height filter)."
                )
                if strategy == STRATEGY_PEAK_EDT
                else HelpButton(
                    "EDT threshold",
                    "Euclidean-distance-transform threshold used to derive seeds inside "
                    "the binary mask (an erosion-style cutoff: seeds = mask voxels with "
                    "EDT ≥ this value).\n"
                    "Higher values erode further, breaking the thin neck between "
                    "touching objects → more splitting. Lower values keep neighbouring "
                    "cores connected → objects stay merged.\n\n"
                    "Note: this is the SAME direction as 'Seed threshold' used by the "
                    "Probability Map + Watershed strategy (there, higher values also "
                    "split more)."
                )
            )
            layout.addLayout(self._pair_row(
                QLabel("EDT threshold:"), self.edt_threshold_spin,
                _edt_help,
                QLabel("Opening px:"), self.opening_nr_pixels_spin,
                HelpButton("Morphological opening",
                    "Number of erosion-then-dilation iterations applied to the mask.\n"
                    "Smooths boundaries and removes small speckles. Set to 0 to disable."),
            ))

            self.segment_size_min_spin = QSpinBox()
            self.segment_size_min_spin.setRange(0, 1000000000)
            self.segment_size_min_spin.setSingleStep(10)
            self.segment_size_min_spin.setValue(
                int(ip.get(f"{cell_type}_segment_size_min", 10))
            )
            self._unit_mgr.register(
                self.segment_size_min_spin, "volume",
                ip.get(f"{cell_type}_segment_size_min", 10),
            )
            self.fill_holes_cb = QCheckBox("Fill holes")
            self.fill_holes_cb.setChecked(
                bool(ip.get(f"{cell_type}_fill_holes", True))
            )
            layout.addLayout(self._solo_row(
                QLabel("Min size:"), self.segment_size_min_spin,
                HelpButton("Minimum segment size",
                    "Segments with fewer voxels than this are discarded after watershed.\n"
                    "Use to filter out noise / debris."),
                self.fill_holes_cb,
                HelpButton("Fill holes",
                    "Fill internal gaps in segmented objects before watershed.\n"
                    "Useful for hollow nuclei or partial-volume effects."),
            ))

            if strategy == STRATEGY_PEAK_EDT:
                self.peak_min_distance_spin = QDoubleSpinBox()
                self.peak_min_distance_spin.setRange(0.0, 1000.0)
                self.peak_min_distance_spin.setSingleStep(0.5)
                self.peak_min_distance_spin.setDecimals(2)
                self.peak_min_distance_spin.setValue(
                    float(ip.get(f"{cell_type}_peak_min_distance", 0.0))
                )
                self._unit_mgr.register(
                    self.peak_min_distance_spin, "distance",
                    ip.get(f"{cell_type}_peak_min_distance", 0.0),
                )
                self.peak_min_distance_spin.setToolTip(
                    "Minimum distance between local EDT peaks used as watershed seeds"
                )

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
                layout.addLayout(self._pair_row(
                    QLabel("Peak min dist:"), self.peak_min_distance_spin,
                    HelpButton("Peak minimum distance",
                        "Minimum distance (µm) between local EDT peaks used as watershed seeds.\n"
                        "Larger values yield fewer, more separated seeds."),
                    QLabel("Peak min ratio:"), self.peak_min_ratio_spin,
                    HelpButton("Peak minimum ratio",
                        "Minimum EDT peak height as a fraction of the local maximum (0–1).\n"
                        "Higher values suppress weak peaks (fewer seeds)."),
                ))

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

    def rebuild_instance_controls(self, strategy=None, initial_params=None):
        """Discard and rebuild the post-processing controls for *strategy*."""
        if initial_params is None:
            initial_params = self._initial_params_cache
        else:
            self._initial_params_cache = dict(initial_params)
        target = _normalize_strategy(strategy or self.strategy)
        if self._is_dead:
            self.strategy = target
            return
        self._build_controls(target)

    def _on_per_tab_strategy_combo_changed(self, new_strategy):
        new_strategy = _normalize_strategy(new_strategy)
        self.rebuild_instance_controls(strategy=new_strategy)
        if callable(self._on_per_tab_strategy_changed):
            self._on_per_tab_strategy_changed(self.cell_type, new_strategy)

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
            params[f"{ct}_edt_threshold"] = float(self.get_native(self.edt_threshold_spin))
        if self.opening_nr_pixels_spin is not None:
            params[f"{ct}_opening_nr_pixels"] = int(self.opening_nr_pixels_spin.value())
        if self.segment_size_min_spin is not None:
            params[f"{ct}_segment_size_min"] = int(round(self.get_native(self.segment_size_min_spin)))
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
                self.get_native(self.peak_min_distance_spin)
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
    """Unified ConvPaint training widget with optional plugin extensions.

    Constructor flags
    -----------------
    ``per_cell_type_strategy`` (default ``False``)
        When ``True``, adds a global "Strategy" combo (with an
        ``Advanced (per cell type)`` option) and a per-tab ``↳ Cell type
        strategy:`` combo on each non-dead tab. The notebook entry point
        never sets this flag.
    ``strategy_resolver`` (default ``None``)
        Optional ``Callable[[str], str]`` that returns the effective strategy
        for a cell type. Takes precedence over the built-in resolution.
    ``extra_toolbar_widgets`` (default ``None``)
        Iterable of ``QWidget`` / ``QLayout`` instances inserted near the
        top of the widget, between the configuration groups and the tabs.
        Used by the plugin to surface session-level controls.
    ``show_device`` (default ``True``)
        When ``False``, the Device combo is built but not added to the
        layout. Pass ``False`` when embedding inside a plugin page that
        already provides its own device selector, to avoid duplication.
    """

    # ── Signals ─────────────────────────────────────────────────────────────
    channels_refreshed = Signal(str)
    training_started = Signal(list)
    training_finished = Signal(list, object)
    instance_preview_started = Signal(str)
    instance_preview_finished = Signal(str)
    strategy_changed = Signal(str)

    STRATEGIES = [STRATEGY_EDT, STRATEGY_PEAK_EDT, STRATEGY_PROB]
    import os
    if os.environ.get("BEHAV3D_DEV_MODE") != "1":
        STRATEGIES.remove(STRATEGY_PEAK_EDT)
    ADVANCED_STRATEGY = "Advanced (per cell type)"

    def __init__(self, viewer, all_images, pixel_class_outdir,
                 all_cell_types, has_death, initial_params=None,
                 on_params_changed=None, convpaint_strategy=None,
                 unified_input_channels=None, death_input_channels=None,
                 parent=None,
                 per_cell_type_strategy=False,
                 strategy_resolver=None,
                 extra_toolbar_widgets=None,
                 show_device=True,
                 external_log=None,
                 show_legend=False,
                 pixel_sizes=None):
        super().__init__(parent)
        self.viewer = viewer
        self.all_images = all_images
        self.pixel_class_outdir = Path(pixel_class_outdir)
        self.all_cell_types = list(all_cell_types)
        self.has_death = has_death
        # Optional resolved pixel sizes ({"xy_um":..., "z_um":...}) threaded
        # down to each per-cell-type tab for the µm/pixel unit toggle on
        # EDT threshold / min size / peak distance.
        self._pixel_sizes = dict(pixel_sizes) if pixel_sizes else {}
        self.unified_input_channels = list(unified_input_channels or [])
        self.death_input_channels = list(death_input_channels or [])
        self._on_params_changed = on_params_changed
        self.ip = initial_params or {}
        self.convpaint_strategy = _normalize_strategy(
            convpaint_strategy or self.ip.get("convpaint_strategy", DEFAULT_STRATEGY)
        )
        self._per_cell_type_strategy = bool(per_cell_type_strategy)
        self._strategy_resolver = strategy_resolver
        self._extra_toolbar_widgets = list(extra_toolbar_widgets or [])
        self._show_device = bool(show_device)
        self._show_legend = bool(show_legend)
        # When provided, log messages are forwarded here instead of the
        # internal log_box (which is then hidden to reclaim vertical space).
        self._external_log = external_log if callable(external_log) else None

        # Active label map (always built fresh from the current cell-type set).
        self.label_map = build_label_map(self.all_cell_types)

        # Per-cell-type segmentation tabs + (optional) Dead tab.
        self._tab_cell_types = list(self.all_cell_types)
        if self.has_death:
            self._tab_cell_types.append("dead")

        self.tabs = {}
        self.legend_tab = None
        self.strategy_combo = None
        self.strategy_help_button = None
        self._bg = BackgroundOperation(self)
        self._build_ui()
        if self._external_log is not None and hasattr(self, "log_box"):
            self.log_box.hide()
        if self.strategy_combo is not None:
            self._apply_strategy_to_tabs(
                self.strategy_combo.currentText(), emit_signal=False
            )

    # ── UI construction ─────────────────────────────────────────────

    def _build_ui(self):
        self.content_widget = QWidget()
        root = QVBoxLayout(self.content_widget)
        root.setContentsMargins(4, 4, 4, 4)
        self._main_layout = root
        root.addWidget(QLabel("<h3>ConvPaint Training (unified)</h3>"))

        legend_intro = QLabel(
            "<i>One classifier predicts every cell type at once. "
            "Annotate every class on the single <b>User Provided Labels</b> layer "
            "(see legend tab below for index ↔ cell type).</i>"
        )
        legend_intro.setWordWrap(True)
        legend_intro.setStyleSheet("color: #444; margin-bottom: 4px;")
        root.addWidget(legend_intro)

        # Optional global Strategy combo (plugin only).
        if self._per_cell_type_strategy:
            strat_row = QHBoxLayout()
            strat_row.addWidget(QLabel("Strategy:"))
            self.strategy_combo = QComboBox()
            options = list(self.STRATEGIES) + [self.ADVANCED_STRATEGY]
            self.strategy_combo.addItems(options)
            initial = self.convpaint_strategy if self.convpaint_strategy in options else options[0]
            self.strategy_combo.setCurrentText(initial)
            _compact_combo(self.strategy_combo, min_chars=12)
            strat_row.addWidget(self.strategy_combo, stretch=1)
            root.addLayout(strat_row)

            strat_desc = QLabel(
                "Strategy determines how the ConvPaint classifier output is"
                " converted to instance segmentation labels. Post-processing"
                " parameters appear inside each cell-type tab."
            )
            strat_desc.setWordWrap(True)
            strat_desc.setStyleSheet("color:#888; font-size:10px; padding: 0 0 4px 0;")
            root.addWidget(strat_desc)

            help_row = QHBoxLayout()
            help_row.setContentsMargins(0, 0, 0, 0)
            help_lbl = QLabel("Instance preview parameters")
            help_lbl.setStyleSheet("color:#888; font-size:10px;")
            self.strategy_help_button = HelpButton(
                *self._strategy_help_content(initial)
            )
            help_row.addWidget(help_lbl)
            help_row.addWidget(self.strategy_help_button)
            help_row.addStretch()
            root.addLayout(help_row)

            self.strategy_combo.currentTextChanged.connect(self._on_global_strategy_changed)

        self.legend_tab = None
        if self._show_legend:
            self.legend_tab = AnnotationLegendTab(self.viewer, self.label_map, has_death=self.has_death)
            root.addWidget(self.legend_tab)

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
        r.addWidget(HelpButton(
            "Feature extractor model",
            "Pre-trained deep network used to extract per-pixel features.\n"
            "  • VGG-16 (default): fast and lightweight; good general baseline.\n"
            "  • DINOv2 / DINOv3: stronger feature representations from "
            "self-supervised vision transformers; slower but often more accurate.\n"
            "Switching models invalidates the trained CatBoost classifier — "
            "retrain after changing this."
        ))
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
        r.addWidget(HelpButton(
            "Channel mode",
            "How the input channels are fed to the feature extractor:\n"
            "  • multi: each channel processed separately and concatenated "
            "(more features, more memory).\n"
            "  • single (mean) / RGB: channels are mixed before extraction.\n"
            "When in doubt keep 'multi' for multi-channel fluorescence."
        ))
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
        r.addWidget(HelpButton(
            "Image downsample",
            "Spatial downsampling factor applied before feature extraction.\n"
            "  • 1 = native resolution.\n"
            "  • >1 = smaller image → faster but coarser predictions.\n"
            "  • <0 = upsample (rarely useful)."
        ))
        r.addWidget(QLabel("Smoothen:"))
        self.smooth_spin = QSpinBox()
        self.smooth_spin.setRange(0, 20)
        self.smooth_spin.setValue(int(self.ip.get("convpaint_seg_smoothening", 1)))
        r.addWidget(self.smooth_spin)
        r.addWidget(HelpButton(
            "Segmentation smoothening",
            "Number of post-classification smoothing iterations applied to the "
            "predicted label map.\n"
            "Reduces single-pixel noise but can erode thin structures."
        ))
        fe_layout.addLayout(r)
        fe_group.setLayout(fe_layout)
        self.fe_group = fe_group
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
        r.addWidget(HelpButton(
            "CatBoost — iterations",
            "Number of boosting iterations (trees) trained by CatBoost.\n"
            "More iterations → richer model but slower training and risk of overfit.\n"
            "Default (100) is usually a good starting point."
        ))
        r.addWidget(QLabel("Depth:"))
        self.depth_spin = QSpinBox()
        self.depth_spin.setRange(1, 16)
        self.depth_spin.setValue(int(self.ip.get("convpaint_clf_depth", 5)))
        r.addWidget(self.depth_spin)
        r.addWidget(HelpButton(
            "CatBoost — tree depth",
            "Maximum depth of each CatBoost decision tree.\n"
            "Shallow trees (3–6) generalise better; deeper trees fit more "
            "complex decision boundaries."
        ))
        clf_layout.addLayout(r)
        r = QHBoxLayout()
        r.addWidget(QLabel("Learning rate:"))
        self.lr_spin = QDoubleSpinBox()
        self.lr_spin.setRange(0.001, 1.0)
        self.lr_spin.setSingleStep(0.01)
        self.lr_spin.setDecimals(3)
        self.lr_spin.setValue(float(self.ip.get("convpaint_clf_learning_rate", 0.1)))
        r.addWidget(self.lr_spin)
        r.addWidget(HelpButton(
            "CatBoost — learning rate",
            "Shrinkage applied to each new tree's contribution.\n"
            "Lower values (e.g. 0.03) need more iterations but often generalise "
            "better; higher values train faster but may overshoot."
        ))
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
        r.addWidget(HelpButton(
            "Class weights",
            "Re-weighting strategy when some classes have far fewer annotations:\n"
            "  • None: train as-is (best when classes are balanced).\n"
            "  • Balanced: inverse-frequency weights (boost rare classes).\n"
            "  • SqrtBalanced: square root of Balanced, milder adjustment."
        ))
        clf_layout.addLayout(r)

        # Tile annotations lives here (not in the prediction-performance
        # group below) because it only affects *training*: it crops feature
        # extraction to bounding boxes around painted regions instead of the
        # whole image. Being inside clf_group means it is correctly disabled
        # together with the rest of the training controls when no training
        # images are loaded (see disable_training_params()).
        r = QHBoxLayout()
        self.tile_annotations_cb = QCheckBox("Tile annotations")
        self.tile_annotations_cb.setToolTip(
            "Only extract features inside annotated regions during training "
            "(faster on sparse annotations)."
        )
        self.tile_annotations_cb.setChecked(
            bool(self.ip.get("convpaint_tile_annotations", False))
        )
        r.addWidget(self.tile_annotations_cb)
        r.addWidget(HelpButton(
            "Tile annotations (training only)",
            "Applies only during 'Train Classifier'.\n\n"
            "Instead of extracting features over the whole image, ConvPaint "
            "crops a bounding box around each painted region and only "
            "extracts features inside those boxes.\n\n"
            "Speeds up training when your annotations are small/sparse "
            "relative to the image size. Has no effect on prediction "
            "(Run instance segmentation / Re-run preview), which is why it "
            "lives with the other training-only controls and is disabled "
            "when no training images are loaded."
        ))
        r.addStretch()
        clf_layout.addLayout(r)

        clf_group.setLayout(clf_layout)
        self.clf_group = clf_group
        root.addWidget(clf_group)

        # Device (global) — build always so _fe_device() is safe, but only
        # add to the layout when the widget is used standalone (not embedded
        # inside a plugin page that already exposes its own device selector).
        self.device_combo = QComboBox()
        for label, dev_str in _detect_torch_devices():
            self.device_combo.addItem(label, dev_str)
        idx = self.device_combo.findData(
            self.ip.get("convpaint_device", _default_device())
        )
        if idx >= 0:
            self.device_combo.setCurrentIndex(idx)
        _compact_combo(self.device_combo, min_chars=10)
        if self._show_device:
            r = QHBoxLayout()
            r.addWidget(QLabel("Device:"))
            r.addWidget(self.device_combo)
            root.addLayout(r)

        # Prediction performance (large images) \u2014 deliberately its own group,
        # separate from Feature Extractor / Classifier: these two controls
        # only affect prediction (Run instance segmentation / Re-run preview
        # / batch processing), so they must stay usable even when no training
        # images are loaded (unlike Tile annotations above, which is
        # training-only and lives in clf_group instead). Never disabled by
        # disable_training_params().
        perf_group = QGroupBox("Prediction Performance (large images)")
        perf_layout = QHBoxLayout()
        self.tile_image_cb = QCheckBox("Tile image")
        self.tile_image_cb.setToolTip(
            "Extract features in tiles during prediction "
            "(reduces peak memory on large images)."
        )
        self.tile_image_cb.setChecked(
            bool(self.ip.get("convpaint_tile_image", False))
        )
        perf_layout.addWidget(self.tile_image_cb)
        perf_layout.addWidget(HelpButton(
            "Tile image",
            "Applies to prediction (Run instance segmentation / Re-run "
            "preview / batch processing) \u2014 available whether or not "
            "training images are currently loaded.\n\n"
            "Splits the image into ~1000\u00d71000 px blocks with a 50 px "
            "overlap margin, extracts features and classifies each block "
            "separately, then stitches the results back together.\n\n"
            "Lower peak memory (only one block's features are held at a "
            "time) at the cost of some recomputation in the overlap "
            "margins \u2014 somewhat slower than processing the whole image "
            "in one pass.\n\n"
            "Turn on when a large image runs out of memory or crashes "
            "during prediction. Leave off for images that comfortably fit "
            "in memory \u2014 it will be faster."
        ))
        self.use_dask_cb = QCheckBox("Use Dask")
        self.use_dask_cb.setToolTip(
            "Parallelize prediction across tiles with Dask. Automatically "
            "enables 'Tile image' \u2014 Dask only takes effect when tiling "
            "is on (per the napari-convpaint API). Beta feature \u2014 may "
            "not be fully optimized yet."
        )
        self.use_dask_cb.setChecked(
            bool(self.ip.get("convpaint_use_dask", False))
        )
        self.use_dask_cb.toggled.connect(self._on_use_dask_toggled)
        perf_layout.addWidget(self.use_dask_cb)
        perf_layout.addWidget(HelpButton(
            "Use Dask (beta)",
            "Only takes effect when 'Tile image' is also on \u2014 checking "
            "this automatically checks 'Tile image' for you.\n\n"
            "Runs the per-block predictions from Tile image in parallel "
            "worker processes via Dask instead of one block at a time. Can "
            "speed up tiled prediction on machines with several free CPU "
            "cores, since blocks are computed concurrently.\n\n"
            "It does not reduce memory further than Tile image alone \u2014 "
            "several blocks are held in memory at once while they compute, "
            "so peak memory can be higher, not lower.\n\n"
            "Beta feature in napari-convpaint \u2014 "
            "may not be fully optimized."
        ))
        perf_layout.addStretch()
        perf_group.setLayout(perf_layout)
        self.perf_group = perf_group
        root.addWidget(perf_group)

        # Extra toolbar widgets (plugin can inject session-level controls
        # before the tabs).
        for item in self._extra_toolbar_widgets:
            if isinstance(item, QWidget):
                root.addWidget(item)
            elif isinstance(item, (QHBoxLayout, QVBoxLayout, QGridLayout)):
                root.addLayout(item)

        # Tabs: per-cell-type segmentation params, then Dead.
        seg_group = QGroupBox("Annotation & Segmentation")
        seg_layout = QVBoxLayout()
        self.tab_widget = QTabWidget()
        self.tab_widget.setElideMode(Qt.ElideRight)
        self.tab_widget.tabBar().setExpanding(False)
        self.tab_widget.tabBar().setUsesScrollButtons(True)

        per_tab_strategies = list(self.STRATEGIES) if self._per_cell_type_strategy else None
        for ct in self._tab_cell_types:
            tab = CellTypeConvPaintTab(
                cell_type=ct,
                strategy=self.convpaint_strategy,
                initial_params=self.ip,
                on_params_changed=self._on_tab_params_changed,
                run_preview_callback=self._run_instance_preview,
                show_strategy_combo=self._per_cell_type_strategy,
                per_tab_strategies=per_tab_strategies,
                on_per_tab_strategy_changed=self._on_per_tab_strategy_changed,
                pixel_sizes=self._pixel_sizes,
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

        scroll = QScrollArea()
        scroll.setWidget(self.content_widget)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        outer_layout = QVBoxLayout(self)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.addWidget(scroll)

    # ── Callbacks ───────────────────────────────────────────────────

    def _on_tab_params_changed(self, params):
        self.ip.update(params)
        if self._on_params_changed:
            self._on_params_changed(params)

    # ── Strategy resolution / per-cell-type plumbing ─────────────────
    def _resolve_strategy(self, ct):
        """Return the effective ConvPaint strategy for *ct*."""
        if callable(self._strategy_resolver):
            try:
                resolved = self._strategy_resolver(ct)
                if resolved:
                    return _normalize_strategy(resolved)
            except Exception:
                pass
        if self.strategy_combo is not None:
            global_choice = self.strategy_combo.currentText()
            if global_choice == self.ADVANCED_STRATEGY:
                tab = self.tabs.get(ct)
                per_combo = getattr(tab, "_per_tab_strategy_combo", None) if tab else None
                if per_combo is not None:
                    return _normalize_strategy(per_combo.currentText())
                return self.STRATEGIES[0]
            return _normalize_strategy(global_choice)
        return self.convpaint_strategy

    def _strategy_help_content(self, strategy_name):
        if strategy_name == self.ADVANCED_STRATEGY:
            return (
                "ConvPaint Preview Parameters (Advanced \u2014 per cell type)",
                "In Advanced mode each cell-type tab has its own strategy selector.\n\n"
                "Select a strategy inside each tab to show the matching post-processing controls:\n"
                "  \u2022 EDT threshold, Opening, Min size, Fill holes  (EDT/Watershed)\n"
                "  \u2022 Peak min dist/ratio                          (Peak EDT/Watershed)\n"
                "  \u2022 Mask threshold, Seed threshold, Opening, Min size  (Probability + Watershed)"
            )
        if strategy_name == STRATEGY_EDT:
            return (
                "ConvPaint Preview Parameters (Mask + EDT/Watershed)",
                "Use the unified classifier mask, then refine instances with EDT + watershed.\n\n"
                "  \u2022 EDT threshold: splits touching objects (lower = more aggressive).\n"
                "  \u2022 Opening px: smooths boundaries and removes noise.\n"
                "  \u2022 Min size: drops tiny segments below this size threshold.\n"
                "  \u2022 Fill holes: fills internal gaps inside segmented objects."
            )
        if strategy_name == STRATEGY_PEAK_EDT:
            return (
                "ConvPaint Preview Parameters (Mask + Peak EDT/Watershed)",
                "Same mask -> EDT pipeline, but watershed is seeded at local EDT peaks.\n\n"
                "  \u2022 Peak min dist: minimum distance between seed peaks.\n"
                "  \u2022 Peak min ratio: peak height vs. local max (0\u20131)."
            )
        if strategy_name == STRATEGY_PROB:
            return (
                "ConvPaint Preview Parameters (Probability + Watershed)",
                "Use the per-cell-type probability map for watershed seeding.\n\n"
                "  \u2022 Mask threshold: foreground cutoff for the watershed mask.\n"
                "  \u2022 Seed threshold: seed cutoff for watershed seeds.\n"
                "  \u2022 Opening px: smooths boundaries and removes noise.\n"
                "  \u2022 Min size: drops tiny segments below this size threshold."
            )
        return (
            "ConvPaint Preview Parameters",
            "Select a strategy to see its post-processing parameters.",
        )

    def _apply_strategy_to_tabs(self, strategy, emit_signal=True):
        strategy = str(strategy)
        if strategy != self.ADVANCED_STRATEGY:
            self.convpaint_strategy = _normalize_strategy(strategy)
        if self.strategy_help_button is not None:
            t, d = self._strategy_help_content(strategy)
            self.strategy_help_button.set_help(t, d)

        is_advanced = strategy == self.ADVANCED_STRATEGY
        for ct, tab in self.tabs.items():
            per_w = getattr(tab, "_per_tab_strategy_widget", None)
            if per_w is not None:
                per_w.setVisible(is_advanced and ct != "dead")
            effective = self._resolve_strategy(ct)
            tab.rebuild_instance_controls(strategy=effective)

        if emit_signal:
            try:
                self.strategy_changed.emit(strategy)
            except Exception:
                pass

    def _on_global_strategy_changed(self, new_strategy):
        self._apply_strategy_to_tabs(new_strategy, emit_signal=True)
        self._persist_all_params()

    def _on_per_tab_strategy_changed(self, cell_type, new_strategy):
        self.ip[f"per_ct_convpaint_strategy_{cell_type}"] = new_strategy
        self._persist_all_params()

    def _on_use_dask_toggled(self, checked):
        # Dask only takes effect when tiling is on (napari-convpaint API),
        # so enabling it here also enables 'Tile image' rather than leaving
        # the user with a checked-but-inert setting.
        if checked and not self.tile_image_cb.isChecked():
            self.tile_image_cb.setChecked(True)

    # ── Enable/disable + teardown helpers ────────────────────────────
    def set_training_enabled(self, enabled):
        enabled = bool(enabled)
        for name in ("btn_train", "btn_rerun_preview", "btn_save_labels"):
            w = getattr(self, name, None)
            if w is not None:
                w.setEnabled(enabled)
        for tab in self.tabs.values():
            btn = getattr(tab, "btn_preview", None)
            if btn is not None:
                btn.setEnabled(enabled)

    def disable_training_params(self):
        """Disable all training-specific controls in the plugin context.

        Leaves segmentation post-processing parameters (EDT threshold,
        probability thresholds, watershed settings, etc.) and their preview
        buttons fully interactive. Only the feature-extractor settings and
        classifier hyper-parameters (including 'Tile annotations', which
        lives inside clf_group since it only affects training) are disabled,
        along with the Train and Save-Labels buttons.

        'Tile image' and 'Use Dask' (perf_group) are deliberately NOT
        disabled here: they control prediction behavior (Run instance
        segmentation / Re-run preview / batch processing), which stays
        usable even when no training images are loaded — e.g. when working
        with an already-trained/loaded classifier.

        Also hides the training log box (it would remain empty since no
        training runs here) to avoid a large blank area in the widget.
        """
        for group_attr in ("fe_group", "clf_group"):
            g = getattr(self, group_attr, None)
            if g is not None:
                g.setEnabled(False)
        for attr in ("btn_train", "btn_save_labels"):
            w = getattr(self, attr, None)
            if w is not None:
                w.setEnabled(False)

    def cleanup(self):
        """No persistent napari event handlers to disconnect for ConvPaint.

        Provided for API symmetry with :class:`APOCTrainingWidget` so the
        plugin can always call ``cleanup()`` before deleting the widget.
        """

    def _on_save_labels(self):
        self.save_user_labels()

    def _log(self, msg):
        if self._external_log is not None:
            self._external_log(str(msg))
        else:
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
        """Train unified (+ death) ConvPaint classifiers on a background thread.

        Phase 1 (Qt thread) — read viewer layers into NumPy arrays, build the
        ConvpaintModel object (which reads Qt widget values), save labels.
        Phase 2 (worker thread) — model.train, segment_per_z, predict_probas.
        Phase 3 (Qt thread, on_done) — write result arrays to viewer layers,
        persist params, refresh legend.
        """
        if self._bg.is_running():
            self._log("⚠️ Training is already in progress.")
            return

        cell_types = list(self._tab_cell_types)

        # ── Phase 1 (Qt thread): emit signal + gather all data ────────────
        try:
            self.training_started.emit(cell_types)
        except Exception:
            pass

        self._log("Auto-saving user labels before training…")
        self.save_user_labels()

        # Read unified annotation layer on the Qt thread.
        unified_annotations = None
        for lyr in self.viewer.layers:
            if lyr.name == UNIFIED_LABELS_LAYER_NAME:
                unified_annotations = np.asarray(lyr.data)
                break
        if unified_annotations is None:
            self._log(f"❌ Label layer '{UNIFIED_LABELS_LAYER_NAME}' not found!")
            try:
                self.training_finished.emit(cell_types, None)
            except Exception:
                pass
            return

        # Read death annotation layer if needed (Qt thread).
        death_annotations = None
        if self.has_death:
            for lyr in self.viewer.layers:
                if lyr.name == DEAD_LABELS_LAYER_NAME:
                    death_annotations = np.asarray(lyr.data)
                    break

        # Build model NOW on Qt thread (reads widget values: fe_combo, etc.).
        try:
            model, fe_device = self._build_model()
        except Exception as e:
            self._log(f"❌ Failed to build model: {e}")
            try:
                self.training_finished.emit(cell_types, None)
            except Exception:
                pass
            return

        # Snapshot all other state needed by the worker.
        unified_input_channels = list(self.unified_input_channels)
        death_input_channels   = list(self.death_input_channels)
        label_map              = dict(self.label_map)
        all_cell_types         = list(self.all_cell_types)
        has_death              = self.has_death
        pixel_class_outdir     = Path(self.pixel_class_outdir)
        all_images             = list(self.all_images)  # references; data is lazy
        convpaint_strategy     = self.convpaint_strategy
        fe_alias               = self.fe_combo.currentData()  # for death model

        # Snapshot per-tab params (Qt thread only). Distance/volume spinboxes
        # (edt_threshold_spin, segment_size_min_spin, peak_min_distance_spin)
        # go through get_native() so the training run always uses the
        # canonical native (px/voxel) value regardless of the tab's current
        # px/µm display setting.
        _native_attrs = {"edt_threshold_spin", "segment_size_min_spin", "peak_min_distance_spin"}
        tabs_params = {}
        for ct, tab in self.tabs.items():
            p = {}
            for attr in (
                "edt_threshold_spin", "opening_nr_pixels_spin",
                "segment_size_min_spin", "peak_min_distance_spin",
                "peak_min_ratio_spin", "prob_mask_threshold_spin",
                "prob_seed_threshold_spin",
            ):
                w = getattr(tab, attr, None)
                if w is not None:
                    p[attr] = tab.get_native(w) if attr in _native_attrs else w.value()
            w2 = getattr(tab, "fill_holes_cb", None)
            if w2 is not None:
                p["fill_holes"] = bool(w2.isChecked())
            tabs_params[ct] = p

        safe_log = ThreadSafeLogger(self._log)

        # ── Phase 2 (worker thread) ───────────────────────────────────────
        def _do_train(progress_cb=None):
            n_steps = 2 if has_death else 1
            result = {}

            # ── Unified classifier ──────────────────────────────────────
            if progress_cb:
                progress_cb(0, n_steps, "Training unified classifier…")

            # Determine active cell types from annotations.
            unique_vals = set(
                int(v) for v in np.unique(unified_annotations) if int(v) > 0
            )
            active_cell_types = [
                ct for ct in all_cell_types
                if int(label_map["celltype_to_label"][ct]) in unique_vals
            ]
            skipped = [ct for ct in all_cell_types if ct not in active_cell_types]
            if skipped:
                safe_log(
                    f"⚠️ No annotations for: {', '.join(skipped)}. "
                    "Those cell types will not be predicted."
                )
            if not active_cell_types:
                safe_log("❌ No annotated cell types — aborting training.")
                result["error"] = "no_active_cell_types"
                return result

            # Remap annotations to dense indices.
            active_label_map = build_label_map(active_cell_types)
            remap = {0: 0, BACKGROUND_LABEL: BACKGROUND_LABEL}
            for ct in active_cell_types:
                old_idx = int(label_map["celltype_to_label"][ct])
                new_idx = int(active_label_map["celltype_to_label"][ct])
                remap[old_idx] = new_idx
            remapped = np.zeros_like(unified_annotations)
            for old_v, new_v in remap.items():
                if old_v != 0:
                    remapped[unified_annotations == old_v] = new_v

            train_images, train_annots = [], []
            for i, img in enumerate(all_images):
                fa = remapped[i]
                if np.any(fa > 0):
                    train_images.append(
                        _slice_image_channels(img, unified_input_channels)
                    )
                    train_annots.append(fa)

            if not train_images:
                safe_log("⚠️ No annotated frames — aborting training.")
                result["error"] = "no_annotated_frames"
                return result

            safe_log(
                f"Training unified classifier on {len(train_images)} frames "
                f"({len(active_cell_types)} cell types: {active_cell_types})…"
            )
            model.train(
                train_images, train_annots,
                fe_use_device=fe_device, clf_use_device=fe_device,
            )

            # Save model + label-map sidecar.
            model_path = unified_model_path(pixel_class_outdir)
            map_path   = unified_label_map_path(pixel_class_outdir)
            model.save(str(model_path))
            save_label_map(map_path, active_label_map)
            save_input_channels(
                unified_input_channels_path(pixel_class_outdir),
                unified_input_channels,
                model_name="unified",
            )
            safe_log(f"✅ Saved model: {model_path.name}")

            # Run inference on the training stack (pure compute).
            safe_log("Running unified inference on the training stack…")
            seg_stack = np.stack(
                [
                    segment_per_z(
                        model,
                        _slice_image_channels(img, unified_input_channels),
                        fe_use_device=fe_device,
                    )
                    for img in all_images
                ],
                axis=0,
            ).astype(np.int16)
            _save_preview_array(
                unified_predicted_labels_path(pixel_class_outdir), seg_stack
            )

            # Per-cell-type probability maps (when needed by strategy).
            needs_prob = convpaint_strategy == STRATEGY_PROB or any(
                tabs_params.get(ct, {}).get("prob_mask_threshold_spin") is not None
                for ct in active_cell_types
            )
            per_ct_prob_stacks = {}
            if needs_prob:
                per_ct_stacks = {ct: [] for ct in active_cell_types}
                for img in all_images:
                    probas = np.asarray(predict_probas_per_z(
                        model,
                        _slice_image_channels(img, unified_input_channels),
                        fe_use_device=fe_device,
                    ))
                    for ct in active_cell_types:
                        k = int(active_label_map["celltype_to_label"][ct])
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
                        _probability_map_path(pixel_class_outdir, ct), stack
                    )
                    per_ct_prob_stacks[ct] = stack
            else:
                # EDT mode: remove stale probability files.
                for ct in active_cell_types:
                    p = _probability_map_path(pixel_class_outdir, ct)
                    if p is not None and p.exists():
                        shutil.rmtree(p)

            # Instance segmentation previews (pure compute, use snapshotted params).
            instance_previews = {}  # ct -> array
            for ct in active_cell_types:
                tp = tabs_params.get(ct, {})
                effective = _normalize_strategy(
                    self._resolve_strategy(ct)  # reads strategy_combo; OK if readonly
                ) if False else convpaint_strategy  # use global strategy as fallback

                k = int(active_label_map["celltype_to_label"][ct])
                mask_stack = (seg_stack == k).astype(np.uint8)

                if effective == STRATEGY_PROB and ct in per_ct_prob_stacks:
                    instances = probability_to_instances(
                        per_ct_prob_stacks[ct],
                        mask_thr=tp.get("prob_mask_threshold_spin", 0.5),
                        seed_thr=tp.get("prob_seed_threshold_spin", 0.5),
                        opening_nr_pixels=int(tp.get("opening_nr_pixels_spin", 0)),
                        segment_size_min=int(tp.get("segment_size_min_spin", 10)),
                    )
                else:
                    instances = mask_to_instances(
                        mask_stack,
                        edt_thr=float(tp.get("edt_threshold_spin", 1.0)),
                        opening_nr_pixels=int(tp.get("opening_nr_pixels_spin", 0)),
                        fill_holes=bool(tp.get("fill_holes", True)),
                        segment_size_min=int(tp.get("segment_size_min_spin", 10)),
                        marker_strategy=(
                            "peak" if effective == STRATEGY_PEAK_EDT else "threshold"
                        ),
                        peak_min_distance=float(tp.get("peak_min_distance_spin", 2.0))
                        if "peak_min_distance_spin" in tp else None,
                        peak_min_ratio=float(tp.get("peak_min_ratio_spin", 0.35))
                        if "peak_min_ratio_spin" in tp else 0.35,
                    )
                # Save to disk (I/O fine in worker).
                pred_path = _predicted_labels_path(pixel_class_outdir, ct)
                if pred_path is not None:
                    _save_preview_array(pred_path, instances)
                instance_previews[ct] = instances

            result["unified"] = {
                "active_label_map": active_label_map,
                "seg_stack": seg_stack,
                "per_ct_prob_stacks": per_ct_prob_stacks,
                "instance_previews": instance_previews,
                "active_cell_types": active_cell_types,
            }

            # ── Death classifier (optional) ─────────────────────────────
            if has_death and death_annotations is not None and np.any(death_annotations):
                if progress_cb:
                    progress_cb(1, n_steps, "Training death classifier…")
                safe_log("Training death classifier…")

                train_d_images, train_d_annots = [], []
                for i, img in enumerate(all_images):
                    if np.any(death_annotations[i] > 0):
                        train_d_images.append(
                            _slice_image_channels(img, death_input_channels)
                        )
                        train_d_annots.append(death_annotations[i])

                if train_d_images:
                    try:
                        # Build a second model with the same params (already
                        # saved the unified model above so model.save is safe).
                        from napari_convpaint import ConvpaintModel
                        d_model = ConvpaintModel(fe_alias)
                        try:
                            d_model.set_params(**model.params)
                        except Exception:
                            pass  # older ConvpaintModel may not support all params
                        d_model.train(
                            train_d_images, train_d_annots,
                            fe_use_device=fe_device, clf_use_device=fe_device,
                        )
                        mp = _death_model_path(pixel_class_outdir)
                        d_model.save(str(mp))
                        save_input_channels(
                            death_input_channels_path(pixel_class_outdir),
                            death_input_channels,
                            model_name="death",
                        )
                        safe_log(f"✅ Saved death model: {mp.name}")

                        # Death preview
                        d_seg_stack = np.stack(
                            [
                                segment_per_z(
                                    d_model,
                                    _slice_image_channels(img, death_input_channels),
                                    fe_use_device=fe_device,
                                )
                                for img in all_images
                            ],
                            axis=0,
                        ).astype(np.int16)
                        death_mask = (d_seg_stack >= 2).astype(np.uint16)
                        pred_path_d = _predicted_labels_path(pixel_class_outdir, "dead")
                        if pred_path_d is not None:
                            _save_preview_array(pred_path_d, death_mask)
                        result["death"] = death_mask
                    except Exception as e:
                        safe_log(f"❌ Death training failed: {e}")
                        import traceback
                        traceback.print_exc()
                else:
                    safe_log("⚠️ No annotations for death — skipping.")
            elif has_death:
                safe_log("⚠️ Death layer missing or empty — skipping death training.")

            return result

        # ── Phase 3 (Qt thread, on_done): write results to viewer ─────────
        def _on_done(result):
            if not result or "error" in result:
                self._log("❌ Training produced no results.")
                try:
                    self.training_finished.emit(cell_types, None)
                except Exception:
                    pass
                return

            # Apply unified results.
            u = result.get("unified", {})
            active_label_map = u.get("active_label_map", {})
            seg_stack        = u.get("seg_stack")
            per_ct_prob      = u.get("per_ct_prob_stacks", {})
            instance_prev    = u.get("instance_previews", {})
            active_cts       = u.get("active_cell_types", [])

            # Update the live label map so downstream tabs decode correctly.
            if active_label_map:
                self.label_map = active_label_map

            if seg_stack is not None:
                self._set_labels_layer(
                    UNIFIED_PREDICTED_LAYER_NAME, seg_stack, visible=False
                )

            for ct, prob_stack in per_ct_prob.items():
                self._set_image_layer(
                    _probability_layer_name(ct), prob_stack, visible=False
                )

            for ct, instances in instance_prev.items():
                self._set_labels_layer(
                    _segments_layer_name(ct), instances
                )
                prob_stack = per_ct_prob.get(ct)
                if prob_stack is not None:
                    self._set_image_layer(
                        _probability_layer_name(ct), prob_stack, visible=False
                    )

            # Hide probability layers for non-active cell types.
            for ct in self.all_cell_types:
                if ct not in active_cts:
                    lname = _probability_layer_name(ct)
                    if lname in self.viewer.layers:
                        self.viewer.layers[lname].visible = False

            # Apply death result.
            death_mask = result.get("death")
            if death_mask is not None:
                self._set_labels_layer(DEAD_PREDICTED_LAYER_NAME, death_mask)
                _set_dead_mask_layer_color(self.viewer.layers[DEAD_PREDICTED_LAYER_NAME])

            _reorder_convpaint_layers(
                self.viewer, self.all_cell_types, has_death=self.has_death
            )
            self._log("✅ Unified inference + previews complete.")

            self._persist_all_params()
            if self.legend_tab is not None:
                self.legend_tab.refresh()
            try:
                self.training_finished.emit(cell_types, None)
            except Exception:
                pass

        def _on_failed(err):
            self._log(f"❌ Training failed: {err}")
            try:
                self.training_finished.emit(cell_types, None)
            except Exception:
                pass

        self._bg.run(
            fn=_do_train,
            desc="ConvPaint training…",
            buttons=[self.btn_train],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
            inject_progress=True,
            indeterminate=False,
        )

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
            self.instance_preview_started.emit(cell_type)
        except Exception:
            pass
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

            effective = self._resolve_strategy(cell_type)
            if effective == STRATEGY_PROB:
                self._preview_probability(cell_type, tab)
            else:
                # Pass the per-cell-type strategy so EDT/Peak-EDT branching
                # is honored when in Advanced mode.
                self._preview_edt(cell_type, tab, strategy=effective)

            _reorder_convpaint_layers(
                self.viewer, self.all_cell_types, has_death=self.has_death,
            )
            self._log(f"\u2705 Preview updated for {cell_type}")
        except Exception as e:
            self._log(f"\u274c Preview failed for {cell_type}: {e}")
            traceback.print_exc()
        finally:
            try:
                self.instance_preview_finished.emit(cell_type)
            except Exception:
                pass

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

    def _preview_edt(self, cell_type, tab, strategy=None):
        unified = self._load_unified_predicted_labels()
        if unified is None:
            self._log(
                "\u26a0\ufe0f No cached unified labels — train the classifier first."
            )
            return
        k = int(self.label_map["celltype_to_label"][cell_type])
        mask_stack = (unified == k).astype(np.uint8)

        effective_strategy = _normalize_strategy(strategy or self.convpaint_strategy)
        instances = mask_to_instances(
            mask_stack,
            edt_thr=float(tab.get_native(tab.edt_threshold_spin)),
            opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
            fill_holes=bool(tab.fill_holes_cb.isChecked())
            if tab.fill_holes_cb is not None else True,
            segment_size_min=int(round(tab.get_native(tab.segment_size_min_spin))),
            marker_strategy=(
                "peak" if effective_strategy == STRATEGY_PEAK_EDT
                else "threshold"
            ),
            peak_min_distance=float(tab.get_native(tab.peak_min_distance_spin)) if tab.peak_min_distance_spin is not None else None,
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
            segment_size_min=int(round(tab.get_native(tab.segment_size_min_spin))),
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
        _set_dead_mask_layer_color(self.viewer.layers[DEAD_PREDICTED_LAYER_NAME])
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

    # Resolve pixel sizes (µm) from metadata so the per-tab µm/pixel unit
    # toggle (EDT threshold / min size / peak distance) is available here
    # too, matching the napari plugin's own ConvPaintWidget.
    from behav3d.core.utils import convert_distance
    if metadata is not None and "pixel_distance_xy" in metadata.columns and "distance_unit" in metadata.columns:
        _unit = str(metadata["distance_unit"].iloc[0])
        _xy_from_md = convert_distance(float(metadata["pixel_distance_xy"].iloc[0]), _unit)
        _z_from_md = convert_distance(float(metadata["pixel_distance_z"].iloc[0]), _unit)
    else:
        _xy_from_md = None
        _z_from_md = None
    pixel_sizes = {
        "xy_um": ip.get("pixel_size_xy") or ip.get("pixel_size_xy_um") or _xy_from_md,
        "z_um": ip.get("pixel_size_z") or ip.get("pixel_size_z_um") or _z_from_md,
    }

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
                    dead_pred_layer = viewer.add_labels(
                        pred, name=DEAD_PREDICTED_LAYER_NAME, opacity=0.8,
                        visible=False,
                    )
                    _set_dead_mask_layer_color(dead_pred_layer)
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
        show_legend=True,
        pixel_sizes=pixel_sizes,
    )
    viewer.window.add_dock_widget(widget, name="ConvPaint", area="right")
    _reorder_convpaint_layers(viewer, all_cell_types, has_death=has_death)
    return viewer
