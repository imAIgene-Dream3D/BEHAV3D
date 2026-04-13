"""
Custom APOC (Accelerated Pixel and Object Classification) widget for BEHAV3D.

This file is self-contained: it handles image loading, the Qt training widget
(with per-cell-type tabs), and APOC train/predict calls.
The original pixel classifier in napari_pixelclassifier.py is left completely untouched.
"""

import os
import gc
import json
import re
import time
import shutil
from pathlib import Path

import numpy as np
import napari
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QLabel, QComboBox,
    QSpinBox, QDoubleSpinBox, QLineEdit, QPushButton as QtPushButton,
    QCheckBox, QGroupBox, QPlainTextEdit, QApplication, QScrollArea,
    QSizePolicy, QTabWidget, QFrame, QMessageBox,
    QDialog, QTableWidget, QTableWidgetItem, QHeaderView,
)
from qtpy.QtCore import Qt
from qtpy.QtGui import QColor

import dask.array as da
import pyclesperanto_prototype as cle

from behav3d.io.images import load_image, load_zarr, save_as_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape

# ---------------------------------------------------------------------------
# Feature grid constants (matching the official napari-apoc widget)
# ---------------------------------------------------------------------------

# Sigma columns (these match the default column headers in the official widget)
APOC_SIGMAS = [0.3, 0.5, 1, 2, 3, 4, 5, 10, 15, 25]

# Feature rows shown in presets
APOC_FEATURES = [
    ("gaussian_blur",               "Gauss"),
    ("difference_of_gaussian",      "DoG"),
    ("laplace_box_of_gaussian_blur","LoG"),
    ("sobel_of_gaussian_blur",      "SoG"),
]
# For custom preset every feature is available (same list here; extend if needed)
APOC_ALL_FEATURES = APOC_FEATURES

# Format sigma values nicely (drop trailing .0)
def _fmt_sigma(s):
    s_float = float(s)
    return str(int(s_float)) if s_float == int(s_float) else str(s_float)


# Each preset is a set of (feature_key, sigma_str) pairs that should be
# pre-checked.  The sigma value is stored as a string matching _fmt_sigma().
FEATURE_PRESETS = {
    "small_quick": {
        "label": "Small / Quick",
        "checked": {
            ("difference_of_gaussian",      "1"),
            ("difference_of_gaussian",      "2"),
            ("laplace_box_of_gaussian_blur","1"),
            ("laplace_box_of_gaussian_blur","2"),
        },
    },
    "medium_quick": {
        "label": "Medium / Quick",
        "checked": {
            ("difference_of_gaussian",      "1"),
            ("difference_of_gaussian",      "2"),
            ("difference_of_gaussian",      "5"),
            ("laplace_box_of_gaussian_blur","1"),
            ("laplace_box_of_gaussian_blur","2"),
            ("laplace_box_of_gaussian_blur","5"),
            ("sobel_of_gaussian_blur",      "1"),
            ("sobel_of_gaussian_blur",      "2"),
            ("sobel_of_gaussian_blur",      "5"),
        },
    },
    "large_quick": {
        "label": "Large / Quick",
        "checked": {
            ("gaussian_blur",               "1"),
            ("gaussian_blur",               "2"),
            ("gaussian_blur",               "5"),
            ("gaussian_blur",               "10"),
            ("gaussian_blur",               "25"),
            ("difference_of_gaussian",      "1"),
            ("difference_of_gaussian",      "2"),
            ("difference_of_gaussian",      "5"),
            ("difference_of_gaussian",      "10"),
            ("difference_of_gaussian",      "25"),
            ("laplace_box_of_gaussian_blur","1"),
            ("laplace_box_of_gaussian_blur","2"),
            ("laplace_box_of_gaussian_blur","5"),
            ("laplace_box_of_gaussian_blur","10"),
            ("laplace_box_of_gaussian_blur","25"),
            ("sobel_of_gaussian_blur",      "1"),
            ("sobel_of_gaussian_blur",      "2"),
            ("sobel_of_gaussian_blur",      "5"),
            ("sobel_of_gaussian_blur",      "10"),
            ("sobel_of_gaussian_blur",      "25"),
        },
    },
    "custom": {
        "label": "Custom",
        "checked": set(),  # user starts with all unchecked
    },
}


def _checked_set_for_preset(preset_name):
    """Return the set of (feature_key, sigma_str) pairs that should be checked for a preset."""
    return set(FEATURE_PRESETS.get(preset_name, FEATURE_PRESETS["medium_quick"])["checked"])


def _build_feature_string_from_checked(checked_set, consider_original=False, current_sigmas=None):
    """Build an APOC feature string from a set of (feature_key, sigma_str) pairs."""
    parts = []
    if consider_original:
        parts.append("original")
    if current_sigmas is None:
        current_sigmas = APOC_SIGMAS
    sigma_strs = [_fmt_sigma(s) for s in current_sigmas]
    # Keep deterministic order: rows in APOC_ALL_FEATURES order, sigmas in current_sigmas order
    for feat_key, _label in APOC_ALL_FEATURES:
        for s_str in sigma_strs:
            if (feat_key, s_str) in checked_set:
                parts.append(f"{feat_key}={s_str}")
    return " ".join(parts)


def _default_grid_sigmas():
    """Return the default APOC sigma grid as display strings."""
    return [_fmt_sigma(s) for s in APOC_SIGMAS]


def _default_grid_sigmas_text():
    """Return the default APOC sigma grid as a text field value."""
    return ", ".join(_default_grid_sigmas())


def _parse_feature_string(feature_string):
    """Convert an APOC feature string into the grid config used by the widget."""
    feature_string = str(feature_string or "").replace(",", " ").replace("\t", " ").strip()
    while "  " in feature_string:
        feature_string = feature_string.replace("  ", " ")

    tokens = [tok for tok in feature_string.split(" ") if tok]
    checked = []
    sigmas = []
    consider_original = False
    valid_features = {feat_key for feat_key, _ in APOC_ALL_FEATURES}

    for token in tokens:
        token = str(token).strip()
        lower = token.lower()
        if lower == "original" or re.fullmatch(r"original_channel\d+", lower):
            consider_original = True
            continue

        feat_key = None
        sigma_text = None
        if "=" in token:
            feat_key, sigma_text = token.split("=", 1)
        else:
            readable = re.sub(r"_channel\d+$", "", token)
            match = re.fullmatch(r"(.+)_sigma([^_]+)", readable)
            if match:
                feat_key, sigma_text = match.groups()

        if feat_key not in valid_features or sigma_text is None:
            continue
        try:
            sigma_text = _fmt_sigma(float(sigma_text))
        except ValueError:
            continue
        checked.append((feat_key, sigma_text))
        if sigma_text not in sigmas:
            sigmas.append(sigma_text)

    checked_set = set(checked)
    feature_preset = "custom"
    if not consider_original:
        for preset_name, preset_cfg in FEATURE_PRESETS.items():
            if preset_name == "custom":
                continue
            if checked_set == set(preset_cfg["checked"]):
                feature_preset = preset_name
                break

    canonical_feature_string = _build_feature_string_from_checked(
        checked_set,
        consider_original=consider_original,
        current_sigmas=sigmas or _default_grid_sigmas(),
    )

    return {
        "feature_string": canonical_feature_string,
        "checked_features": [list(pair) for pair in checked],
        "sigmas": ",".join(sigmas),
        "grid_sigmas": ", ".join(sigmas),
        "consider_original": consider_original,
        "feature_preset": feature_preset,
    }


def _read_classifier_header_value(opencl_path, key, default=None):
    """Read a single APOC metadata value from the header of a trained .cl file."""
    path = Path(opencl_path)
    if not path.exists():
        return default

    prefix = f"{key} = "
    with path.open() as f:
        line = ""
        count = 0
        while line != "*/" and count < 50:
            line = f.readline()
            if not line:
                break
            count += 1
            line = line.strip()
            if line.startswith(prefix):
                return line[len(prefix):].strip()
    return default


def _classifier_path(pixel_class_outdir, cell_type):
    """Return the expected APOC classifier path for a cell type."""
    if not pixel_class_outdir:
        return None
    if cell_type == "dead":
        fname = "PixelClassifier_Death.cl"
    else:
        fname = f"PixelClassifier_{cell_type.capitalize()}.cl"
    return Path(pixel_class_outdir) / fname


def _predicted_labels_path(pixel_class_outdir, cell_type):
    """Return the expected predicted-label zarr path for a cell type."""
    if not pixel_class_outdir:
        return None
    if cell_type == "dead":
        fname = "PixelClassifier_Death_PredictedLabels.zarr"
    else:
        fname = f"PixelClassifier_{cell_type.capitalize()}_PredictedLabels.zarr"
    return Path(pixel_class_outdir) / fname


def _probability_map_path(pixel_class_outdir, cell_type):
    """Return the expected probability-map zarr path for a cell type."""
    if not pixel_class_outdir:
        return None
    if cell_type == "dead":
        fname = "PixelClassifier_Death_ProbabilityMap.zarr"
    else:
        fname = f"PixelClassifier_{cell_type.capitalize()}_ProbabilityMap.zarr"
    return Path(pixel_class_outdir) / fname


def _reorder_apoc_layers(viewer, all_cell_types, has_death=False):
    """Group APOC layers by cell type with segments above probability above labels."""
    layer_names = [layer.name for layer in viewer.layers]
    channel_names = [name for name in layer_names if name.startswith("Channel ")]
    ordered_names = list(channel_names)

    for cell_type in all_cell_types:
        ordered_names.extend([
            f"User Provided Labels ({cell_type.capitalize()})",
            f"Probability Map ({cell_type.capitalize()})",
            f"{cell_type.capitalize()} Segments",
        ])

    if has_death:
        ordered_names.extend([
            "User Provided Labels (Dead)",
            "Probability Map (Dead)",
            "Pixel Classification (Dead)",
        ])

    target_names = [name for name in ordered_names if name in layer_names]
    trailing_names = [name for name in layer_names if name not in target_names]
    final_order = target_names + trailing_names

    for target_idx, name in enumerate(final_order):
        current_idx = viewer.layers.index(name)
        if current_idx != target_idx:
            viewer.layers.move(current_idx, target_idx)


def _load_classifier_restore_config(pixel_class_outdir, cell_type):
    """Recover APOC widget defaults from a trained classifier file when available."""
    clf_path = _classifier_path(pixel_class_outdir, cell_type)
    if clf_path is None or not clf_path.exists():
        return {}

    cfg = {}
    feature_string = _read_classifier_header_value(clf_path, "feature_specification")
    if feature_string:
        cfg.update(_parse_feature_string(feature_string))

    max_depth = _read_classifier_header_value(clf_path, "max_depth")
    if max_depth is not None:
        try:
            cfg["max_depth"] = int(max_depth)
        except (TypeError, ValueError):
            pass

    num_trees = _read_classifier_header_value(clf_path, "num_trees")
    if num_trees is not None:
        try:
            cfg["num_ensembles"] = int(num_trees)
        except (TypeError, ValueError):
            pass

    return cfg


def _extract_cell_type_config(initial_params, cell_type):
    """Gather APOC settings for one cell type from flat initial params."""
    cfg = {}
    prefix = f"apoc_{cell_type}_"
    for key, value in (initial_params or {}).items():
        if key.startswith(prefix):
            cfg[key[len(prefix):]] = value
    return cfg


def _finalize_segments(segments, segment_size_min):
    """Apply shared cleanup to a labeled segmentation volume."""
    from behav3d.preprocessing.segmentation import segment_size_filter, segment_2d_filter

    segments = np.asarray(segments).copy()
    segments = segment_size_filter(segments, size_min=segment_size_min)
    if segments.ndim == 3:
        segments = segment_2d_filter(segments)
    return np.asarray(segments).astype(np.uint16)


def _probability_volume_to_segments(prob_map, mask_thr, seed_thr, opening_nr_pixels, segment_size_min):
    """Convert a single probability-map volume to instances."""
    from scipy import ndimage as ndi
    from skimage.measure import label as sk_label
    from skimage.segmentation import watershed
    from behav3d.preprocessing.segmentation.segmentation_utils import postprocess_mask

    prob_map = np.asarray(prob_map).astype(np.float32)
    mask_out = postprocess_mask(
        prob_map > mask_thr,
        fill_holes=False,
        opening_nr_pixels=opening_nr_pixels,
    ).astype(bool)

    cc_labels = sk_label(mask_out)
    seed_mask = (prob_map > seed_thr) & mask_out
    segments = np.zeros_like(cc_labels, dtype=np.uint16)
    next_id = 0

    for comp_idx, slc in enumerate(ndi.find_objects(cc_labels), start=1):
        if slc is None:
            continue
        comp_mask = cc_labels[slc] == comp_idx
        sub_seeds = sk_label(seed_mask[slc] & comp_mask)
        n_seeds = int(sub_seeds.max())
        if n_seeds <= 1:
            next_id += 1
            segments[slc][comp_mask] = next_id
            continue

        sub_result = watershed(-prob_map[slc], markers=sub_seeds, mask=comp_mask)
        for seed_id in range(1, int(sub_result.max()) + 1):
            next_id += 1
            segments[slc][sub_result == seed_id] = next_id

    return _finalize_segments(segments, segment_size_min)


def _probability_array_to_segments(prob_map, mask_thr, seed_thr, opening_nr_pixels, segment_size_min):
    """Convert a 3D/4D probability map array to instances."""
    prob_map = np.asarray(prob_map)
    if prob_map.ndim == 4:
        return np.stack([
            _probability_volume_to_segments(prob_map[t], mask_thr, seed_thr, opening_nr_pixels, segment_size_min)
            for t in range(prob_map.shape[0])
        ], axis=0)
    return _probability_volume_to_segments(prob_map, mask_thr, seed_thr, opening_nr_pixels, segment_size_min)


def _mask_array_to_segments(mask, edt_thr, opening_nr_pixels, segment_size_min, fill_holes):
    """Convert a binary mask array to instances using the EDT workflow."""
    from behav3d.preprocessing.segmentation.segmentation_utils import postprocess_mask, segment_mask

    mask = np.asarray(mask)

    def _segment_one(mask_volume):
        proc_mask = postprocess_mask(
            mask_volume.astype(bool),
            fill_holes=fill_holes,
            opening_nr_pixels=opening_nr_pixels,
        )
        return np.asarray(
            segment_mask(
                proc_mask,
                edt_thr=edt_thr,
                edt_thr_refined=None,
                segment_size_min=segment_size_min,
                use_dims=3,
                n_workers=1,
            )
        ).astype(np.uint16)

    if mask.ndim == 4:
        return np.stack([_segment_one(mask[t]) for t in range(mask.shape[0])], axis=0)
    return _segment_one(mask)


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
            
            cached = load_zarr(image_outpath)
            # The cached data is (C, T, Z, Y, X)
            # We need to return a list of (C, Z, Y, X) arrays, one per timepoint
            all_images = [np.asarray(cached[:, t, :, :, :]) for t in range(cached.shape[1])]
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
        # load_image already normalizes to TCZYX (default_axis_order in images.py)
        img = load_image(raw_image_path, axis_order=axis_order)

        # Since load_image returned TCZYX, T is always axis 0
        t_axis = 0
            
        n_timepoints = img.shape[t_axis]
        n_to_select = min(examples_per_sample, n_timepoints)

        # Select equally spaced timepoints (first, last, and middle)
        if n_to_select <= 1:
            t_indices = [0]
        else:
            t_indices = np.round(np.linspace(0, n_timepoints - 1, n_to_select)).astype(int)
            t_indices = sorted(list(set(t_indices)))
            
        t_indices_list = [int(t) for t in t_indices]
        print(f"  {sample_name}: selected {len(t_indices)} equidistant timepoints {t_indices_list}")

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
    # Transpose to (C, T, Z, Y, X) to match original classifier
    all_images_stack = all_images_stack.transpose(1, 0, 2, 3, 4)
    save_as_zarr(all_images_stack, image_outpath)
    del all_images_stack
    gc.collect()

    # Re-load to ensure we have the right shape
    cached = load_zarr(image_outpath)
    # The cached data is (C, T, Z, Y, X)
    # We need to return a list of (C, Z, Y, X) arrays, one per timepoint
    all_images = [np.asarray(cached[:, t, :, :, :]) for t in range(cached.shape[1])]

    return all_images, pixel_class_outdir, has_death, all_cell_types


# ---------------------------------------------------------------------------
# Per-cell-type tab panel
# ---------------------------------------------------------------------------

class CellTypeTab(QWidget):
    """
    Widget for a single cell type tab.
    Contains:
      - channel selection group
      - preset dropdown (kept)
      - collapsible "Tune Features" QGroupBox with a sigma × feature checkbox grid
      - "Consider original image as well" standalone checkbox
      - "Show classifier statistics" button
      - max_depth / num_ensembles RF parameters
    """

    def __init__(
        self,
        cell_type,
        viewer,
        pixel_class_outdir=None,
        initial_params=None,
        apoc_strategy="APOC (Direct Instance Segmentation)",
        on_params_changed=None,
        run_instance_callback=None,
        parent=None,
    ):
        super().__init__(parent)
        self.cell_type = cell_type
        self.viewer = viewer
        self._pixel_class_outdir = pixel_class_outdir
        self._on_params_changed = on_params_changed
        self._run_instance_callback = run_instance_callback
        self._apoc_strategy = str(apoc_strategy)
        root_params = initial_params or {}
        saved_cfg = _extract_cell_type_config(initial_params, cell_type)
        trained_cfg = _load_classifier_restore_config(pixel_class_outdir, cell_type)
        cfg = dict(saved_cfg)
        cfg.update({key: value for key, value in trained_cfg.items() if key != "grid_sigmas"})
        cfg["grid_sigmas"] = (
            saved_cfg.get("grid_sigmas")
            or trained_cfg.get("grid_sigmas")
            or _default_grid_sigmas_text()
        )
        ip = {f"apoc_{cell_type}_{key}": value for key, value in cfg.items()}
        self._default_channel_names = list(saved_cfg.get("channels", []))

        layout = QVBoxLayout()
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # ── Channel selection ────────────────────────────────────────────────
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

        # ── Preset dropdown ──────────────────────────────────────────────────
        feat_row = QHBoxLayout()
        feat_row.addWidget(QLabel("Preset:"))
        self.feature_combo = QComboBox()
        self.feature_combo.addItems(list(FEATURE_PRESETS.keys()))
        saved_preset = ip.get(f"apoc_{cell_type}_feature_preset", "medium_quick")
        if saved_preset not in FEATURE_PRESETS:
            saved_preset = "medium_quick"
        self.feature_combo.setCurrentText(saved_preset)
        feat_row.addWidget(self.feature_combo)
        feat_row.addStretch()
        layout.addLayout(feat_row)

        # ── Tune Features collapsible group ─────────────────────────────────
        self.tune_group = QGroupBox("Tune Features")
        self.tune_group.setCheckable(True)
        self.tune_group.setChecked(False)  # collapsed by default
        self.tune_layout = QVBoxLayout()
        self.tune_layout.setContentsMargins(4, 4, 4, 4)
        self.tune_layout.setSpacing(2)

        # ── Sigmas Input Row ──
        sigma_row = QHBoxLayout()
        sigma_row.addWidget(QLabel("Custom Sigmas:"))
        self.sigma_input = QLineEdit(", ".join(_fmt_sigma(s) for s in APOC_SIGMAS))
        self.sigma_input.setToolTip("Comma or space-separated list of sigmas to use for APOC filtering.")
        sigma_row.addWidget(self.sigma_input)
        self.update_grid_btn = QtPushButton("Update Grid")
        self.update_grid_btn.clicked.connect(self._on_update_grid)
        sigma_row.addWidget(self.update_grid_btn)
        self.tune_layout.addLayout(sigma_row)

        self.current_sigmas = list(APOC_SIGMAS)
        self._feat_sigma_checks = {}
        self._grid_widget = None

        self.tune_group.setLayout(self.tune_layout)
        layout.addWidget(self.tune_group)

        # ── Consider original image checkbox ─────────────────────────────────
        self.consider_original_cb = QCheckBox("Consider original image as well")
        saved_orig = bool(ip.get(f"apoc_{cell_type}_consider_original", False))
        self.consider_original_cb.setChecked(saved_orig)
        self.consider_original_cb.stateChanged.connect(self._update_preview)
        layout.addWidget(self.consider_original_cb)

        # ── Feature string preview ───────────────────────────────────────────
        self.preview_label = QLabel("")
        self.preview_label.setWordWrap(True)
        self.preview_label.setStyleSheet(
            "color: #888; font-size: 10px; padding: 2px 4px; "
            "background: rgba(0,0,0,0.05); border-radius: 3px;"
        )
        layout.addWidget(self.preview_label)

        # ── Show classifier statistics button ────────────────────────────────
        self.stats_btn = QtPushButton("Show classifier statistics")
        self.stats_btn.setToolTip(
            "Display the feature importance distribution of the trained classifier."
        )
        self.stats_btn.clicked.connect(self._on_show_statistics)
        layout.addWidget(self.stats_btn)

        # ── RF parameters ────────────────────────────────────────────────────
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
        rf_row.addStretch()
        layout.addLayout(rf_row)

        self.run_instance_btn = None
        self.instance_group = None
        self.prob_mask_threshold_spin = None
        self.prob_seed_threshold_spin = None
        self.edt_threshold_spin = None
        self.segment_size_min_spin = None
        self.opening_nr_pixels_spin = None
        self.fill_holes_cb = None
        self._build_instance_controls(root_params)

        layout.addStretch()
        self.setLayout(layout)

        # Wire up preset changes LAST so initial_params can be applied first
        # (we do not auto-reset when the user changes the preset manually —
        #  we only reset if the preset is changed *and* the user hasn't
        #  customised the grid yet; simplest UX: always reset on preset change)
        self.feature_combo.currentTextChanged.connect(self._on_preset_changed)
        self.tune_group.toggled.connect(self._on_tune_toggled)

        # Restore saved grid configuration
        saved_sigmas = ip.get(f"apoc_{cell_type}_grid_sigmas")
        if saved_sigmas:
            self.sigma_input.setText(saved_sigmas)
            self._on_update_grid() # This parses sigmas, builds the grid, etc.
        else:
            self._build_grid()

        # Restore saved checkbox state from config (if available)
        saved_checked = ip.get(f"apoc_{cell_type}_checked_features")
        if saved_checked:
            # saved_checked is a list of [feat_key, sigma_str] pairs
            saved_set = {tuple(pair) for pair in saved_checked}
            self._apply_checked_set(saved_set)
        else:
            # Apply preset defaults
            self._apply_preset_defaults(saved_preset)

        self._update_preview()

    # ────────────────────────────────────────────────────────────────────────
    # Internal helpers
    # ────────────────────────────────────────────────────────────────────────

    def _apply_checked_set(self, checked_set):
        """Check/uncheck grid cells according to the given set of (feat_key, sigma_str)."""
        for key, cb in self._feat_sigma_checks.items():
            cb.blockSignals(True)
            cb.setChecked(key in checked_set)
            cb.blockSignals(False)

    def _apply_preset_defaults(self, preset_name):
        """Reset the grid to the preset's default checked cells."""
        self._apply_checked_set(_checked_set_for_preset(preset_name))

    def _get_current_checked_set(self):
        """Return the set of (feat_key, sigma_str) currently checked in the grid."""
        return {key for key, cb in self._feat_sigma_checks.items() if cb.isChecked()}

    def _build_instance_controls(self, initial_params):
        """Create per-tab instance-segmentation preview controls."""
        if self.cell_type == "dead" or self._apoc_strategy == "APOC (Direct Instance Segmentation)":
            return

        group = QGroupBox("Instance Segmentation Preview")
        group_layout = QVBoxLayout()
        group_layout.setContentsMargins(4, 4, 4, 4)
        group_layout.setSpacing(4)
        group_layout.addWidget(QLabel(f"Uses notebook strategy: {self._apoc_strategy}"))

        if self._apoc_strategy == "APOC Probability Map + Watershed":
            row1 = QHBoxLayout()
            row2 = QHBoxLayout()
            self.prob_mask_threshold_spin = QDoubleSpinBox()
            self.prob_mask_threshold_spin.setRange(0.0, 1.0)
            self.prob_mask_threshold_spin.setSingleStep(0.05)
            self.prob_mask_threshold_spin.setValue(float(initial_params.get(f"{self.cell_type}_prob_mask_threshold", 0.5)))
            row1.addWidget(QLabel("Mask threshold:"))
            row1.addWidget(self.prob_mask_threshold_spin)

            self.prob_seed_threshold_spin = QDoubleSpinBox()
            self.prob_seed_threshold_spin.setRange(0.0, 1.0)
            self.prob_seed_threshold_spin.setSingleStep(0.05)
            self.prob_seed_threshold_spin.setValue(float(initial_params.get(f"{self.cell_type}_prob_seed_threshold", 0.8)))
            row1.addWidget(QLabel("Seed threshold:"))
            row1.addWidget(self.prob_seed_threshold_spin)
            group_layout.addLayout(row1)

            self.opening_nr_pixels_spin = QSpinBox()
            self.opening_nr_pixels_spin.setRange(0, 10)
            self.opening_nr_pixels_spin.setValue(int(initial_params.get(f"{self.cell_type}_opening_nr_pixels", 0)))
            row2.addWidget(QLabel("Opening px:"))
            row2.addWidget(self.opening_nr_pixels_spin)

            self.segment_size_min_spin = QSpinBox()
            self.segment_size_min_spin.setRange(0, 100000)
            self.segment_size_min_spin.setSingleStep(10)
            self.segment_size_min_spin.setValue(int(initial_params.get(f"{self.cell_type}_segment_size_min", 10)))
            row2.addWidget(QLabel("Min size:"))
            row2.addWidget(self.segment_size_min_spin)
            group_layout.addLayout(row2)

        elif self._apoc_strategy == "APOC Mask + EDT/Watershed Resegmentation":
            row1 = QHBoxLayout()
            row2 = QHBoxLayout()
            self.edt_threshold_spin = QDoubleSpinBox()
            self.edt_threshold_spin.setRange(0.0, 50.0)
            self.edt_threshold_spin.setSingleStep(0.5)
            self.edt_threshold_spin.setValue(float(initial_params.get(f"{self.cell_type}_edt_threshold", 1.0)))
            row1.addWidget(QLabel("EDT threshold:"))
            row1.addWidget(self.edt_threshold_spin)

            self.opening_nr_pixels_spin = QSpinBox()
            self.opening_nr_pixels_spin.setRange(0, 10)
            self.opening_nr_pixels_spin.setValue(int(initial_params.get(f"{self.cell_type}_opening_nr_pixels", 0)))
            row1.addWidget(QLabel("Opening px:"))
            row1.addWidget(self.opening_nr_pixels_spin)
            group_layout.addLayout(row1)

            self.segment_size_min_spin = QSpinBox()
            self.segment_size_min_spin.setRange(0, 100000)
            self.segment_size_min_spin.setSingleStep(10)
            self.segment_size_min_spin.setValue(int(initial_params.get(f"{self.cell_type}_segment_size_min", 10)))
            row2.addWidget(QLabel("Min size:"))
            row2.addWidget(self.segment_size_min_spin)

            self.fill_holes_cb = QCheckBox("Fill holes")
            self.fill_holes_cb.setChecked(bool(initial_params.get(f"{self.cell_type}_fill_holes", True)))
            row2.addWidget(self.fill_holes_cb)
            group_layout.addLayout(row2)

        self.run_instance_btn = QtPushButton("Run instance segmentation")
        self.run_instance_btn.clicked.connect(self._on_run_instance_segmentation)
        group_layout.addWidget(self.run_instance_btn)
        group.setLayout(group_layout)
        self.instance_group = group

        for widget in [
            self.prob_mask_threshold_spin,
            self.prob_seed_threshold_spin,
            self.edt_threshold_spin,
            self.segment_size_min_spin,
            self.opening_nr_pixels_spin,
        ]:
            if widget is not None:
                widget.valueChanged.connect(self._emit_params_changed)
        if self.fill_holes_cb is not None:
            self.fill_holes_cb.stateChanged.connect(self._emit_params_changed)

    def _emit_params_changed(self, *_args):
        """Notify the notebook widget that Napari-side params changed."""
        if callable(self._on_params_changed):
            self._on_params_changed()

    def _on_run_instance_segmentation(self):
        """Run instance-segmentation preview for this tab."""
        if callable(self._run_instance_callback):
            self._run_instance_callback(self.cell_type)

    def _on_preset_changed(self, preset_name):
        """Reset grid checkboxes to the new preset's defaults."""
        if preset_name != "custom":
            self.sigma_input.setText(", ".join(_fmt_sigma(s) for s in APOC_SIGMAS))
            self.current_sigmas = list(APOC_SIGMAS)
            self._build_grid()
            self._apply_preset_defaults(preset_name)
        self._update_preview()

    def _on_tune_toggled(self, checked):
        """Show/hide the inner grid widget when the group box is toggled."""
        if self._grid_widget:
            self._grid_widget.setVisible(checked)

    def _on_update_grid(self):
        """Parse custom sigmas from input and rebuild the feature grid."""
        text = self.sigma_input.text()
        try:
            parts = text.replace(",", " ").split()
            new_sigmas = [float(p) for p in parts if p.strip()]
            if not new_sigmas:
                raise ValueError("Empty list")
            self.current_sigmas = new_sigmas
            self.feature_combo.setCurrentText("custom")
            self._build_grid()
        except ValueError:
            QMessageBox.warning(self, "Invalid Sigmas", "Could not parse sigma values. Please enter numbers separated by spaces or commas.")

    def _build_grid(self):
        """Rebuild the checkbox grid based on current_sigmas."""
        old_checked = set()
        if self._feat_sigma_checks:
            old_checked = self._get_current_checked_set()

        if self._grid_widget is not None:
            self.tune_layout.removeWidget(self._grid_widget)
            self._grid_widget.deleteLater()

        self._grid_widget = QWidget()
        grid = QGridLayout()
        grid.setSpacing(2)
        grid.setContentsMargins(0, 0, 0, 0)
        self._feat_sigma_checks = {}

        sigma_header = QLabel("sigma")
        sigma_header.setStyleSheet("font-weight: bold; font-size: 11px;")
        grid.addWidget(sigma_header, 0, 0)

        for col_idx, s in enumerate(self.current_sigmas):
            s_str = _fmt_sigma(s)
            lbl = QLabel(s_str)
            lbl.setStyleSheet("font-size: 10px;")
            lbl.setAlignment(Qt.AlignCenter)
            grid.addWidget(lbl, 0, col_idx + 1)

        for row_idx, (feat_key, feat_label) in enumerate(APOC_ALL_FEATURES):
            feat_lbl = QLabel(feat_label)
            feat_lbl.setStyleSheet("font-size: 11px;")
            grid.addWidget(feat_lbl, row_idx + 1, 0)
            for col_idx, s in enumerate(self.current_sigmas):
                s_str = _fmt_sigma(s)
                cb = QCheckBox()
                cb.setChecked(False)
                cb.setStyleSheet("QCheckBox { margin: 0px; padding: 0px; }")
                cb.stateChanged.connect(self._update_preview)
                grid.addWidget(cb, row_idx + 1, col_idx + 1, alignment=Qt.AlignCenter)
                self._feat_sigma_checks[(feat_key, s_str)] = cb

        self._grid_widget.setLayout(grid)
        self.tune_layout.insertWidget(1, self._grid_widget)

        # Retain checked state for checkboxes that still exist in the new grid
        self._apply_checked_set(old_checked)
        self._update_preview()

    def _update_preview(self):
        """Rebuild the feature string preview from the current grid state."""
        feat_str = self.get_feature_string()
        # Show a truncated preview if very long
        if len(feat_str) > 120:
            display = feat_str[:117] + "…"
        else:
            display = feat_str
        self.preview_label.setText(f"<b>Features:</b> {display}")
        self.preview_label.setToolTip(feat_str)

    def _get_clf_path(self):
        """Return the expected .cl file path for this cell type."""
        if not hasattr(self, '_pixel_class_outdir') or not self._pixel_class_outdir:
            return None
        ct = self.cell_type
        if ct == "dead":
            fname = "PixelClassifier_Death.cl"
        else:
            fname = f"PixelClassifier_{ct.capitalize()}.cl"
        return Path(self._pixel_class_outdir) / fname

    def _on_show_statistics(self):
        """Show a window with the classifier feature statistics."""
        clf_path = self._get_clf_path()
        if clf_path is None or not Path(clf_path).exists():
            QMessageBox.information(
                self,
                "Classifier Statistics",
                f"No trained classifier found for '{self.cell_type}'.\n"
                "Train the classifier first, then click this button."
            )
            return
        try:
            import apoc
            clf = apoc.ObjectSegmenter(opencl_filename=str(clf_path))
            feature_importances = clf.feature_importances()  # dict feature -> importance

            dlg = QDialog(self)
            dlg.setWindowTitle(f"Classifier Statistics — {self.cell_type.capitalize()}")
            dlg_layout = QVBoxLayout(dlg)

            table = QTableWidget(len(feature_importances), 2)
            table.setHorizontalHeaderLabels(["Feature", "Importance"])
            table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
            for i, (feat, imp) in enumerate(sorted(
                feature_importances.items(), key=lambda x: -x[1]
            )):
                table.setItem(i, 0, QTableWidgetItem(str(feat)))
                item = QTableWidgetItem(f"{imp:.4f}")
                item.setTextAlignment(Qt.AlignCenter)
                # Colour code: green for high importance
                r = int(max(0, 255 - imp * 1200))
                g = int(min(255, 80 + imp * 1200))
                item.setBackground(QColor(r, g, 80))
                table.setItem(i, 1, item)

            dlg_layout.addWidget(table)
            dlg.setMinimumSize(500, 400)
            dlg.exec_()
        except Exception as exc:
            QMessageBox.critical(
                self,
                "Classifier Statistics Error",
                f"Could not load statistics for '{self.cell_type}':\n{exc}"
            )

    # ────────────────────────────────────────────────────────────────────────
    # Public API
    # ────────────────────────────────────────────────────────────────────────

    def get_feature_string(self):
        """Return the full APOC feature string based on current grid state."""
        checked = self._get_current_checked_set()
        consider_orig = self.consider_original_cb.isChecked()
        return _build_feature_string_from_checked(checked, consider_original=consider_orig, current_sigmas=self.current_sigmas)

    def refresh_channel_checkboxes(self):
        """Rebuild channel checkboxes from current Napari image layers."""
        existing_names = {cb.text() for cb in self.channel_checkboxes}
        checked_names = {cb.text() for cb in self.channel_checkboxes if cb.isChecked()}
        default_names = set(self._default_channel_names)
        use_defaults = not self.channel_checkboxes

        for cb in self.channel_checkboxes:
            self.chan_checkbox_layout.removeWidget(cb)
            cb.deleteLater()
        self.channel_checkboxes = []

        for layer in self.viewer.layers:
            if (
                isinstance(layer, napari.layers.Image)
                and not layer.name.startswith("Pixel Classification")
                and not layer.name.startswith("Probability Map")
                and not layer.name.startswith("Instance Segmentation")
            ):
                cb = QCheckBox(layer.name)
                if use_defaults and default_names:
                    checked = layer.name in default_names
                elif layer.name in existing_names:
                    checked = layer.name in checked_names
                elif default_names:
                    checked = layer.name in default_names
                else:
                    checked = True
                cb.setChecked(checked)
                self.chan_checkbox_layout.addWidget(cb)
                self.channel_checkboxes.append(cb)
        self._default_channel_names = [
            cb.text() for cb in self.channel_checkboxes if cb.isChecked()
        ]

    def get_config(self):
        """Return a dict with all current widget values."""
        checked_set = self._get_current_checked_set()
        return {
            "feature_preset":     self.feature_combo.currentText(),
            "grid_sigmas":        self.sigma_input.text().strip(),
            # Legacy keys kept for backward compat:
            "sigmas":             ",".join(
                sorted({s for _, s in checked_set},
                       key=lambda x: APOC_SIGMAS.index(float(x)) if float(x) in APOC_SIGMAS else 999)
            ),
            "custom_feature_string": "",
            "feature_string":     self.get_feature_string(),
            "consider_original":  self.consider_original_cb.isChecked(),
            # Rich state for round-tripping
            "checked_features":   [list(pair) for pair in checked_set],
            "max_depth":          self.max_depth_spin.value(),
            "num_ensembles":      self.num_ensembles_spin.value(),
            "channels":           [
                cb.text() for cb in self.channel_checkboxes if cb.isChecked()
            ],
            "prob_mask_threshold": (
                float(self.prob_mask_threshold_spin.value()) if self.prob_mask_threshold_spin is not None else None
            ),
            "prob_seed_threshold": (
                float(self.prob_seed_threshold_spin.value()) if self.prob_seed_threshold_spin is not None else None
            ),
            "edt_threshold": (
                float(self.edt_threshold_spin.value()) if self.edt_threshold_spin is not None else None
            ),
            "segment_size_min": (
                int(self.segment_size_min_spin.value()) if self.segment_size_min_spin is not None else None
            ),
            "opening_nr_pixels": (
                int(self.opening_nr_pixels_spin.value()) if self.opening_nr_pixels_spin is not None else None
            ),
            "fill_holes": (
                bool(self.fill_holes_cb.isChecked()) if self.fill_holes_cb is not None else None
            ),
        }

    def apply_config(self, cfg):
        """Restore widget values from a config dict (used by 'Apply to all tabs')."""
        if "feature_preset" in cfg:
            # Block the auto-reset signal while we apply config
            self.feature_combo.blockSignals(True)
            self.feature_combo.setCurrentText(cfg["feature_preset"])
            self.feature_combo.blockSignals(False)
        if "grid_sigmas" in cfg:
            self.sigma_input.setText(cfg["grid_sigmas"])
            self._on_update_grid()
        if "checked_features" in cfg and cfg["checked_features"]:
            saved_set = {tuple(pair) for pair in cfg["checked_features"]}
            self._apply_checked_set(saved_set)
        elif "feature_preset" in cfg:
            self._apply_preset_defaults(cfg["feature_preset"])
        if "consider_original" in cfg:
            self.consider_original_cb.setChecked(bool(cfg["consider_original"]))
        if "max_depth" in cfg:
            self.max_depth_spin.setValue(int(cfg["max_depth"]))
        if "num_ensembles" in cfg:
            self.num_ensembles_spin.setValue(int(cfg["num_ensembles"]))
        if "channels" in cfg:
            for cb in self.channel_checkboxes:
                cb.setChecked(cb.text() in cfg["channels"])
            self._default_channel_names = list(cfg["channels"])
        if "prob_mask_threshold" in cfg and self.prob_mask_threshold_spin is not None and cfg["prob_mask_threshold"] is not None:
            self.prob_mask_threshold_spin.setValue(float(cfg["prob_mask_threshold"]))
        if "prob_seed_threshold" in cfg and self.prob_seed_threshold_spin is not None and cfg["prob_seed_threshold"] is not None:
            self.prob_seed_threshold_spin.setValue(float(cfg["prob_seed_threshold"]))
        if "edt_threshold" in cfg and self.edt_threshold_spin is not None and cfg["edt_threshold"] is not None:
            self.edt_threshold_spin.setValue(float(cfg["edt_threshold"]))
        if "segment_size_min" in cfg and self.segment_size_min_spin is not None and cfg["segment_size_min"] is not None:
            self.segment_size_min_spin.setValue(int(cfg["segment_size_min"]))
        if "opening_nr_pixels" in cfg and self.opening_nr_pixels_spin is not None and cfg["opening_nr_pixels"] is not None:
            self.opening_nr_pixels_spin.setValue(int(cfg["opening_nr_pixels"]))
        if "fill_holes" in cfg and self.fill_holes_cb is not None and cfg["fill_holes"] is not None:
            self.fill_holes_cb.setChecked(bool(cfg["fill_holes"]))
        self._update_preview()


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
        self._apoc_strategy = str(self._initial_params.get("apoc_strategy", "APOC (Direct Instance Segmentation)"))
        self._log_fn = print

        # All tab labels = cell types + optional Death
        self._tab_cell_types = list(all_cell_types)
        if has_death:
            self._tab_cell_types.append("dead")

        self._build_ui()
        self._connect_signals()
        self._refresh_instance_controls()
        self._refresh_all_channels()

        # Listen for layer changes
        self.viewer.layers.events.inserted.connect(lambda _: self._refresh_all_channels())
        self.viewer.layers.events.removed.connect(lambda _: self._refresh_all_channels())

    def set_log_fn(self, log_fn):
        """Attach an optional GUI log sink."""
        self._log_fn = log_fn or print

    def _log(self, message):
        """Write a message to stdout and, when available, to the GUI log."""
        text = str(message)
        print(text)
        if callable(self._log_fn) and self._log_fn is not print:
            self._log_fn(text)

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
            tab = CellTypeTab(
                ct,
                self.viewer,
                pixel_class_outdir=self.pixel_class_outdir,
                initial_params=self._initial_params,
                apoc_strategy=self._apoc_strategy,
                on_params_changed=self._persist_params,
                run_instance_callback=self._run_instance_preview,
            )
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

        self.instance_controls_widget = QWidget()
        self.instance_controls_layout = QVBoxLayout()
        self.instance_controls_layout.setContentsMargins(0, 0, 0, 0)
        self.instance_controls_layout.setSpacing(6)
        self.instance_controls_widget.setLayout(self.instance_controls_layout)
        self.instance_controls_widget.setVisible(False)
        layout.addWidget(self.instance_controls_widget)

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
        self.tab_widget.currentChanged.connect(lambda _idx: self._refresh_instance_controls())

    def _refresh_instance_controls(self):
        """Show the current tab's instance-segmentation controls in the main dock."""
        while self.instance_controls_layout.count():
            item = self.instance_controls_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)

        current_ct = self._tab_cell_types[self.tab_widget.currentIndex()]
        current_tab = self.tabs[current_ct]
        group = getattr(current_tab, "instance_group", None)

        if group is None:
            self.instance_controls_widget.setVisible(False)
            return

        self.instance_controls_layout.addWidget(group)
        self.instance_controls_widget.setVisible(True)

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
        self._persist_params()
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

    def _get_selected_channel_names(self, ct):
        """Return selected channel layer names in tab order."""
        tab = self.tabs[ct]
        return [cb.text() for cb in tab.channel_checkboxes if cb.isChecked()]

    def _normalize_feature_string(self, feature_string):
        feature_string = str(feature_string or "").replace(",", " ").replace("\t", " ").strip()
        while "  " in feature_string:
            feature_string = feature_string.replace("  ", " ")
        return feature_string

    def _feature_tokens(self, feature_string):
        normalized = self._normalize_feature_string(feature_string)
        return [tok for tok in normalized.split(" ") if tok]

    def _format_feature_token_name(self, token):
        """Convert an APOC feature token into a readable stable name."""
        token = str(token or "").strip()
        if not token:
            return ""
        if token.lower() == "original":
            return "original"
        if "=" not in token:
            return token.replace(" ", "_")
        feature_name, sigma_text = token.split("=", 1)
        sigma_text = sigma_text.strip()
        return f"{feature_name.strip()}_sigma{sigma_text}"

    def _expanded_feature_names(self, ct, feature_string):
        """Return one readable feature name per cached feature plane."""
        channel_names = self._get_selected_channel_names(ct)
        feature_tokens = self._feature_tokens(feature_string)
        expanded = []
        for channel_idx, _channel_name in enumerate(channel_names):
            for token in feature_tokens:
                feature_name = self._format_feature_token_name(token)
                if feature_name:
                    expanded.append(f"{feature_name}_channel{channel_idx}")
        return expanded

    def _generate_feature_list_for_timepoint(self, images, feature_string, t_idx=0, *, push_to_device=False):
        """Generate APOC features for one timepoint without caching them on disk."""
        import apoc

        feature_tokens = self._feature_tokens(feature_string)
        if not feature_tokens:
            raise RuntimeError("No APOC features selected. Choose at least one feature before training.")

        features = []
        for img in images:
            img_t = np.asarray(img[t_idx]) if getattr(img, "ndim", 0) == 4 else np.asarray(img)
            generated = apoc.generate_feature_stack(img_t, feature_string)
            if push_to_device:
                features.extend(generated)
            else:
                features.extend(np.asarray(feat) for feat in generated)
        return features

    def _set_labels_layer(self, name, data, visible=True, opacity=0.8):
        """Create or update a Napari labels layer."""
        if name in [layer.name for layer in self.viewer.layers]:
            self.viewer.layers[name].data = data
            self.viewer.layers[name].visible = visible
        else:
            self.viewer.add_labels(data, name=name, opacity=opacity, visible=visible)

    def _set_image_layer(self, name, data, visible=False, opacity=0.6):
        """Create or update a Napari image layer."""
        if name in [layer.name for layer in self.viewer.layers]:
            self.viewer.layers[name].data = data
            self.viewer.layers[name].visible = visible
        else:
            self.viewer.add_image(
                data,
                name=name,
                opacity=opacity,
                blending="additive",
                colormap="magma",
                contrast_limits=(0.0, 1.0),
                visible=visible,
            )

    def _save_preview_array(self, path, data):
        """Persist a preview array to zarr, replacing the old preview if needed."""
        path = Path(path)
        if path.exists():
            shutil.rmtree(path)
        save_as_zarr(data, path)

    def _predict_classifier_outputs(self, ct, clf=None):
        """Predict raw labels and probability maps for the selected tab images."""
        import apoc
        from behav3d.preprocessing.segmentation.apoc_segment import _get_probability_classifier

        images = self._get_images_for_tab(ct)
        if not images:
            raise RuntimeError(f"No image layers selected for '{ct}'.")
        feature_string = self.tabs[ct].get_config()["feature_string"]

        if clf is None:
            clf_path = self.tabs[ct]._get_clf_path()
            if clf_path is None or not Path(clf_path).exists():
                raise FileNotFoundError(f"No trained classifier found for '{ct}'.")
            clf = apoc.ObjectSegmenter(opencl_filename=str(clf_path))

        prob_clf = _get_probability_classifier(clf)
        if getattr(images[0], "ndim", 0) == 4:
            results = []
            prob_results = []
            for t in range(images[0].shape[0]):
                feats_t = self._generate_feature_list_for_timepoint(
                    images,
                    feature_string,
                    t,
                    push_to_device=True,
                )
                results.append(np.asarray(clf.predict(features=feats_t)).astype(np.int16))
                prob_results.append(np.asarray(prob_clf.predict(features=feats_t)).astype(np.float32))
            return np.stack(results, axis=0), np.stack(prob_results, axis=0)

        features = self._generate_feature_list_for_timepoint(
            images,
            feature_string,
            0,
            push_to_device=True,
        )
        return (
            np.asarray(clf.predict(features=features)).astype(np.int16),
            np.asarray(prob_clf.predict(features=features)).astype(np.float32),
        )

    def _update_prediction_layers(self, ct, segments_result, prob_result):
        """Push classifier preview outputs into Napari layers."""
        seg_layer_name = "Pixel Classification (Dead)" if ct == "dead" else f"{ct.capitalize()} Segments"
        self._set_labels_layer(seg_layer_name, segments_result, visible=True, opacity=0.8)

        prob_layer_name = "Probability Map (Dead)" if ct == "dead" else f"Probability Map ({ct.capitalize()})"
        self._set_image_layer(prob_layer_name, prob_result, visible=False, opacity=0.6)
        _reorder_apoc_layers(self.viewer, self.all_cell_types, has_death=self.has_death)

    def _build_display_segments(self, ct, raw_prediction, prob_prediction):
        """Convert raw classifier outputs into the segment layer shown in Napari."""
        if ct == "dead" or self._apoc_strategy == "APOC (Direct Instance Segmentation)":
            return np.asarray(raw_prediction).astype(np.int16)

        tab = self.tabs[ct]
        if self._apoc_strategy == "APOC Probability Map + Watershed":
            return _probability_array_to_segments(
                prob_prediction,
                mask_thr=float(tab.prob_mask_threshold_spin.value()),
                seed_thr=float(tab.prob_seed_threshold_spin.value()),
                opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
                segment_size_min=int(tab.segment_size_min_spin.value()),
            )

        return _mask_array_to_segments(
            raw_prediction > 0,
            edt_thr=float(tab.edt_threshold_spin.value()),
            opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
            segment_size_min=int(tab.segment_size_min_spin.value()),
            fill_holes=bool(tab.fill_holes_cb.isChecked()),
        )

    def _run_instance_preview(self, ct):
        """Run instance segmentation preview for a non-dead APOC tab."""
        if ct == "dead":
            return

        tab = self.tabs[ct]
        strategy = self._apoc_strategy
        if strategy == "APOC (Direct Instance Segmentation)":
            self.status_label.setText("Direct APOC does not use instance resegmentation preview.")
            return

        try:
            self.status_label.setText(f"Running instance preview for {ct}...")
            QApplication.processEvents()

            prob_layer_name = f"Probability Map ({ct.capitalize()})"
            prob_path = _probability_map_path(self.pixel_class_outdir, ct)
            prob_prediction = None

            if prob_layer_name in [layer.name for layer in self.viewer.layers]:
                prob_prediction = np.asarray(self.viewer.layers[prob_layer_name].data)

            if strategy == "APOC Probability Map + Watershed":
                if prob_prediction is None or prob_path is None or not Path(prob_path).exists():
                    _raw_prediction, prob_prediction = self._predict_classifier_outputs(ct)
                instance_preview = _probability_array_to_segments(
                    prob_prediction,
                    mask_thr=float(tab.prob_mask_threshold_spin.value()),
                    seed_thr=float(tab.prob_seed_threshold_spin.value()),
                    opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
                    segment_size_min=int(tab.segment_size_min_spin.value()),
                )
            else:
                raw_prediction, prob_prediction = self._predict_classifier_outputs(ct)
                instance_preview = _mask_array_to_segments(
                    raw_prediction > 0,
                    edt_thr=float(tab.edt_threshold_spin.value()),
                    opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
                    segment_size_min=int(tab.segment_size_min_spin.value()),
                    fill_holes=bool(tab.fill_holes_cb.isChecked()),
                )

            self._update_prediction_layers(ct, instance_preview, prob_prediction)

            pred_path = _predicted_labels_path(self.pixel_class_outdir, ct)
            if pred_path is not None:
                self._save_preview_array(pred_path, instance_preview)
            if prob_path is not None and prob_prediction is not None:
                self._save_preview_array(prob_path, prob_prediction)

            self.status_label.setText(f"✅ Instance preview updated for {ct}")
        except Exception as exc:
            self.status_label.setText(f"❌ Instance preview failed for {ct}: {exc}")

    def save_user_labels(self, log=print):
        """Save all user-provided labels for all cell types and Dead (if present)."""
        for cell_type in self.all_cell_types:
            lname = f"User Provided Labels ({cell_type.capitalize()})"
            if lname not in [l.name for l in self.viewer.layers]:
                continue
            label_layer = self.viewer.layers[lname]
            label_data = np.asarray(label_layer.data)
            outpath = Path(self.pixel_class_outdir, f"PixelClassifier_User{cell_type.capitalize()}Labels.zarr")
            if outpath.exists():
                shutil.rmtree(outpath)
            save_as_zarr(label_data, outpath)
            log(f"Saved {cell_type} labels → {outpath}")

        if self.has_death and "User Provided Labels (Dead)" in [l.name for l in self.viewer.layers]:
            dead_layer = self.viewer.layers["User Provided Labels (Dead)"]
            dead_label_data = np.asarray(dead_layer.data)
            dead_outpath = Path(self.pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
            if dead_outpath.exists():
                shutil.rmtree(dead_outpath)
            save_as_zarr(dead_label_data, dead_outpath)
            log(f"Saved Death labels → {dead_outpath}")

        log("✅ All user labels saved!")

    def _run_training(self, cell_types_to_train):
        """Train (and apply) APOC classifiers for the given list of cell types."""
        import apoc
        training_start = time.time()

        # Auto-save labels before training
        self._log("Auto-saving user labels before training...")
        self.save_user_labels(log=self._log)

        successes = []
        try:
            for ct in cell_types_to_train:
                tab = self.tabs[ct]
                cfg = tab.get_config()
                feature_string = cfg["feature_string"]
                max_depth = cfg["max_depth"]
                num_ensembles = cfg["num_ensembles"]
                expanded_feature_names = self._expanded_feature_names(ct, feature_string)
                expanded_feature_spec = " ".join(expanded_feature_names)

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

                # Old APOC incremental-style training path kept for reference:
                # if annotation.ndim == 4:
                #     for t in range(n_timepoints):
                #         # Load only the current timepoint slice into memory
                #         ann_t = np.asarray(annotation[t])
                #         if not np.any(ann_t):
                #             continue
                #         feats_t = self._generate_feature_list_for_timepoint(images, feature_string, t)
                #         clf.feature_specification = expanded_feature_spec
                #         clf.train(feats_t, ann_t, continue_training=has_trained)
                #         clf.feature_specification = expanded_feature_spec
                #         has_trained = True
                # else:
                #     # Single timepoint
                #     ann_np = np.asarray(annotation)
                #     feats_np = self._generate_feature_list_for_timepoint(images, feature_string, 0)
                #     clf.feature_specification = expanded_feature_spec
                #     clf.train(feats_np, ann_np)
                #     clf.feature_specification = expanded_feature_spec
                #     has_trained = True

                from sklearn.ensemble import RandomForestClassifier

                X_parts = []
                y_parts = []
                gt_ndim = None

                if annotation.ndim == 4:
                    for t in range(n_timepoints):
                        ann_t = np.asarray(annotation[t])
                        if not np.any(ann_t):
                            continue
                        feats_t = self._generate_feature_list_for_timepoint(images, feature_string, t)
                        X_t, y_t = clf._to_np(feats_t, ann_t)
                        if X_t.size == 0 or y_t.size == 0:
                            continue
                        X_parts.append(X_t)
                        y_parts.append(y_t)
                        gt_ndim = ann_t.ndim
                else:
                    ann_np = np.asarray(annotation)
                    feats_np = self._generate_feature_list_for_timepoint(images, feature_string, 0)
                    X_t, y_t = clf._to_np(feats_np, ann_np)
                    if X_t.size > 0 and y_t.size > 0:
                        X_parts.append(X_t)
                        y_parts.append(y_t)
                        gt_ndim = ann_np.ndim

                if X_parts:
                    X = np.concatenate(X_parts, axis=0)
                    y = np.concatenate(y_parts, axis=0)
                    fitted_rf = RandomForestClassifier(
                        max_depth=max_depth,
                        n_estimators=num_ensembles,
                        random_state=0,
                    )
                    fitted_rf.fit(X, y)
                    clf.classifier = fitted_rf
                    clf._feature_importances = fitted_rf.feature_importances_
                    clf._X = X
                    clf._y = y
                    clf.num_features = X.shape[1]
                    clf.num_ground_truth_dimensions = gt_ndim
                    clf.feature_specification = expanded_feature_spec
                    has_trained = True

                if not has_trained:
                    continue

                clf.feature_specification = expanded_feature_spec
                clf.to_opencl_file(clf_path)

                # Predict (visual confirmation)
                raw_prediction, prob_result = self._predict_classifier_outputs(ct, clf=clf)
                display_segments = self._build_display_segments(ct, raw_prediction, prob_result)
                self._update_prediction_layers(ct, display_segments, prob_result)

                pred_path = _predicted_labels_path(self.pixel_class_outdir, ct)
                if pred_path is not None:
                    self._save_preview_array(pred_path, display_segments)
                prob_path = _probability_map_path(self.pixel_class_outdir, ct)
                if prob_path is not None:
                    self._save_preview_array(prob_path, prob_result)

                successes.append(ct)

            if successes:
                elapsed_s = time.time() - training_start
                elapsed_txt = f"{elapsed_s:.1f}s"
                self.status_label.setText(f"✅ Trained: {', '.join(successes)} ({elapsed_txt})")
                self._log(f"✅ Training finished in {elapsed_txt}")
                # Auto-show statistics for each successfully trained classifier
                for ct in successes:
                    if ct in self.tabs:
                        self.tabs[ct]._on_show_statistics()
            else:
                elapsed_s = time.time() - training_start
                elapsed_txt = f"{elapsed_s:.1f}s"
                self.status_label.setText(f"⚠️ No cell types were trained (check labels). ({elapsed_txt})")
                self._log(f"⚠️ No cell types were trained. Elapsed time: {elapsed_txt}")

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
        params["apoc_strategy"] = self._apoc_strategy
        for ct, tab in self.tabs.items():
            cfg = tab.get_config()
            params[f"apoc_{ct}_feature_preset"]        = cfg["feature_preset"]
            params[f"apoc_{ct}_grid_sigmas"]           = cfg["grid_sigmas"]
            params[f"apoc_{ct}_sigmas"]                = cfg["sigmas"]
            params[f"apoc_{ct}_custom_feature_string"] = cfg["custom_feature_string"]
            params[f"apoc_{ct}_feature_string"]        = cfg["feature_string"]
            params[f"apoc_{ct}_consider_original"]     = cfg["consider_original"]
            params[f"apoc_{ct}_checked_features"]      = cfg["checked_features"]
            params[f"apoc_{ct}_max_depth"]             = cfg["max_depth"]
            params[f"apoc_{ct}_num_ensembles"]         = cfg["num_ensembles"]
            params[f"apoc_{ct}_channels"]              = cfg["channels"]
            if cfg["prob_mask_threshold"] is not None:
                params[f"{ct}_prob_mask_threshold"] = cfg["prob_mask_threshold"]
            if cfg["prob_seed_threshold"] is not None:
                params[f"{ct}_prob_seed_threshold"] = cfg["prob_seed_threshold"]
            if cfg["edt_threshold"] is not None:
                params[f"{ct}_edt_threshold"] = cfg["edt_threshold"]
            if cfg["segment_size_min"] is not None:
                params[f"{ct}_segment_size_min"] = cfg["segment_size_min"]
            if cfg["opening_nr_pixels"] is not None:
                params[f"{ct}_opening_nr_pixels"] = cfg["opening_nr_pixels"]
            if cfg["fill_holes"] is not None:
                params[f"{ct}_fill_holes"] = cfg["fill_holes"]
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

    # image_list is a list of timepoints: [t1(C,Z,Y,X), t2(C,Z,Y,X), ...]
    # Stacking along axis 0 gives (T_total, C, Z, Y, X)
    stacked = np.stack(image_list, axis=0)

    # Since stacked was created by stacking (C,Z,Y,X) frames along axis 0,
    # it is ALWAYS (T_total, C, Z, Y, X).
    c_axis_in_stacked = 1

    # --- Create Napari viewer ---
    viewer = napari.Viewer()

    # stacked is (T_total, C, Z, Y, X).
    # Axis 0 = timepoints, axis 1 = channels.
    T_total = stacked.shape[0]
    n_channels = stacked.shape[1]   # ← was incorrectly stacked.shape[0]
    channel_colors = [
        "cyan", "yellow", "red", "green", "magenta", "blue",
        "gray", "turbo", "viridis", "plasma", "inferno", "twilight",
    ]

    for ch in range(n_channels):
        channel_data = stacked[:, ch, :, :, :]  # (T_total, Z, Y, X)
        nonzero = channel_data[channel_data > 0]
        clim = (0, float(np.percentile(nonzero, 99.8))) if nonzero.size > 0 else (0, 1e-3)
        img_layer = viewer.add_image(
            channel_data,
            name=f"Channel {ch}",
            contrast_limits=clim,
            colormap=channel_colors[ch % len(channel_colors)],
            blending="additive",
            opacity=0.8,
        )
        img_layer.contrast_limits_range = (0, float(channel_data.max()))

    # --- Annotation label layers per cell type ---
    # Labels must be (T_total, Z, Y, X) — NOT (C, Z, Y, X)
    label_shape = (T_total,) + stacked.shape[2:]  # (T_total, Z, Y, X)

    ip = initial_params or {}

    def _restore_saved_labels(saved_path, expected_shape, label_name):
        if not saved_path.exists():
            return np.zeros(expected_shape, dtype=np.int16)
        try:
            existing = np.asarray(load_zarr(saved_path))
        except Exception as exc:
            print(f"  ⚠️ Could not restore saved labels for '{label_name}' from {saved_path}: {exc}")
            return np.zeros(expected_shape, dtype=np.int16)
        if existing.shape == expected_shape:
            print(f"  ↩ Restored saved labels for '{label_name}' ({expected_shape})")
            return existing
        print(f"  ⚠️ Saved labels shape {existing.shape} ≠ expected {expected_shape} — starting fresh")
        return np.zeros(expected_shape, dtype=np.int16)

    # 1. User Provided Labels (all cell types)
    for cell_type in all_cell_types:
        saved_path = Path(pixel_class_outdir, f"PixelClassifier_User{cell_type.capitalize()}Labels.zarr")
        user_labels = _restore_saved_labels(saved_path, label_shape, cell_type)

        viewer.add_labels(
            user_labels,
            name=f"User Provided Labels ({cell_type.capitalize()})",
            opacity=0.5,
        )

    if has_death:
        dead_path = Path(pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
        dead_labels = _restore_saved_labels(dead_path, label_shape, "dead")
        viewer.add_labels(dead_labels, name="User Provided Labels (Dead)", opacity=0.5)

    # 2. Results (Pixel Classification / Segments) on top
    # Dead result
    if has_death:
        pred_death_path = _predicted_labels_path(pixel_class_outdir, "dead")
        if pred_death_path.exists():
            pred_death = np.asarray(load_zarr(pred_death_path))
            if pred_death.shape != label_shape:
                pred_death = np.zeros(label_shape, dtype=np.int16)
        else:
            pred_death = np.zeros(label_shape, dtype=np.int16)
        viewer.add_labels(pred_death, name="Pixel Classification (Dead)", opacity=0.8, visible=False)

        prob_death_path = _probability_map_path(pixel_class_outdir, "dead")
        if prob_death_path.exists():
            prob_death = np.asarray(load_zarr(prob_death_path))
            if prob_death.shape != label_shape:
                prob_death = np.zeros(label_shape, dtype=np.float32)
        else:
            prob_death = np.zeros(label_shape, dtype=np.float32)
        viewer.add_image(
            prob_death,
            name="Probability Map (Dead)",
            opacity=0.6,
            blending="additive",
            colormap="magma",
            contrast_limits=(0.0, 1.0),
            visible=False,
        )

    # Cell type segments
    for cell_type in all_cell_types:
        pred_path = _predicted_labels_path(pixel_class_outdir, cell_type)
        if pred_path.exists():
            pred_data = np.asarray(load_zarr(pred_path))
            if pred_data.shape != label_shape:
                pred_data = np.zeros(label_shape, dtype=np.int16)
        else:
            pred_data = np.zeros(label_shape, dtype=np.int16)

        viewer.add_labels(
            pred_data,
            name=f"{cell_type.capitalize()} Segments",
            opacity=0.8,
            visible=False,
        )

        prob_path = _probability_map_path(pixel_class_outdir, cell_type)
        if prob_path.exists():
            prob_data = np.asarray(load_zarr(prob_path))
            if prob_data.shape != label_shape:
                prob_data = np.zeros(label_shape, dtype=np.float32)
        else:
            prob_data = np.zeros(label_shape, dtype=np.float32)

        viewer.add_image(
            prob_data,
            name=f"Probability Map ({cell_type.capitalize()})",
            opacity=0.6,
            blending="additive",
            colormap="magma",
            contrast_limits=(0.0, 1.0),
            visible=False,
        )

    _reorder_apoc_layers(viewer, all_cell_types, has_death=has_death)

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
    apoc_widget.set_log_fn(log_output.appendPlainText)

    save_button = QtPushButton("💾 Save User Labels")
    save_button.setStyleSheet("background-color: #FF9800; color: white; font-weight: bold; padding: 6px;")
    save_button.clicked.connect(lambda: apoc_widget.save_user_labels(log=log_output.appendPlainText))
    apoc_widget.layout().addWidget(save_button)

    return viewer
