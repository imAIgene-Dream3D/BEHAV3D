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
from qtpy.QtCore import Qt, Signal, QTimer

from behav3d.napari._background_runner import BackgroundOperation, ThreadSafeLogger
from qtpy.QtGui import QColor

import dask.array as da
import pyclesperanto_prototype as cle

from behav3d.io.images import load_image, load_zarr, save_as_zarr
from behav3d.preprocessing import zeropad_image_to_match_shape
from behav3d.core.qt_help import HelpButton, make_help_row
from behav3d.core.utils import ignore_missing_rmtree_error

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
    "small_preset": {
        "label": "Small Preset",
        "original": True,
        "checked": {
            ("gaussian_blur",               "1"),
            ("gaussian_blur",               "2"),
            ("gaussian_blur",               "5"),
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
    "medium_preset": {
        "label": "Medium Preset",
        "original": True,
        "checked": {
            ("gaussian_blur",               "1"),
            ("gaussian_blur",               "2"),
            ("gaussian_blur",               "5"),
            ("gaussian_blur",               "15"),
            ("difference_of_gaussian",      "1"),
            ("difference_of_gaussian",      "2"),
            ("difference_of_gaussian",      "5"),
            ("difference_of_gaussian",      "15"),
            ("laplace_box_of_gaussian_blur","1"),
            ("laplace_box_of_gaussian_blur","2"),
            ("laplace_box_of_gaussian_blur","5"),
            ("laplace_box_of_gaussian_blur","15"),
            ("sobel_of_gaussian_blur",      "1"),
            ("sobel_of_gaussian_blur",      "2"),
            ("sobel_of_gaussian_blur",      "5"),
            ("sobel_of_gaussian_blur",      "15"),
        },
    },
    "large_preset": {
        "label": "Large Preset",
        "original": True,
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
        "original": False,
        "checked": set(),
    },
}


def _checked_set_for_preset(preset_name):
    """Return the set of (feature_key, sigma_str) pairs that should be checked for a preset."""
    return set(FEATURE_PRESETS.get(preset_name, FEATURE_PRESETS["large_preset"])["checked"])


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
        while line != "*/" and count < 80:
            line = f.readline()
            if not line:
                break
            count += 1
            line = line.strip()
            if line.startswith(prefix):
                return line[len(prefix):].strip()
    return default


def _parse_channel_indices(value):
    """Parse a comma-separated list of channel indices from .cl metadata."""
    if value is None:
        return []
    indices = []
    for token in str(value).split(","):
        token = token.strip()
        if not token:
            continue
        try:
            idx = int(token)
        except (TypeError, ValueError):
            continue
        if idx >= 0 and idx not in indices:
            indices.append(idx)
    return indices


def _parse_channel_names(value):
    """Parse a pipe-separated list of channel names from .cl metadata."""
    if value is None:
        return []
    names = []
    for token in str(value).split("|"):
        name = token.strip()
        if name and name not in names:
            names.append(name)
    return names


def _channel_index_from_name(name):
    """Extract the integer index from a canonical Napari layer name."""
    match = re.fullmatch(r"Channel\s+(\d+)", str(name or "").strip())
    if not match:
        return None
    try:
        return int(match.group(1))
    except (TypeError, ValueError):
        return None


def _channel_name_from_index(index):
    """Format a channel index as the canonical Napari training layer name."""
    return f"Channel {int(index)}"


def _normalize_channel_names(channel_names=None, channel_indices=None):
    """Normalize restored channel metadata to stable Napari layer names."""
    normalized = []
    for name in channel_names or []:
        clean_name = str(name).strip()
        if clean_name and clean_name not in normalized:
            normalized.append(clean_name)

    for index in channel_indices or []:
        try:
            canonical = _channel_name_from_index(index)
        except (TypeError, ValueError):
            continue
        if canonical not in normalized:
            normalized.append(canonical)

    return normalized


def _read_classifier_channel_metadata(opencl_path):
    """Return channel names and indices embedded in a trained classifier."""
    channel_indices = _parse_channel_indices(
        _read_classifier_header_value(opencl_path, "input_channel_indices")
    )
    channel_names = _parse_channel_names(
        _read_classifier_header_value(opencl_path, "input_channel_names")
    )
    normalized_names = _normalize_channel_names(
        channel_names=channel_names,
        channel_indices=channel_indices,
    )
    return {
        "channel_indices": channel_indices,
        "channel_names": normalized_names,
    }


def _upsert_classifier_header_values(opencl_path, values):
    """Insert or replace metadata lines inside the header comment of a .cl file."""
    path = Path(opencl_path)
    if not path.exists() or not values:
        return

    content = path.read_text()
    newline = "\r\n" if "\r\n" in content else "\n"
    lines = content.splitlines(keepends=True)

    header_end = None
    for idx, line in enumerate(lines[:80]):
        if line.strip() == "*/":
            header_end = idx
            break

    if header_end is None:
        header_lines = ["/*" + newline]
        for key, value in values.items():
            header_lines.append(f"{key} = {value}{newline}")
        header_lines.append("*/" + newline)
        path.write_text("".join(header_lines) + content)
        return

    new_lines = []
    replaced = set()
    for idx, line in enumerate(lines):
        if idx < header_end:
            stripped = line.strip()
            replaced_line = False
            for key, value in values.items():
                prefix = f"{key} = "
                if stripped.startswith(prefix):
                    new_lines.append(f"{key} = {value}{newline}")
                    replaced.add(key)
                    replaced_line = True
                    break
            if replaced_line:
                continue

        if idx == header_end:
            for key, value in values.items():
                if key not in replaced:
                    new_lines.append(f"{key} = {value}{newline}")
            new_lines.append(line)
            continue

        new_lines.append(line)

    path.write_text("".join(new_lines))


def _write_classifier_channel_metadata(opencl_path, channel_names):
    """Persist selected input channels into a trained APOC classifier header."""
    normalized_names = _normalize_channel_names(channel_names=channel_names)
    channel_indices = []
    for name in normalized_names:
        index = _channel_index_from_name(name)
        if index is not None and index not in channel_indices:
            channel_indices.append(index)

    values = {
        "input_channel_indices": ",".join(str(idx) for idx in channel_indices),
        "input_channel_names": "|".join(normalized_names),
    }
    _upsert_classifier_header_values(opencl_path, values)


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

    channel_meta = _read_classifier_channel_metadata(clf_path)
    if channel_meta["channel_names"]:
        cfg["channels"] = channel_meta["channel_names"]

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


def _mask_array_to_segments(
    mask,
    edt_thr,
    opening_nr_pixels,
    segment_size_min,
    fill_holes,
    marker_strategy="threshold",
    peak_min_distance=None,
    peak_min_ratio=0.35,
):
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
                marker_strategy=marker_strategy,
                peak_min_distance=peak_min_distance,
                peak_min_ratio=peak_min_ratio,
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
# Training data persistence helpers
# ---------------------------------------------------------------------------
# Training data is stored as a single zip bundle: PixelClassifier_TrainingData.zip
# containing training_metadata.yml and per-cell-type {Stem}_X.npy / {Stem}_y.npy.
# The bundle is written with a read-modify-write pattern so each save call
# (per cell type, then metadata) updates only its own members.
# Old-style TrainingData/ folders are still readable for backward compat.
# ---------------------------------------------------------------------------

def _training_bundle_path(pixel_class_outdir):
    """Return the path to the single-file training bundle (.zip), or None."""
    if not pixel_class_outdir:
        return None
    return Path(pixel_class_outdir) / "PixelClassifier_TrainingData.zip"


def _training_data_dir(pixel_class_outdir):
    """Return the legacy TrainingData/ subdirectory path (backward compat)."""
    if not pixel_class_outdir:
        return None
    return Path(pixel_class_outdir) / "TrainingData"


def _celltype_file_stem(cell_type):
    """Return the file stem used for saving a cell type's data ('Death' / capitalized)."""
    return "Death" if cell_type == "dead" else cell_type.capitalize()


def _zip_read_all_members(zip_path):
    """Return a dict {member_name: bytes} for an existing zip, or {} if absent."""
    import zipfile
    if not zip_path.exists():
        return {}
    with zipfile.ZipFile(zip_path, "r") as zf:
        return {name: zf.read(name) for name in zf.namelist()}


def _zip_write_members(zip_path, members):
    """Write *members* dict {name: bytes} to *zip_path*, replacing the file."""
    import zipfile
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for name, data in members.items():
            zf.writestr(name, data)


def _save_training_data(pixel_class_outdir, cell_type, X, y):
    """Save flattened feature matrix X and label vector y for one cell type."""
    import io
    bundle = _training_bundle_path(pixel_class_outdir)
    if bundle is None:
        return
    bundle.parent.mkdir(parents=True, exist_ok=True)
    stem = _celltype_file_stem(cell_type)

    members = _zip_read_all_members(bundle)

    buf_X = io.BytesIO()
    np.save(buf_X, X.astype(np.float32))
    members[f"{stem}_X.npy"] = buf_X.getvalue()

    buf_y = io.BytesIO()
    np.save(buf_y, y.astype(np.uint8))
    members[f"{stem}_y.npy"] = buf_y.getvalue()

    _zip_write_members(bundle, members)


def _load_training_data_for_celltype(path_or_dir, cell_type):
    """Load (X, y) arrays for one cell type; returns (None, None) if not found.

    Accepts either a .zip bundle path or a legacy TrainingData/ directory.
    """
    import io, zipfile
    stem = _celltype_file_stem(cell_type)
    p = Path(path_or_dir)

    if p.suffix == ".zip":
        if not p.exists():
            return None, None
        with zipfile.ZipFile(p, "r") as zf:
            names = zf.namelist()
            x_name, y_name = f"{stem}_X.npy", f"{stem}_y.npy"
            if x_name not in names or y_name not in names:
                return None, None
            X = np.load(io.BytesIO(zf.read(x_name)))
            y = np.load(io.BytesIO(zf.read(y_name)))
        return X, y

    # Backward compat: old TrainingData/ folder
    x_path = p / f"{stem}_X.npz"
    y_path = p / f"{stem}_y.npy"
    if not x_path.exists() or not y_path.exists():
        return None, None
    X = np.load(str(x_path))["X"]
    y = np.load(str(y_path))
    return X, y


def _save_training_metadata(pixel_class_outdir, metadata_dict):
    """Write training_metadata.yml into the training bundle zip."""
    import yaml
    bundle = _training_bundle_path(pixel_class_outdir)
    if bundle is None:
        return
    bundle.parent.mkdir(parents=True, exist_ok=True)

    members = _zip_read_all_members(bundle)
    members["training_metadata.yml"] = yaml.safe_dump(
        metadata_dict, default_flow_style=False, sort_keys=False
    ).encode("utf-8")
    _zip_write_members(bundle, members)


def _load_training_metadata(path_to_zip_yml_or_dir):
    """Load training metadata; accepts a .zip bundle, a .yml path, or a directory.

    Returns the parsed dict, or None if the file does not exist or fails to parse.
    """
    import yaml, zipfile
    if path_to_zip_yml_or_dir is None:
        return None
    p = Path(path_to_zip_yml_or_dir)

    if p.suffix == ".zip":
        if not p.exists():
            return None
        try:
            with zipfile.ZipFile(p, "r") as zf:
                return yaml.safe_load(zf.read("training_metadata.yml"))
        except Exception:
            return None

    # Backward compat: .yml path or directory
    if p.is_dir():
        p = p / "training_metadata.yml"
    if not p.exists():
        return None
    try:
        with open(p, "r") as fh:
            return yaml.safe_load(fh)
    except Exception:
        return None


def _load_training_bundle(zip_path):
    """Load a training bundle zip; returns (metadata_dict, {cell_type: (X, y)}).

    Raises on unreadable files so callers can show a meaningful error dialog.
    """
    import io, yaml, zipfile
    p = Path(zip_path)
    with zipfile.ZipFile(p, "r") as zf:
        meta = yaml.safe_load(zf.read("training_metadata.yml"))
        cell_types = list(meta.get("cell_types", []))
        if meta.get("has_death"):
            cell_types.append("dead")
        names = zf.namelist()
        result = {}
        for ct in cell_types:
            stem = _celltype_file_stem(ct)
            x_name, y_name = f"{stem}_X.npy", f"{stem}_y.npy"
            if x_name in names and y_name in names:
                X = np.load(io.BytesIO(zf.read(x_name)))
                y = np.load(io.BytesIO(zf.read(y_name)))
                result[ct] = (X, y)
    return meta, result


def _build_training_metadata(
    cell_types,
    has_death,
    params_per_ct,
    n_input_channels,
    channel_names,
    pixel_size_xy_um,
    pixel_size_z_um,
    pixel_counts_by_ct,
    imported_from=None,
    imported_counts=None,
    new_counts=None,
):
    """Build the metadata dict that is written to training_metadata.yml."""
    import datetime

    all_cts = list(cell_types) + (["dead"] if has_death else [])
    features_section = {}
    for ct in all_cts:
        cfg = params_per_ct.get(ct, {})
        raw_sigmas = cfg.get("grid_sigmas", "") or ""
        sigma_list = []
        for tok in raw_sigmas.replace(",", " ").split():
            try:
                sigma_list.append(float(tok))
            except ValueError:
                pass
        features_section[ct] = {
            "feature_string":    cfg.get("feature_string", ""),
            "sigmas":            sigma_list,
            "consider_original": bool(cfg.get("consider_original", False)),
            "channels_used":     list(cfg.get("channels", [])),
            "max_depth":         int(cfg.get("max_depth", 5)),
            "num_ensembles":     int(cfg.get("num_ensembles", 100)),
        }

    training_data_section = {}
    for ct, counts in (pixel_counts_by_ct or {}).items():
        training_data_section[ct] = {
            "n_pixels":   int(counts.get("total", 0)),
            "n_positive": int(counts.get("positive", 0)),
            "n_negative": int(counts.get("negative", 0)),
        }

    return {
        "version":    1,
        "created_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "cell_types": [str(ct) for ct in cell_types],
        "has_death":  bool(has_death),
        "image_metadata": {
            "n_input_channels": int(n_input_channels) if n_input_channels is not None else None,
            "channel_names":    list(channel_names),
            "pixel_size_xy_um": float(pixel_size_xy_um) if pixel_size_xy_um is not None else None,
            "pixel_size_z_um":  float(pixel_size_z_um)  if pixel_size_z_um  is not None else None,
        },
        "features":      features_section,
        "training_data": training_data_section,
        "provenance": {
            "imported_from":         str(imported_from) if imported_from else None,
            "imported_pixel_counts": {
                str(k): int(v) for k, v in (imported_counts or {}).items()
            },
            "new_pixel_counts": {
                str(k): int(v) for k, v in (new_counts or {}).items()
            },
        },
    }


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
        instance_controls_inline=False,
        show_strategy_combo=False,
        per_tab_strategies=None,
        on_per_tab_strategy_changed=None,
    ):
        super().__init__(parent)
        self.cell_type = cell_type
        self.viewer = viewer
        self._pixel_class_outdir = pixel_class_outdir
        self._on_params_changed = on_params_changed
        self._run_instance_callback = run_instance_callback
        self._apoc_strategy = str(apoc_strategy)
        # Inline mode: keep the instance-controls group inside this tab instead
        # of letting the parent training widget reparent it into a shared dock.
        self._instance_controls_inline = bool(instance_controls_inline)
        self._show_strategy_combo = bool(show_strategy_combo) and cell_type != "dead"
        self._per_tab_strategies = list(per_tab_strategies or [])
        self._on_per_tab_strategy_changed = on_per_tab_strategy_changed
        self._per_tab_strategy_combo = None
        self._per_tab_strategy_widget = None
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
        self._default_channel_names = list(cfg.get("channels", []))

        layout = QVBoxLayout()
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        # ── Channel selection ────────────────────────────────────────────────
        self.chan_group = QGroupBox("Image Channel Inputs")
        chan_layout = QVBoxLayout()
        chan_layout.setSpacing(2)
        self.channel_checkboxes = []
        self.chan_checkbox_container = QWidget()
        self.chan_checkbox_layout = QVBoxLayout()
        self.chan_checkbox_layout.setContentsMargins(0, 0, 0, 0)
        self.chan_checkbox_container.setLayout(self.chan_checkbox_layout)
        chan_layout.addWidget(self.chan_checkbox_container)
        self.chan_group.setLayout(chan_layout)
        layout.addWidget(self.chan_group)

        # ── Optional per-tab strategy combo ──────────────────────────────────
        # Only shown in the napari GUI plugin (per_cell_type_strategy=True). The
        # notebook entry point does not set this flag, so the combo never
        # appears there and the notebook behaviour is identical to before.
        if self._show_strategy_combo and self._per_tab_strategies:
            saved_per_tab = root_params.get(
                f"per_ct_strategy_{cell_type}", self._apoc_strategy
            )
            if saved_per_tab not in self._per_tab_strategies:
                saved_per_tab = self._per_tab_strategies[0]
            self._per_tab_strategy_widget = QWidget()
            per_tab_lay = QHBoxLayout(self._per_tab_strategy_widget)
            per_tab_lay.setContentsMargins(0, 2, 0, 4)
            per_tab_lay.addWidget(QLabel("\u21b3 Cell type strategy:"))
            self._per_tab_strategy_combo = QComboBox()
            self._per_tab_strategy_combo.addItems(self._per_tab_strategies)
            self._per_tab_strategy_combo.setCurrentText(saved_per_tab)
            self._per_tab_strategy_combo.currentTextChanged.connect(
                self._on_per_tab_strategy_combo_changed
            )
            per_tab_lay.addWidget(self._per_tab_strategy_combo, stretch=1)
            layout.addWidget(self._per_tab_strategy_widget)

        # Labeling hint
        self._hint_label = QLabel("Labels: <b>1</b> = background&nbsp;&nbsp; <b>2</b> = foreground")
        self._hint_label.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(self._hint_label)

        # ── Preset dropdown ──────────────────────────────────────────────────
        self._feat_row_widget = QWidget()
        feat_row = QHBoxLayout(self._feat_row_widget)
        feat_row.setContentsMargins(0, 0, 0, 0)
        feat_row.addWidget(QLabel("Preset:"))
        self.feature_combo = QComboBox()
        self.feature_combo.addItems(list(FEATURE_PRESETS.keys()))
        saved_preset = ip.get(f"apoc_{cell_type}_feature_preset", "large_preset")
        if saved_preset not in FEATURE_PRESETS:
            saved_preset = "large_preset"
        self.feature_combo.setCurrentText(saved_preset)
        feat_row.addWidget(self.feature_combo)
        feat_row.addWidget(HelpButton(
            "Feature preset",
            "Pre-configured combination of feature filters and sigmas.\n"
            "  • Small Preset: all features × sigmas 1, 2, 5 + original image.\n"
            "  • Medium Preset: all features × sigmas 1, 2, 5, 15 + original image.\n"
            "  • Large Preset: all features × sigmas 1, 2, 5, 10, 25 + original image.\n"
            "  • Custom: manually select features and sigmas.\n"
            "Open 'Tune Features' below to customise the grid."
        ))
        feat_row.addStretch()
        layout.addWidget(self._feat_row_widget)

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
        sigma_row.addWidget(HelpButton(
            "Custom sigmas",
            "Sigma values (in pixels) used to build the multiscale feature grid.\n"
            "Each sigma is combined with each enabled feature filter (Gaussian, DoG, "
            "Laplacian, ...) to form one column of the feature matrix.\n"
            "Common sets: '1, 2, 4' (small structures), '2, 4, 8, 16' (medium-large)."
        ))
        self.tune_layout.addLayout(sigma_row)

        self.current_sigmas = list(APOC_SIGMAS)
        self._feat_sigma_checks = {}
        self._grid_widget = None

        self.tune_group.setLayout(self.tune_layout)
        layout.addWidget(self.tune_group)

        # ── Consider original image checkbox ─────────────────────────────────
        self._orig_row_widget = QWidget()
        orig_row = QHBoxLayout(self._orig_row_widget)
        orig_row.setContentsMargins(0, 0, 0, 0)
        self.consider_original_cb = QCheckBox("Consider original image as well")
        preset_orig_default = FEATURE_PRESETS.get(saved_preset, {}).get("original", False)
        saved_orig = bool(ip.get(f"apoc_{cell_type}_consider_original", preset_orig_default))
        self.consider_original_cb.setChecked(saved_orig)
        self.consider_original_cb.stateChanged.connect(self._update_preview)
        orig_row.addWidget(self.consider_original_cb)
        orig_row.addWidget(HelpButton(
            "Consider original image",
            "When enabled, raw pixel intensity is added as an extra feature column "
            "alongside the multiscale filters.\n"
            "Useful when intensity alone is a strong cue (e.g. very bright nuclei)."
        ))
        orig_row.addStretch()
        layout.addWidget(self._orig_row_widget)

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
        self._rf_row_widget = QWidget()
        rf_row = QHBoxLayout(self._rf_row_widget)
        rf_row.setContentsMargins(0, 0, 0, 0)
        rf_row.addWidget(QLabel("Max depth:"))
        self.max_depth_spin = QSpinBox()
        self.max_depth_spin.setRange(1, 20)
        self.max_depth_spin.setValue(int(ip.get(f"apoc_{cell_type}_max_depth", 5)))
        rf_row.addWidget(self.max_depth_spin)
        rf_row.addWidget(HelpButton(
            "Random Forest — max depth",
            "Maximum depth of each decision tree in the APOC random forest.\n"
            "Shallow trees (2–5) train fast and generalise well; deeper trees "
            "(10+) can overfit on small label sets."
        ))
        rf_row.addWidget(QLabel("Trees:"))
        self.num_ensembles_spin = QSpinBox()
        self.num_ensembles_spin.setRange(10, 1000)
        self.num_ensembles_spin.setSingleStep(10)
        self.num_ensembles_spin.setValue(int(ip.get(f"apoc_{cell_type}_num_ensembles", 100)))
        rf_row.addWidget(self.num_ensembles_spin)
        rf_row.addWidget(HelpButton(
            "Random Forest — number of trees",
            "Number of decision trees in the random forest ensemble.\n"
            "More trees → smoother predictions but slower training and inference.\n"
            "Default (100) is a good starting point."
        ))
        rf_row.addStretch()
        layout.addWidget(self._rf_row_widget)

        self.run_instance_btn = None
        self.instance_group = None
        self.prob_mask_threshold_spin = None
        self.prob_seed_threshold_spin = None
        self.edt_threshold_spin = None
        self.segment_size_min_spin = None
        self.opening_nr_pixels_spin = None
        self.fill_holes_cb = None
        self.peak_min_distance_spin = None
        self.peak_min_ratio_spin = None
        # Cached so future rebuilds (e.g. on strategy change) can pick up the
        # same persisted post-processing parameters.
        self._initial_params_cache = dict(root_params)
        # Placeholder where the instance group is inserted in inline mode.
        # In docked mode this placeholder is unused; the group is reparented
        # into the parent training widget's shared dock area.
        self._instance_group_anchor = QWidget()
        self._instance_group_anchor.setVisible(False)
        self._instance_group_anchor_layout = QVBoxLayout(self._instance_group_anchor)
        self._instance_group_anchor_layout.setContentsMargins(0, 0, 0, 0)
        self._instance_group_anchor_layout.setSpacing(0)
        self._build_instance_controls(root_params)
        if self._instance_controls_inline and self.instance_group is not None:
            self._instance_group_anchor_layout.addWidget(self.instance_group)
            self._instance_group_anchor.setVisible(True)
        layout.addWidget(self._instance_group_anchor)

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

        self._locked_feature_string = None
        self._update_preview()

    def set_classifier_params_visible(self, visible):
        """Toggle the visibility of all classifier-specific parameters."""
        self.chan_group.setVisible(visible)
        self._hint_label.setVisible(visible)
        self._feat_row_widget.setVisible(visible)
        self.tune_group.setVisible(visible)
        self._orig_row_widget.setVisible(visible)
        self.preview_label.setVisible(visible)
        self.stats_btn.setVisible(visible)
        self._rf_row_widget.setVisible(visible)

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

    def _build_instance_controls(self, initial_params, strategy=None):
        """Create per-tab instance-segmentation preview controls.

        When *strategy* is provided it overrides ``self._apoc_strategy`` for
        this build (used by :meth:`rebuild_instance_controls` after a strategy
        change). Spinbox/checkbox attributes are reset to ``None`` upfront so
        callers can always rely on ``getattr(tab, "<attr>_spin", None)``.
        """
        # Reset references to controls that depend on the chosen strategy.
        self.run_instance_btn = None
        self.instance_group = None
        self.prob_mask_threshold_spin = None
        self.prob_seed_threshold_spin = None
        self.edt_threshold_spin = None
        self.segment_size_min_spin = None
        self.opening_nr_pixels_spin = None
        self.fill_holes_cb = None
        self.peak_min_distance_spin = None
        self.peak_min_ratio_spin = None

        effective_strategy = str(strategy or self._apoc_strategy)

        if self.cell_type == "dead" or effective_strategy == "APOC (Direct Instance Segmentation)":
            return

        group = QGroupBox("Instance Segmentation Preview")
        group_layout = QVBoxLayout()
        group_layout.setContentsMargins(4, 4, 4, 4)
        group_layout.setSpacing(4)

        if effective_strategy == "APOC Probability Map + Watershed":
            row1 = QHBoxLayout()
            row2 = QHBoxLayout()
            self.prob_mask_threshold_spin = QDoubleSpinBox()
            self.prob_mask_threshold_spin.setRange(0.0, 1.0)
            self.prob_mask_threshold_spin.setSingleStep(0.05)
            self.prob_mask_threshold_spin.setValue(float(initial_params.get(f"{self.cell_type}_prob_mask_threshold", 0.5)))
            row1.addWidget(QLabel("Mask threshold:"))
            row1.addWidget(self.prob_mask_threshold_spin)
            row1.addWidget(HelpButton(
                "Mask threshold",
                "Foreground cutoff applied to the probability map.\n"
                "Pixels with probability above this value are kept as foreground."
            ))

            self.prob_seed_threshold_spin = QDoubleSpinBox()
            self.prob_seed_threshold_spin.setRange(0.0, 1.0)
            self.prob_seed_threshold_spin.setSingleStep(0.05)
            self.prob_seed_threshold_spin.setValue(float(initial_params.get(f"{self.cell_type}_prob_seed_threshold", 0.8)))
            row1.addWidget(QLabel("Seed threshold:"))
            row1.addWidget(self.prob_seed_threshold_spin)
            row1.addWidget(HelpButton(
                "Seed threshold",
                "Higher cutoff used to place watershed seeds.\n"
                "Should be ≥ Mask threshold (typical: 0.8).\n"
                "Lower values produce more seeds (split more touching objects)."
            ))
            group_layout.addLayout(row1)

            self.opening_nr_pixels_spin = QSpinBox()
            self.opening_nr_pixels_spin.setRange(0, 10)
            self.opening_nr_pixels_spin.setValue(int(initial_params.get(f"{self.cell_type}_opening_nr_pixels", 0)))
            row2.addWidget(QLabel("Opening px:"))
            row2.addWidget(self.opening_nr_pixels_spin)
            row2.addWidget(HelpButton(
                "Morphological opening",
                "Number of erosion-then-dilation iterations applied to the mask.\n"
                "Smooths boundaries and removes small speckles. Set to 0 to disable."
            ))

            self.segment_size_min_spin = QSpinBox()
            self.segment_size_min_spin.setRange(0, 100000)
            self.segment_size_min_spin.setSingleStep(10)
            self.segment_size_min_spin.setValue(int(initial_params.get(f"{self.cell_type}_segment_size_min", 0)))
            row2.addWidget(QLabel("Min size:"))
            row2.addWidget(self.segment_size_min_spin)
            row2.addWidget(HelpButton(
                "Minimum segment size",
                "Segments with fewer voxels than this are discarded after watershed.\n"
                "Use to filter out noise / debris."
            ))
            group_layout.addLayout(row2)

        elif effective_strategy in {
            "APOC Mask + EDT/Watershed Resegmentation",
            "APOC Mask + Peak EDT/Watershed Resegmentation",
        }:
            row1 = QHBoxLayout()
            row2 = QHBoxLayout()
            self.edt_threshold_spin = QDoubleSpinBox()
            self.edt_threshold_spin.setRange(0.0, 50.0)
            self.edt_threshold_spin.setSingleStep(0.5)
            self.edt_threshold_spin.setValue(float(initial_params.get(f"{self.cell_type}_edt_threshold", 1.0)))
            row1.addWidget(QLabel("EDT threshold:"))
            row1.addWidget(self.edt_threshold_spin)
            row1.addWidget(HelpButton(
                "EDT threshold",
                "Euclidean-distance-transform threshold used to derive seeds inside "
                "the binary mask.\n"
                "Lower values give more aggressive splitting of touching objects."
            ))

            self.opening_nr_pixels_spin = QSpinBox()
            self.opening_nr_pixels_spin.setRange(0, 10)
            self.opening_nr_pixels_spin.setValue(int(initial_params.get(f"{self.cell_type}_opening_nr_pixels", 0)))
            row1.addWidget(QLabel("Opening px:"))
            row1.addWidget(self.opening_nr_pixels_spin)
            row1.addWidget(HelpButton(
                "Morphological opening",
                "Number of erosion-then-dilation iterations applied to the mask.\n"
                "Smooths boundaries and removes small speckles. Set to 0 to disable."
            ))
            group_layout.addLayout(row1)

            self.segment_size_min_spin = QSpinBox()
            self.segment_size_min_spin.setRange(0, 100000)
            self.segment_size_min_spin.setSingleStep(10)
            self.segment_size_min_spin.setValue(int(initial_params.get(f"{self.cell_type}_segment_size_min", 0)))
            row2.addWidget(QLabel("Min size:"))
            row2.addWidget(self.segment_size_min_spin)
            row2.addWidget(HelpButton(
                "Minimum segment size",
                "Segments with fewer voxels than this are discarded after watershed.\n"
                "Use to filter out noise / debris."
            ))

            self.fill_holes_cb = QCheckBox("Fill holes")
            self.fill_holes_cb.setChecked(bool(initial_params.get(f"{self.cell_type}_fill_holes", True)))
            row2.addWidget(self.fill_holes_cb)
            row2.addWidget(HelpButton(
                "Fill holes",
                "Fill internal gaps in segmented objects before watershed.\n"
                "Useful for hollow nuclei or partial-volume effects."
            ))
            group_layout.addLayout(row2)

            if effective_strategy == "APOC Mask + Peak EDT/Watershed Resegmentation":
                row3 = QHBoxLayout()
                self.peak_min_distance_spin = QDoubleSpinBox()
                self.peak_min_distance_spin.setRange(0.0, 50.0)
                self.peak_min_distance_spin.setSingleStep(0.5)
                self.peak_min_distance_spin.setValue(float(initial_params.get(f"{self.cell_type}_peak_min_distance", 0.0)))
                self.peak_min_distance_spin.setToolTip("Minimum distance (µm) between local EDT peaks used as watershed seeds")
                row3.addWidget(QLabel("Peak min dist:"))
                row3.addWidget(self.peak_min_distance_spin)
                row3.addWidget(HelpButton(
                    "Peak minimum distance",
                    "Minimum distance (µm) between local EDT peaks used as watershed seeds.\n"
                    "Larger values yield fewer, more separated seeds."
                ))

                self.peak_min_ratio_spin = QDoubleSpinBox()
                self.peak_min_ratio_spin.setRange(0.0, 1.0)
                self.peak_min_ratio_spin.setSingleStep(0.05)
                self.peak_min_ratio_spin.setValue(float(initial_params.get(f"{self.cell_type}_peak_min_ratio", 0.35)))
                self.peak_min_ratio_spin.setToolTip("Minimum EDT peak height as a fraction of the local maximum (0–1)")
                row3.addWidget(QLabel("Peak min ratio:"))
                row3.addWidget(self.peak_min_ratio_spin)
                row3.addWidget(HelpButton(
                    "Peak minimum ratio",
                    "Minimum EDT peak height as a fraction of the local maximum (0–1).\n"
                    "Higher values suppress weak peaks (fewer seeds)."
                ))
                group_layout.addLayout(row3)

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
            self.peak_min_distance_spin,
            self.peak_min_ratio_spin,
        ]:
            if widget is not None:
                widget.valueChanged.connect(self._emit_params_changed)
        if self.fill_holes_cb is not None:
            self.fill_holes_cb.stateChanged.connect(self._emit_params_changed)

    def rebuild_instance_controls(self, strategy=None, initial_params=None):
        """Discard the current ``instance_group`` and rebuild it for *strategy*.

        Used when the per-tab strategy changes (plugin Advanced mode) or when
        the global strategy combo switches and the tab is in inline mode.
        """
        if initial_params is None:
            initial_params = self._initial_params_cache
        else:
            self._initial_params_cache = dict(initial_params)

        # Remove the old group widget from the anchor layout (if attached).
        old_group = self.instance_group
        if old_group is not None:
            try:
                self._instance_group_anchor_layout.removeWidget(old_group)
            except Exception:
                pass
            old_group.setParent(None)
            old_group.deleteLater()

        target_strategy = str(strategy or self._apoc_strategy)
        # Persist the new strategy on the tab so backend helpers that read
        # ``tab._apoc_strategy`` see the up-to-date value.
        self._apoc_strategy = target_strategy
        self._build_instance_controls(initial_params, strategy=target_strategy)

        if self._instance_controls_inline:
            if self.instance_group is not None:
                self._instance_group_anchor_layout.addWidget(self.instance_group)
                self._instance_group_anchor.setVisible(True)
            else:
                self._instance_group_anchor.setVisible(False)

    def _on_per_tab_strategy_combo_changed(self, new_strategy):
        """Handle the user changing the per-tab strategy combo."""
        if self._instance_controls_inline:
            self.rebuild_instance_controls(strategy=str(new_strategy))
        if callable(self._on_per_tab_strategy_changed):
            self._on_per_tab_strategy_changed(self.cell_type, str(new_strategy))

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
            preset_data = FEATURE_PRESETS.get(preset_name, {})
            if "original" in preset_data:
                self.consider_original_cb.blockSignals(True)
                self.consider_original_cb.setChecked(preset_data["original"])
                self.consider_original_cb.blockSignals(False)
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
        """Return the full APOC feature string based on current grid state.

        When features are locked (imported mode) the locked string is returned
        unchanged so the training loop always uses the imported classifier's
        feature specification.
        """
        if getattr(self, "_locked_feature_string", None):
            return self._locked_feature_string
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
        if self.channel_checkboxes:
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
            "peak_min_distance": (
                float(self.peak_min_distance_spin.value()) if self.peak_min_distance_spin is not None else None
            ),
            "peak_min_ratio": (
                float(self.peak_min_ratio_spin.value()) if self.peak_min_ratio_spin is not None else None
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
        if "peak_min_distance" in cfg and self.peak_min_distance_spin is not None and cfg["peak_min_distance"] is not None:
            self.peak_min_distance_spin.setValue(float(cfg["peak_min_distance"]))
        if "peak_min_ratio" in cfg and self.peak_min_ratio_spin is not None and cfg["peak_min_ratio"] is not None:
            self.peak_min_ratio_spin.setValue(float(cfg["peak_min_ratio"]))
        self._update_preview()

    # ── Import / feature-locking API ──────────────────────────────────────────

    def lock_features(self, feature_string, sigmas=None, consider_original=None, channels_used=None):
        """Disable all feature controls and fix the feature string (imported mode).

        Applies the imported settings visually before locking so greyed-out
        controls show the imported values rather than the last user-edited state.
        """
        self._locked_feature_string = feature_string

        # --- Apply imported values visually before disabling ---
        parsed = _parse_feature_string(feature_string)

        # Sigmas: explicit list from metadata takes priority over parsed
        if sigmas:
            sigma_str = ", ".join(_fmt_sigma(s) for s in sigmas)
        else:
            sigma_str = parsed.get("grid_sigmas", "")
        if sigma_str:
            self.sigma_input.setText(sigma_str)
            try:
                parts = sigma_str.replace(",", " ").split()
                new_sigmas = [float(p) for p in parts if p.strip()]
                if new_sigmas:
                    self.current_sigmas = new_sigmas
                    self._build_grid()
            except ValueError:
                pass

        # Feature grid checkboxes
        checked_set = {tuple(pair) for pair in parsed.get("checked_features", [])}
        self._apply_checked_set(checked_set)

        # Consider-original: explicit value takes priority over parsed
        co_value = parsed.get("consider_original", False) if consider_original is None else consider_original
        self.consider_original_cb.blockSignals(True)
        self.consider_original_cb.setChecked(co_value)
        self.consider_original_cb.blockSignals(False)

        # Channel checkboxes
        if channels_used is not None:
            self._default_channel_names = list(channels_used)
            channels_used_set = set(channels_used)
            for cb in self.channel_checkboxes:
                cb.setChecked(cb.text() in channels_used_set)

        # --- Disable all controls ---
        self.feature_combo.setEnabled(False)
        self.sigma_input.setEnabled(False)
        self.update_grid_btn.setEnabled(False)
        self.consider_original_cb.setEnabled(False)
        for cb in self._feat_sigma_checks.values():
            cb.setEnabled(False)
        for ch_cb in self.channel_checkboxes:
            ch_cb.setEnabled(False)
        self.chan_group.setEnabled(False)
        self.tune_group.setEnabled(False)

        short = feature_string if len(feature_string) <= 110 else feature_string[:107] + "…"
        self.preview_label.setText(
            f"<b>🔒 Features locked (imported):</b><br/>"
            f"<span style='font-family:monospace;'>{short}</span>"
        )
        self.preview_label.setToolTip(feature_string)
        self.preview_label.setStyleSheet(
            "color: #c06000; font-size: 10px; padding: 4px 6px; "
            "background: rgba(220,120,0,0.08); border-radius: 3px; "
            "border: 1px solid #c06000;"
        )

    def unlock_features(self):
        """Re-enable all feature controls and clear the import lock."""
        self._locked_feature_string = None

        self.feature_combo.setEnabled(True)
        self.sigma_input.setEnabled(True)
        self.update_grid_btn.setEnabled(True)
        self.consider_original_cb.setEnabled(True)
        for cb in self._feat_sigma_checks.values():
            cb.setEnabled(True)
        for ch_cb in self.channel_checkboxes:
            ch_cb.setEnabled(True)
        self.chan_group.setEnabled(True)
        self.tune_group.setEnabled(True)

        self.preview_label.setStyleSheet(
            "color: #888; font-size: 10px; padding: 2px 4px; "
            "background: rgba(0,0,0,0.05); border-radius: 3px;"
        )
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

    The widget supports two modes selected via constructor flags:

    * ``instance_controls_mode='docked'`` (default, used by the notebook):
      a single shared instance-controls dock below the tab list is shown for
      the currently active tab, identical to the original notebook layout.
    * ``instance_controls_mode='inline'`` (used by the napari plugin):
      each tab owns its own instance-controls group inline, so switching
      tabs does not reparent widgets.

    ``per_cell_type_strategy=True`` (plugin only) adds a global "Strategy"
    combo (including an ``"Advanced (per cell type)"`` option) plus a
    ``\u21b3 Cell type strategy:`` combo on each non-dead tab. In ``"Advanced"``
    mode each tab uses its own strategy; otherwise every tab follows the
    global selection.
    """

    # ── Signals ─────────────────────────────────────────────────────────────
    # Emitted on every successful refresh of one tab's channel checkboxes.
    channels_refreshed = Signal(str)
    # Emitted when a training run starts/finishes (cell types involved).
    training_started = Signal(list)
    training_finished = Signal(list, object)
    # Emitted around per-tab instance preview runs.
    instance_preview_started = Signal(str)
    instance_preview_finished = Signal(str)
    # Emitted when the global strategy combo value changes (plugin mode).
    strategy_changed = Signal(str)
    # Emitted when imported training data is applied or cleared.
    import_applied = Signal()
    import_cleared = Signal()

    STRATEGIES = [
        "APOC (Direct Instance Segmentation)",
        "APOC Mask + EDT/Watershed Resegmentation",
        "APOC Mask + Peak EDT/Watershed Resegmentation",
        "APOC Probability Map + Watershed",
    ]
    import os
    if os.environ.get("BEHAV3D_DEV_MODE") != "1":
        STRATEGIES.remove("APOC Mask + Peak EDT/Watershed Resegmentation")
    ADVANCED_STRATEGY = "Advanced (per cell type)"

    def __init__(
        self,
        viewer,
        pixel_class_outdir,
        all_cell_types,
        has_death,
        initial_params=None,
        on_params_changed=None,
        parent=None,
        instance_controls_mode="docked",
        per_cell_type_strategy=False,
        strategy_resolver=None,
        extra_toolbar_widgets=None,
        pixel_sizes=None,
    ):
        super().__init__(parent)
        self._bg = BackgroundOperation(self)
        self.viewer = viewer
        self.pixel_class_outdir = pixel_class_outdir
        self.all_cell_types = all_cell_types
        self.has_death = has_death
        self._initial_params = initial_params or {}
        self._on_params_changed = on_params_changed
        self._apoc_strategy = str(self._initial_params.get("apoc_strategy", "APOC Probability Map + Watershed"))
        self._log_fn = print
        # Optional pixel sizes dict: {"xy_um": float, "z_um": float}
        self._pixel_sizes = dict(pixel_sizes) if pixel_sizes else {}
        # Imported training data state
        self._imported_training_data: dict = {}
        self._imported_metadata: dict = {}
        self._cell_type_mapping: dict = {}
        self._imported_metadata_path: str = ""

        if instance_controls_mode not in ("docked", "inline"):
            raise ValueError(
                "instance_controls_mode must be 'docked' or 'inline', "
                f"got {instance_controls_mode!r}"
            )
        self._instance_controls_mode = instance_controls_mode
        self._per_cell_type_strategy = bool(per_cell_type_strategy)
        # Optional caller-provided resolver; takes precedence over the
        # built-in Advanced mode logic when set.
        self._strategy_resolver = strategy_resolver
        self._extra_toolbar_widgets = list(extra_toolbar_widgets or [])
        # Recorded so ``cleanup()`` can disconnect them deterministically.
        self._layer_event_handles = []

        # All tab labels = cell types + optional Death
        self._tab_cell_types = list(all_cell_types)
        if has_death:
            self._tab_cell_types.append("dead")

        self._build_ui()
        self._connect_signals()
        if self._instance_controls_mode == "docked":
            self._refresh_instance_controls()
        else:
            # Inline mode: never reparent; keep the shared dock hidden.
            self.instance_controls_widget.setVisible(False)
        self._refresh_all_channels()
        # Apply the initial strategy combo state (plugin Advanced mode).
        if self.strategy_combo is not None:
            self._apply_strategy_to_tabs(
                self.strategy_combo.currentText(),
                emit_signal=False,
            )

        # Listen for layer changes (recorded so cleanup() can disconnect).
        ins_cb = lambda _evt=None: self._refresh_all_channels()
        rem_cb = lambda _evt=None: self._refresh_all_channels()
        self.viewer.layers.events.inserted.connect(ins_cb)
        self.viewer.layers.events.removed.connect(rem_cb)
        self._layer_event_handles.append(
            (self.viewer.layers.events.inserted, ins_cb)
        )
        self._layer_event_handles.append(
            (self.viewer.layers.events.removed, rem_cb)
        )

    def set_log_fn(self, log_fn):
        """Attach an optional GUI log sink."""
        self._log_fn = log_fn or print

    def _log(self, message):
        """Write a message to stdout and, when available, to the GUI log."""
        text = str(message)
        print(text)
        if callable(self._log_fn) and self._log_fn is not print:
            self._log_fn(text)

    # ── Imported training data API ────────────────────────────────────────────

    def apply_import(self, imported_metadata, data_by_celltype, cell_type_mapping, source_path=""):
        """Apply imported training data: store X/y pairs and lock per-tab features.

        Parameters
        ----------
        imported_metadata : dict
            Parsed training_metadata.yml from the source experiment.
        data_by_celltype : dict
            {imported_cell_type: (X, y)} numpy arrays.
        cell_type_mapping : dict
            {imported_cell_type: local_cell_type | None}.  None means skip.
        source_path : str
            Path to the source training_metadata.yml (stored for provenance).
        """
        self._imported_metadata = dict(imported_metadata)
        self._cell_type_mapping = dict(cell_type_mapping)
        self._imported_training_data = {}
        self._imported_metadata_path = str(source_path)

        features_info = imported_metadata.get("features", {})
        for imp_ct, local_ct in cell_type_mapping.items():
            if local_ct is None:
                continue
            if imp_ct in data_by_celltype:
                imp_X, imp_y = data_by_celltype[imp_ct]
                self._imported_training_data[local_ct] = (
                    imp_X.astype(np.float32),
                    imp_y.astype(np.int32),
                )
            if local_ct in self.tabs:
                ct_features = features_info.get(imp_ct, {})
                feature_string = ct_features.get("feature_string", "")
                if feature_string:
                    self.tabs[local_ct].lock_features(
                        feature_string,
                        sigmas=ct_features.get("sigmas"),
                        consider_original=ct_features.get("consider_original"),
                        channels_used=ct_features.get("channels_used"),
                    )

        try:
            self.import_applied.emit()
        except Exception:
            pass

    def clear_import(self):
        """Remove imported training data and unlock all feature controls."""
        self._imported_training_data = {}
        self._imported_metadata = {}
        self._cell_type_mapping = {}
        self._imported_metadata_path = ""
        for tab in self.tabs.values():
            tab.unlock_features()
        try:
            self.import_cleared.emit()
        except Exception:
            pass

    def get_import_summary(self):
        """Return a human-readable summary of the active import (or empty string)."""
        if not self._imported_metadata:
            return ""
        lines = []
        img_meta = self._imported_metadata.get("image_metadata", {})
        xy = img_meta.get("pixel_size_xy_um")
        z  = img_meta.get("pixel_size_z_um")
        if xy is not None:
            lines.append(f"Pixel size: xy={xy} µm, z={z} µm")
        ch_names = img_meta.get("channel_names", [])
        if ch_names:
            lines.append(f"Channels: {', '.join(ch_names)}")
        td_info = self._imported_metadata.get("training_data", {})
        for imp_ct, local_ct in self._cell_type_mapping.items():
            if local_ct is None:
                continue
            counts = td_info.get(imp_ct, {})
            n = counts.get("n_pixels", 0)
            pos = counts.get("n_positive", 0)
            neg = counts.get("n_negative", 0)
            lines.append(f"  {imp_ct} → {local_ct}: {n} px ({pos} fg / {neg} bg)")
        return "\n".join(lines)

    def _build_ui(self):
        # Content widget (scrollable)
        content_widget = QWidget()
        layout = QVBoxLayout()
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(6)
        self._main_layout = layout

        title = QLabel("<h3>APOC Pixel Classification</h3>")
        layout.addWidget(title)

        # Strategy combo (only when per-cell-type strategy is enabled). Lives
        # above the tab widget so users see it before drilling into a tab.
        self.strategy_combo = None
        self.strategy_help_button = None
        if self._per_cell_type_strategy:
            strat_row = QHBoxLayout()
            strat_row.addWidget(QLabel("Strategy:"))
            self.strategy_combo = QComboBox()
            strategies = list(self.STRATEGIES) + [self.ADVANCED_STRATEGY]
            self.strategy_combo.addItems(strategies)
            initial_strategy = self._apoc_strategy if self._apoc_strategy in strategies else strategies[0]
            self.strategy_combo.setCurrentText(initial_strategy)
            strat_row.addWidget(self.strategy_combo, stretch=1)
            layout.addLayout(strat_row)

            strat_desc = QLabel(
                "Strategy determines how the APOC classifier output is converted to"
                " instance segmentation labels. Post-processing parameters appear in"
                " each cell-type tab."
            )
            strat_desc.setWordWrap(True)
            strat_desc.setStyleSheet("color:#888; font-size:10px; padding: 0 0 4px 0;")
            layout.addWidget(strat_desc)

            help_row = QHBoxLayout()
            help_row.setContentsMargins(0, 0, 0, 0)
            help_lbl = QLabel("Instance preview parameters")
            help_lbl.setStyleSheet("color:#888; font-size:10px;")
            self.strategy_help_button = HelpButton(
                *self._strategy_help_content(initial_strategy)
            )
            help_row.addWidget(help_lbl)
            help_row.addWidget(self.strategy_help_button)
            help_row.addStretch()
            layout.addLayout(help_row)

        # --- Tab widget ---
        self.tab_widget = QTabWidget()
        self.tabs = {}  # cell_type -> CellTypeTab
        per_tab_strategies = list(self.STRATEGIES) if self._per_cell_type_strategy else None

        for ct in self._tab_cell_types:
            tab = CellTypeTab(
                ct,
                self.viewer,
                pixel_class_outdir=self.pixel_class_outdir,
                initial_params=self._initial_params,
                apoc_strategy=self._apoc_strategy,
                on_params_changed=self._persist_params,
                run_instance_callback=self._run_instance_preview,
                instance_controls_inline=(self._instance_controls_mode == "inline"),
                show_strategy_combo=self._per_cell_type_strategy,
                per_tab_strategies=per_tab_strategies,
                on_per_tab_strategy_changed=self._on_per_tab_strategy_changed,
            )
            self.tabs[ct] = tab
            self.tab_widget.addTab(tab, ct.capitalize())

        layout.addWidget(self.tab_widget)

        # --- Separator ---
        line = QFrame()
        line.setFrameShape(QFrame.HLine)
        line.setFrameShadow(QFrame.Sunken)
        layout.addWidget(line)

        # Extra toolbar widgets/layouts from the host plugin (e.g. GPU device,
        # workers, run-batch buttons). They appear above the global button row.
        for item in self._extra_toolbar_widgets:
            if isinstance(item, QWidget):
                layout.addWidget(item)
            elif isinstance(item, (QHBoxLayout, QVBoxLayout, QGridLayout)):
                layout.addLayout(item)

        # --- Global buttons ---
        global_row1 = QHBoxLayout()

        self.apply_all_btn = QtPushButton("⬇ Apply config to all tabs")
        self.apply_all_btn.setStyleSheet("padding: 5px 10px;")
        self.apply_all_btn.setToolTip("Copy the current tab's preset, sigmas, depth and trees to ALL other tabs")
        global_row1.addWidget(self.apply_all_btn)

        self.save_labels_btn = QtPushButton("Save User Labels")
        self.save_labels_btn.setStyleSheet(
            "background-color: #fd7e14; color: white; font-weight: bold; padding: 5px 10px;"
        )
        self.save_labels_btn.setToolTip("Save all user-provided label layers to disk")
        self.save_labels_btn.clicked.connect(lambda: self.save_user_labels(log=self._log))
        global_row1.addWidget(self.save_labels_btn)
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
        content_widget.setLayout(layout)

        # Wrap everything in a scroll area so the full widget is reachable
        scroll = QScrollArea()
        scroll.setWidget(content_widget)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        outer_layout = QVBoxLayout()
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.addWidget(scroll)
        self.setLayout(outer_layout)

    def _connect_signals(self):
        self.apply_all_btn.clicked.connect(self._on_apply_to_all)
        self.train_current_btn.clicked.connect(self._on_train_current)
        self.train_all_btn.clicked.connect(self._on_train_all)
        # Forward background progress label to the status_label for live updates.
        self._bg.progress.connect(
            lambda cur, tot, lbl: self.status_label.setText(lbl),
            Qt.QueuedConnection,
        )
        # In docked mode, switching tabs reparents the instance group; in
        # inline mode the group already lives in the tab so we do nothing.
        if self._instance_controls_mode == "docked":
            self.tab_widget.currentChanged.connect(lambda _idx: self._refresh_instance_controls())
        self.tab_widget.currentChanged.connect(self._update_train_current_btn_text)
        if self.strategy_combo is not None:
            self.strategy_combo.currentTextChanged.connect(self._on_global_strategy_changed)
        self._update_train_current_btn_text()

    def _update_train_current_btn_text(self, *_args):
        if hasattr(self, "train_current_btn") and hasattr(self, "tab_widget") and self._tab_cell_types:
            try:
                current_ct = self._tab_cell_types[self.tab_widget.currentIndex()]
                self.train_current_btn.setText(f"▶ Train {current_ct.capitalize()} classifier")
            except IndexError:
                pass


    def _refresh_instance_controls(self):
        """Show the current tab's instance-segmentation controls in the main dock.

        In inline mode this is a no-op because each tab owns its own group.
        """
        if self._instance_controls_mode != "docked":
            return
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
        for ct, tab in self.tabs.items():
            tab.refresh_channel_checkboxes()
            try:
                self.channels_refreshed.emit(ct)
            except Exception:
                pass

    # ------------------------------------------------------------------
    # Strategy resolution / per-cell-type plumbing
    # ------------------------------------------------------------------

    def _resolve_strategy(self, ct):
        """Return the effective APOC strategy for *ct*.

        Resolution order:

        1. Caller-provided ``strategy_resolver(ct)`` (highest precedence).
        2. ``Advanced (per cell type)`` -> per-tab combo value.
        3. The global combo value when the combo exists.
        4. The constructor-supplied ``apoc_strategy`` value.
        """
        if callable(self._strategy_resolver):
            try:
                resolved = self._strategy_resolver(ct)
                if resolved:
                    return str(resolved)
            except Exception:
                pass

        if self.strategy_combo is not None:
            global_choice = self.strategy_combo.currentText()
            if global_choice == self.ADVANCED_STRATEGY:
                tab = self.tabs.get(ct)
                per_combo = getattr(tab, "_per_tab_strategy_combo", None) if tab else None
                if per_combo is not None:
                    return per_combo.currentText()
                return self.STRATEGIES[0]
            return global_choice

        return self._apoc_strategy

    def _strategy_help_content(self, strategy_name):
        """Return ``(title, description)`` for the strategy help popup."""
        if strategy_name == self.ADVANCED_STRATEGY:
            return (
                "APOC Preview Parameters (Advanced \u2014 per cell type)",
                "In Advanced mode each cell-type tab has its own strategy selector.\n\n"
                "Select a strategy inside each tab to show the matching post-processing controls:\n"
                "  \u2022 EDT threshold, Opening, Min size, Fill holes  (EDT/Watershed)\n"
                "  \u2022 Mask threshold, Seed threshold, Opening, Min size  (Probability Map)\n"
                "  \u2022 No extra controls  (Direct)"
            )
        if strategy_name == "APOC Mask + EDT/Watershed Resegmentation":
            return (
                "APOC Preview Parameters (Mask + EDT/Watershed)",
                "Use the classifier mask, then refine instances with EDT + watershed.\n\n"
                "  \u2022 EDT threshold: controls split sensitivity for touching objects.\n"
                "  \u2022 Opening px: smooths boundaries and removes noise.\n"
                "  \u2022 Min size: removes tiny segments below this size.\n"
                "  \u2022 Fill holes: fills internal gaps inside segmented objects."
            )
        if strategy_name == "APOC Mask + Peak EDT/Watershed Resegmentation":
            return (
                "APOC Preview Parameters (Mask + Peak EDT/Watershed)",
                "Same as Mask + EDT/Watershed, but seeds the watershed at local EDT peaks.\n\n"
                "  \u2022 Peak min dist: minimum distance between seed peaks.\n"
                "  \u2022 Peak min ratio: minimum peak height vs. local max (0\u20131)."
            )
        if strategy_name == "APOC Probability Map + Watershed":
            return (
                "APOC Preview Parameters (Probability Map + Watershed)",
                "Use the probability map and watershed for instance splitting.\n\n"
                "  \u2022 Mask threshold: foreground cutoff for the watershed mask.\n"
                "  \u2022 Seed threshold: seed cutoff for watershed seeds.\n"
                "  \u2022 Opening px: smooths boundaries and removes noise.\n"
                "  \u2022 Min size: removes tiny segments below this size."
            )
        return (
            "APOC Preview Parameters (Direct Instance Segmentation)",
            "Direct strategy predicts instances directly from the classifier output."
            " No post-processing parameters are available in preview mode."
        )

    def _apply_strategy_to_tabs(self, strategy, emit_signal=True):
        """Propagate a strategy change to all tabs (inline mode rebuilds)."""
        strategy = str(strategy)
        self._apoc_strategy = (
            strategy
            if strategy != self.ADVANCED_STRATEGY
            else self._apoc_strategy
        )

        if self.strategy_help_button is not None:
            title, desc = self._strategy_help_content(strategy)
            self.strategy_help_button.set_help(title, desc)

        is_advanced = strategy == self.ADVANCED_STRATEGY

        for ct, tab in self.tabs.items():
            per_widget = getattr(tab, "_per_tab_strategy_widget", None)
            if per_widget is not None:
                per_widget.setVisible(is_advanced and ct != "dead")
            if self._instance_controls_mode == "inline":
                effective = self._resolve_strategy(ct)
                tab.rebuild_instance_controls(strategy=effective)

        if self._instance_controls_mode == "docked":
            self._refresh_instance_controls()

        if emit_signal:
            try:
                self.strategy_changed.emit(strategy)
            except Exception:
                pass

    def _on_global_strategy_changed(self, new_strategy):
        """Handle the global strategy combo changing."""
        self._apply_strategy_to_tabs(new_strategy, emit_signal=True)
        self._persist_params()

    def _on_per_tab_strategy_changed(self, cell_type, new_strategy):
        """Handle a per-tab strategy combo changing.

        Triggered by ``CellTypeTab`` when ``per_cell_type_strategy=True``.
        Persists the per-tab choice so it survives a reload.
        """
        self._initial_params[f"per_ct_strategy_{cell_type}"] = new_strategy
        self._persist_params()

    # ------------------------------------------------------------------
    # Enable/disable + teardown helpers
    # ------------------------------------------------------------------

    def set_training_enabled(self, enabled):
        """Enable or disable the training-related controls in one call.

        The plugin uses this to lock the widget until training data has
        been loaded.
        """
        enabled = bool(enabled)
        self.apply_all_btn.setEnabled(enabled)
        self.train_current_btn.setEnabled(enabled)
        self.train_all_btn.setEnabled(enabled)
        for tab in self.tabs.values():
            for cb in getattr(tab, "channel_checkboxes", []):
                cb.setEnabled(enabled)
            btn = getattr(tab, "run_instance_btn", None)
            if btn is not None:
                btn.setEnabled(enabled)

    def set_classifier_params_visible(self, visible):
        """Show or hide the classifier-specific parameters while keeping segmentation parameters active."""
        self.apply_all_btn.setVisible(visible)
        self.train_current_btn.setVisible(visible)
        self.train_all_btn.setVisible(visible)
        self.save_labels_btn.setVisible(visible)
        for tab in self.tabs.values():
            if hasattr(tab, "set_classifier_params_visible"):
                tab.set_classifier_params_visible(visible)

    def cleanup(self):
        """Disconnect viewer-layer callbacks before the widget is destroyed.

        The widget connects to ``viewer.layers.events.inserted`` and
        ``removed`` for live channel refresh; if those callbacks remain
        connected after the widget is hidden/deleted, napari can crash when
        layers change later. Call this from the host before ``deleteLater``.
        """
        for emitter, callback in list(self._layer_event_handles):
            try:
                emitter.disconnect(callback)
            except Exception:
                pass
        self._layer_event_handles = []

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
        self._run_training_bg([current_ct])

    def _on_train_all(self):
        """Train classifiers for ALL cell types using each tab's individual config."""
        self._run_training_bg(self._tab_cell_types)

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
        strategy = self._resolve_strategy(ct)
        if ct == "dead" or strategy == "APOC (Direct Instance Segmentation)":
            return np.asarray(raw_prediction).astype(np.int16)

        tab = self.tabs[ct]
        if strategy == "APOC Probability Map + Watershed":
            return _probability_array_to_segments(
                prob_prediction,
                mask_thr=float(tab.prob_mask_threshold_spin.value()),
                seed_thr=float(tab.prob_seed_threshold_spin.value()),
                opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
                segment_size_min=int(tab.segment_size_min_spin.value()),
            )

        marker_strategy = (
            "peak"
            if strategy == "APOC Mask + Peak EDT/Watershed Resegmentation"
            else "threshold"
        )
        return _mask_array_to_segments(
            raw_prediction > 0,
            edt_thr=float(tab.edt_threshold_spin.value()),
            opening_nr_pixels=int(tab.opening_nr_pixels_spin.value()),
            segment_size_min=int(tab.segment_size_min_spin.value()),
            fill_holes=bool(tab.fill_holes_cb.isChecked()),
            marker_strategy=marker_strategy,
            peak_min_distance=float(tab.peak_min_distance_spin.value()) if tab.peak_min_distance_spin is not None else None,
            peak_min_ratio=float(tab.peak_min_ratio_spin.value()) if tab.peak_min_ratio_spin is not None else 0.35,
        )

    def _run_instance_preview(self, ct):
        """Run instance segmentation preview for a non-dead APOC tab."""
        if ct == "dead":
            return

        tab = self.tabs[ct]
        strategy = self._resolve_strategy(ct)
        if strategy == "APOC (Direct Instance Segmentation)":
            self.status_label.setText("Direct APOC does not use instance resegmentation preview.")
            return

        try:
            self.instance_preview_started.emit(ct)
        except Exception:
            pass

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
                mask_thr = float(tab.prob_mask_threshold_spin.value())
                seed_thr = float(tab.prob_seed_threshold_spin.value())
                opening_nr_pixels = int(tab.opening_nr_pixels_spin.value())
                segment_size_min = int(tab.segment_size_min_spin.value())
                print(f"  ⚙ {ct} preview params: mask_thr={mask_thr}, seed_thr={seed_thr}, min_size={segment_size_min}, opening_px={opening_nr_pixels}")
                instance_preview = _probability_array_to_segments(
                    prob_prediction,
                    mask_thr=mask_thr,
                    seed_thr=seed_thr,
                    opening_nr_pixels=opening_nr_pixels,
                    segment_size_min=segment_size_min,
                )
            else:
                raw_prediction, prob_prediction = self._predict_classifier_outputs(ct)
                marker_strategy = (
                    "peak"
                    if strategy == "APOC Mask + Peak EDT/Watershed Resegmentation"
                    else "threshold"
                )
                edt_thr = float(tab.edt_threshold_spin.value())
                opening_nr_pixels = int(tab.opening_nr_pixels_spin.value())
                segment_size_min = int(tab.segment_size_min_spin.value())
                fill_holes = bool(tab.fill_holes_cb.isChecked())
                peak_min_distance = float(tab.peak_min_distance_spin.value()) if tab.peak_min_distance_spin is not None else None
                peak_min_ratio = float(tab.peak_min_ratio_spin.value()) if tab.peak_min_ratio_spin is not None else 0.35
                print(f"  ⚙ {ct} preview params: edt_thr={edt_thr}, min_size={segment_size_min}, fill_holes={fill_holes}, opening_px={opening_nr_pixels}, peak_min_distance={peak_min_distance}, peak_min_ratio={peak_min_ratio}")
                instance_preview = _mask_array_to_segments(
                    raw_prediction > 0,
                    edt_thr=edt_thr,
                    opening_nr_pixels=opening_nr_pixels,
                    segment_size_min=segment_size_min,
                    fill_holes=fill_holes,
                    marker_strategy=marker_strategy,
                    peak_min_distance=peak_min_distance,
                    peak_min_ratio=peak_min_ratio,
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
        finally:
            try:
                self.instance_preview_finished.emit(ct)
            except Exception:
                pass

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
                shutil.rmtree(outpath, onexc=ignore_missing_rmtree_error)
            save_as_zarr(label_data, outpath)
            log(f"Saved {cell_type} labels → {outpath}")

        if self.has_death and "User Provided Labels (Dead)" in [l.name for l in self.viewer.layers]:
            dead_layer = self.viewer.layers["User Provided Labels (Dead)"]
            dead_label_data = np.asarray(dead_layer.data)
            dead_outpath = Path(self.pixel_class_outdir, "PixelClassifier_UserDeadLabels.zarr")
            if dead_outpath.exists():
                shutil.rmtree(dead_outpath, onexc=ignore_missing_rmtree_error)
            save_as_zarr(dead_label_data, dead_outpath)
            log(f"Saved Death labels → {dead_outpath}")

        log("✅ All user labels saved!")

    def _run_training_bg(self, cell_types_to_train):
        """Background-threaded version of :meth:`_run_training`.

        Phase 1 (Qt thread): read viewer layers & widget values into plain
        NumPy arrays and APOC config dicts.
        Phase 2 (worker thread): fit classifiers, run predictions.
        Phase 3 (Qt thread, via on_done): write result arrays to viewer layers,
        persist params, update status label.
        """
        import apoc

        if self._bg.is_running():
            self._log("⚠️ Training is already in progress.")
            return

        cell_types_to_train = list(cell_types_to_train)

        # ── Phase 1 (Qt thread): emit signal + gather all data ────────────
        try:
            self.training_started.emit(cell_types_to_train)
        except Exception:
            pass

        self._log("Auto-saving user labels before training…")
        self.save_user_labels(log=self._log)

        # Collect per-cell-type config, images, and annotation arrays.
        ct_inputs = {}
        for ct in cell_types_to_train:
            tab = self.tabs[ct]
            cfg = tab.get_config()
            feature_string = cfg["feature_string"]

            # Images: convert to numpy so the worker never touches viewer layers.
            images_np = [
                np.asarray(img) for img in self._get_images_for_tab(ct)
            ]
            if not images_np:
                self._log(f"  ⚠️ No image layers selected for '{ct}' — skipping.")
                continue

            # Annotation layer data.
            if ct == "dead":
                layer_name = "User Provided Labels (Dead)"
            else:
                layer_name = f"User Provided Labels ({ct.capitalize()})"
            try:
                annotation = np.asarray(self.viewer.layers[layer_name].data)
            except KeyError:
                self._log(f"  ⚠️ Annotation layer '{layer_name}' not found — skipping.")
                continue
            if not np.any(annotation):
                self._log(f"  ⚠️ No labels drawn for '{ct}' — skipping.")
                continue

            # Snapshot per-tab segmentation params (must happen on Qt thread).
            strategy = self._resolve_strategy(ct)
            tab_params = {}
            if hasattr(tab, "edt_threshold_spin") and tab.edt_threshold_spin is not None:
                tab_params["edt_thr"] = float(tab.edt_threshold_spin.value())
            if hasattr(tab, "opening_nr_pixels_spin") and tab.opening_nr_pixels_spin is not None:
                tab_params["opening"] = int(tab.opening_nr_pixels_spin.value())
            if hasattr(tab, "segment_size_min_spin") and tab.segment_size_min_spin is not None:
                tab_params["min_size"] = int(tab.segment_size_min_spin.value())
            if hasattr(tab, "fill_holes_cb") and tab.fill_holes_cb is not None:
                tab_params["fill_holes"] = bool(tab.fill_holes_cb.isChecked())
            if hasattr(tab, "prob_mask_threshold_spin") and tab.prob_mask_threshold_spin is not None:
                tab_params["prob_mask_thr"] = float(tab.prob_mask_threshold_spin.value())
            if hasattr(tab, "prob_seed_threshold_spin") and tab.prob_seed_threshold_spin is not None:
                tab_params["prob_seed_thr"] = float(tab.prob_seed_threshold_spin.value())
            if hasattr(tab, "peak_min_distance_spin") and tab.peak_min_distance_spin is not None:
                tab_params["peak_min_dist"] = float(tab.peak_min_distance_spin.value())
            if hasattr(tab, "peak_min_ratio_spin") and tab.peak_min_ratio_spin is not None:
                tab_params["peak_min_ratio"] = float(tab.peak_min_ratio_spin.value())

            if ct == "dead":
                clf_name = "PixelClassifier_Death.cl"
            else:
                clf_name = f"PixelClassifier_{ct.capitalize()}.cl"
            clf_path = str(Path(self.pixel_class_outdir, clf_name))

            ct_inputs[ct] = {
                "cfg": cfg,
                "feature_string": feature_string,
                "images_np": images_np,
                "annotation": annotation,
                "clf_path": clf_path,
                "strategy": strategy,
                "tab_params": tab_params,
            }

        if not ct_inputs:
            self._log("⚠️ No cell types with valid data to train.")
            try:
                self.training_finished.emit(cell_types_to_train, [])
            except Exception:
                pass
            return

        pixel_class_outdir = Path(self.pixel_class_outdir)
        safe_log = ThreadSafeLogger(self._log)
        training_start = time.time()
        new_pixel_counts: dict = {}     # {ct: n_new_pixels} populated in _do_train
        pixel_counts_by_ct: dict = {}   # {ct: {"total", "positive", "negative"}}

        # ── Phase 2 (worker thread) ───────────────────────────────────────
        def _do_train(progress_cb=None):
            import apoc as _apoc
            from sklearn.ensemble import RandomForestClassifier

            cts = list(ct_inputs.keys())
            n = len(cts)
            results = {}  # ct -> {"display_segments", "prob_result"}

            for i, ct in enumerate(cts):
                if progress_cb:
                    progress_cb(i, n, f"Training {ct}…")
                safe_log(f"  Training {ct}…")

                inp = ct_inputs[ct]
                cfg            = inp["cfg"]
                feature_string = inp["feature_string"]
                images_np      = inp["images_np"]
                annotation     = inp["annotation"]
                clf_path       = inp["clf_path"]
                strategy       = inp["strategy"]
                tab_params     = inp["tab_params"]
                max_depth      = cfg["max_depth"]
                num_ensembles  = cfg["num_ensembles"]

                # Erase existing classifier
                if Path(clf_path).exists():
                    _apoc.erase_classifier(clf_path)

                clf = _apoc.ObjectSegmenter(
                    opencl_filename=clf_path,
                    max_depth=max_depth,
                    num_ensembles=num_ensembles,
                    positive_class_identifier=2,
                )

                # Build X/y arrays
                X_parts, y_parts = [], []
                gt_ndim = None
                n_timepoints = annotation.shape[0] if annotation.ndim == 4 else 1

                if annotation.ndim == 4:
                    for t in range(n_timepoints):
                        ann_t = annotation[t]
                        if not np.any(ann_t):
                            continue
                        feats_t = self._generate_feature_list_for_timepoint(
                            images_np, feature_string, t
                        )
                        X_t, y_t = clf._to_np(feats_t, ann_t)
                        if X_t.size == 0 or y_t.size == 0:
                            continue
                        X_parts.append(X_t)
                        y_parts.append(y_t)
                        gt_ndim = ann_t.ndim
                else:
                    feats_np = self._generate_feature_list_for_timepoint(
                        images_np, feature_string, 0
                    )
                    X_t, y_t = clf._to_np(feats_np, annotation)
                    if X_t.size > 0 and y_t.size > 0:
                        X_parts.append(X_t)
                        y_parts.append(y_t)
                        gt_ndim = annotation.ndim

                if not X_parts:
                    safe_log(f"  ⚠️ No labeled pixels for '{ct}' — skipping.")
                    continue

                X = np.concatenate(X_parts, axis=0)
                y = np.concatenate(y_parts, axis=0)

                # Track new-label pixel counts before import prepend
                new_pixel_counts[ct] = int(len(y))

                # Prepend imported training data when available for this type
                imp_pair = self._imported_training_data.get(ct)
                if imp_pair is not None:
                    imp_X, imp_y = imp_pair
                    safe_log(
                        f"  {ct}: combining {len(imp_y)} imported + "
                        f"{len(y)} new = {len(imp_y) + len(y)} total pixels"
                    )
                    X = np.concatenate([imp_X, X], axis=0)
                    y = np.concatenate([imp_y, y], axis=0)

                fitted_rf = RandomForestClassifier(
                    max_depth=max_depth, n_estimators=num_ensembles, random_state=0
                )
                fitted_rf.fit(X, y)
                clf.classifier              = fitted_rf
                clf._feature_importances   = fitted_rf.feature_importances_
                clf._X                     = X
                clf._y                     = y
                clf.num_features           = X.shape[1]
                clf.num_ground_truth_dimensions = gt_ndim
                clf.feature_specification  = feature_string
                clf.to_opencl_file(clf_path)
                _write_classifier_channel_metadata(clf_path, cfg["channels"])
                safe_log(f"  Saved classifier: {Path(clf_path).name}")

                # Save flattened training data (combined) for future import
                try:
                    _save_training_data(pixel_class_outdir, ct, X, y)
                    n_pos = int(np.sum(y == 2))
                    n_neg = int(np.sum(y == 1))
                    pixel_counts_by_ct[ct] = {
                        "total":    len(y),
                        "positive": n_pos,
                        "negative": n_neg,
                    }
                except Exception as _exc:
                    safe_log(f"  ⚠️ Could not save training data for '{ct}': {_exc}")

                # Predict (for visual confirmation in viewer).
                if progress_cb:
                    progress_cb(i, n, f"Predicting {ct}…")
                raw_prediction, prob_result = self._predict_classifier_outputs(
                    ct, clf=clf
                )

                # Build display segments using snapshotted tab params
                # (mirrors _build_display_segments but without touching Qt widgets).
                if ct == "dead" or strategy == "APOC (Direct Instance Segmentation)":
                    display_segments = np.asarray(raw_prediction).astype(np.int16)
                elif strategy == "APOC Probability Map + Watershed":
                    display_segments = _probability_array_to_segments(
                        prob_result,
                        mask_thr=tab_params.get("prob_mask_thr", 0.5),
                        seed_thr=tab_params.get("prob_seed_thr", 0.5),
                        opening_nr_pixels=tab_params.get("opening", 0),
                        segment_size_min=tab_params.get("min_size", 10),
                    )
                else:
                    marker_strategy = (
                        "peak"
                        if strategy == "APOC Mask + Peak EDT/Watershed Resegmentation"
                        else "threshold"
                    )
                    display_segments = _mask_array_to_segments(
                        raw_prediction > 0,
                        edt_thr=tab_params.get("edt_thr", 1.0),
                        opening_nr_pixels=tab_params.get("opening", 0),
                        segment_size_min=tab_params.get("min_size", 10),
                        fill_holes=tab_params.get("fill_holes", True),
                        marker_strategy=marker_strategy,
                        peak_min_distance=tab_params.get("peak_min_dist"),
                        peak_min_ratio=tab_params.get("peak_min_ratio", 0.35),
                    )

                # Save preview arrays to disk (file I/O is fine in the worker).
                pred_path = _predicted_labels_path(pixel_class_outdir, ct)
                if pred_path is not None:
                    self._save_preview_array(pred_path, display_segments)
                prob_path = _probability_map_path(pixel_class_outdir, ct)
                if prob_path is not None:
                    self._save_preview_array(prob_path, prob_result)

                results[ct] = {
                    "display_segments": display_segments,
                    "prob_result": prob_result,
                }

            return results

        # ── Phase 3 (Qt thread, on_done): write results to viewer ─────────
        def _on_done(results):
            successes = list(results.keys())
            for ct, r in results.items():
                self._update_prediction_layers(ct, r["display_segments"], r["prob_result"])
                # Auto-show statistics
                if ct in self.tabs:
                    self.tabs[ct]._on_show_statistics()

            elapsed_s   = time.time() - training_start
            elapsed_txt = f"{elapsed_s:.1f}s"
            if successes:
                self.status_label.setText(
                    f"✅ Trained: {', '.join(successes)} ({elapsed_txt})"
                )
                self._log(f"✅ Training finished in {elapsed_txt}")
            else:
                self.status_label.setText(
                    f"⚠️ No cell types were trained (check labels). ({elapsed_txt})"
                )
                self._log(f"⚠️ No cell types were trained. Elapsed: {elapsed_txt}")

            self._persist_params()

            # Save training metadata YAML when at least one cell type was trained
            if successes:
                try:
                    params_per_ct = {}
                    for ct in list(self.all_cell_types) + (["dead"] if self.has_death else []):
                        if ct in self.tabs:
                            params_per_ct[ct] = self.tabs[ct].get_config()
                    image_layers = [
                        la for la in self.viewer.layers
                        if hasattr(la, "data") and la.name.startswith("Channel ")
                    ]
                    n_channels = len(image_layers) if image_layers else None
                    ch_names = [la.name for la in image_layers]
                    imported_counts = {}
                    for imp_ct, local_ct in self._cell_type_mapping.items():
                        if local_ct is None:
                            continue
                        imp_pair = self._imported_training_data.get(local_ct)
                        if imp_pair is not None:
                            imported_counts[imp_ct] = int(len(imp_pair[1]))
                    meta = _build_training_metadata(
                        cell_types=self.all_cell_types,
                        has_death=self.has_death,
                        params_per_ct=params_per_ct,
                        n_input_channels=n_channels,
                        channel_names=ch_names,
                        pixel_size_xy_um=self._pixel_sizes.get("xy_um"),
                        pixel_size_z_um=self._pixel_sizes.get("z_um"),
                        pixel_counts_by_ct=pixel_counts_by_ct,
                        imported_from=self._imported_metadata_path or None,
                        imported_counts=imported_counts if imported_counts else None,
                        new_counts=new_pixel_counts if new_pixel_counts else None,
                    )
                    _save_training_metadata(self.pixel_class_outdir, meta)
                    self._log("📄 Training data bundle saved to PixelClassifier_TrainingData.zip")
                except Exception as _exc:
                    self._log(f"  ⚠️ Could not save training metadata: {_exc}")

            try:
                self.training_finished.emit(cell_types_to_train, successes)
            except Exception:
                pass

        def _on_failed(err):
            self.status_label.setText(f"❌ Error: {err}")
            self._log(f"❌ Training error: {err}")
            try:
                self.training_finished.emit(cell_types_to_train, [])
            except Exception:
                pass

        self._bg.run(
            fn=_do_train,
            desc="APOC training…",
            buttons=[self.train_current_btn, self.train_all_btn],
            viewer=self.viewer,
            on_done=_on_done,
            on_failed=_on_failed,
            inject_progress=True,
            indeterminate=False,
        )

    def _run_training(self, cell_types_to_train):
        """Train (and apply) APOC classifiers for the given list of cell types."""
        import apoc
        training_start = time.time()

        try:
            self.training_started.emit(list(cell_types_to_train))
        except Exception:
            pass

        # Auto-save labels before training
        self._log("Auto-saving user labels before training...")
        self.save_user_labels(log=self._log)

        successes = []
        pixel_counts_by_ct: dict = {}
        new_pixel_counts: dict = {}
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
                    self._log(f"Skipping '{ct}': annotation layer not found")
                    continue

                if not np.any(annotation):
                    self._log(f"Skipping '{ct}': no labels drawn")
                    continue

                if ct == "dead":
                    clf_name = "PixelClassifier_Death.cl"
                else:
                    clf_name = f"PixelClassifier_{ct.capitalize()}.cl"

                clf_path = str(Path(self.pixel_class_outdir, clf_name))

                # Erase existing classifier
                if Path(clf_path).exists():
                    apoc.erase_classifier(clf_path)

                clf = apoc.ObjectSegmenter(
                    opencl_filename=clf_path,
                    max_depth=max_depth,
                    num_ensembles=num_ensembles,
                    positive_class_identifier=2,
                )

                has_trained = False
                n_timepoints = annotation.shape[0] if annotation.ndim == 4 else 1

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
                    X_new = np.concatenate(X_parts, axis=0)
                    y_new = np.concatenate(y_parts, axis=0)

                    # Track new-label pixel counts before adding imported data
                    new_pixel_counts[ct] = int(len(y_new))

                    # Prepend imported training data when available for this type
                    imp_pair = self._imported_training_data.get(ct)
                    if imp_pair is not None:
                        imp_X, imp_y = imp_pair
                        X = np.concatenate([imp_X, X_new], axis=0)
                        y = np.concatenate([imp_y, y_new], axis=0)
                        self._log(
                            f"  {ct}: combined {len(imp_y)} imported + "
                            f"{len(y_new)} new = {len(y)} total pixels"
                        )
                    else:
                        X, y = X_new, y_new

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
                    clf.feature_specification = feature_string
                    has_trained = True

                    # Save flattened training data (combined) for future import
                    try:
                        _save_training_data(self.pixel_class_outdir, ct, X, y)
                        n_pos = int(np.sum(y == 2))
                        n_neg = int(np.sum(y == 1))
                        pixel_counts_by_ct[ct] = {
                            "total":    len(y),
                            "positive": n_pos,
                            "negative": n_neg,
                        }
                    except Exception as _exc:
                        self._log(f"  ⚠️ Could not save training data for '{ct}': {_exc}")

                if not has_trained:
                    continue

                clf.feature_specification = feature_string
                clf.to_opencl_file(clf_path)
                _write_classifier_channel_metadata(clf_path, cfg["channels"])

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
                for ct in successes:
                    if ct in self.tabs:
                        self.tabs[ct]._on_show_statistics()
            else:
                elapsed_s = time.time() - training_start
                elapsed_txt = f"{elapsed_s:.1f}s"
                self.status_label.setText(f"⚠️ No cell types were trained (check labels). ({elapsed_txt})")
                self._log(f"⚠️ No cell types were trained. Elapsed time: {elapsed_txt}")

            self._persist_params()

            # Save training metadata YAML when at least one cell type was trained
            if successes:
                try:
                    params_per_ct = {}
                    for ct in list(self.all_cell_types) + (["dead"] if self.has_death else []):
                        if ct in self.tabs:
                            cfg = self.tabs[ct].get_config()
                            params_per_ct[ct] = cfg
                    image_layers = [
                        la for la in self.viewer.layers
                        if hasattr(la, "data") and la.name.startswith("Channel ")
                    ]
                    n_channels = len(image_layers) if image_layers else None
                    ch_names = [la.name for la in image_layers]
                    imported_counts = {}
                    for imp_ct, local_ct in self._cell_type_mapping.items():
                        if local_ct is None:
                            continue
                        imp_pair = self._imported_training_data.get(local_ct)
                        if imp_pair is not None:
                            imported_counts[imp_ct] = int(len(imp_pair[1]))
                    meta = _build_training_metadata(
                        cell_types=self.all_cell_types,
                        has_death=self.has_death,
                        params_per_ct=params_per_ct,
                        n_input_channels=n_channels,
                        channel_names=ch_names,
                        pixel_size_xy_um=self._pixel_sizes.get("xy_um"),
                        pixel_size_z_um=self._pixel_sizes.get("z_um"),
                        pixel_counts_by_ct=pixel_counts_by_ct,
                        imported_from=self._imported_metadata_path or None,
                        imported_counts=imported_counts if imported_counts else None,
                        new_counts=new_pixel_counts if new_pixel_counts else None,
                    )
                    _save_training_metadata(self.pixel_class_outdir, meta)
                    self._log("📄 Training data bundle saved to PixelClassifier_TrainingData.zip")
                except Exception as _exc:
                    self._log(f"  ⚠️ Could not save training metadata: {_exc}")

        except Exception as e:
            self.status_label.setText(f"❌ Error: {e}")
            import traceback
            traceback.print_exc()
        finally:
            try:
                self.training_finished.emit(list(cell_types_to_train), list(successes))
            except Exception:
                pass

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
            if cfg["peak_min_distance"] is not None:
                params[f"{ct}_peak_min_distance"] = cfg["peak_min_distance"]
            if cfg["peak_min_ratio"] is not None:
                params[f"{ct}_peak_min_ratio"] = cfg["peak_min_ratio"]
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
