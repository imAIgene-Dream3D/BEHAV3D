"""Cellpose prediction utilities.

This module provides a single convenience function `run_cellpose_prediction` that
wraps the ad-hoc logic from the original `cellpose_TKI_tcell_run.py` demo script
into a re-usable function that can be imported and called from notebooks (e.g.
`run_behav3d.ipynb`).

The function
1. Loads the image (OME-Tiff, regular Tiff, or Zarr) using `behav3d.io.images`.
2. If the input is a Tiff, it is automatically converted to Zarr (chunked per
   time point) next to the original file so that subsequent calls can re-use
   the Zarr version directly.
3. Runs a pretrained 3-D Cellpose model **per time point** and returns the masks
   as a single NumPy array with shape `(T, Z, Y, X)`.
4. Optionally writes the masks to disk (Zarr).

The implementation keeps all heavy-weight dependencies local to the function so
that importing the module is cheap (important for lightweight environments such
as documentation builds).
"""
from __future__ import annotations

import shutil
import yaml
from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import torch
import pandas as pd
from skimage.filters import threshold_otsu

from behav3d.io.images import save_as_zarr, load_zarr, load_image, append_to_zarr
from behav3d.preprocessing.segmentation import segment_2d_filter
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    has_dead_channel,
    resolve_metadata_csv_path,
)

# Must match CELLPOSE_CHANNEL_UNUSED in behav3d.widgets.segmentation
_CELLPOSE_CHANNEL_UNUSED_LABEL = "(unused)"


def _label_to_channel_from_stored_map(stored: dict | None) -> dict[str, int]:
    """Invert {ch_idx: label} to {label: ch_idx}; skip '(unused)' / empty slots."""
    if not stored:
        return {}
    out: dict[str, int] = {}
    for ch_idx, label in stored.items():
        if label is None:
            continue
        s = str(label).strip()
        if not s or s.lower() == "none" or s == _CELLPOSE_CHANNEL_UNUSED_LABEL:
            continue
        out[str(label)] = int(ch_idx)
    return out

# Cellpose import is relatively expensive - only load when we actually need it.

try:
    from cellpose import models  # noqa: WPS433 (allow external import)
except ImportError as err:  # pragma: no cover - handled at runtime
    raise ImportError(
        "cellpose is required for `behav3d.preprocessing.segmentation.cellpose_prediction`. "
        "Please install it via `pip install cellpose==3.1.1.2` (or compatible).",
    ) from err

def run_cellpose_prediction(  # noqa: WPS231 (complexity) - unavoidable for pipeline function
    image: np.ndarray,  #already loaded image
    pretrained_model_dir: str | Path,
    raw_image_path: str | Path, 
    channels: Optional[Sequence[int]] = None,
    diameter: Optional[int | float] = None,
    min_size: int = 30,
    nchan: int = 2,
    anisotropy: Optional[float] = 4.655 / 1.17365196872,
    timepoint_range: Optional[Tuple[int, int]] = None,
    device: Optional[str] = None,
    convert_tiff_to_zarr: bool = True,
    save_masks: bool = True,
    masks_suffix: str = "_cellpose_masks.zarr",
    masks_outpath: Optional[str | Path] = None,
    verbose: bool = True,
) -> Tuple[np.ndarray, Path]:
    """Run 3-D Cellpose prediction on *image_path*.

    Parameters
    ----------
    image_path:
        Path to input image. Supported extensions: `.tif`, `.tiff`, `.zarr`.
    pretrained_model_dir:
        Directory that contains the pretrained Cellpose weights (same argument
        as *pretrained_model* in ``cellpose.models.CellposeModel``).
    channels, diameter, anisotropy, min_size:
        Passed through to ``CellposeModel.eval``.
    nchan:
        Passed when building the model through ``CellposeModel``.
    device:
        CUDA device string (e.g. ``"cuda:0"``); when ``None`` the function will
        automatically pick GPU if available otherwise CPU.
    convert_tiff_to_zarr:
        If *image_path* is a Tiff and this flag is ``True`` the image is stored
        as Zarr next to the original file for future reuse.
    save_masks:
        When ``True`` the resulting mask array is saved as Zarr (path is the
        original image stem + *masks_suffix*, unless *masks_outpath* is provided).
    masks_suffix:
        Suffix to use when writing the masks to disk (ignored when
        *save_masks* is ``False`` or *masks_outpath* is provided).
    masks_outpath:
        Optional explicit output path for the masks. If provided, overrides the
        default naming convention (raw_image_stem + masks_suffix).
    verbose:
        Emit progress information via ``print``.

    Returns
    -------
    masks:
        NumPy array with shape ``(T, Z, Y, X)`` storing the label masks.
    masks_path:
        Path to the saved mask Zarr file when *save_masks* is ``True``;
        otherwise the path is ``Path()`` (empty path).
    """
    # ---------------------------------------------------------------------
    # 1. Prepare input & model
    # ---------------------------------------------------------------------
    if verbose:
        print("Loaded image with shape (T, C, Z, Y, X):", image.shape)
    
    T_total = image.shape[0]
    # Optional time sub-setting - we no longer slice the image here to maintain global indexing
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, T_total - 1))
        end_t = max(start_t, min(end_t, T_total - 1))
        if verbose:
            print(f"Range requested: {start_t}-{end_t} (Global total: {T_total})")
            
    T, C, Z, Y, X = image.shape

    # Cellpose works with at least 2 channels. If we have one, include an additional empty one
    if (C==1):
        image=np.concatenate((image,np.zeros_like(image)),axis=1)
        nchan=2

    # Determine computation device.
    torch_device = (
        torch.device(device) if device is not None else torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    )

    # Initialise Cellpose model.
    model = models.CellposeModel(
        pretrained_model=str(pretrained_model_dir),
        gpu=True,
        nchan=nchan,
        device=torch_device
    )

    # ---------------------------------------------------------------------
    # 2. Run prediction per time point
    # ---------------------------------------------------------------------
    # Global indexing: Initialize masks array with the FULL length of the input data
    # (T_total if we have it, else T from the current image view)
    masks = np.zeros((T_total, Z, Y, X), dtype=np.uint16)
    
    # Try to load existing results to support incremental updates
    existing_masks_loaded = False
    if masks_outpath is not None:
        target_path = Path(masks_outpath)
        if target_path.exists():
            try:
                existing_arr = load_image(target_path)
                if existing_arr.shape == masks.shape:
                    print(f"  [INCREMENTAL] Loading existing data from {target_path.name}...")
                    masks[:] = existing_arr.compute()
                    existing_masks_loaded = True
                else:
                    print(f"  [OVERWRITE] Existing array shape {existing_arr.shape} != {masks.shape}. Re-initializing.")
            except Exception as e:
                print(f"  [WARNING] Could not load existing masks: {e}")

    # If we are using a range, we iterate through the requested range.
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        print(f"  Global indexing: processing T {start_t} to {end_t} (Total: {T_total})")
        indices_to_process = range(start_t, end_t + 1)
    else:
        indices_to_process = range(T_total)

    if verbose:
        print("Starting Cellpose prediction...")

    for t in indices_to_process:
        # Extract data for one time point - shape (C, Z, Y, X).
        img_t = image[t]

        if verbose:
            print(f"  - Timepoint {t+1}/{T_total} -", end="", flush=True)
        start_ts = torch.cuda.Event(enable_timing=True) if torch_device.type == "cuda" else None
        end_ts = torch.cuda.Event(enable_timing=True) if torch_device.type == "cuda" else None
        if start_ts is not None:
            start_ts.record()

        # Cellpose expects channels first (C, Z, Y, X) with channel_axis=0, z_axis=1.
        mask_t, *_ = model.eval(
            img_t,
            diameter=diameter,
            channels=channels,
            min_size=min_size,
            do_3D=True,
            normalize=True,
            channel_axis=0,
            z_axis=1,
            anisotropy=anisotropy,
        )
        mask_t=segment_2d_filter(mask_t)

        masks[t] = mask_t.astype(np.uint16)

        if start_ts is not None and end_ts is not None:
            end_ts.record()
            torch.cuda.synchronize()
            dur = start_ts.elapsed_time(end_ts) / 1000  # s
            if verbose:
                print(f" done in {dur:.2f} s.")
        elif verbose:
            print(" done.")

    # 3. Save masks
    masks_path = Path()
    if save_masks:
        if masks_outpath is not None:
            # Use provided path directly
            masks_path = Path(masks_outpath)
        else:
            # Use default naming convention
            masks_path = Path(raw_image_path).with_name(Path(raw_image_path).stem + masks_suffix).with_suffix(".zarr")
        
        save_as_zarr(masks, masks_path)  # type: ignore[arg-type]
        if verbose:
            print("Saved masks to:", masks_path)

    return masks, masks_path


def run_cellpose_segmentation(
    output_dir: str | Path,
    metadata,
    pretrained_model_dir: str | Path,
    input_channels: list,
    label_name: str = "tcell", 
    timepoint_range: Optional[Tuple[int, int]] = None,
    channel_labels_config: Optional[dict] = None,
    progress_cb=None,
    **cellpose_kwargs,
):
    """Batch run Cellpose segmentation for all samples in *metadata*.

    The function mirrors the API of ``run_pixel_classifier_segmentation`` so it
    can be dropped into existing notebooks with minimal changes.

    Parameters
    ----------
    output_dir:
        Root output directory (same one you pass to other behav3d helpers).
    metadata:
        Pandas ``DataFrame`` returned by ``behav3d.utils.load_behav3d_metadata``.
        The frame **must** contain the columns ``sample_name`` and
        ``raw_image_path``.
    pretrained_model_dir:
        Folder that contains the pretrained Cellpose weights.
    input_channels:
        List of channel labels to use for Cellpose evaluation (e.g., ['organoid1', 'tcell1']).
        The order determines which label maps to which model channel.
    label_name:
        Name of the cell type to segment (e.g., "tcell", "organoid", "macro").
        The category (organoid/immune/other) is auto-detected from metadata
        to determine the correct prefix and min_size.
    timepoint_range:
        Optional ``(start_t, end_t)`` tuple - if given the returned *mask*
        Zarrs are cropped to that time-window. ``None`` means all time points.
    channel_labels_config:
        Optional dict with channel label configuration. Should have:
        - 'labels_mode': 'same_for_all' or 'per_sample'
        - 'channel_labels': dict mapping channel index to label (for same_for_all)
        - 'per_sample_channel_labels': dict mapping sample_name -> {channel_idx: label}
        If None, the function will try to read from behav3d_parameters.yml in output_dir.
    **cellpose_kwargs:
        Extra keyword arguments passed straight to
        :pyfunc:`run_cellpose_prediction` (e.g. ``anisotropy=...``).

    Returns
    -------
    metadata:
        The input ``DataFrame`` with an additional (or overwritten) column
        ``{prefix}_{label_name}_segments_image_path`` that points at the generated segments Zarr.
        Where prefix is 'or' for organoids, 'im' for immune cells, 'ot' for other cells.
    """
    
    assert isinstance(metadata, pd.DataFrame), "metadata must be a pandas DataFrame"

    output_dir = Path(output_dir)
    
    # Load channel labels config if not provided
    if channel_labels_config is None:
        config_path = output_dir / "behav3d_parameters.yml"
        if config_path.exists():
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f) or {}
            channel_labels_config = config.get("cellpose", {})
        else:
            channel_labels_config = {}
    
    labels_mode = channel_labels_config.get("labels_mode", "same_for_all")
    global_channel_labels = channel_labels_config.get("channel_labels", {})
    per_sample_channel_labels = channel_labels_config.get("per_sample_channel_labels", {})
    
    # Auto-detect cell type categories from metadata
    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)
    
    # Determine prefix and min_size based on which category label_name belongs to
    if label_name in organoid_types:
        prefix = 'or'
        min_size = 1000  # Large objects (organoids)
        category_name = "organoid"
    elif label_name in immune_types:
        prefix = 'im'
        min_size = 10    # Small objects (immune cells)
        category_name = "immune"
    elif label_name in other_types:
        prefix = 'ot'
        min_size = 10    # Small objects (other cells)
        category_name = "other"
    else:
        raise ValueError(
            f"Cell type '{label_name}' not found in metadata. "
            f"Available types - Organoids: {organoid_types}, Immune: {immune_types}, Other: {other_types}"
        )
    
    print(f"Detected '{label_name}' as {category_name} type (prefix: {prefix}_, min_size: {min_size})")

    summary = {"processed": [], "skipped": []}
    
    # Ensure the path column exists with correct dtype
    path_col = f'{prefix}_{label_name}_segments_image_path'
    if path_col not in metadata.columns:
        metadata[path_col] = pd.NA
    metadata[path_col] = metadata[path_col].astype('object')

    _total_samples = len(metadata)
    for _i, (idx, sample) in enumerate(metadata.iterrows()):
        sample_name = sample['sample_name']
        if progress_cb is not None:
            try:
                progress_cb(_i, _total_samples, f"cellpose {label_name} / {sample_name}")
            except Exception:
                pass

        print(f"Segmenting '{label_name}' ({category_name}) with Cellpose for: {sample_name}")
        cellpose_dir = output_dir / "images" / f"{sample_name}"
        cellpose_dir.mkdir(parents=True, exist_ok=True)
        raw_image_path = Path(sample['raw_image_path'])
        raw_image_zarr = Path(output_dir, "images", sample_name, f"{sample_name}.zarr")
        
        if not raw_image_zarr.exists():
            print(f"  [SKIP] Zarr file missing for {sample_name}, run image preprocessing first.")
            summary["skipped"].append(sample_name)
            continue

        images = load_image(raw_image_zarr)
        

        masks_outpath = cellpose_dir / f"{sample_name}_{label_name}_segments.zarr"

        # Always remove existing file if present (running cell = want new segmentation)
        if masks_outpath.exists():
            print(f"  [OVERWRITE] Removing existing: {masks_outpath}")
            shutil.rmtree(masks_outpath)

        # Build label-to-channel mapping for this sample
        label_to_channel = {}
        
        if labels_mode == "per_sample" and sample_name in per_sample_channel_labels:
            # Per-sample mode: use sample-specific channel labels from config
            sample_channel_labels = per_sample_channel_labels[sample_name]
            label_to_channel = _label_to_channel_from_stored_map(sample_channel_labels)
            print(f"  Using per-sample channel config for '{sample_name}'")
        elif labels_mode == "per_sample" and sample_name not in per_sample_channel_labels:
            # Per-sample mode but this sample is missing
            raise ValueError(
                f"Sample '{sample_name}' not found in per_sample_channel_labels config. "
                f"Available samples: {list(per_sample_channel_labels.keys())}. "
                f"Please configure channel labels for this sample in the 'Configure Channel Labels' panel."
            )
        elif global_channel_labels:
            # Same for all mode: use global channel labels from config
            label_to_channel = _label_to_channel_from_stored_map(global_channel_labels)
            print(f"  Using global channel config (same for all samples)")
        else:
            # No channel labels configured at all
            raise ValueError(
                f"No channel labels configured. "
                f"Please configure channel labels in the 'Configure Channel Labels' panel before running segmentation."
            )
        
        # Extract the channel indices based on input_channels (labels user selected)
        channels = []
        for i, label in enumerate(input_channels):
            if label in label_to_channel:
                # Use the channel index from the config
                channel_idx = label_to_channel[label]
                channels.append(channel_idx)
                print(f"    Model channel {i} <- Image channel {channel_idx} (label: {label})")
            else:
                raise ValueError(
                    f"Channel label '{label}' not found for sample '{sample_name}'. "
                    f"Available labels: {list(label_to_channel.keys())}. "
                    f"Make sure this label is configured in the 'Configure Channel Labels' panel."
                )

        # Extract the anisotropy
        anisotropy = sample['pixel_distance_z'] / sample['pixel_distance_xy']

        # Run Cellpose with explicit output path for consistent naming
        masks, _ = run_cellpose_prediction(
            image=images[:, channels, ...],
            pretrained_model_dir=pretrained_model_dir,
            raw_image_path=raw_image_path,
            save_masks=True,
            masks_outpath=masks_outpath, 
            verbose=True,
            timepoint_range=timepoint_range,
            channels=None,
            anisotropy=anisotropy,
            min_size=min_size,
            nchan=len(channels),
            **cellpose_kwargs,
        )
        print(f"  [SAVED] Segments at {masks_outpath}")

        metadata.at[idx, path_col] = str(masks_outpath)
        summary["processed"].append(sample_name)

    if progress_cb is not None:
        try:
            progress_cb(_total_samples, _total_samples, f"cellpose {label_name} done")
        except Exception:
            pass
    return metadata, summary


def run_cellpose_and_sync_metadata(
    output_dir,
    metadata_loader,
    pretrained_model_dir,
    input_channels,
    label_name,
    timepoint_range=None,
    progress_cb=None,
    **cellpose_kwargs,
):
    """
    Wrapper around run_cellpose_segmentation that automatically handles
    metadata CSV reading, updating, and saving.
    """
    # 1. Run the segmentation
    updated_metadata, summary = run_cellpose_segmentation(
        output_dir=output_dir,
        metadata=metadata_loader.metadata,
        pretrained_model_dir=pretrained_model_dir,
        input_channels=input_channels,
        label_name=label_name,
        timepoint_range=timepoint_range,
        progress_cb=progress_cb,
        **cellpose_kwargs,
    )

    # 2. Update the in-memory metadata
    metadata_loader.metadata = updated_metadata

    # 3. Save
    csv_path = resolve_metadata_csv_path(metadata_loader)
    if not csv_path:
        raise RuntimeError(
            "Cellpose finished but could not find a metadata CSV path to save to "
            "(no metadata_csv_path/_loaded_csv_path, and no paths.metadata_csv in "
            "behav3d_parameters). The updated metadata is still in memory - save it "
            "manually to avoid losing this run's results."
        )
    updated_metadata.to_csv(csv_path, index=False)

    return updated_metadata, summary


def run_otsu_threshold_segmentation_from_zarr(
    output_dir: str | Path,
    metadata,
    mask_suffix: str = "_mask_dead",
    timepoint_range: tuple[int, int] | None = None,
    progress_cb=None,
):
    """Run Otsu thresholding segmentation on the dead channel for all samples.
    
    Parameters
    ----------
    output_dir:
        Root output directory.
    metadata:
        Pandas DataFrame with sample information.
    mask_suffix:
        Suffix for the output mask file (default: "_mask_dead").
    timepoint_range:
        Optional (start_t, end_t) tuple for time cropping.
    
    Returns
    -------
    metadata:
        Updated DataFrame with 'dead_mask_path' column populated.
    """
    assert isinstance(metadata, pd.DataFrame), "metadata must be a pandas DataFrame"
    
    # Check if dead channel exists in metadata
    if not has_dead_channel(metadata):
        print("[SKIP] No dead channel detected in metadata. Skipping dead mask segmentation.")
        return metadata, {"processed": [], "skipped": []}
    
    output_dir = Path(output_dir)
    
    # Load channel labels config from behav3d_parameters.yml
    config_path = output_dir / "behav3d_parameters.yml"
    channel_labels_config = {}
    if config_path.exists():
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f) or {}
        channel_labels_config = config.get("cellpose", {})
    
    labels_mode = channel_labels_config.get("labels_mode", "same_for_all")
    global_channel_labels = channel_labels_config.get("channel_labels", {})
    per_sample_channel_labels = channel_labels_config.get("per_sample_channel_labels", {})
    
    # Ensure dead_mask_path column exists with correct dtype
    if 'dead_mask_path' not in metadata.columns:
        metadata['dead_mask_path'] = pd.NA
    metadata['dead_mask_path'] = metadata['dead_mask_path'].astype('object')

    summary = {"processed": [], "skipped": []}

    _total_samples = len(metadata)
    for _i, (idx, sample) in enumerate(metadata.iterrows()):
        sample_name = sample['sample_name']
        if progress_cb is not None:
            try:
                progress_cb(_i, _total_samples, f"dead mask / {sample_name}")
            except Exception:
                pass
        print(f"\nProcessing dead mask for {sample_name}...")
        
        # Find the death channel from config file
        death_channel = None
        
        if labels_mode == "per_sample" and sample_name in per_sample_channel_labels:
            # Per-sample mode: look for 'dead' label in this sample's config
            sample_channel_labels = per_sample_channel_labels[sample_name]
            for ch_idx, label in sample_channel_labels.items():
                if label == 'dead':
                    death_channel = int(ch_idx)
                    print(f"  Using per-sample channel config: dead channel = {death_channel}")
                    break
        elif global_channel_labels:
            # Same for all mode: look for 'dead' label in global config
            for ch_idx, label in global_channel_labels.items():
                if label == 'dead':
                    death_channel = int(ch_idx)
                    print(f"  Using global channel config: dead channel = {death_channel}")
                    break
        
        if death_channel is None:
            print(f"[SKIP] No 'dead' channel label found in config for {sample_name}")
            print(f"       Please configure 'dead' label in the 'Configure Channel Labels' panel.")
            summary["skipped"].append(sample_name)
            continue
        
        mask_dir = output_dir / "images" / sample_name
        mask_dir.mkdir(parents=True, exist_ok=True)

        raw_image_zarr = mask_dir / f"{sample_name}.zarr"
        if not raw_image_zarr.exists():
            print(f"[SKIP] Zarr file missing for {sample_name}, run image preprocessing first.")
            summary["skipped"].append(sample_name)
            continue

        masks_outpath = mask_dir / f"{sample_name}{mask_suffix}.zarr"
        
        # Always remove existing file if present (running cell = want new segmentation)
        if masks_outpath.exists():
            print(f"[OVERWRITE] Removing existing: {masks_outpath}")
            shutil.rmtree(masks_outpath)

        # Load image in (T, C, Z, Y, X)
        images = load_image(raw_image_zarr)
        T_total = images.shape[0]
        spatial_shape = images.shape[2:]

        # Select time indices for the loop. For dask, the actual array slice must use slice()
        # or NumPy-style :, not a bare Python range — dask treats range as one fancy-index
        # axis and fails inside normalize_index (TypeError: range vs int).
        if timepoint_range is not None:
            start_t, end_t = timepoint_range
            start_t = max(0, int(start_t))
            end_t = min(T_total - 1, int(end_t))
            indices_to_process = list(range(start_t, end_t + 1))
            print(f"  Otsu global indexing: processing T {start_t} to {end_t} (Total: {T_total})")
            t_slice = slice(start_t, end_t + 1)
        else:
            indices_to_process = list(range(T_total))
            t_slice = slice(None)  # all T — same as images[:, death_channel]

        # Initialize full-length masks array
        masks = np.zeros((T_total,) + spatial_shape, dtype=np.uint8)

        # Apply thresholding
        # To be efficient and follow user 1-point-at-a-time pattern, we'll process timepoints
        # but compute a global threshold if possible or use the requested range.
        
        # Select single channel for relevant timepoints: (T_subset, Z, Y, X)
        channel_img_subset = images[t_slice, death_channel]

        # Flatten to compute threshold on the requested range
        flat_vals = channel_img_subset.ravel()
        global_thresh = threshold_otsu(flat_vals.compute())
        print(f"  [INFO] Otsu threshold (computed on requested range): {global_thresh}")

        # Apply threshold and write to the correct absolute indices in 'masks'
        # Dask array comparison is lazy, compute() it for use.
        for i, t in enumerate(indices_to_process):
            masks[t] = (channel_img_subset[i] > global_thresh).compute().astype(np.uint8)

        save_as_zarr(masks, masks_outpath)
        print(f"[SAVED] Dead mask at {masks_outpath}")
        
        # Update metadata with the path
        metadata.at[idx, 'dead_mask_path'] = str(masks_outpath)
        summary["processed"].append(sample_name)

    if progress_cb is not None:
        try:
            progress_cb(_total_samples, _total_samples, "dead mask done")
        except Exception:
            pass
    return metadata, summary

def run_otsu_and_sync_metadata(
    output_dir,
    metadata_loader,
    mask_suffix="_mask_dead",
    timepoint_range=None,
    progress_cb=None,
):
    """
    Wrapper around run_otsu_threshold_segmentation_from_zarr that automatically handles
    metadata CSV reading, updating, and saving.
    """
    # 1. Run the segmentation (The Worker)
    updated_metadata, summary = run_otsu_threshold_segmentation_from_zarr(
        output_dir=output_dir,
        metadata=metadata_loader.metadata,
        mask_suffix=mask_suffix,
        timepoint_range=timepoint_range,
        progress_cb=progress_cb,
    )

    # 2. Update the in-memory metadata
    metadata_loader.metadata = updated_metadata

    # 3. Save
    csv_path = resolve_metadata_csv_path(metadata_loader)
    if not csv_path:
        raise RuntimeError(
            "Otsu thresholding finished but could not find a metadata CSV path to save "
            "to (no metadata_csv_path/_loaded_csv_path, and no paths.metadata_csv in "
            "behav3d_parameters). The updated metadata is still in memory - save it "
            "manually to avoid losing this run's results."
        )
    updated_metadata.to_csv(csv_path, index=False)

    return updated_metadata, summary