"""Cellpose prediction utilities.

This module provides a single convenience function `run_cellpose_prediction` that
wraps the ad-hoc logic from the original `cellpose_TKI_tcell_run.py` demo script
into a re-usable function that can be imported and called from notebooks (e.g.
`run_behav3d.ipynb`).

The function
1. Loads the image (OME-Tiff, regular Tiff, or Zarr) using `behav3d.utils.fileio`.
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
import napari
from skimage.filters import threshold_otsu

from behav3d.io.images import save_as_zarr, load_zarr, load_image, append_to_zarr
from skimage.filters import threshold_otsu
from behav3d.preprocessing.segmentation import segment_2d_filter
from behav3d.core.metadata import (
    detect_organoid_types_from_metadata,
    detect_immune_cell_types_from_metadata,
    detect_other_cell_types_from_metadata,
    has_dead_channel,
)

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
    
    # Optional time sub-setting
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, image.shape[0] - 1))
        end_t = max(start_t, min(end_t, image.shape[0] - 1))
        image = image[start_t : end_t + 1]
        if verbose:
            print(f"Using timepoints {start_t}-{end_t} (total {image.shape[0]})")
            
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
    masks = np.zeros((T, Z, Y, X), dtype=np.uint16)

    if verbose:
        print("Starting Cellpose prediction...")

    for t in range(T):
        # Extract data for one time point - shape (C, Z, Y, X).
        img_t = image[t]

        if verbose:
            print(f"  - Timepoint {t+1}/{T} -", end="", flush=True)
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

        if end_ts is not None:
            end_ts.record()
            torch.cuda.synchronize()
            dur = start_ts.elapsed_time(end_ts) / 1000  # ms ? s
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

    for idx, sample in metadata.iterrows():
        sample_name = sample['sample_name']

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
            label_to_channel = {str(label): int(ch_idx) for ch_idx, label in sample_channel_labels.items()}
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
            label_to_channel = {str(label): int(ch_idx) for ch_idx, label in global_channel_labels.items()}
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

    return metadata, summary


def run_cellpose_and_sync_metadata(
    output_dir,
    metadata_loader,
    pretrained_model_dir,
    input_channels,
    label_name,
    timepoint_range=None,
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
        **cellpose_kwargs,
    )

    # 2. Update the in-memory metadata
    metadata_loader.metadata = updated_metadata

    # 3. Save
    updated_metadata.to_csv(metadata_loader.metadata_csv_path, index=False)
    
    return updated_metadata, summary


def run_otsu_threshold_segmentation_from_zarr(
    output_dir: str | Path,
    metadata,
    mask_suffix: str = "_mask_dead",
    timepoint_range: tuple[int, int] | None = None,
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

    for idx, sample in metadata.iterrows():
        sample_name = sample['sample_name']
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

        # Time cropping
        if timepoint_range is not None:
            start_t, end_t = timepoint_range
            images = images[start_t:end_t + 1]

        # Select single channel: (T, Z, Y, X)
        channel_img = images[:, death_channel]

        # Flatten to compute global threshold
        flat_vals = channel_img.ravel()
        global_thresh = threshold_otsu(flat_vals.compute())
        print(f"[INFO] Global Otsu threshold: {global_thresh}")

        # Apply threshold to the entire 4D array (T, Z, Y, X)
        masks = (channel_img > global_thresh).astype(np.uint8)

        save_as_zarr(masks, masks_outpath)
        print(f"[SAVED] Dead mask at {masks_outpath}")
        
        # Update metadata with the path
        metadata.at[idx, 'dead_mask_path'] = str(masks_outpath)
        summary["processed"].append(sample_name)
    
    return metadata, summary

def run_otsu_and_sync_metadata(
    output_dir,
    metadata_loader,
    mask_suffix="_mask_dead",
    timepoint_range=None,
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
    )

    # 2. Update the in-memory metadata
    metadata_loader.metadata = updated_metadata

    # 3. Save
    updated_metadata.to_csv(metadata_loader.metadata_csv_path, index=False)
    
    return updated_metadata, summary

# -----------------------------------------------------------------------------
# 3. Quick visualization helper (Napari)
# -----------------------------------------------------------------------------
def visualize_cellpose_sample(
    output_dir: str | Path,
    sample_name: str,
    metadata=None,
    timepoint_range: Optional[Tuple[int, int]] = None,
    channel_colors: Sequence[str] = (
    "cyan", "yellow", "red", "green", "magenta", "blue",
    "gray", "turbo", "viridis", "plasma", "inferno", "twilight"),
) -> None:
    """Visualize a sample with all detected cell type segments in Napari.
    
    Parameters
    ----------
    output_dir:
        Root output directory.
    sample_name:
        Name of the sample to visualize.
    metadata:
        Optional metadata DataFrame. If provided, cell types are auto-detected
        and all available segment masks are loaded dynamically.
    timepoint_range:
        Optional (start_t, end_t) tuple for time slicing.
    channel_colors:
        Colors for image channels.
    """
    print(f"Sample selected: {sample_name}")
    print(f"Timepoint range: {timepoint_range}")

    output_dir = Path(output_dir)
    raw_image_zarr = Path(output_dir, "images", sample_name, f"{sample_name}.zarr")

    # Build mask paths dynamically based on metadata or fall back to common names
    mask_files = {}
    
    if metadata is not None:
        # Auto-detect cell types from metadata
        organoid_types = detect_organoid_types_from_metadata(metadata)
        immune_types = detect_immune_cell_types_from_metadata(metadata)
        other_types = detect_other_cell_types_from_metadata(metadata)
        has_death = has_dead_channel(metadata)
        
        print(f"Detected cell types - Organoids: {organoid_types}, Immune: {immune_types}, Other: {other_types}")
        
        # Get the sample row from metadata
        sample_row = metadata[metadata['sample_name'] == sample_name]
        if sample_row.empty:
            print(f"[Warning] Sample '{sample_name}' not found in metadata")
        else:
            sample_row = sample_row.iloc[0]
            
            # Add segment paths if they exist in metadata
            for cell_type in organoid_types:
                path_col = f'or_{cell_type}_segments_image_path'
                if path_col in metadata.columns:
                    path_str = sample_row.get(path_col)
                    if path_str and not pd.isna(path_str):
                        mask_files[f"{cell_type}_segments"] = Path(path_str)
                    # Skip if path is empty, segmentation not done yet
            
            for cell_type in immune_types:
                path_col = f'im_{cell_type}_segments_image_path'
                if path_col in metadata.columns:
                    path_str = sample_row.get(path_col)
                    if path_str and not pd.isna(path_str):
                        mask_files[f"{cell_type}_segments"] = Path(path_str)
                    # Skip if path is empty, segmentation not done yet
            
            for cell_type in other_types:
                path_col = f'ot_{cell_type}_segments_image_path'
                if path_col in metadata.columns:
                    path_str = sample_row.get(path_col)
                    if path_str and not pd.isna(path_str):
                        mask_files[f"{cell_type}_segments"] = Path(path_str)
                    # Skip if path is empty, segmentation not done yet
            
            # Add dead mask from metadata
            if has_death and 'dead_mask_path' in metadata.columns:
                dead_path_str = sample_row.get('dead_mask_path')
                if dead_path_str and not pd.isna(dead_path_str):
                    mask_files["mask_dead"] = Path(dead_path_str)
                # Skip if path is empty, segmentation not done yet
    else:
        # Fallback: scan directory for *_segments.zarr files
        sample_dir = output_dir / "images" / sample_name
        if sample_dir.exists():
            for zarr_path in sample_dir.glob("*_segments.zarr"):
                # Extract cell type name from filename: {sample}_{celltype}_segments.zarr
                name = zarr_path.stem.replace(f"{sample_name}_", "").replace("_segments", "")
                mask_files[f"{name}_segments"] = zarr_path
            
            # Check for dead mask
            dead_mask_path = sample_dir / f"{sample_name}_mask_dead.zarr"
            if dead_mask_path.exists():
                mask_files["mask_dead"] = dead_mask_path

    for name, path in mask_files.items():
        print(f"Mask '{name}' path: {path}")

    # Load masks that exist
    loaded_masks = {}
    for name, path in mask_files.items():
        if path.exists():
            mask_arr = load_zarr(path)
            loaded_masks[name] = mask_arr
        else:
            print(f"[Warning] Mask '{name}' file not found: {path}")

    # Load raw image
    img_arr = load_zarr(raw_image_zarr)
    print(f"Original image shape: {img_arr.shape}")

    # Slice timepoints
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, img_arr.shape[0] - 1))
        end_t = max(start_t, min(end_t, img_arr.shape[0] - 1))
        print(f"Slicing timepoints from {start_t} to {end_t}")
        img_arr = img_arr[start_t:end_t + 1]
        for name in loaded_masks:
            loaded_masks[name] = loaded_masks[name][start_t:end_t + 1]

    # Convert to NumPy
    img_np = np.asarray(img_arr)

    for name in loaded_masks:
        loaded_masks[name] = np.asarray(loaded_masks[name])

    # Launch viewer
    viewer = napari.Viewer()

    if img_np.ndim != 5:
        raise ValueError(f"Expected image shape (T, C, Z, Y, X), got {img_np.shape}")

    T, C, Z, Y, X = img_np.shape
    for ch in range(C):
        channel_data = img_np[:, ch]  # (T, Z, Y, X)
        viewer.add_image(
            channel_data,
            name=f"channel_{ch}",
            colormap=channel_colors[ch] if ch < len(channel_colors) else 'gray',
            scale=(1, 1, 1, 1),  # T, Z, Y, X
            blending="additive",
            channel_axis=None,
        )

    # Add all masks as separate label layers
    if loaded_masks:
        for name, mask_np in loaded_masks.items():
            print(f"Adding mask layer '{name}' to viewer...")
            viewer.add_labels(mask_np, name=name, scale=(1, 1, 1, 1))
    else:
        print("No mask layers added.")

    print("Launching Napari viewer...")
    napari.run()