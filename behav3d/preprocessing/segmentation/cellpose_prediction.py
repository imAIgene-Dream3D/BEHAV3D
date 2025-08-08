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

from pathlib import Path
from typing import Optional, Sequence, Tuple

import numpy as np
import torch
import pandas as pd

from behav3d.utils.fileio import save_as_zarr, load_zarr, load_image, append_to_zarr
from skimage.filters import threshold_otsu
from behav3d.utils.segmentation import segment_2d_filter

# Cellpose import is relatively expensive – only load when we actually need it.
try:
    from cellpose import models  # noqa: WPS433 (allow external import)
except ImportError as err:  # pragma: no cover – handled at runtime
    raise ImportError(
        "cellpose is required for `behav3d.preprocessing.segmentation.cellpose_prediction`. "
        "Please install it via `pip install cellpose>=4.0` (or compatible).",
    ) from err



def run_cellpose_prediction(  # noqa: WPS231 (complexity) – unavoidable for pipeline function
    image: np.ndarray,  #already loaded image
    pretrained_model_dir: str | Path,
    *,
    raw_image_path: str | Path, 
    nchan: int = 3,
    channels: Optional[Sequence[int]] = None,
    diameter: Optional[int | float] = None,
    min_size: int = 25,
    anisotropy: Optional[float] = 4.655 / 1.17365196872,
    manual_dim_order: Optional[str | Sequence[str]] = None,
    timepoint_range: Optional[Tuple[int, int]] = None,
    device: Optional[str] = None,
    convert_tiff_to_zarr: bool = True,
    save_masks: bool = True,
    masks_suffix: str = "_cellpose_masks.zarr",
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
    nchan, channels, diameter, anisotropy, min_size:
        Passed through to ``CellposeModel.eval``.
    device:
        CUDA device string (e.g. ``"cuda:0"``); when ``None`` the function will
        automatically pick GPU if available otherwise CPU.
    convert_tiff_to_zarr:
        If *image_path* is a Tiff and this flag is ``True`` the image is stored
        as Zarr next to the original file for future reuse.
    save_masks:
        When ``True`` the resulting mask array is saved as Zarr (path is the
        original image stem + *masks_suffix*).
    masks_suffix:
        Suffix to use when writing the masks to disk (ignored when
        *save_masks* is ``False``).
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
    image_np = image

    if verbose:
        print("Loaded image with shape (T, C, Z, Y, X):", image_np.shape)

    # Optional time sub-setting
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, image_np.shape[0] - 1))
        end_t = max(start_t, min(end_t, image_np.shape[0] - 1))
        image_np = image_np[start_t : end_t + 1]
        if verbose:
            print(f"Using timepoints {start_t}–{end_t} (total {image_np.shape[0]})")

    # Determine computation device.
    torch_device = (
        torch.device(device) if device is not None else torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    )

    # Initialise Cellpose model.
    model = models.CellposeModel(
        pretrained_model=str(pretrained_model_dir),
        gpu=torch_device.type == "cuda",
        nchan=nchan,
        device=torch_device,
    )

    # ---------------------------------------------------------------------
    # 2. Run prediction per time point
    # ---------------------------------------------------------------------
    T, C, Z, Y, X = image_np.shape
    masks = np.zeros((T, Z, Y, X), dtype=np.uint16)

    if verbose:
        print("Starting Cellpose prediction...")

    for t in range(T):
        # Extract data for one time point – shape (C, Z, Y, X).
        img_t = image_np[t]

        if verbose:
            print(f"  • Time-point {t+1}/{T} …", end="", flush=True)
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
        mask_t = segment_2d_filter(mask_t)
        masks[t] = mask_t.astype(np.uint16)
        
        if end_ts is not None:
            end_ts.record()
            torch.cuda.synchronize()
            dur = start_ts.elapsed_time(end_ts) / 1000  # ms → s
            if verbose:
                print(f" done in {dur:.2f} s.")
        elif verbose:
            print(" done.")

    # ---------------------------------------------------------------------
    # 3. Save masks (optional)
    # ---------------------------------------------------------------------
    masks_path = Path()
    if save_masks:
        masks_path = Path(raw_image_path).with_name(Path(raw_image_path).stem + masks_suffix).with_suffix(".zarr")
        save_as_zarr(masks, masks_path)  # type: ignore[arg-type]
        if verbose:
            print("Saved masks to:", masks_path)

        return masks, masks_path


def run_cellpose_segmentation(
    output_dir: str | Path,
    metadata,
    pretrained_model_dir: str | Path,
    *,
    label_name: str = "tcell", 
    manual_dim_order=None, 
    timepoint_range: Optional[Tuple[int, int]] = None,
    overwrite: bool = False,
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
    timepoint_range:
        Optional ``(start_t, end_t)`` tuple – if given the returned *mask*
        Zarrs are cropped to that time-window. ``None`` means all time points.
    overwrite:
        If ``True`` existing mask Zarr files are recomputed.
    **cellpose_kwargs:
        Extra keyword arguments passed straight to
        :pyfunc:`run_cellpose_prediction` (e.g. ``anisotropy=...``).

    Returns
    -------
    metadata:
        The input ``DataFrame`` with an additional (or overwritten) column
        ``{label_name}_segments_image_path`` that points at the generated mask Zarr.
    """
    import pandas as pd  # local import to avoid heavy dependency at module import

    assert isinstance(metadata, pd.DataFrame), "metadata must be a pandas DataFrame"

    output_dir = Path(output_dir)
    

    for idx, sample in metadata.iterrows():
        sample_name = sample['sample_name']

        print(f"Calculating features for: {sample_name}")
        cellpose_dir = output_dir / "images" / f"{sample_name}"
        cellpose_dir.mkdir(parents=True, exist_ok=True)
        raw_image_path = Path(sample['raw_image_path'])
        raw_image_zarr = Path(output_dir, "images", sample_name, f"{sample_name}.zarr")

        if not raw_image_zarr.exists():
            print(f"- Converting raw image to .zarr for memory efficiency...")
            images = load_image(raw_image_path)
            print(f"Original image shape: {images.shape}")
            default_order = ("T", "C", "Z", "Y", "X")  # T, C, Z, Y, X
            if manual_dim_order is not None:
                # Convert to tuple of characters if string
                if isinstance(manual_dim_order, str):
                    manual_dim_order = tuple(manual_dim_order)
                if manual_dim_order != default_order:
                    # Compute permutation: for each axis in default_order, find its index in manual_dim_order
                    perm = [manual_dim_order.index(ax) for ax in default_order]
                    print(f"Transposing image from {manual_dim_order} to {default_order} using permutation {perm}")
                    images = images.transpose(perm)
                    print(f"New image shape: {images.shape}")
                else:
                    print("Image is already in default order (T, C, Z, Y, X).")
            else:
                print("No manual_dim_order provided, assuming image is already in default order.")
            chunksize = (1,) + images.shape[1:]
            save_as_zarr(
                img=images,
                path=raw_image_zarr,
                chunks=chunksize
            )

        images = load_image(raw_image_zarr)
        max_t = images.shape[0] - 1
        print(images.shape)

        masks_outpath = cellpose_dir / f"{sample_name}_{label_name}_segments.zarr"

        if masks_outpath.exists() and not overwrite:
            # Skip computation; just record the path.
            metadata.at[idx, f"{label_name}_segments_image_path"] = str(masks_outpath)
            continue

        # Run Cellpose – we force saving directly to *masks_outpath* by setting
        # save_masks=False and writing ourselves to ensure consistent naming.
        masks, _ = run_cellpose_prediction(
            image=images,
            pretrained_model_dir=pretrained_model_dir,
            raw_image_path=raw_image_path,
            save_masks=True,
            verbose=True,
            timepoint_range=timepoint_range,
            **cellpose_kwargs,
        )

        # Write to disk
        save_as_zarr(masks, masks_outpath)

        metadata.at[idx, f"{label_name}_segments_image_path"] = str(masks_outpath)

    return metadata



def run_otsu_threshold_segmentation_from_zarr(
    output_dir: str | Path,
    metadata,
    *,
    channel: int = 2,
    mask_suffix: str = "_mask_dead",
    timepoint_range: tuple[int, int] | None = None,
    overwrite: bool = False,
):
    from pathlib import Path
    import numpy as np
    from skimage.filters import threshold_otsu
    from behav3d.utils.fileio import load_image, save_as_zarr

    output_dir = Path(output_dir)

    for idx, sample in metadata.iterrows():
        sample_name = sample['sample_name']
        print(f"\nProcessing {sample_name}...")
        
        mask_dir = output_dir / "images" / sample_name
        mask_dir.mkdir(parents=True, exist_ok=True)

        raw_image_zarr = mask_dir / f"{sample_name}.zarr"
        if not raw_image_zarr.exists():
            print(f"[SKIP] Zarr file missing for {sample_name}, run image preprocessing first.")
            continue

        # Load image in (T, C, Z, Y, X)
        images = load_image(raw_image_zarr)

        # Time cropping
        if timepoint_range is not None:
            start_t, end_t = timepoint_range
            images = images[start_t:end_t + 1]

        # Select single channel: (T, Z, Y, X)
        channel_img = images[:, channel]

        # Flatten to compute global threshold
        flat_vals = channel_img.ravel()
        global_thresh = threshold_otsu(flat_vals)
        print(f"[INFO] Global Otsu threshold: {global_thresh}")

        # Apply threshold to the entire 4D array (T, Z, Y, X)
        masks = (channel_img > global_thresh).astype(np.uint8)

        masks_outpath = mask_dir / f"{sample_name}{mask_suffix}.zarr"
        if masks_outpath.exists() and not overwrite:
            print(f"[SKIP] Mask already exists: {masks_outpath}")
            continue

        save_as_zarr(masks, masks_outpath)
        print(f"[SAVED] 3D mask at {masks_outpath}")

# -----------------------------------------------------------------------------
# 3. Quick visualization helper (Napari)
# -----------------------------------------------------------------------------
def visualize_cellpose_sample(
    output_dir: str | Path,
    sample_name: str,
    *,
    timepoint_range: Optional[Tuple[int, int]] = None,
    channel_colors: Sequence[str] = ("red", "green", "blue", "cyan", "magenta", "yellow"),
) -> None:
    import napari
    import numpy as np
    from pathlib import Path
    from behav3d.utils.fileio import load_zarr as _lz

    print(f"Sample selected: {sample_name}")
    print(f"Timepoint range: {timepoint_range}")

    output_dir = Path(output_dir)
    raw_image_zarr = Path(output_dir, "images", sample_name, f"{sample_name}.zarr")

    # Define all possible masks with descriptive names
    mask_files = {
        "tcell_segments": Path(output_dir, "images", sample_name, f"{sample_name}_tcell_segments.zarr"),
        "organoid_segments": Path(output_dir, "images", sample_name, f"{sample_name}_organoid_segments.zarr"),
        "mask_dead": Path(output_dir, "images", sample_name, f"{sample_name}_mask_dead.zarr"),
    }

    for name, path in mask_files.items():
        print(f"Mask '{name}' path: {path}")

    # Load masks that exist
    loaded_masks = {}
    for name, path in mask_files.items():
        if path.exists():
            mask_arr = _lz(path)
            loaded_masks[name] = mask_arr
        else:
            print(f"[Warning] Mask '{name}' file not found: {path}")

    # Load raw image
    img_arr = _lz(raw_image_zarr)
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
            name=f"channel_{ch+1}",
            colormap=channel_colors[ch % len(channel_colors)],
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