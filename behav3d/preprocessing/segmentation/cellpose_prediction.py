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

from behav3d.utils.fileio import save_as_zarr, load_zarr, load_image, append_to_zarr

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
    min_size: int = 10,
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
        ``tcell_segments_image_path`` that points at the generated mask Zarr.
    """
    import pandas as pd  # local import to avoid heavy dependency at module import

    assert isinstance(metadata, pd.DataFrame), "metadata must be a pandas DataFrame"

    output_dir = Path(output_dir)
    cellpose_dir = output_dir / "images" / "CellPose"
    cellpose_dir.mkdir(parents=True, exist_ok=True)

    for idx, sample in metadata.iterrows():
        sample_name = sample['sample_name']

        print(f"Calculating features for: {sample_name}")

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

        masks_outpath = cellpose_dir / f"{sample_name}_tcell_segments.zarr"

        if masks_outpath.exists() and not overwrite:
            # Skip computation; just record the path.
            metadata.at[idx, "tcell_segments_image_path"] = str(masks_outpath)
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

        metadata.at[idx, "tcell_segments_image_path"] = str(masks_outpath)

    return metadata


# -----------------------------------------------------------------------------
# 3. Quick visualization helper (Napari)
# -----------------------------------------------------------------------------

def visualize_cellpose_sample(
    metadata,
    sample_name: str,
    *,
    timepoint_range: Optional[Tuple[int, int]] = None,
    channel_colors: Sequence[str] = ("red", "green", "blue", "cyan", "magenta", "yellow"),
) -> None:
    """Open a Napari viewer showing raw channels + Cellpose masks.

    Parameters
    ----------
    image_path:
        Path to the raw image file.
    masks:
        Either a path to the mask *zarr* (or tiff) produced by
        :pyfunc:`run_cellpose_prediction` **or** an in-memory :class:`numpy.ndarray`.
        Expected shape ``(T,Z,Y,X)``.
    manual_dim_order:
        Same meaning as in :pyfunc:`run_cellpose_prediction` – supply when the
        file is stored in a non-standard axis order so it can be brought to
        default TCZYX.
    timepoint_range:
        Optional ``(start_t, end_t)`` tuple to restrict the viewer to a subset of
        time points.
    channel_colors:
        Colors for the first six channels; repeat or extend as needed.
    """
    # Lazy import so that the module can be used in headless environments.
    import napari  # noqa: WPS433

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    import pandas as pd  # local import

    assert isinstance(metadata, pd.DataFrame), "metadata must be a DataFrame"
    match = metadata[metadata["sample_name"] == sample_name]
    if match.empty:
        raise ValueError(f"Sample '{sample_name}' not found in metadata.")
    if len(match) > 1:
        raise ValueError(f"Multiple rows found for sample '{sample_name}'. Provide unique sample_name.")
    raw_image_col = "raw_image_path"
    masks_col = "tcell_segments_image_path"
    row = match.iloc[0]
    raw_image_path = Path(row[raw_image_col])
    mask_path      = Path(row[masks_col]) if not pd.isna(row[masks_col]) else None

    if mask_path is None or not mask_path.exists():
        raise ValueError(f"Mask file '{mask_path}' does not exist for sample '{sample_name}'.")

    # ------------------------------------------------------------------
    # Slice and load ONLY requested timepoints (lazy for Zarr)
    # ------------------------------------------------------------------
    # Load raw image and masks lazily (dask arrays) so we can slice before
    # materialising anything into memory. We skip `load_image` here because
    # axis-order is already TCZYX inside our own Zarr files.
    from behav3d.utils.fileio import load_zarr

    # Use Zarr loader only if the file is a Zarr store, otherwise fall back to load_image
    from behav3d.utils.fileio import load_zarr as _lz, load_image as _li

    if raw_image_path.suffix == ".zarr" or str(raw_image_path).endswith(".zarr.zip"):
        img_arr = _lz(raw_image_path)
    else:
        img_arr = _li(raw_image_path)   # may return dask array depending on backend

    mask_arr = _lz(mask_path)

    # ------------------------------------------------------------------
    # Slice ONLY requested timepoints (lazy) before computing to NumPy.
    # ------------------------------------------------------------------
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, img_arr.shape[0]-1))
        end_t   = max(start_t, min(end_t,   img_arr.shape[0]-1))
        img_arr  = img_arr[start_t:end_t+1]
        mask_arr = mask_arr[start_t:end_t+1]

    # Bring slices into memory for Napari
    img_np   = np.asarray(img_arr)
    masks_np = np.asarray(mask_arr)

    # ------------------------------------------------------------------
    # Optional cropping in time
    # ------------------------------------------------------------------
    if timepoint_range is not None:
        start_t, end_t = timepoint_range
        start_t = max(0, min(start_t, img_np.shape[0] - 1))
        end_t = max(start_t, min(end_t, img_np.shape[0] - 1))
        img_np = img_np[start_t : end_t + 1]
        masks_np = masks_np[start_t : end_t + 1]

    # ------------------------------------------------------------------
    # Launch viewer
    # ------------------------------------------------------------------
    viewer = napari.Viewer()

    T, C, Z, Y, X = img_np.shape
    for ch in range(C):
        channel_data = img_np[:, ch]  # shape (T,Z,Y,X)
        # Dynamic contrast per channel
        perc_99 = float(np.percentile(channel_data, 99))
        viewer.add_image(
            channel_data,
            name=f"channel_{ch+1}",
            contrast_limits=(0, perc_99),
            colormap=channel_colors[ch % len(channel_colors)],
            scale=(1, 1, 1, 1),  # T, Z, Y, X
            blending="additive",
        )

    viewer.add_labels(masks_np, name="cellpose_masks", scale=(1, 1, 1, 1))

    napari.run()



