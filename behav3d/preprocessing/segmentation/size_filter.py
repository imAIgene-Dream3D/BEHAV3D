"""
Size-based filtering of labelled segmentation volumes stored as Zarr.

The filtering is performed **per-timepoint** so that the full 5-D image
never needs to be loaded into memory at once.
"""

import numpy as np
import zarr
from pathlib import Path
from scipy import ndimage
from tqdm import tqdm
import sys


def filter_segments_by_size(
    segments_zarr_path,
    mask_zarr_path,
    min_size_voxels,
):
    """Remove connected components smaller than *min_size_voxels* from a
    segmentation zarr and regenerate the corresponding binary mask zarr.

    Parameters
    ----------
    segments_zarr_path : str | Path
        Path to the ``*_segments.zarr`` array (T, Z, Y, X) with uint16 labels.
    mask_zarr_path : str | Path
        Path to the matching ``*_mask.zarr`` array that will be regenerated.
    min_size_voxels : int
        Minimum 3-D volume (in voxels) a segment must have to be kept.
    """
    segments_zarr_path = Path(segments_zarr_path)
    mask_zarr_path = Path(mask_zarr_path)

    if not segments_zarr_path.exists():
        print(f"  ⚠️ Segments file not found: {segments_zarr_path}")
        return

    seg = zarr.open(str(segments_zarr_path), mode="r+")
    n_timepoints = seg.shape[0]

    # Open / create the mask with the same shape
    mask = zarr.open(
        str(mask_zarr_path),
        mode="w",
        shape=seg.shape,
        chunks=seg.chunks,
        dtype="uint16",
    )

    removed_total = 0

    pbar = tqdm(
        range(n_timepoints),
        desc="    Filtering",
        unit="tp",
        file=sys.stdout,
        dynamic_ncols=True,
    )

    for t in pbar:
        vol = np.asarray(seg[t])  # single 3-D volume (Z, Y, X)

        # Label connected components in the 3-D volume
        labelled, n_labels = ndimage.label(vol > 0)

        if n_labels == 0:
            mask[t] = np.zeros_like(vol, dtype=np.uint16)
            continue

        # Count voxels per label (index 0 = background)
        sizes = np.bincount(labelled.ravel())

        # Find labels that are too small (skip background at index 0)
        small_labels = np.where(sizes[1:] < min_size_voxels)[0] + 1
        removed_total += len(small_labels)

        if len(small_labels) > 0:
            # Build a boolean mask of voxels belonging to small labels
            remove_mask = np.isin(labelled, small_labels)
            vol[remove_mask] = 0
            seg[t] = vol

        # Regenerate binary mask
        mask[t] = (vol > 0).astype(np.uint16)

        pbar.set_postfix(removed=removed_total)

    return removed_total
