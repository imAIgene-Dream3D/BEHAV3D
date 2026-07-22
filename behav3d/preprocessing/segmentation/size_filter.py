"""Size filtering for segmentation label images.

All sizes are **native voxel counts** (3-D volumes), matching how they are
persisted in ``behav3d_parameters.yml``. GUI panels may display µm³ via the unit
manager, but convert back to voxels before calling in here.

The per-volume work lives in
:func:`behav3d.preprocessing.segmentation.segment_size_filter` so that the
interactive preview, the Cellpose-SAM run-time filter and the on-disk filter
below all behave identically.
"""
import sys
from pathlib import Path

import numpy as np
import zarr
from tqdm import tqdm

from behav3d.preprocessing.segmentation import segment_size_filter


def compute_label_sizes(volume):
    """Return ``(labels, sizes)`` in voxels for a labelled volume, excluding background.

    Uses ``np.bincount``, which is what makes the interactive size-filter preview
    responsive enough to re-filter on every slider move without re-running
    segmentation.
    """
    volume = np.asarray(volume)
    if not np.any(volume):
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    counts = np.bincount(volume.ravel())
    labels = np.nonzero(counts)[0]
    labels = labels[labels != 0]
    return labels, counts[labels]


def filter_labels_by_size(volume, size_min=None, size_max=None, fill_holes=True):
    """Filter one in-memory volume. Thin alias kept for call-site readability."""
    return segment_size_filter(
        np.asarray(volume).copy(),
        size_min=size_min or None,
        size_max=size_max or None,
        fill_holes=fill_holes,
    )


def _open_segments_array(segments_zarr_path):
    """Open a segments zarr ``r+``, unwrapping a single-array group if needed."""
    seg = zarr.open(str(segments_zarr_path), mode="r+")
    if isinstance(seg, zarr.Group):
        keys = list(seg.array_keys())
        if not keys:
            print(f"  ⚠️ No arrays found in Zarr group: {segments_zarr_path}")
            return None
        print(f"  📂 Zarr is a group. Using internal array: {keys[0]}")
        return seg[keys[0]]
    return seg


def filter_segments_by_size(
    segments_zarr_path,
    min_size_voxels=None,
    fill_holes=True,
    max_size_voxels=None,
):
    """Remove segments outside ``[min_size_voxels, max_size_voxels]`` from a segments zarr.

    Rewrites the array **in place**, one timepoint at a time. Either bound may be
    ``None``/0 to leave that side unbounded.

    Parameters
    ----------
    segments_zarr_path : str | Path
        Path to a ``*_segments.zarr`` array, ``(T, Z, Y, X)`` uint16.
    min_size_voxels : int | None
        Minimum 3-D volume in voxels a segment must have to be kept.
    fill_holes : bool
        Fill holes left by removed segments that were embedded in a larger one.
    max_size_voxels : int | None
        Maximum 3-D volume in voxels; larger segments are removed.

    Returns
    -------
    int
        Number of segments removed across all timepoints.
    """
    segments_zarr_path = Path(segments_zarr_path)
    if not segments_zarr_path.exists():
        print(f"  ⚠️ Segments file not found: {segments_zarr_path}")
        return 0

    size_min = int(min_size_voxels) if min_size_voxels else None
    size_max = int(max_size_voxels) if max_size_voxels else None
    if size_min is None and size_max is None:
        return 0

    seg = _open_segments_array(segments_zarr_path)
    if seg is None:
        return 0

    removed_total = 0
    pbar = tqdm(
        range(seg.shape[0]),
        desc="    Filtering",
        unit="tp",
        file=sys.stdout,
        dynamic_ncols=True,
    )

    for t in pbar:
        vol = np.asarray(seg[t])  # single 3-D volume (Z, Y, X)
        if not np.any(vol):
            continue

        n_before = len(np.unique(vol))
        filtered = segment_size_filter(
            vol, size_min=size_min, size_max=size_max, fill_holes=fill_holes
        )
        n_after = len(np.unique(filtered))

        removed = max(0, n_before - n_after)
        if removed:
            bounds = f">= {size_min}" if size_max is None else (
                f"<= {size_max}" if size_min is None else f"in [{size_min}, {size_max}]"
            )
            print(f"    [T{t}] removed {removed} segments not {bounds} voxels "
                  f"({n_before - 1} -> {n_after - 1} objects)")
            removed_total += removed

        seg[t] = filtered
        pbar.set_postfix(removed=removed_total)

    return removed_total
