import numpy as np
import zarr
from pathlib import Path
from scipy import ndimage
from tqdm import tqdm
import sys


def filter_segments_by_size(
    segments_zarr_path,
    min_size_voxels,
    fill_holes=True,
):
    """Remove segments smaller than *min_size_voxels* from a segmentation zarr.
    If *fill_holes* is True, holes left by removed segments that were embedded
    inside larger segments are filled with the surrounding label.

    Parameters
    ----------
    segments_zarr_path : str | Path
        Path to the ``*_segments.zarr`` array (T, Z, Y, X) with uint16 labels.
    min_size_voxels : int
        Minimum 3-D volume (in voxels) a segment must have to be kept.
    fill_holes : bool
        Whether to fill holes left by removed "embedded" segments.
    """
    segments_zarr_path = Path(segments_zarr_path)

    if not segments_zarr_path.exists():
        print(f"  ⚠️ Segments file not found: {segments_zarr_path}")
        return 0

    seg = zarr.open(str(segments_zarr_path), mode="r+")
    
    # If the user opened a group instead of an array, pick the first array
    if isinstance(seg, zarr.Group):
        if len(list(seg.array_keys())) > 0:
            dataset_name = list(seg.array_keys())[0]
            print(f"  📂 Zarr is a group. Using internal array: {dataset_name}")
            seg = seg[dataset_name]
        else:
            print(f"  ⚠️ No arrays found in Zarr group: {segments_zarr_path}")
            return 0

    n_timepoints = seg.shape[0]
    print(f"  📦 Processing {n_timepoints} timepoints in {segments_zarr_path.name}")
    print(f"    📝 Path: {segments_zarr_path.absolute()}")

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
        n_labels_before = len(np.unique(vol))
        
        if not np.any(vol):
            continue

        # Count voxels per label directly (no re-labeling)
        sizes = np.bincount(vol.ravel())
        if len(sizes) <= 1:
            if n_labels_before > 1: # only print if there's data
                 print(f"    [T{t}] 0 segments < {min_size_voxels} voxels (Total labels: {n_labels_before})")
            continue

        # Find labels that are too small (skip background at index 0 and missing labels)
        labels = np.arange(len(sizes))
        small_mask = (sizes < min_size_voxels) & (sizes > 0)
        small_mask[0] = False # background never small
        small_labels = labels[small_mask]

        if len(small_labels) == 0:
            # We still might want to check for holes if fill_holes is True
            # but if no segments were removed and no segments are small, 
            # there are probably no new holes. 
            # For now, let's skip if no small labels AND no holes to fill.
            if not fill_holes:
                if n_labels_before > 1:
                     print(f"    [T{t}] 0 segments < {min_size_voxels} voxels (Total labels: {n_labels_before})")
                continue
        else:
             print(f"    [T{t}] Found {len(small_labels)} segments < {min_size_voxels} voxels (Total labels: {n_labels_before})")

        removed_total += len(small_labels)

        # Create mask of voxels to remove
        to_remove_mask = np.isin(vol, small_labels)
        
        if not fill_holes:
            vol[to_remove_mask] = 0
            seg[t] = vol
            pbar.set_postfix(removed=removed_total)
            continue

        # --- Hole Filling Logic ---
        # 1. Mask of kept voxels
        kept_mask = (vol > 0) & ~to_remove_mask
        
        # 2. Find holes in the kept mask
        filled_mask = ndimage.binary_fill_holes(kept_mask)
        holes_to_fill = filled_mask & ~kept_mask
        
        # We only want to fill holes that were specifically created/occupied by 
        # the removed segments (though binary_fill_holes might fill other 
        # pre-existing background holes too).
        
        # 3. Assign labels to new voxels in holes
        # For each voxel in holes_to_fill, assign it to the nearest kept label
        if np.any(holes_to_fill):
            # Distance transform to find nearest non-zero labels
            # We use the kept volume as source of labels
            kept_vol = vol.copy()
            kept_vol[to_remove_mask] = 0
            
            # Label connected components of holes to fill them efficiently?
            # Or just use a nearest-neighbor approach
            indices = ndimage.distance_transform_edt(
                kept_vol == 0, 
                return_distances=False, 
                return_indices=True
            )
            
            # Map hole voxels to their nearest kept neighbor's value
            # indices is (ndim, Z, Y, X). We want the values at those indices
            new_vol = kept_vol[tuple(indices)]
            
            # Only update the voxels that were part of holes
            vol[holes_to_fill] = new_vol[holes_to_fill]
            
            # Also ensure to_remove_mask that were NOT holes are set to 0
            # (If a small segment was on the edge, it's not a hole)
            not_hole = to_remove_mask & ~holes_to_fill
            vol[not_hole] = 0
        else:
            # No holes to fill, just clear small labels
            vol[to_remove_mask] = 0

        seg[t] = vol
        
        # VERIFICATION: Read back and check
        check_vol = np.asarray(seg[t])
        n_labels_after = len(np.unique(check_vol))
        
        if n_labels_before != n_labels_after:
             print(f"    ✅ T{t} updated: {n_labels_before} → {n_labels_after} labels.")
        
        if not np.array_equal(vol, check_vol):
            print(f"    ⚠️ T{t} written volume does not match read volume!")

    print(f"  ✨ Finished filtering {segments_zarr_path.name}. Total removed: {removed_total}")
    return removed_total
