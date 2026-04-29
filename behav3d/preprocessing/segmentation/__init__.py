import warnings
from collections import defaultdict
import numpy as np
from skimage.measure import label
from scipy import ndimage
from .multicolor_segment_processing import (
    calculate_multicolor_overlap,
    process_multicolor_segments,
    visualize_multicolor_overlap,
    apply_multicolor_segment_correction_for_base,
)

def keep_largest_connected_components(segments):
    """
    Select all separate connected components of the same segment ID
    Only keep the largest one, set other subsegments to 0
    """
    relabeled_segments = label(segments>0)
    for segment in np.unique(segments):
        if segment == 0: 
            continue
        mask = segments == segment

        # Calculate the size of the connected component and only keep the largest
        subsegments, subsizes = np.unique(relabeled_segments[mask], return_counts=True)
        biggest_segment = subsegments[np.argmax(subsizes)]
        segments[(mask) & (relabeled_segments != biggest_segment)]=0
    return segments

def segment_size_filter(segments, size_min=None, size_max=None):
    def _filter_one_volume(vol):
        vol = np.asarray(vol).copy()
        if not np.any(vol):
            return vol

        labels, counts = np.unique(vol, return_counts=True)
        remove_labels = []
        for lbl, count in zip(labels, counts):
            if lbl == 0:
                continue
            if size_min is not None and count < size_min:
                remove_labels.append(lbl)
                continue
            if size_max is not None and count > size_max:
                remove_labels.append(lbl)

        if not remove_labels:
            return vol

        to_remove_mask = np.isin(vol, remove_labels)
        kept_mask = (vol > 0) & ~to_remove_mask
        if not np.any(kept_mask):
            vol[to_remove_mask] = 0
            return vol

        filled_mask = ndimage.binary_fill_holes(kept_mask)
        holes_to_fill = filled_mask & ~kept_mask

        if np.any(holes_to_fill):
            kept_vol = vol.copy()
            kept_vol[to_remove_mask] = 0
            indices = ndimage.distance_transform_edt(
                kept_vol == 0,
                return_distances=False,
                return_indices=True,
            )
            new_vol = kept_vol[tuple(indices)]
            vol[holes_to_fill] = new_vol[holes_to_fill]
            vol[to_remove_mask & ~holes_to_fill] = 0
        else:
            vol[to_remove_mask] = 0
        return vol

    segments = np.asarray(segments)
    if segments.ndim == 4:
        out = np.zeros_like(segments)
        for t in range(segments.shape[0]):
            out[t] = _filter_one_volume(segments[t])
        segments[...] = out
        return segments

    filtered = _filter_one_volume(segments)
    segments[...] = filtered
    return segments

def segment_2d_filter(segments):
    """
    Remove segments from a 3D labeled array that are 2D in any axis (Z, Y, or X).
    
    Parameters:
    - segments (ndarray): 3D labeled array (z, y, x)
    
    Returns:
    - ndarray: 3D labeled array with 2D segments removed (set to 0)
    """
    segments = segments.copy()
    coords = np.argwhere(segments > 0)  # [z, y, x]
    if coords.size == 0:
        return segments

    labels = segments[segments > 0]
    if labels.size == 0:
        return segments

    # Build axis-wise mappings
    label_to_z = defaultdict(set)
    label_to_y = defaultdict(set)
    label_to_x = defaultdict(set)

    for (z, y, x), label in zip(coords, labels):
        label_to_z[label].add(z)
        label_to_y[label].add(y)
        label_to_x[label].add(x)

    # Identify flat (2D) labels: only 1 slice along any axis
    flat_labels = {
        label for label in set(labels)
        if len(label_to_z[label]) == 1 or len(label_to_y[label]) == 1 or len(label_to_x[label]) == 1
    }

    # Zero out flat segments
    mask = np.isin(segments, list(flat_labels))
    segments[mask] = 0

    return segments

def get_border_segments(segments):
    """
    For segmentation result, return a list of the IDs of the segments that touch the borders of the image
    """
    shape = segments.shape
    ndim = segments.ndim
    border_values = set()
    # Iterate over each dimension to collect the border values
    for dim in range(ndim):
        slices = [slice(None)] * len(shape)
        # Collect the values from the minimum edge of this dimension
        slices[dim] = 0
        border_values.update(segments[tuple(slices)].flatten())
        # Collect the values from the maximum edge of this dimension
        slices[dim] = -1
        border_values.update(segments[tuple(slices)].flatten())
    border_values = list(border_values)
    border_values = [int(i) for i in border_values]
    return(border_values)
    
def remove_boundary_segments(segments, border_segments=None):
    """
    Based on the list returned by 'get_border_segments',
    remove these IDs from the segmentation image
    """
    if border_segments is None:          
        border_segments = get_border_segments(segments)
    border_mask=~np.isin(segments, border_segments)
    segments = segments * border_mask
    return(segments)
