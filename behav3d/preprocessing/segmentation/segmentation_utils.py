import numpy as np


def postprocess_mask(mask, fill_holes=True, opening_nr_pixels=1):
    """
    Postprocess a binary mask by filling holes and/or applying morphological opening.

    Parameters
    ----------
    mask : np.ndarray
        Binary mask (bool or 0/1). Can be 2D (Y, X), 3D (Z, Y, X).
    fill_holes : bool
        Whether to fill holes in the mask.
    opening_nr_pixels : int
        Number of pixels for morphological opening (0 to disable).

    Returns
    -------
    np.ndarray
        Postprocessed binary mask (bool).
    """
    from scipy.ndimage import binary_fill_holes
    from behav3d.preprocessing import open_mask

    mask = mask.astype(bool)

    if fill_holes:
        if mask.ndim == 2:
            mask = binary_fill_holes(mask)
        elif mask.ndim == 3:
            filled_mask = np.zeros_like(mask)
            for i in range(mask.shape[0]):
                filled_mask[i] = binary_fill_holes(mask[i])
            mask = filled_mask
        elif mask.ndim == 4:
            filled_mask = np.zeros_like(mask)
            for t in range(mask.shape[0]):
                for z in range(mask.shape[1]):
                    filled_mask[t, z] = binary_fill_holes(mask[t, z])
            mask = filled_mask

    if opening_nr_pixels and opening_nr_pixels > 0:
        mask = open_mask(mask, nr_pixels=opening_nr_pixels)

    return mask.astype(bool)


def segment_mask(mask, edt_thr=1.5, edt_thr_refined=None, segment_size_min=15, use_dims=2, n_workers=1):
    """
    Instance-segment a binary mask using EDT + watershed refinement.

    Delegates to the proven original implementation in napari_pixelclassifier
    which uses per-label EDT refinement (never discards entire objects).

    Parameters
    ----------
    mask : np.ndarray
        Binary mask (bool or 0/1).
    edt_thr : float
        EDT threshold for initial seed detection.
    edt_thr_refined : float or list of float, optional
        Additional refinement threshold(s).
    segment_size_min : int
        Minimum segment size in voxels.
    use_dims : int
        Dimensionality for EDT (2 or 3).
    n_workers : int
        Number of threads for parallel refinement.

    Returns
    -------
    np.ndarray
        Labeled segments (uint16).
    """
    from behav3d.preprocessing.segmentation.napari_pixelclassifier import segment_mask as _original_segment_mask
    return _original_segment_mask(
        mask, edt_thr=edt_thr, edt_thr_refined=edt_thr_refined,
        segment_size_min=segment_size_min, use_dims=use_dims, n_workers=n_workers,
    )
