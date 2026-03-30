from tifffile import imwrite, imread, TiffFile
from pathlib import Path
import numpy as np


def get_tiff_shape(path):
    """Get the shape of a TIFF image without loading the full data."""
    path = Path(path)
    if path.is_dir():
        tiff_paths = sorted(
            p for p in path.iterdir()
            if p.is_file() and p.suffix.lower() in {".tif", ".tiff"}
        )
        if not tiff_paths:
            return ()
        with TiffFile(str(tiff_paths[0])) as tf:
            slice_shape = tf.series[0].shape if tf.series else tf.pages[0].shape
        return (len(tiff_paths),) + tuple(slice_shape)
    with TiffFile(str(path)) as tf:
        if tf.series:
            return tuple(tf.series[0].shape)
        return tuple(tf.pages[0].shape)


def get_tiff_dimension_order(path):
    """
    Try to detect axis order from TIFF metadata (OME-TIFF / ImageJ TIFF).
    Returns a 5-char axis string if detection succeeds, None otherwise.
    """
    path = Path(path)
    if path.is_dir():
        return None
    try:
        with TiffFile(str(path)) as tf:
            if not tf.series:
                return None
            axes = tf.series[0].axes.upper()
            axes = axes.replace("S", "C").replace("I", "Z").replace("Q", "")
            known = set("TCZYX")
            axes = "".join(ax for ax in axes if ax in known)
            if len(axes) == 5 and set(axes) == known:
                return axes
    except Exception:
        pass
    return None


def load_tiff(path):
    """
    Load TIFF data from either:
    - a single .tif/.tiff file, or
    - a directory containing 2D .tif/.tiff files (stacked along axis 0)
    """
    path = Path(path)

    if path.is_dir():
        tiff_paths = sorted(
            [
                p for p in path.iterdir()
                if p.is_file() and p.suffix.lower() in {".tif", ".tiff"}
            ]
        )
        if len(tiff_paths) == 0:
            raise FileNotFoundError(f"No .tif/.tiff files found in folder:\n{path}")

        slices = [imread(p) for p in tiff_paths]
        ref_shape = slices[0].shape
        for p, arr in zip(tiff_paths, slices):
            if arr.shape != ref_shape:
                raise ValueError(
                    f"All TIFF slices must have same shape for stacking. "
                    f"Expected {ref_shape}, got {arr.shape} for:\n{p}"
                )

        return np.stack(slices, axis=0)

    img = imread(path)
    return img

def save_as_tiff(img, path):
    path = Path(path)
    if not path.parent.exists():
         raise FileNotFoundError(f"Parent folder of supplied .tiff path does not exist:\n{path}")
    imwrite(path, img)
