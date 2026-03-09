from tifffile import imwrite, imread
from pathlib import Path
import numpy as np

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
