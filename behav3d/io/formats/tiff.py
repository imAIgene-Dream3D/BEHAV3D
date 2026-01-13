from tifffile import imwrite, imread
from pathlib import Path

def load_tiff(path):
    """
    Loading .tif/.tiff images
    """
    img = imread(path)
    return(img)

def save_as_tiff(img, path):
    path = Path(path)
    if not path.parent.exists():
         raise FileNotFoundError(f"Parent folder of supplied .tiff path does not exist:\n{path}")
    imwrite(path, img)
