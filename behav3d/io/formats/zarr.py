import zarr
import dask.array as da
import shutil
from pathlib import Path

def load_zarr(path, group=None, mode="r"):
    """
    Loading .zarr images
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Path to zarr file does not exist:\n{path}")
    
    if path.suffix==".zip":
        zarr_store = zarr.storage.ZipStore(path)
    else:
        zarr_store = zarr.storage.LocalStore(path)
        
    if group:
        dask_img = da.from_zarr(zarr_store, component=group)
    else:
        dask_img = da.from_zarr(zarr_store)
    return(dask_img)

def save_as_zarr(
    img, 
    path, 
    chunks=None, 
    group=None
    ):
    path=Path(path)
    zipping=False
    if path.suffix == ".zip":
        zipping=True
        path = Path(path.parent, path.stem)
        
    if path.suffix==".zarr":
        if chunks is None:
            chunks = (1,) + img.shape[1:]
        if not isinstance(img, da.Array):
            img = da.from_array(img, chunks=chunks)
 
        zarr_store = zarr.storage.LocalStore(path)
        # If group is specified, create or open that group
        if group:
            da.to_zarr(img, zarr_store, component=group, compute=True, overwrite=True)
        else:
            da.to_zarr(img, zarr_store, compute=True, overwrite=True)
        
    img=None
    if zipping:
        shutil.make_archive(path, "zip", path)
        shutil.rmtree(path)

def append_to_zarr(img, outpath, chunks=None):
    """
    Append a timepoint to an existing .zarr array
    If non-existent, create the .zarr array.
    chunks excludes the leading time axis.
    """
    outpath = Path(outpath)
    if chunks is None:
        chunks = img.shape[1:]

    if not outpath.exists():
        zarr_file = zarr.open(
            outpath, 
            mode='w', 
            shape=(0,) + img.shape[1:], 
            chunks=(1,) + tuple(chunks), 
            dtype=img.dtype
            )
    else:
        zarr_file = zarr.open(outpath, mode='a')
    zarr_file.append(img)

def zip_zarr(path):
    outpath = Path(f"{path}.zip")
    if outpath.exists():
        outpath.unlink()
    shutil.make_archive(path, "zip", path)
    shutil.rmtree(path)
    return(outpath)
