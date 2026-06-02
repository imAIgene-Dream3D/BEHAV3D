import zarr
import dask.array as da
import numpy as np
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

def write_zarr_parallel(
    outpath,
    index=None,
    data=None,
    *,
    group=None,
    shape=None,
    dtype=None,
    chunks=None,
    overwrite=False,
    attrs=None,
):
    """
    Initialize or write to a pre-allocated .zarr array by fixed index.
    """
    outpath = Path(outpath)
    if outpath.suffix != ".zarr":
        raise ValueError("Parallel zarr writing only supports directory-backed .zarr outputs")

    def _normalize_chunks(shape, chunks):
        if chunks is None:
            return (1,) + tuple(shape[1:])
        chunks = tuple(chunks)
        if len(chunks) == len(shape) - 1:
            return (1,) + chunks
        if len(chunks) != len(shape):
            raise ValueError(
                f"Chunks must have {len(shape) - 1} or {len(shape)} dimensions, got {len(chunks)}"
            )
        return chunks

    if shape is not None:
        shape = tuple(shape)
        if len(shape) < 1:
            raise ValueError("Shape must include the leading indexed axis")

    if data is not None:
        data = np.asarray(data)
        if shape is None and index is None:
            shape = (1,) + tuple(data.shape)
        if dtype is None:
            dtype = data.dtype

    dtype_obj = np.dtype(dtype) if dtype is not None else None

    normalized_chunks = None if shape is None else _normalize_chunks(shape, chunks)

    if group is not None:
        root = zarr.open_group(str(outpath), mode="a")
        if shape is not None and dtype_obj is not None and (overwrite or group not in root):
            if group in root:
                del root[group]
            arr = root.create_dataset(
                str(group),
                shape=shape,
                chunks=normalized_chunks,
                dtype=dtype_obj,
                overwrite=True,
            )
            if attrs:
                arr.attrs.update(attrs)
        else:
            if group not in root:
                raise ValueError(f"Zarr group '{group}' must be initialized before fixed-index writes")
            arr = root[str(group)]
            if shape is not None and tuple(arr.shape) != shape:
                raise ValueError(f"Existing zarr shape {tuple(arr.shape)} does not match requested shape {shape}")
            if dtype_obj is not None and np.dtype(arr.dtype) != dtype_obj:
                raise ValueError(f"Existing zarr dtype {arr.dtype} does not match requested dtype {dtype_obj}")
            if attrs:
                arr.attrs.update(attrs)
    else:
        if shape is not None and dtype_obj is not None and (not outpath.exists() or overwrite):
            if outpath.exists():
                if outpath.is_dir():
                    shutil.rmtree(outpath)
                else:
                    outpath.unlink()
            arr = zarr.open(
                str(outpath),
                mode="w",
                shape=shape,
                chunks=normalized_chunks,
                dtype=dtype_obj,
            )
            if attrs:
                arr.attrs.update(attrs)
        else:
            if not outpath.exists():
                raise ValueError("Zarr array must be initialized before fixed-index writes")
            arr = zarr.open(str(outpath), mode="r+")
            if shape is not None and tuple(arr.shape) != shape:
                raise ValueError(f"Existing zarr shape {tuple(arr.shape)} does not match requested shape {shape}")
            if dtype_obj is not None and np.dtype(arr.dtype) != dtype_obj:
                raise ValueError(f"Existing zarr dtype {arr.dtype} does not match requested dtype {dtype_obj}")
            if attrs:
                arr.attrs.update(attrs)

    if index is None:
        return arr

    if data is None:
        raise ValueError("Data must be provided when writing by index")
    if not 0 <= int(index) < arr.shape[0]:
        raise IndexError(f"Index {index} is out of bounds for zarr with length {arr.shape[0]}")
    if tuple(data.shape) != tuple(arr.shape[1:]):
        raise ValueError(
            f"Data shape {tuple(data.shape)} does not match target item shape {tuple(arr.shape[1:])}"
        )

    arr[int(index)] = data
    return arr

def zip_zarr(path):
    outpath = Path(f"{path}.zip")
    if outpath.exists():
        outpath.unlink()
    shutil.make_archive(path, "zip", path)
    shutil.rmtree(path)
    return(outpath)
