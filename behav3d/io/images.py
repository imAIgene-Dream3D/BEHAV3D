from concurrent.futures import ThreadPoolExecutor
import numpy as np
import zarr
from pathlib import Path
from .formats.czi import load_czi, load_czi_metadata, load_elsizes_czi, get_czi_shape, load_czi_timepoint_czyx
from .formats.h5 import load_h5
from .formats.ims import load_ims, load_ims_metadata, load_ims_timepoint_czyx
from .formats.zarr import load_zarr, save_as_zarr, append_to_zarr
from .formats.tiff import load_tiff, get_tiff_shape, get_tiff_dimension_order
from .formats.liff import load_liff, get_liff_shape, load_liff_metadata, load_elsizes_liff, load_liff_timepoint_czyx

TARGET_AXIS_ORDER = "TCZYX"
import shutil
from tqdm import tqdm
from .formats.czi import load_czi, load_czi_metadata, load_elsizes_czi, get_czi_shape
from .formats.h5 import load_h5
from .formats.ims import load_ims, load_ims_metadata
from .formats.zarr import load_zarr, save_as_zarr, append_to_zarr, write_zarr_parallel
from .formats.tiff import load_tiff
from .formats.liff import load_liff, get_liff_shape, load_liff_metadata, load_elsizes_liff

def get_filepath_stem(path):
    path=Path(path)
    if path.suffix == ".zip" or path.suffix ==".bz2":
        path = Path(path.stem)
    return(path.stem)

def reorder_axes(img, source_order, target_order=TARGET_AXIS_ORDER):
    """
    Reorder image axes from source_order to target_order via np.transpose.

    Both orders must be strings of the same 5 axis labels (e.g. "TZCYX" → "TCZYX").
    Returns the image unchanged if the orders already match.
    """
    if source_order == target_order:
        return img
    perm = [source_order.index(ax) for ax in target_order]
    return np.transpose(img, axes=perm)

def load_image(path, group=None, mode="r", **kwargs):
    """
    Load an image from any supported format.
    Returns the array in its native/constructed axis order (no reordering).
    """
    path = Path(path)
    if path.suffix==".czi":
        img = load_czi(path, **kwargs)
    elif path.suffix==".h5":
        img = load_h5(path, substruct=group)
    elif path.suffix==".ims":
        img = load_ims(path)
    elif path.suffix==".zarr" or str(path).endswith(".zarr.zip"):
        img = load_zarr(path, group=group, mode=mode)
    elif path.suffix==".tif" or path.suffix==".tiff":
        img = load_tiff(path)
    elif path.suffix==".lif" or path.suffix==".liff":
        img = load_liff(path)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    return(img)

def load_image_metadata(path):
    path = Path(path)
    if path.suffix==".czi":
        metadata = load_czi_metadata(path)
    elif path.suffix==".h5":
        pass
    elif path.suffix==".ims":
        metadata = load_ims_metadata(path)
    elif path.suffix==".zarr" or str(path).endswith(".zarr.zip"):
        pass
    elif path.suffix==".tif" or path.suffix==".tiff":
        pass
    elif path.suffix==".lif" or path.suffix==".liff":
        metadata = load_liff_metadata(path)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    return(metadata)

def _load_czyx_for_timepoint(path, t, source_order):
    """
    Load a single timepoint from any supported format.
    Returns a numpy array in CZYX order (C, Z, Y, X).

    For CZI, LIF and IMS, only the requested timepoint is read from disk.
    For TIFF, H5 and Zarr, the full image is loaded and sliced (fallback).

    Parameters
    ----------
    path : Path
        Source image path.
    t : int
        Timepoint index in the source image.
    source_order : str
        5-char axis order of the source data (e.g. "TCZYX", "TZCYX").
    """
    path = Path(path)

    if path.suffix == ".czi":
        return load_czi_timepoint_czyx(path, t, source_order)

    elif path.suffix in (".lif", ".liff"):
        return load_liff_timepoint_czyx(path, t)

    elif path.suffix == ".ims":
        return load_ims_timepoint_czyx(path, t)

    elif path.suffix == ".zarr" or str(path).endswith(".zarr.zip"):
        # Zarr is lazy: slice without loading other timepoints
        dask_arr = load_zarr(path)
        if source_order != TARGET_AXIS_ORDER:
            import dask.array as da
            axes = [source_order.index(ax) for ax in TARGET_AXIS_ORDER]
            dask_arr = da.transpose(dask_arr, axes)
        return np.asarray(dask_arr[t])  # (C, Z, Y, X)

    else:
        # TIFF / H5 fallback: load full image, reorder to TCZYX, then slice.
        # This requires source_order to be correctly set (no default assumed).
        img = load_image(path)
        if source_order != TARGET_AXIS_ORDER:
            img = reorder_axes(img, source_order=source_order, target_order=TARGET_AXIS_ORDER)
        t_idx = list(TARGET_AXIS_ORDER).index("T")
        return np.asarray(np.take(img, t, axis=t_idx))  # (C, Z, Y, X)

def _normalize_n_workers(n_workers):
    return max(1, int(n_workers or 1))

def _resolve_timepoints(n_timepoints, t_start=None, t_end=None):
    timepoints = range(n_timepoints)
    if t_start is None or t_end is None:
        return timepoints, None

    t_start = max(0, t_start)
    t_end = min(n_timepoints - 1, t_end)
    if t_start > t_end:
        return range(n_timepoints), f"Warning: Invalid range {t_start} to {t_end}. Defaulting to full image."
    return range(t_start, t_end + 1), f"  Clipping timepoints: {t_start} to {t_end} (of {n_timepoints} total)"

def _slice_timepoints(img, timepoints):
    if isinstance(timepoints, range):
        return img[timepoints.start:timepoints.stop:timepoints.step]
    return img[timepoints]

def convert_file_to_zarr(
    path,
    outpath=None,
    axis_order=None,
    chunks=None,
    overwrite=False,
    t_start=None,
    t_end=None,
    n_workers=1,
    ):
    """
    Convert a single image file to .zarr in TCZYX axis order.

    Memory-efficient: reads and writes one timepoint at a time. Only one (C, Z, Y, X)
    volume is held in RAM at any moment (except for TIFF/H5 which require a full load).

    Parameters
    ----------
    path : Path
        Source image path.
    outpath : Path
        Destination .zarr path.
    axis_order : str
        The axis order of the source image (e.g. "TZCYX").
        Must be explicitly provided — no default is assumed.
    """
    if (
        (path.suffix == ".zarr" and path.exists()) or
        (outpath.exists() and not overwrite)
    ):
        print("Skipping conversion to zarr, as file already exists")
        return

    print(f"Converting {path.name} → {outpath}")

    if axis_order is None:
        raise ValueError(
            f"axis_order must be specified for '{path.name}'.\n"
            "Cannot assume a default axis order for any format. "
            "Ensure dimension_order is set in the metadata (e.g. via the Preprocessing widget)."
        )
    source_order = axis_order

    # Get source shape without loading data
    source_shape = get_image_shape(path)

    # Compute the TCZYX output shape by mapping source dims
    dim_size = dict(zip(source_order, source_shape))
    tczyx_shape = tuple(dim_size[ax] for ax in TARGET_AXIS_ORDER)
    n_t_total = tczyx_shape[0]
    czyx_shape = tczyx_shape[1:]

    # Apply time range clipping
    t_s = max(0, t_start) if t_start is not None else 0
    t_e = min(n_t_total - 1, t_end) if t_end is not None else n_t_total - 1
    if t_s > t_e:
        print(f"  Warning: invalid time range [{t_s}, {t_e}], using full range.")
        t_s, t_e = 0, n_t_total - 1
    n_t_out = t_e - t_s + 1

    if t_start is not None or t_end is not None:
        print(f"  Clipping timepoints: {t_s} to {t_e} (of {n_t_total} total)")

    out_shape = (n_t_out,) + czyx_shape
    if chunks is None:
        chunks = (1,) + czyx_shape

    # Load the first timepoint to get dtype (minimal I/O)
    first = _load_czyx_for_timepoint(path, t_s, source_order)
    dtype = first.dtype

    # Create the zarr store with the final TCZYX shape
    outpath.parent.mkdir(parents=True, exist_ok=True)
    store = zarr.open(
        zarr.storage.LocalStore(str(outpath)),
        mode="w",
        shape=out_shape,
        chunks=chunks,
        dtype=dtype,
    )

    print(f"  Writing {n_t_out} timepoints as TCZYX {out_shape}, dtype={dtype} ...")
    store[0] = np.asarray(first)
    del first

    for t_out, t_in in enumerate(range(t_s + 1, t_e + 1), start=1):
        czyx = _load_czyx_for_timepoint(path, t_in, source_order)
        store[t_out] = np.asarray(czyx)
        print(f"  Timepoint {t_in}/{t_e}", end="\r")

    print(f"\n  Done: {outpath.name}")

def _ensure_zarr(path, label="Image"):
    """Raise a clear error if a path does not point to a .zarr file."""
    path = Path(path)
    if path.suffix != ".zarr" and not str(path).endswith(".zarr.zip"):
        raise ValueError(
            f"{label} '{path.name}' is not in .zarr format ({path.suffix}).\n"
            "Please convert all images to Zarr in the Preprocessing step first."
        )
    return path
    path = Path(path)
    outpath = Path(outpath)
    n_workers = _normalize_n_workers(n_workers)
    if  (
        (path.suffix == ".zarr" and path.exists()) or 
        outpath.exists() and not overwrite
        ):
        print("Skipping conversion to zarr, as file already exists")
    
    else:
        if overwrite and outpath.exists():
            if outpath.is_dir():
                shutil.rmtree(outpath)
            else:
                outpath.unlink()
        print(f"Converting {path} to {outpath}")
        if path.suffix == ".czi" and outpath.suffix == ".zarr":
            img_shape = tuple(get_czi_shape(path))
            n_timepoints = img_shape[0]
            timepoints, timepoint_msg = _resolve_timepoints(n_timepoints, t_start=t_start, t_end=t_end)
            if timepoint_msg is not None:
                print(timepoint_msg)

            if n_workers == 1:
                for t in tqdm(timepoints, desc="Converting timepoints"):
                    img = np.expand_dims(load_czi(path, t=t), axis=0)
                    append_to_zarr(
                        img=img,
                        outpath=outpath,
                        chunks=chunks,
                    )
            else:
                timepoints = list(timepoints)
                if not timepoints:
                    raise ValueError("No timepoints selected for conversion")
                first_img = np.asarray(load_czi(path, t=timepoints[0]))
                write_zarr_parallel(
                    outpath=outpath,
                    shape=(len(timepoints),) + tuple(first_img.shape),
                    dtype=first_img.dtype,
                    chunks=chunks,
                    overwrite=overwrite,
                )
                write_zarr_parallel(outpath=outpath, index=0, data=first_img)

                def _write_czi_timepoint(args):
                    write_idx, t = args
                    img_t = np.asarray(load_czi(path, t=t))
                    write_zarr_parallel(outpath=outpath, index=write_idx, data=img_t)
                    return write_idx

                with ThreadPoolExecutor(max_workers=n_workers) as executor:
                    remaining = enumerate(timepoints[1:], start=1)
                    for _ in tqdm(executor.map(_write_czi_timepoint, remaining), total=max(0, len(timepoints) - 1), desc="Converting timepoints"):
                        pass
        else:
            img = load_image(path)
            timepoints, timepoint_msg = _resolve_timepoints(img.shape[0], t_start=t_start, t_end=t_end)
            if timepoint_msg is not None:
                print(timepoint_msg)
            img = _slice_timepoints(img, timepoints)

            if n_workers > 1 and outpath.suffix == ".zarr":
                write_zarr_parallel(
                    outpath=outpath,
                    shape=tuple(img.shape),
                    dtype=img.dtype,
                    chunks=chunks,
                    overwrite=overwrite,
                )

                def _write_loaded_timepoint(t):
                    write_zarr_parallel(outpath=outpath, index=t, data=np.asarray(img[t]))
                    return t

                with ThreadPoolExecutor(max_workers=n_workers) as executor:
                    for _ in tqdm(executor.map(_write_loaded_timepoint, range(img.shape[0])), total=img.shape[0], desc="Converting timepoints"):
                        pass
            else:
                if n_workers > 1 and outpath.suffix != ".zarr":
                    print("Parallel zarr writing only supports directory-backed .zarr outputs. Falling back to sequential conversion.")
                if chunks is None:
                    chunks = (1,) + img.shape[1:]
                save_as_zarr(
                    img=img, 
                    path=outpath, 
                    chunks=chunks
                    )

def convert_input_files_to_zarr(
    sample_name,
    current_cell_segments_path,
    raw_image_path,
    output_dir=None,
    chunks=None,
    overwrite=False,
    axis_order=None,
    ):
    """
    Convert cell segments and raw image to .zarr format for memory efficiency.

    Parameters
    ----------
    axis_order : str
        Source axis order (e.g. "TZCYX"). Passed through to convert_file_to_zarr
        for reordering to TCZYX. Must be explicitly provided.
    """
    current_cell_segments_path = Path(current_cell_segments_path)
    raw_image_path = Path(raw_image_path)
    
    filename = current_cell_segments_path.stem
    
    if output_dir is None:
        current_cell_zarr_out_path = current_cell_segments_path if current_cell_segments_path.suffix == ".zarr" or str(current_cell_segments_path).endswith(".zarr") else Path(current_cell_segments_path.parent, f"{filename}.zarr")
        raw_image_zarr_out_path = raw_image_path if raw_image_path.suffix == ".zarr" or str(raw_image_path).endswith(".zarr") else Path(raw_image_path.parent, f"{sample_name}.zarr")
    else:
        current_cell_zarr_out_path = Path(output_dir, f"{filename}.zarr")
        raw_image_zarr_out_path = Path(output_dir, f"{sample_name}.zarr")
    
    convert_file_to_zarr(
        path=current_cell_segments_path, 
        outpath=current_cell_zarr_out_path, 
        chunks=chunks,
        overwrite=overwrite,
        axis_order=axis_order,
    )
    
    convert_file_to_zarr(
        path=raw_image_path, 
        outpath=raw_image_zarr_out_path, 
        chunks=chunks,
        overwrite=overwrite,
        axis_order=axis_order,
    )
    
    return (current_cell_zarr_out_path, raw_image_zarr_out_path)

def load_elsizes(path):
    path = Path(path)
    if path.suffix==".czi":
        elsizes = load_elsizes_czi(path)
    elif path.suffix==".h5":
        pass
    elif path.suffix==".ims":
        pass
    elif path.suffix==".zarr" or str(path).endswith(".zarr.zip"):
        pass
    elif path.suffix==".tif" or path.suffix==".tiff":
        pass
    elif path.suffix==".lif" or path.suffix==".liff":
        elsizes = load_elsizes_liff(path)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    
    return(elsizes)

def get_image_shape(path):
    """
    Get the shape of an image without fully loading it into memory.
    Supports: .zarr, .zarr.zip, .czi, .lif/.liff, .tif/.tiff, .ims, .h5
    """
    path = Path(path)
    if path.suffix == ".zarr" or str(path).endswith(".zarr.zip"):
        img = load_zarr(path)
        shape = img.shape
    elif path.suffix == ".czi":
        shape = get_czi_shape(path)
    elif path.suffix == ".lif" or path.suffix == ".liff":
        shape = get_liff_shape(path)
    elif path.suffix == ".tif" or path.suffix == ".tiff":
        shape = get_tiff_shape(path)
    elif path.suffix == ".ims":
        shape = load_ims_metadata(path)
    elif path.suffix == ".h5":
        import h5py
        with h5py.File(str(path), "r") as f:
            keys = list(f.keys())
            if keys:
                shape = f[keys[0]].shape
            else:
                shape = ()
    else:
        raise ValueError(f"get_image_shape: unsupported format {path.suffix}")
    return tuple(shape)

def get_image_dimension_order(path):
    """
    Try to read the axis order directly from the **image file's** metadata.

    Only returns a value when the axis order is explicitly stored in the file
    (e.g. CZI dimension tags, OME-TIFF XML axes). For every other format the
    axis order is unknown and the user must select it manually.

    Returns
    -------
    str or None
        5-char axis order string (e.g. "TCZYX") read from the image file,
        or None if the file format does not contain axis order metadata.
    """
    path = Path(path)
    try:
        if path.suffix == ".czi":
            from aicspylibczi import CziFile
            czifile = CziFile(path)
            raw_dims = czifile.dims
            dim_order = "".join(ax for ax in raw_dims if ax in "TCZYX")
            if len(dim_order) == 5:
                return dim_order
            return None

        elif path.suffix == ".tif" or path.suffix == ".tiff":
            return get_tiff_dimension_order(path)  # None if not OME/ImageJ

        else:
            # LIF, IMS, H5, Zarr — no axis order metadata in the file.
            return None

    except Exception:
        return None
