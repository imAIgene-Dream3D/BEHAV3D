from concurrent.futures import ThreadPoolExecutor
import numpy as np
from pathlib import Path
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

def convert_axis_order(img, default_axis_order, user_axis_order):
    """
    Convert axis order of image
    """
    if user_axis_order==default_axis_order:
        return img
    else:
        user_axes = list(user_axis_order)
        default_axes = list(default_axis_order)
        axis_map = {user_axes[i]: default_axes[i] for i in range(len(user_axes))}
        new_order = [axis_map[axis] for axis in default_axes]
        return np.transpose(img, axes=new_order)

def load_image(path, axis_order="TCZYX", group=None, mode="r", **kwargs):
    path = Path(path)
    default_axis_order = "TCZYX"
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
    
    img = convert_axis_order(
        img, 
        default_axis_order=default_axis_order, 
        user_axis_order = axis_order
        )
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
    chunks=None,
    overwrite=False,
    t_start=None,
    t_end=None,
    n_workers=1,
    ):
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
    overwrite=False
    ):
    """
    Convert cell segments and raw image to .zarr format for memory efficiency.
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
        overwrite=overwrite
    )
    
    convert_file_to_zarr(
        path=raw_image_path, 
        outpath=raw_image_zarr_out_path, 
        chunks=chunks,
        overwrite=overwrite
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
    Loads an image and prints its shape and axis sizes.
    """
    path = Path(path)
    if path.suffix == ".zarr":
        img = load_zarr(path)
        shape = img.shape
    elif path.suffix == ".czi":
        shape = get_czi_shape(path)
    elif path.suffix==".lif" or path.suffix==".liff":
        shape = get_liff_shape(path)
    shape = tuple(shape)
    return shape

def get_image_dimension_order(path):
    """
    Loads an image and prints its shape and axis sizes.
    """
    path = Path(path)
    if path.suffix == ".czi":
        from aicspylibczi import CziFile
        czifile = CziFile(path)
        dim_order = czifile.dims[-5:]
    elif path.suffix==".lif" or path.suffix==".liff":
        # readlif doesn't expose dims directly; assume standard TCZYX
        dim_order = "TCZYX"
    else:
        dim_order = "TCZYX"
    return dim_order
