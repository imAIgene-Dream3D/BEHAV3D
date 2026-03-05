import numpy as np
from pathlib import Path
from .formats.czi import load_czi, load_czi_metadata, load_elsizes_czi, get_czi_shape
from .formats.h5 import load_h5
from .formats.ims import load_ims, load_ims_metadata
from .formats.zarr import load_zarr, save_as_zarr, append_to_zarr
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

def convert_file_to_zarr(
    path,
    outpath=None,
    chunks=None,
    overwrite=False
    ):
    if  (
        (path.suffix == ".zarr" and path.exists()) or 
        outpath.exists() and not overwrite
        ):
        print("Skipping conversion to zarr, as file already exists")
    
    else:
        print(f"Converting {path} to {outpath}")
        img = load_image(path)
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
