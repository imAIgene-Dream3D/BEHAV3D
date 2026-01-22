import time
import warnings
from pathlib import Path
import numpy as np
from scipy.ndimage import distance_transform_edt
from skimage.draw import ellipsoid
from skimage.filters import median, threshold_sauvola
from skimage.morphology import (
    binary_closing,
    binary_dilation,
    binary_erosion,
    binary_opening,
    disk,
)
from behav3d.io.images import convert_file_to_zarr, get_image_shape
from behav3d.core.utils import rel_elsize

def calculate_edt(image, use_dims=3, elsize=None):
    """
    Perform median filtering for each pixel/voxel in the image based on a certain shape and radius
    """
    img_dim = image.ndim
    if elsize is None:
        elsize = [1]*img_dim
    elsize=rel_elsize(elsize)
    
    edt_result = np.zeros_like(image, dtype=np.float32)

    if use_dims == 2:
        # Loop through all leading dimensions except the last two
        for index in np.ndindex(image.shape[:-2]):  
            edt_result[index] = distance_transform_edt(image[index], sampling=elsize[-2:])
    
    elif use_dims == 3:
        # Loop through all leading dimensions except the last three
        for index in np.ndindex(image.shape[:-3]):  
            edt_result[index] = distance_transform_edt(image[index], sampling=elsize[-3:])
    
    return edt_result


def draw_ellipsoid(shape=[], radius=1, elsize=[1,1,1], remove_border_zeros=True):
    """
    Create a numpy array containing an elipsoid of ones
    """
    if shape:
        z,y,x = shape
    elif radius:
        z,y,x = [radius]*3

    elsize_scaled=rel_elsize(elsize)
    # print(elsize_scaled)
    z_temp, y_temp, x_temp = z,y,x
    if z==0:
        z_temp=1
    if y==0:
        y_temp=1
    if x==0:
        x_temp=1

    shape=ellipsoid(z_temp, y_temp, x_temp, spacing=elsize_scaled)

    if z==0:
        shape=shape[[2],:,:]
    if y==0:
        shape=shape[:,[2],:]
    if x==0:
        shape=shape[:,:,[2]]
    
    if remove_border_zeros:
        shape=trim_zeros(shape)
    return(shape)

def create_footprint(img_dims, use_dimensions, nr_pixels):
    """
    Create the footprint (disk or ellipse) that is often used for dilation, closing, etc. operation
    It matches the footprint to then target image nr of dimensions
    """
    def expand_selem(selem, dims):
        while selem.ndim != dims:
            selem = np.expand_dims(selem, axis=0)
        return(selem)
    
    if use_dimensions == 2:
        k = disk(nr_pixels)
    elif use_dimensions == 3:
        filter_shape=[nr_pixels]*3
        k = draw_ellipsoid(shape=filter_shape)
        
    while k.ndim < img_dims:
        k = expand_selem(k, img_dims)
    return(k)

def dilate_mask(image, use_dimensions=2, nr_pixels=1):
    """
    Perform mask dilation, increasing the size of the mask
    """
    img_dim = image.ndim
    k = create_footprint(img_dim, use_dimensions=use_dimensions, nr_pixels=nr_pixels)
    data_expanded = binary_dilation(image=image, footprint=k)
    return(data_expanded)

def erode_mask(image, use_dimensions=2, nr_pixels=1):
    """
    Perform mask erosion, decreasing the size of the mask
    """
    img_dim = image.ndim
    k = create_footprint(img_dim, use_dimensions=use_dimensions, nr_pixels=nr_pixels)
    # print(f"Performing median smoothing with a ellipsoid of shape {k.shape}")
    data_expanded = binary_erosion(image=image, footprint=k)
    return(data_expanded)

def open_mask(image, use_dimensions=2, nr_pixels=1):
    """
    Perform mask opening (Smooths the mask), doing subsequent erosion then dilation 
    to remove small separated parts of the mask while keeping shape of larger parts. 
    """
    img_dim = image.ndim
    k = create_footprint(img_dim, use_dimensions=use_dimensions, nr_pixels=nr_pixels)
    data_expanded = binary_opening(image=image, footprint=k)
    return(data_expanded)

def close_mask(image, use_dimensions=2, nr_pixels=1):
    """
    Perform mask opening (Tries to fill holes or dents in the mask), doing subsequent dilation then erosion 
    to fill small holes or connect unconnect parts of the mask while not generating new mask
    """
    img_dim = image.ndim
    k = create_footprint(img_dim, use_dimensions=use_dimensions, nr_pixels=nr_pixels)
    data_expanded = binary_closing(image=image, footprint=k)
    return(data_expanded)

def filter_median(image, use_dimensions=3, radius=None, filter_shape=None, elsize=None):
    """
    Perform median filtering for each pixel/voxel in the image based on a certain shape and radius
    """
    assert radius is not None or filter_shape is not None, "Supply either 'radius' or 'filter_shape'"
    img_dim = image.ndim
    if elsize is None:
        elsize = [1]*img_dim
    elsize=rel_elsize(elsize)
    
    def expand_selem(selem, dims):
        while selem.ndim != dims:
            selem = np.expand_dims(selem, axis=0)
        return(selem)
    
    if use_dimensions == 2:
        k = disk(radius)
    elif use_dimensions == 3:
        if filter_shape is None and radius:
            filter_shape=[radius]*3
        k = draw_ellipsoid(shape=filter_shape)
    # !!!! Change shape here!!!
    while k.ndim < image.ndim:
        k = expand_selem(k, img_dim)
    # print(f"Performing median smoothing with a ellipsoid of shape {k.shape}")
    data_smooth = median(image=image, footprint=k)
    return(data_smooth)

def sauvola_thresholding(
    image,
    use_dimensions=3,
    window_size=5,
    sd_weight=0.2, #k, sd=standard deviation
    max_sd=None, #r
    elsize=[1,1,1],
    absmin=0,
    override_thr=None,
    in_plane=False
):
    """
    Perform sauvola thresholding. Creates a mask by keeping local pixel/voxel intensities into account

    k: For low contrast images, make k as small as possible, otherwise high k will do

    https://hal.archives-ouvertes.fr/hal-02181880/document
    """
    
    mask=image>absmin

    if isinstance(window_size, int):
        window_size_tmp=[window_size]*use_dimensions
        while len(window_size_tmp) < image.ndim:
            window_size_tmp = [1] + window_size_tmp
        while len(elsize) < image.ndim:
            elsize = [1] + elsize
            
        window_size_tmp=[int(round(x/y)) for x,y in zip(window_size_tmp, rel_elsize(elsize))]

        window_size=[]
        for el in window_size_tmp:
            if el%2 != 1:
                el=int(el)+1
            window_size.append(el)
        if in_plane:
            window_size=window_size[1:]
                
    if any([x%2 != 1 for x in list(window_size)]):
        raise ValueError("One of the window_size values is not odd. Sauvola requires odd window_size")

    # while len(window_size) != image.ndim:
    #     window_size = [1] + window_size
            
    thr = threshold_sauvola(image, window_size=window_size, k=sd_weight, r=max_sd)

    mask_sauvola = image > thr
    mask &= mask_sauvola

    if override_thr:
        override_mask = image > override_thr
        mask |= override_mask

    return(mask)


def trim_zeros(arr):
    """
    Trim the zeros off the borders of an array
    """
    slices = tuple(slice(idx.min(), idx.max()+1) for idx in np.nonzero(arr))
    return(arr[slices])

def zeropad_image_to_match_shape(img, target_shape, axes=[-3, -2, -1]):
    """
    Zero-pad an image array only along specified axes to match the target shape.

    Parameters
    ----------
    img : np.ndarray
        Input image of shape (e.g. (C, T, Z, Y, X)).
    target_shape : tuple of int
        Desired final shape.
    axes : list of int
        Axes along which to pad (default: last three spatial axes [-3, -2, -1]).
        
    Returns
    -------
    np.ndarray
        Padded image matching target_shape along the specified axes.
    """
    # Ensure positive axis indices
    ndim = img.ndim
    axes = [a if a >= 0 else ndim + a for a in axes]
    
    # Initialize output same as input
    padded = img

    # Pad each specified axis if needed
    pad_width = [(0, 0)] * ndim
    for ax in axes:
        diff = target_shape[ax] - img.shape[ax]
        if diff > 0:
            pad_width[ax] = (0, diff)  # pad zeros at the end

    # Apply padding only where needed
    padded = np.pad(padded, pad_width=pad_width, mode='constant', constant_values=0)

    return padded

def calc_z_projection(im, z_axis=-3, projection='max'):
    """
    Calculate the z projection of a 3D czi_utils

    Args:
        im (array-like): 3D czi_utils, dimension order is ZYX by default
        z_axis (int): The index of Z-plane
        projection (string, or tuple of strings): The method to compute z-projection. by default ('max', 'mean', 'std')

    Returns:
        results (list of tuples): [(method, value), ...], the projection method and the total intensity of the composite czi_utils after z-projection
    """

    if projection == 'max':
        return np.max(im, axis=z_axis).astype(im.dtype)
    elif projection == 'mean':
        return np.mean(im, axis=z_axis).astype(im.dtype)
    elif projection == 'std':
        return np.std(im, axis=z_axis).astype(im.dtype)
    else:
        warnings.warn(f'Projection {projection} is not defined')
    return -1

def convert_input_files_to_zarr(
    output_dir,
    metadata,
    ):
    
    for idx, sample in metadata.iterrows():
        print(f"Processing sample: {sample['sample_name']}") 
        start_time = time.time()
        
        sample_name = sample['sample_name']
        raw_image_path = Path(sample['raw_image_path'])
        raw_image_zarr =  Path(output_dir, "images", sample_name, f"{sample_name}.zarr")

        convert_file_to_zarr(
            path=raw_image_path, 
            outpath=raw_image_zarr, 
            overwrite=False
        )
                
        metadata.at[idx, "raw_image_path"] = str(raw_image_zarr)
    
    return(metadata)

def get_all_image_shapes(metadata):
    """
    Get the shapes of all images in the metadata.
    """
    shapes = []
    for _, row in metadata.iterrows():
        path = Path(row['raw_image_path'])
        if path.exists():
            shape = get_image_shape(path)
            shapes.append((row['sample_name'], shape))
        else:
            shapes.append((row['sample_name'], "File not found"))
    return shapes