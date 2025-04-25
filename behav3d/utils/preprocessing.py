import numpy as np

from skimage.filters import median, threshold_sauvola
from skimage.morphology import binary_dilation, binary_erosion, binary_opening, binary_closing

from skimage.draw import ellipsoid
from skimage.morphology import disk

from behav3d.utils import rel_elsize
import warnings

def trim_zeros(arr):
    """
    Trim the zeros off the borders of an array
    """
    slices = tuple(slice(idx.min(), idx.max()+1) for idx in np.nonzero(arr))
    return(arr[slices])


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

def filter_mean(image, use_dimensions=3, radius=None, filter_shape=None, elsize=None):
    """
    Perform mean filtering for each pixel/voxel in the image based on a certain shape and radius
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
    print(f"Performing median smoothing with a ellipsoid of shape {k.shape}")
    
    data_smooth = skmean(image=image, footprint=k)
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