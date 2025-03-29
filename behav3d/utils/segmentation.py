import numpy as np
from skimage.measure import label, find_contours
from behav3d.utils import rel_elsize
from scipy.ndimage import distance_transform_edt

def keep_largest_connected_components(segments):
    """
    Select all separate connected components of the same segment ID
    Only keep the largest one, set other subsegments to 0
    """
    relabeled_segments = label(segments>0)
    for segment in np.unique(segments):
        if segment == 0: 
            continue
        mask = segments == segment

        # Calculate the size of the connected component and only keep the largest
        subsegments, subsizes = np.unique(relabeled_segments[mask], return_counts=True)
        biggest_segment = subsegments[np.argmax(subsizes)]
        segments[(mask) & (relabeled_segments != biggest_segment)]=0
    return segments

def segment_size_filter(segments, size_min=None, size_max=None):
    labels, counts = np.unique(segments, return_counts=True)
    for label, count in zip(labels, counts):
        if size_min != None:
            if count < size_min:
                segments[segments == label] = 0
        if size_max!= None:
            if count > size_max:
                segments[segments == label] = 0
    return(segments)

def get_border_segments(segments):
    """
    For segmentation result, return a list of the IDs of the segments that touch the borders of the image
    """
    shape = segments.shape
    ndim = segments.ndim
    border_values = set()
    # Iterate over each dimension to collect the border values
    for dim in range(ndim):
        slices = [slice(None)] * len(shape)
        # Collect the values from the minimum edge of this dimension
        slices[dim] = 0
        border_values.update(segments[tuple(slices)].flatten())
        # Collect the values from the maximum edge of this dimension
        slices[dim] = -1
        border_values.update(segments[tuple(slices)].flatten())
    border_values = list(border_values)
    border_values = [int(i) for i in border_values]
    return(border_values)
    
def remove_boundary_segments(segments, border_segments=None):
    """
    Based on the list returned by 'get_border_segments',
    remove these IDs from the segmentation image
    """
    if border_segments is None:          
        border_segments = get_border_segments(segments)
    border_mask=~np.isin(segments, border_segments)
    segments = segments * border_mask
    return(segments)

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