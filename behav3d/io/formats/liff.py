from pathlib import Path
import numpy as np

# pip install readlif
from readlif.reader import LifFile


def _get_liff_series_dims(img_series):
    """Extract dimension sizes from a readlif image series object."""
    n_channels = img_series.channels
    n_z = img_series.nz if hasattr(img_series, 'nz') else (
        img_series.dims.z if hasattr(img_series, 'dims') and hasattr(img_series.dims, 'z') else 1
    )
    n_t = img_series.nt if hasattr(img_series, 'nt') else 1
    height = img_series.dims.y if hasattr(img_series, 'dims') and hasattr(img_series.dims, 'y') else img_series.dims[1]
    width = img_series.dims.x if hasattr(img_series, 'dims') and hasattr(img_series.dims, 'x') else img_series.dims[0]
    return n_t, n_channels, n_z, height, width


def _get_liff_frame(img_series, z, t, c):
    """Get a single frame, handling different readlif API versions."""
    try:
        return np.asarray(img_series.get_frame(z=z, t=t, c=c))
    except TypeError:
        return np.asarray(img_series.get_frame(z=z, t=t, channel=c))


def load_liff_timepoint_czyx(path, t):
    """
    Load a single timepoint from a LIF file.
    Returns array in CZYX order (C, Z, Y, X) without loading other timepoints.
    """
    path = Path(path)
    lf = LifFile(str(path))
    if len(lf.image_list) == 0:
        raise ValueError(f"No images found in LIF file: {path}")
    img_series = lf.get_image(0)
    _, n_channels, n_z, _, _ = _get_liff_series_dims(img_series)

    z_slices = []
    for z in range(n_z):
        c_slices = [_get_liff_frame(img_series, z=z, t=t, c=c) for c in range(n_channels)]
        z_slices.append(np.stack(c_slices, axis=0))  # (C, Y, X)
    return np.stack(z_slices, axis=1)  # (C, Z, Y, X)


def _load_liff_raw(path):
    """
    Internal helper: open a .lif file and return (array, dims_string).
    Returns array in (T, C, Z, Y, X) order.
    """
    path = Path(path)
    lf = LifFile(str(path))
    if len(lf.image_list) == 0:
        raise ValueError(f"No images found in LIF file: {path}")
    img_series = lf.get_image(0)
    n_t, n_channels, n_z, _, _ = _get_liff_series_dims(img_series)

    frames = []
    for t in range(n_t):
        time_frames = []
        for z in range(n_z):
            channel_frames = [_get_liff_frame(img_series, z=z, t=t, c=c) for c in range(n_channels)]
            time_frames.append(np.stack(channel_frames, axis=0))  # (C, Y, X)
        frames.append(np.stack(time_frames, axis=0))  # (Z, C, Y, X)

    arr = np.stack(frames, axis=0)              # (T, Z, C, Y, X)
    arr = np.transpose(arr, axes=(0, 2, 1, 3, 4))  # → (T, C, Z, Y, X)
    dims = "TCZYX"
    return arr, dims


def load_liff(path):
    """
    Load a .lif/.liff file and return data as (T, C, Z, Y, X).
    
    The readlif library always addresses frames by (t, z, c), so the
    returned array is always 5D TCZYX — including singletons at size 1
    (e.g. T=1, C=1, or Z=1).
    """
    arr, _dims = _load_liff_raw(path)
    return arr


def get_liff_dimension_order(path=None):
    """
    Return the axis order for LIF files.

    LIF always has all 5 dimensions (T, C, Z, Y, X) — the readlif API addresses
    frames by (t, z, c), so the order is always TCZYX even when T=1, C=1, or Z=1.
    """
    return "TCZYX"


def get_liff_shape(path, take_dims="TCZYX"):
    """
    Return the shape of a .lif/.liff file for the requested axes (default: TCZYX).
    
    Gets shape without loading the full image data into memory.
    """
    path = Path(path)
    lf = LifFile(str(path))
    
    if len(lf.image_list) == 0:
        return [1] * len(take_dims)
    
    img_series = lf.get_image(0)
    
    # Extract dimensions
    n_channels = img_series.channels
    n_z = img_series.nz if hasattr(img_series, 'nz') else (img_series.dims.z if hasattr(img_series, 'dims') and hasattr(img_series.dims, 'z') else 1)
    n_t = img_series.nt if hasattr(img_series, 'nt') else 1
    height = img_series.dims.y if hasattr(img_series, 'dims') and hasattr(img_series.dims, 'y') else img_series.dims[1]
    width = img_series.dims.x if hasattr(img_series, 'dims') and hasattr(img_series.dims, 'x') else img_series.dims[0]
    
    size_by_dim = {
        'T': n_t,
        'C': n_channels,
        'Z': n_z,
        'Y': height,
        'X': width
    }
    
    out_shape = [int(size_by_dim.get(ax, 1)) for ax in take_dims]
    return out_shape


def load_liff_metadata(path):
    """
    Load metadata from a .lif/.liff file (pixel sizes, time interval, etc).
    
    Returns a dict with keys: elsize_x, elsize_y, elsize_z, time_interval
    (in micrometers and seconds, respectively).
    """
    path = Path(path)
    lf = LifFile(str(path))
    
    metadata = {}
    
    if len(lf.image_list) == 0:
        return {
            'elsize_x': None,
            'elsize_y': None,
            'elsize_z': None,
            'time_interval': None
        }
    
    img_series = lf.get_image(0)
    
    # Try to extract pixel sizes from scale attributes
    # readlif may expose these as scale_x, scale_y, scale_z (in meters) or similar
    if hasattr(img_series, 'scale'):
        scale = img_series.scale
        if isinstance(scale, (list, tuple)) and len(scale) >= 3:
            # Convert from meters to micrometers
            metadata['elsize_x'] = float(scale[0]) * 1e6 if scale[0] else None
            metadata['elsize_y'] = float(scale[1]) * 1e6 if scale[1] else None
            metadata['elsize_z'] = float(scale[2]) * 1e6 if scale[2] else None
    
    # Try individual scale attributes
    for attr, key, factor in [
        ('scale_x', 'elsize_x', 1e6),  # meters to µm
        ('scale_y', 'elsize_y', 1e6),
        ('scale_z', 'elsize_z', 1e6),
        ('scale_t', 'time_interval', 1),  # already in seconds
    ]:
        if key not in metadata and hasattr(img_series, attr):
            val = getattr(img_series, attr)
            if val is not None and val != 0:
                metadata[key] = float(val) * factor
    
    # Try info dict if available
    if hasattr(img_series, 'info') and isinstance(img_series.info, dict):
        info = img_series.info
        for src_key, dst_key in [
            ('scale_x', 'elsize_x'),
            ('scale_y', 'elsize_y'),
            ('scale_z', 'elsize_z'),
        ]:
            if dst_key not in metadata and src_key in info:
                val = info[src_key]
                if val is not None and val != 0:
                    metadata[dst_key] = float(val) * 1e6
    
    # Set defaults for missing values (user can override in CSV)
    metadata.setdefault('elsize_x', None)
    metadata.setdefault('elsize_y', None)
    metadata.setdefault('elsize_z', None)
    metadata.setdefault('time_interval', None)
    
    return metadata


def load_elsizes_liff(path):
    """
    Load element sizes (pixel distances) from a .lif/.liff file.
    
    Returns a dict with keys: x, y, z (in micrometers), t (time interval), 
    elsize_unit, time_interval_unit.
    
    Note: If metadata cannot be extracted from the file, values will be None
    and should be specified manually in the metadata CSV.
    """
    metadata = load_liff_metadata(path)
    
    elsizes = {
        "x": metadata.get('elsize_x'),
        "y": metadata.get('elsize_y'),
        "z": metadata.get('elsize_z'),
        "elsize_unit": "µm",
        "t": metadata.get('time_interval'),
        "time_interval_unit": "s"
    }
    return elsizes
