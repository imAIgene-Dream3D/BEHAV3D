import h5py
import numpy as np
import dask.array as da
from datetime import datetime, timedelta

def load_ims_as_numpy(path):
    """
    Loading .ims images.
    Always returns a 5-D array in TCZYX order, even when T=1 or C=1.
    """
    imsfile = h5py.File(name=path, mode="r")

    nr_timepoints = len([x for x in imsfile['/DataSet/ResolutionLevel 0/']])
    nr_channels = len([x for x in imsfile['/DataSet/ResolutionLevel 0/TimePoint 0']])

    image_channels = []
    for channel in range(nr_channels):
        image_timepoints = []
        for timepoint in range(nr_timepoints):
            time_img = imsfile[f'/DataSet/ResolutionLevel 0/TimePoint {timepoint}/Channel {channel}/Data'][:]
            image_timepoints.append(time_img)
        image_channels.append(np.stack(image_timepoints, axis=0))  # (T, Z, Y, X)

    img = np.stack(image_channels, axis=0)  # (C, T, Z, Y, X)
    img = np.transpose(img, (1, 0, 2, 3, 4))  # → (T, C, Z, Y, X)
    return img

def load_ims_timepoint_czyx(path, t):
    """
    Load a single timepoint from an IMS file.
    Returns array in CZYX order (C, Z, Y, X) without loading other timepoints.
    """
    with h5py.File(str(path), "r") as f:
        g0 = f["/DataSet/ResolutionLevel 0"]
        chan_keys = sorted(
            [k for k in g0[f"TimePoint {t}"].keys() if k.startswith("Channel ")],
            key=lambda s: int(s.split()[-1])
        )
        c_arrays = [g0[f"TimePoint {t}/{ck}/Data"][:] for ck in chan_keys]
    return np.stack(c_arrays, axis=0)  # (C, Z, Y, X)


def load_ims(path, as_numpy=False):
    f = h5py.File(path, "r")  # keep open while dask graph exists
    g0 = f["/DataSet/ResolutionLevel 0"]

    # Find all timepoints and channels
    time_keys = sorted(
        [k for k in g0.keys() if k.startswith("TimePoint ")],
        key=lambda s: int(s.split()[-1])
    )
    chan_keys = sorted(
        [k for k in g0["TimePoint 0"].keys() if k.startswith("Channel ")],
        key=lambda s: int(s.split()[-1])
    )

    # Look at first dataset to get full ZYX shape
    sample = g0[f"{time_keys[0]}/{chan_keys[0]}/Data"]
    chunks_spatial = sample.shape  # (Z, Y, X)

    t_blocks = []
    for t in time_keys:
        c_blocks = []
        for c in chan_keys:
            ds = g0[f"{t}/{c}/Data"]  # h5py.Dataset
            a = da.from_array(
                ds,
                chunks=chunks_spatial,  # entire ZYX volume in one chunk
                lock=True,
                asarray=False,
            )
            a = a[None, None, ...]  # -> (1,1,Z,Y,X)
            c_blocks.append(a)
        ccat = da.concatenate(c_blocks, axis=1)  # (1,C,Z,Y,X)
        t_blocks.append(ccat)

    darr = da.concatenate(t_blocks, axis=0)      # (T,C,Z,Y,X)
    darr = darr.rechunk((1,1) + chunks_spatial)  # enforce chunking
    
    if as_numpy:
        darr = darr.compute()
    return darr

def load_ims_metadata(path):
    with h5py.File(path, "r") as f:
        g0 = f["/DataSet/ResolutionLevel 0"]

        time_keys = sorted([k for k in g0.keys() if k.startswith("TimePoint ")],
                        key=lambda s: int(s.split()[-1]))
        chan_keys = sorted([k for k in g0[time_keys[0]].keys() if k.startswith("Channel ")],
                        key=lambda s: int(s.split()[-1]))

        ds = g0[f"{time_keys[0]}/{chan_keys[0]}/Data"]
        z, y, x = ds.shape

        t = len(time_keys)
        c = len(chan_keys)

        return (t, c, z, y, x)


def get_ims_dimension_order(path=None):
    """
    Return the axis order for IMS files.

    IMS (Imaris) HDF5 always stores data as TimePoint N / Channel N / Data,
    so the order is always TCZYX — including singletons at size 1.
    """
    return "TCZYX"


def get_ims_shape(path):
    """
    Return the shape of an IMS file as (T, C, Z, Y, X) without loading data.
    """
    return load_ims_metadata(path)

def save_as_imaris(
    img,
    outpath,
    voxel_size_x,
    voxel_size_y,
    voxel_size_z,
    time_interval_seconds,
    ImarisFileConverter_dll_path="C:/Program Files/Bitplane/ImarisFileConverter 10.0.0/bpImarisWriter100.dll",
    number_of_threads=1
    ):
    # Note: PW (PyImarisWriter) needs to be available in the environment
    try:
        from PyImarisWriter import PyImarisWriter as PW
    except ImportError:
        raise ImportError("PyImarisWriter is required for save_as_imaris but not found.")

    class ImageConverter(PW.ImageConverter):
        def __init__(
            self,
            datatype : str,
            image_size : PW.ImageSize,
            sample_size : PW.ImageSize,
            dimension_sequence : PW.DimensionSequence,
            block_size: PW.ImageSize,
            output_filename : str,
            options : PW.Options,
            application_name : str,
            application_version : str,
            progress_callback_class: PW.CallbackClass,
            ImarisFileConverter_dll_path = "/Applications/ImarisFileConverter 9.7.2.app/Contents/Frameworks/libbpImarisWriter.9.7.dylib",
            ):

            self.ImarisFileConverter_dll_path = ImarisFileConverter_dll_path
            super().__init__(
                output_filename = output_filename, 
                image_size = image_size,
                sample_size = sample_size,
                dimension_sequence = dimension_sequence,
                block_size = block_size,
                application_name = application_name,
                application_version=application_version,
                options=options,
                datatype=datatype,
                progress_callback_class = progress_callback_class
                )
        

        def _get_dll_filename(self):
            return(self.ImarisFileConverter_dll_path)
                
    import os
    from pathlib import Path

    outdir = outpath.parent
    os.chdir(outdir)
    outname = outpath.name
    
    img = np.ascontiguousarray(img).copy()
    image_shape = img.shape
    image_size = PW.ImageSize(
        x=image_shape[-1], 
        y=image_shape[-2], 
        z=image_shape[-3], 
        c=image_shape[-4], 
        t=image_shape[-5]
        )


    dimension_sequence = PW.DimensionSequence('x', 'y', 'z', 'c', 't')
    block_size = image_size
    sample_size = PW.ImageSize(x=1, y=1, z=1, c=1, t=1)
        
    application_name = 'PyImarisWriter'
    application_version = '1.0.0'
    callback_class = PW.CallbackClass()

    options = PW.Options()
    options.mNumberOfThreads = number_of_threads
            
    converter = ImageConverter(
        ImarisFileConverter_dll_path = ImarisFileConverter_dll_path,
        output_filename = outname, 
        image_size = image_size,
        sample_size = sample_size,
        dimension_sequence = dimension_sequence,
        block_size = block_size,
        application_name = application_name,
        application_version=application_version,
        options=options,
        datatype='uint16',
        progress_callback_class = callback_class
        )

    num_blocks = image_size / block_size
    block_index = PW.ImageSize()
    for c in range(num_blocks.c):
        block_index.c = c
        for t in range(num_blocks.t):
            block_index.t = t
            for z in range(num_blocks.z):
                block_index.z = z
                for y in range(num_blocks.y):
                    block_index.y = y
                    for x in range(num_blocks.x):
                        block_index.x = x
                        if converter.NeedCopyBlock(block_index):
                            converter.CopyBlock(img, block_index)
                            
    adjust_color_range = True
    extent_x = image_size.x * voxel_size_x
    extent_y = image_size.y * voxel_size_y
    extent_z = image_size.z * voxel_size_z
    image_extents = PW.ImageExtents(0, 0, 0, extent_x, extent_y, extent_z)
    parameters = PW.Parameters()
    now = datetime.today()
    time_infos = [now + timedelta(seconds=t*time_interval_seconds) for t in range(img.shape[0])]
    color_infos = [PW.ColorInfo() for _ in range(image_size.c)]
    converter.Finish(image_extents, parameters, time_infos, color_infos, adjust_color_range)
    converter.Destroy()
