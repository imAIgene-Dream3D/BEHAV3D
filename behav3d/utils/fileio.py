
# from PyImarisWriter import PyImarisWriter as PW
import numpy as np
from datetime import datetime, timedelta
import os
from pathlib import Path
from tifffile import imwrite, imread
from aicspylibczi import CziFile
import h5py
import zarr
import dask.array as da
import shutil
from time import sleep
from behav3d.utils import element_to_dict

def get_filepath_stem(path):
    path=Path(path)
    if path.suffix == ".zip" or path.suffix ==".bz2":
        path = Path(path.stem)
    return(path.stem)

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
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    return(metadata)

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

def load_tiff(path):
    """
    Loading .tif/.tiff images
    """
    img = imread(path)
    return(img)

def save_as_tiff(img, path):
    path = Path(path)
    assert path.parent.exists(), f"Parent folder of supplied .tiff path does not exist:\n{path}"
    imwrite(path, img)

def load_czi(path, t=None, z=None, c=None):
    """
    Loading .czi images
    """
    czifile=CziFile(path)
    # Load specific timepoints, z-slice, etc.
    img, _ = czifile.read_image()
    # img, _ = czifile.read_image(T=t, Z=z, C=c)

    # The shape of img contain a lot of singular dimensions
    # img.shape = (1, 1, 1, 1, 1, 1, 350, 3, 36, 200, 200)
    # np.squeeze removes preceding or trailing "empty" dimensions
    img = np.squeeze(img) # img.shape = (350, 3, 36, 200, 200)
    return(img)


def load_czi_metadata(path):
    """
    Loading .czi metadata
    """
    czifile=CziFile(path)
    metadata = element_to_dict(czifile.meta)
        
    try:
        elsize_x = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingX']["ScalingX"])*(10**6)
        elsize_y = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingY']["ScalingY"])*(10**6)
        elsize_z = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingZ']["ScalingZ"])*(10**6)
    except:
        elsize_x = float(metadata["Metadata"]["Scaling"]["Items"]["Distance"][0]["Value"]["Value"])*(10**6)
        elsize_y = float(metadata["Metadata"]["Scaling"]["Items"]["Distance"][1]["Value"]["Value"])*(10**6)
        elsize_z = float(metadata["Metadata"]["Scaling"]["Items"]["Distance"][2]["Value"]["Value"])*(10**6)
    
    #LSM880
    try:
        time_interval = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["TimeSeriesSetup"]["Switches"]["Switch"]["SwitchAction"]["SetIntervalAction"]["Interval"]["TimeSpan"]["Value"]["Value"])
        interval_unit = metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["TimeSeriesSetup"]["Switches"]["Switch"]["SwitchAction"]["SetIntervalAction"]["Interval"]["TimeSpan"]["DefaultUnitFormat"]["DefaultUnitFormat"]
        if interval_unit == "ms":
            time_interval = time_interval / 1000
    except:
        time_interval = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["SubDimensionSetups"]["TimeSeriesSetup"]["Interval"]["TimeSpan"]["Value"]["Value"])
        interval_unit = metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["SubDimensionSetups"]["TimeSeriesSetup"]["Interval"]["TimeSpan"]["DefaultUnitFormat"]["DefaultUnitFormat"]
    
    if interval_unit == "ms":
        time_interval = time_interval / 1000
    elif interval_unit == "min":
        time_interval = time_interval * 60
        
    metadata = {
        "elsize_x": elsize_x,
        "elsize_y": elsize_y,
        "elsize_z": elsize_z,
        "time_interval": time_interval,
    }
    return(metadata)

def load_h5(path, substruct):
    """
    Loading .h5 images
    """
    h5file = h5py.File(name=path, mode="r")
    img = h5file[substruct][:]
    return(img)

def load_ims_as_numpy(path):
    """
    Loading .ims images
    """
    imsfile = h5py.File(name=path, mode="r")

    nr_timepoints = len([x for x in imsfile['/DataSet/ResolutionLevel 0/']])
    nr_channels = len([x for x in imsfile['/DataSet/ResolutionLevel 0/TimePoint 0']])

    image_channels = []
    for channel in range(nr_channels):
        # print("Loading Channel {}".format(channel))
        if nr_timepoints>1:
            image_timepoints = []
            for timepoint in range(nr_timepoints):
                # Loads in single timepoint and single channel
                # print("Loading TimePoint {}".format(timepoint))
                time_img=imsfile[f'/DataSet/ResolutionLevel 0/TimePoint {timepoint}/Channel {channel}/Data'][:]
                image_timepoints.append(time_img)
                
            #Stack all timepoints to create single channel image over all timepoints
            image_timepoints = np.stack(image_timepoints, axis=0)
        else:
            image_timepoints = imsfile[f'/DataSet/ResolutionLevel 0/TimePoint 0/Channel {channel}/Data'][:]
            
        image_channels.append(image_timepoints)
    
    #Stack all channels to create full image
    img=np.stack(image_channels, axis=0)
    return(img)

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
    with h5py.File("image.ims", "r") as f:
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
    
def append_to_zarr(img, outpath):
        """
        Append a timepoint to an existing .zarr array
        If non-existent, create the .zarr array
        """
        outpath = Path(outpath)
        if not outpath.exists():
            zarr_file = zarr.open(
                outpath, 
                mode='w', 
                shape=(0,) + img.shape[1:], 
                chunks=(1,) + img.shape[1:], 
                dtype=img.dtype
                )
        else:
            zarr_file = zarr.open(outpath, mode='a')
        zarr_file.append(img)
        
def load_zarr(path, group=None, mode="r"):
    """
    Loading .zarr images
    """
    path = Path(path)
    assert path.exists(), f"Path to zarr file does not exist:\n{path}"
    if path.suffix==".zip":
        zarr_store = zarr.storage.ZipStore(path)
    else:
        zarr_store = zarr.storage.LocalStore(path)
        
    
    if group:
        # zarr_obj = zarr.open(zarr_store, mode=mode)[group]
        dask_img = da.from_zarr(zarr_store, component=group)
    else:
        dask_img = da.from_zarr(zarr_store)
    # dask_img = da.from_zarr(zarr_obj)
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

def zip_zarr(path):
    outpath = Path(f"{path}.zip")
    if outpath.exists():
        outpath.unlink()
    shutil.make_archive(path, "zip", path)
    shutil.rmtree(path)
    return(outpath)
    
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
    ImarisFileConverter_dll_path=str(ImarisFileConverter_dll_path)
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
                
    # Define the input and output file paths
    # input_czi_path = "D:/Sharing/personal/Sam/FUNC/FUNC_BV1_Exp020/FUNC_BV1_Exp020_Img1.czi"
    # input_czi_path = Path(input_czi_path)

    outdir = outpath.parent
    os.chdir(outdir)
    outname = outpath.name
    
    # Step 1: Read the CZI file
    # czi = CziFile(input_czi_path)
    # image_data, _ = czi.read_image() 
    # metadata = element_to_dict(czi.meta)
    # Extract img data as a NumPy array
    img = np.ascontiguousarray(img).copy()
    # image_data = np.squeeze(image_data)

    # czi.get_dims_shape()
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
    # options.mCompressionAlgorithmType = PW.eCompressionAlgorithmGzipLevel2
    # options.mEnableLogProgress = True
            
    # Step 3: Create an IMS file with ImarisWriter
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
        )  # Adjust data type as needed

    # voxel_size_x = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingX']["ScalingX"])*(10**6)
    # voxel_size_y = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingY']["ScalingY"])*(10**6)
    # voxel_size_z = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingZ']["ScalingZ"])*(10**6)

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
    tcell_segments_path,
    organoid_segments_path,
    raw_image_path,
    output_dir = None,
    chunks=None,
    overwrite=False
    ):
    
    tcell_segments_path = Path(tcell_segments_path)
    organoid_segments_path = Path(organoid_segments_path)
    raw_image_path = Path(raw_image_path)
    
    if output_dir is None:
        tcell_zarr_out_path = Path(tcell_segments_path, f"{sample_name}_tcell_tracked.zarr")
        organoid_zarr_out_path = Path(organoid_segments_path, f"{sample_name}_organoid_tracked.zarr")
        raw_image_zarr_out_path = Path(raw_image_path, f"{sample_name}.zarr")
    else:   
        tcell_zarr_out_path = Path(output_dir, f"{sample_name}_tcell_tracked.zarr")
        organoid_zarr_out_path = Path(output_dir, f"{sample_name}_organoid_tracked.zarr")
        raw_image_zarr_out_path = Path(output_dir, f"{sample_name}.zarr")
    
    convert_file_to_zarr(
            path=tcell_segments_path, 
            outpath=tcell_zarr_out_path, 
            chunks=chunks,
            overwrite=overwrite
            )

    convert_file_to_zarr(
            path=organoid_segments_path, 
            outpath=organoid_zarr_out_path, 
            chunks=chunks,
            overwrite=overwrite
            )
    
    convert_file_to_zarr(
            path=raw_image_path, 
            outpath=raw_image_zarr_out_path, 
            chunks=chunks,
            overwrite=overwrite
            )
    
    return(
        tcell_zarr_out_path,
        organoid_zarr_out_path,
        raw_image_zarr_out_path
    )
    
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
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
    
    return(elsizes)

def load_elsizes_czi(path):
    czi = CziFile(path)
    metadata = element_to_dict(czi.meta)
    try:
        elsize_x = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingX']["ScalingX"])*(10**6)
        elsize_y = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingY']["ScalingY"])*(10**6)
        elsize_z = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]['AcquisitionBlock']['AcquisitionModeSetup']['ScalingZ']["ScalingZ"])*(10**6)
    except:
        elsize_x = float(metadata["Metadata"]["Scaling"]["Items"]["Distance"][0]["Value"]["Value"])*(10**6)
        elsize_y = float(metadata["Metadata"]["Scaling"]["Items"]["Distance"][1]["Value"]["Value"])*(10**6)
        elsize_z = float(metadata["Metadata"]["Scaling"]["Items"]["Distance"][2]["Value"]["Value"])*(10**6)
    
    ### Interval unit does not make sense as it seems to be in seconds but "defaultunit is minutes"
    try:
        time_interval = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["TimeSeriesSetup"]["Switches"]["Switch"]["SwitchAction"]["SetIntervalAction"]["Interval"]["TimeSpan"]["Value"]["Value"])
        interval_unit = metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["TimeSeriesSetup"]["Switches"]["Switch"]["SwitchAction"]["SetIntervalAction"]["Interval"]["TimeSpan"]["DefaultUnitFormat"]["DefaultUnitFormat"]
    except:
        time_interval = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["SubDimensionSetups"]["TimeSeriesSetup"]["Interval"]["TimeSpan"]["Value"]["Value"])
        interval_unit = metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["SubDimensionSetups"]["TimeSeriesSetup"]["Interval"]["TimeSpan"]["DefaultUnitFormat"]["DefaultUnitFormat"]
    
    elsizes = {
        "x": elsize_x,
        "y": elsize_y,
        "z": elsize_z,
        "elsize_unit": "µm",
        "t": time_interval,
        # For some reason there doesnt seem 
        "time_interval_unit": "s"   
    }
    return(elsizes)


def get_image_shape(path):
    """
    Loads an image and prints its shape and axis sizes.
    Args:
        image_path: Path to the image file (zarr, tif, etc.)
    """
    path = Path(path)
    if path.suffix == ".zarr":
        # Load zarr image
        img = load_zarr(path)
        shape = img.shape
    elif path.suffix == ".czi":
        shape = get_czi_shape(path)
    shape = tuple(shape)
    return shape

def get_czi_shape(path, take_dims="TCZYX"):
    """
    Get the shape of a CZI image.
    Args:
    """
    path = Path(path)
    if path.suffix == ".czi":
        czifile = CziFile(path)
        dim_order = czifile.dims
        shape = czifile.get_dims_shape()[0]
        size_by_dim = {ax: int(sz[1]) for ax, sz in shape.items()}

        shape = [size_by_dim.get(ax, 1) for ax in dim_order if ax in take_dims]
        return(shape)
    else:
        print(f"{path} is not a .czi")
        
def get_image_dimension_order(path):
    """
    Loads an image and prints its shape and axis sizes.
    Args:
        image_path: Path to the image file (zarr, tif, etc.)
    """
    path = Path(path)
    if path.suffix == ".czi":
        czifile = CziFile(path)
        dim_order = czifile.dims[-5:]
    else:
        dim_order = "TCZYX"
    return dim_order
