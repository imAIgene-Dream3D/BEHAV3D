from aicspylibczi import CziFile
from behav3d.core.utils import element_to_dict
import numpy as np

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
    except (KeyError, TypeError, IndexError):
        try:
            time_interval = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["SubDimensionSetups"]["TimeSeriesSetup"]["Interval"]["TimeSpan"]["Value"]["Value"])
            interval_unit = metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"]["SubDimensionSetups"]["TimeSeriesSetup"]["Interval"]["TimeSpan"]["DefaultUnitFormat"]["DefaultUnitFormat"]

        except (KeyError, TypeError, IndexError):
            time_interval = float(metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"][0]["SubDimensionSetups"]["TimeSeriesSetup"]["Interval"]["TimeSpan"]["Value"]["Value"])
            interval_unit = metadata["Metadata"]["Experiment"]["ExperimentBlocks"]["AcquisitionBlock"][0]["SubDimensionSetups"]["TimeSeriesSetup"]["Interval"]["TimeSpan"]["DefaultUnitFormat"]["DefaultUnitFormat"]

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

def get_czi_shape(path, take_dims="TCZYX"):
    """
    Get the shape of a CZI image.
    Args:
    """
    czifile = CziFile(path)
    dim_order = czifile.dims
    shape = czifile.get_dims_shape()[0]
    size_by_dim = {ax: int(sz[1]) for ax, sz in shape.items()}

    shape = [size_by_dim.get(ax, 1) for ax in dim_order if ax in take_dims]
    return(shape)
