from aicspylibczi import CziFile
from behav3d.core.utils import element_to_dict
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET

_DATA_AXES = set("TCZYX")


def _squeeze_wrapper_dims(img, czi_dims, keep="TCZYX"):
    """
    Remove only CZI wrapper dimensions (Scene, Block, Mosaic, etc.)
    while preserving all data axes in *keep*, even when their size is 1.
    """
    keep_set = set(keep)
    axes_to_squeeze = tuple(
        i for i, ax in enumerate(czi_dims) if ax not in keep_set
    )
    if axes_to_squeeze:
        img = np.squeeze(img, axis=axes_to_squeeze)
    return img


def load_czi(path, t=None, z=None, c=None):
    """
    Load a CZI image as-is.
    Wrapper dimensions (Scene, Block, etc.) are removed but all TCZYX data
    axes present in the file are preserved — even when their size is 1.
    The axis order matches the file's native data dimension order.
    """
    czifile = CziFile(path)
    read_kwargs = {}
    if t is not None:
        read_kwargs["T"] = t
    if z is not None:
        read_kwargs["Z"] = z
    if c is not None:
        read_kwargs["C"] = c

    if read_kwargs:
        img, _ = czifile.read_image(**read_kwargs)
    else:
        img, _ = czifile.read_image()

    img = _squeeze_wrapper_dims(img, czifile.dims, keep="TCZYX")
    return img


def load_czi_timepoint_czyx(path, t, source_order):
    """
    Load a single timepoint from a CZI file as-is.
    Wrapper dims are removed and the T axis is collapsed (since we loaded
    a specific timepoint).  The returned array keeps whatever CZYX subset
    the file actually has, in the file's native order.
    Singleton expansion and CZYX reorder happen in _load_czyx_for_timepoint.
    """
    czifile = CziFile(path)
    img, _ = czifile.read_image(T=t)

    img = _squeeze_wrapper_dims(img, czifile.dims, keep="TCZYX")

    remaining_order = "".join(ax for ax in czifile.dims if ax in _DATA_AXES)

    if "T" in remaining_order:
        t_pos = remaining_order.index("T")
        img = np.squeeze(img, axis=t_pos)

    return img

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

    # Prefer the actual measured interval (from real per-timepoint
    # AcquisitionTime stamps) over the nominal/configured SetIntervalAction
    # value: the configured interval is only a target and, especially for
    # multi-scene/multi-Z time-lapses, real cadence can run slower (or
    # otherwise drift) than requested. See get_actual_time_intervals_czi.
    try:
        actual = get_actual_time_intervals_czi(path)
        time_interval = actual["median_interval"]
    except Exception:
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


def get_actual_time_intervals_czi(path, s=0, z=0, c=0):
    """
    Measure the actual per-timepoint interval (in seconds) from the real
    AcquisitionTime recorded in each subblock's own metadata, as opposed to
    the nominal/configured interval reported by load_czi_metadata /
    load_elsizes_czi (SetIntervalAction). The configured interval is only a
    target: if a full T-cycle (all scenes/Z-planes/channels) takes longer to
    actually scan than the configured interval, ZEN does not wait for it, so
    the real cadence can run slower than requested, or otherwise drift.

    Fixes a single scene/Z-plane/channel (s, z, c) so exactly one real
    timestamp exists per timepoint, then diffs consecutive timepoints'
    timestamps to get the actual interval between two subsequent timepoints.

    Args:
        path: Path to the .czi file.
        s, z, c: Scene / Z-plane / Channel index to fix while T varies.
            Defaults to the first of each.
    Returns:
        dict with:
            "timepoints": sorted T indices (np.ndarray[int])
            "timestamps": actual acquisition datetime per timepoint (pd.Series)
            "intervals": actual interval in seconds between each pair of
                subsequent timepoints (np.ndarray[float], length len(timepoints)-1)
            "mean_interval", "median_interval": summary stats (float)
    """
    czifile = CziFile(path)
    subblocks = czifile.read_subblock_metadata(S=s, Z=z, C=c, unified_xml=False)

    records = []
    for dims, xml_str in subblocks:
        t = dims.get("T", 0)
        root = ET.fromstring(xml_str)
        acquisition_time = root.find(".//AcquisitionTime")
        if acquisition_time is None or acquisition_time.text is None:
            continue
        records.append((t, acquisition_time.text))

    if not records:
        raise ValueError(
            f"No AcquisitionTime found in subblock metadata for S={s}, Z={z}, C={c} in {path}."
        )

    timestamps = pd.DataFrame(records, columns=["t", "timestamp"]).drop_duplicates("t")
    timestamps["timestamp"] = pd.to_datetime(timestamps["timestamp"])
    timestamps = timestamps.sort_values("t").reset_index(drop=True)

    intervals = timestamps["timestamp"].diff().dt.total_seconds().to_numpy()[1:]

    return {
        "timepoints": timestamps["t"].to_numpy(),
        "timestamps": timestamps["timestamp"],
        "intervals": intervals,
        "mean_interval": float(np.mean(intervals)) if len(intervals) else float("nan"),
        "median_interval": float(np.median(intervals)) if len(intervals) else float("nan"),
    }

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
