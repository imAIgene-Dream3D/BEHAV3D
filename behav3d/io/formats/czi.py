from aicspylibczi import CziFile
# aicspylibczi has no attachment-reading API, so the "TimeStamps" attachment
# (see get_actual_time_intervals_czi) is read with czifile instead. Imported
# under a different name to avoid clashing with the `czifile` variable name
# used throughout this module for aicspylibczi.CziFile instances.
import czifile as _czifile_pkg
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
        time_interval = actual["mean_interval"]
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


def _actual_time_intervals_from_timestamps_attachment(path):
    """
    Read the file's "TimeStamps" attachment (CZTIMS segment), if present:
    a dedicated array with one real acquisition time (seconds elapsed since
    some fixed reference, not necessarily timepoint 0) per T-index. ZEN
    writes this for essentially all time-lapse experiments, and it is
    present even when per-subblock AcquisitionTime metadata is not --
    which is the common case for large multi-scene time-lapses, where ZEN
    omits per-subblock metadata entirely to save space.

    ZEN pre-allocates the attachment for the originally planned number of
    timepoints and zero-fills any that were never reached (e.g. an
    acquisition stopped early), so this keeps only the leading run of
    strictly increasing values and discards the rest as padding.

    Returns None (rather than raising) when there is no TimeStamps
    attachment, or fewer than 1 valid entry in it, so callers can fall back
    to another method.
    """
    try:
        with _czifile_pkg.CziFile(path) as reader:
            arr = None
            for attachment in reader.attachments():
                if attachment.attachment_entry.name == "TimeStamps":
                    arr = np.asarray(attachment.data(), dtype=float)
                    break
    except Exception:
        return None

    if arr is None or arr.size == 0:
        return None

    diffs = np.diff(arr)
    non_increasing = np.where(diffs <= 0)[0]
    n_valid = int(non_increasing[0] + 1) if len(non_increasing) else len(arr)
    if n_valid == 0:
        return None
    seconds = arr[:n_valid]
    seconds = seconds - seconds[0]

    intervals = np.diff(seconds)
    return {
        "timepoints": np.arange(len(seconds)),
        "timestamps": pd.Series(seconds),
        "intervals": intervals,
        "mean_interval": float(np.mean(intervals)) if len(intervals) else float("nan"),
        "median_interval": float(np.median(intervals)) if len(intervals) else float("nan"),
    }


def _actual_time_intervals_from_subblock_metadata(path, s=None, z=None, c=None):
    """
    Fallback for files with no TimeStamps attachment: measure the actual
    per-timepoint interval from the real AcquisitionTime recorded in each
    subblock's own metadata.

    Fixes a single scene/Z-plane/channel (s, z, c) so exactly one real
    timestamp exists per timepoint, then diffs consecutive timepoints'
    timestamps to get the actual interval between two subsequent timepoints.

    Args:
        path: Path to the .czi file.
        s, z, c: Scene / Z-plane / Channel index to fix while T varies.
            Defaults to the first valid index reported by the file for that
            axis (via get_dims_shape()), which is not necessarily 0 -- some
            multi-scene CZI files have dimension ranges that don't start at
            0. An axis absent from the file is left unconstrained.
    """
    czifile = CziFile(path)
    dim_shape = czifile.get_dims_shape()[0]

    def _resolve(axis, value):
        if value is not None:
            return value
        return int(dim_shape[axis][0]) if axis in dim_shape else None

    s = _resolve("S", s)
    z = _resolve("Z", z)
    c = _resolve("C", c)

    fix_kwargs = {k: v for k, v in (("S", s), ("Z", z), ("C", c)) if v is not None}
    subblocks = czifile.read_subblock_metadata(unified_xml=False, **fix_kwargs)

    records = []
    for dims, xml_str in subblocks:
        t = dims.get("T", 0)
        if not xml_str:
            continue
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


def get_actual_time_intervals_czi(path, s=None, z=None, c=None):
    """
    Measure the actual per-timepoint interval (in seconds) from real
    acquisition timing recorded in the file, as opposed to the
    nominal/configured interval reported by load_czi_metadata /
    load_elsizes_czi (SetIntervalAction). The configured interval is only a
    target: if a full T-cycle (all scenes/Z-planes/channels) takes longer to
    actually scan than the configured interval, ZEN does not wait for it, so
    the real cadence can run slower than requested, or otherwise drift.

    Tries two sources, in order:
        1. The file's "TimeStamps" attachment -- one real timestamp per
           T-index, independent of scene/Z/channel. Present for essentially
           all time-lapse acquisitions, including large multi-scene ones
           where per-subblock metadata is stripped (see
           _actual_time_intervals_from_timestamps_attachment).
        2. Per-subblock AcquisitionTime metadata for a fixed (s, z, c), as a
           fallback for files with no TimeStamps attachment (see
           _actual_time_intervals_from_subblock_metadata).

    Args:
        path: Path to the .czi file.
        s, z, c: Scene / Z-plane / Channel index to fix if falling back to
            per-subblock metadata (ignored when a TimeStamps attachment is
            found). See _actual_time_intervals_from_subblock_metadata.
    Returns:
        dict with:
            "timepoints": sorted T indices (np.ndarray[int])
            "timestamps": acquisition time per timepoint (pd.Series). From
                the TimeStamps attachment this is seconds elapsed since the
                first timepoint (float); from subblock metadata it is an
                actual datetime.
            "intervals": actual interval in seconds between each pair of
                subsequent timepoints (np.ndarray[float], length len(timepoints)-1)
            "mean_interval", "median_interval": summary stats (float)
    Raises:
        ValueError: if neither source yields any real per-timepoint timing.
    """
    result = _actual_time_intervals_from_timestamps_attachment(path)
    if result is not None:
        return result
    return _actual_time_intervals_from_subblock_metadata(path, s, z, c)

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
