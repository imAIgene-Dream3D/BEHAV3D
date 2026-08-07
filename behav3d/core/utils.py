from datetime import datetime
import fnmatch
import inspect
import shutil
import numpy as np

def ignore_missing_rmtree_error(func, path, exc):
    err = exc[1] if isinstance(exc, tuple) else exc
    if not isinstance(err, FileNotFoundError):
        raise err

# shutil.rmtree() only gained the `onexc` keyword in Python 3.12; on older
# versions (e.g. the pinned 3.11 conda env) it must be called as `onerror`
# instead. ignore_missing_rmtree_error() already handles both calling
# conventions, so pick whichever keyword this interpreter supports.
_RMTREE_SUPPORTS_ONEXC = "onexc" in inspect.signature(shutil.rmtree).parameters

def rmtree_ignore_missing(path):
    """Remove a directory tree, ignoring races where a file is already gone."""
    if _RMTREE_SUPPORTS_ONEXC:
        shutil.rmtree(path, onexc=ignore_missing_rmtree_error)
    else:
        shutil.rmtree(path, onerror=ignore_missing_rmtree_error)

def format_time(
    start_time,
    end_time
):
    elapsed_time = end_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)
    return(hours, minutes, seconds)

def convert_time(time_interval, time_unit, convert_to="h"):
        time_conversions = {
            "s": 1,
            "m": 60,
            "h": 3600,
        }
        assert time_unit in time_conversions.keys(), f"time unit needs to be one of: {time_conversions.keys()}, is {time_unit}"
        assert convert_to in time_conversions.keys(), f"convert_to needs to be one of: {time_conversions.keys()}, is {convert_to}"
        time_interval = time_interval * time_conversions[time_unit]
        time_interval = time_interval / time_conversions[convert_to]
        return(time_interval)

def get_current_time():
    return(datetime.now().strftime("%H:%M:%S"))

def convert_distance(distance, distance_unit):
        distance_conversions={
            "nm":1000,
            "μm":1,
            "um":1,
            "mm":0.001
        }
        assert distance_unit in list(distance_conversions.keys()), f"distance unit needs to be one of: {list(distance_conversions.keys())}, is {distance_unit}"
        distance = distance/distance_conversions[distance_unit]
        return(distance)
      
def unit_conversion_factor(kind, pixel_distance_xy, pixel_distance_z=None):
    """Native→physical multiplier for a parameter of the given ``kind``.

    Segmentation parameters are stored in native pixel/voxel units; the GUI
    can display them in physical units (µm). This returns the factor ``f``
    such that ``physical = native * f``:

    - ``"distance"`` (lateral distance, e.g. EDT threshold, search radius):
      ``f = pixel_distance_xy`` (µm per pixel, XY).
    - ``"volume"`` (e.g. minimum segment size in voxels):
      ``f = pixel_distance_xy**2 * pixel_distance_z`` (µm³ per voxel).
    - anything else (dimensionless / fixed-by-definition): ``f = 1.0``.
    """
    xy = float(pixel_distance_xy)
    if kind == "distance":
        return xy
    if kind == "volume":
        z = float(pixel_distance_z if pixel_distance_z is not None else xy)
        return xy * xy * z
    return 1.0


def native_to_physical(native, kind, pixel_distance_xy, pixel_distance_z=None):
    """Convert a native pixel/voxel value to physical units (µm / µm³)."""
    return float(native) * unit_conversion_factor(kind, pixel_distance_xy, pixel_distance_z)


def physical_to_native(physical, kind, pixel_distance_xy, pixel_distance_z=None):
    """Convert a physical-unit value (µm / µm³) back to native pixels/voxels."""
    f = unit_conversion_factor(kind, pixel_distance_xy, pixel_distance_z)
    if f == 0:
        return float(physical)
    return float(physical) / f


def resolution_from_metadata(metadata):
    """Return ``(xy, z, valid)`` pixel sizes (µm) from a metadata DataFrame.

    Uses the first sample row and assumes uniform resolution across samples.
    ``valid`` is False when the columns are missing/unusable, in which case
    callers should fall back to pixel units and disable the physical toggle.
    """
    try:
        if metadata is None or len(metadata) == 0:
            return 1.0, 1.0, False
        if "pixel_distance_xy" not in metadata.columns or "pixel_distance_z" not in metadata.columns:
            return 1.0, 1.0, False
        xy = float(metadata["pixel_distance_xy"].iloc[0])
        z = float(metadata["pixel_distance_z"].iloc[0])
        if not (np.isfinite(xy) and np.isfinite(z)) or xy <= 0 or z <= 0:
            return 1.0, 1.0, False
        return xy, z, True
    except Exception:
        return 1.0, 1.0, False


def hours_per_frame_from_metadata(metadata):
    """Return ``(hours_per_frame, valid)`` from a metadata DataFrame.

    Uses the first sample row's ``time_interval``/``time_unit`` columns and
    assumes a uniform frame interval across samples, mirroring
    :func:`resolution_from_metadata`. ``valid`` is False when the columns
    are missing/unusable, in which case callers should fall back to frame
    units and disable any hours display toggle.
    """
    try:
        if metadata is None or len(metadata) == 0:
            return 1.0, False
        if "time_interval" not in metadata.columns or "time_unit" not in metadata.columns:
            return 1.0, False
        interval = float(metadata["time_interval"].iloc[0])
        unit = str(metadata["time_unit"].iloc[0])
        if not np.isfinite(interval) or interval <= 0:
            return 1.0, False
        return float(convert_time(interval, unit, convert_to="h")), True
    except Exception:
        return 1.0, False


def minutes_per_frame_from_metadata(metadata):
    """Return ``(minutes_per_frame, valid)`` from a metadata DataFrame.

    Mirrors :func:`hours_per_frame_from_metadata` (same "first sample row, uniform interval"
    assumption) but converts to minutes instead of hours, since per-track contact-bout lengths are
    usually a handful of timepoints and hours would round most of them down to 0.
    """
    hours, valid = hours_per_frame_from_metadata(metadata)
    return (float(hours) * 60.0, True) if valid else (1.0, False)


_TIME_UNIT_DISPLAY_NAMES = {"s": "seconds", "m": "minutes", "h": "hours"}


def format_timepoints_as_time(n_timepoints, metadata):
    """Return ``"= <value> <unit>"`` for ``n_timepoints`` converted to real time.

    Uses the first sample row's ``time_interval``/``time_unit`` metadata
    columns (mirroring :func:`hours_per_frame_from_metadata`), without
    forcing a conversion to hours: the value stays in whatever unit the
    metadata specifies (e.g. 3 timepoints at time_interval=2, time_unit="m"
    -> "= 6 minutes"). Returns "" when metadata is missing/unusable.
    """
    try:
        if metadata is None or len(metadata) == 0:
            return ""
        if "time_interval" not in metadata.columns or "time_unit" not in metadata.columns:
            return ""
        interval = float(metadata["time_interval"].iloc[0])
        unit = str(metadata["time_unit"].iloc[0])
        if not np.isfinite(interval) or interval <= 0 or unit not in _TIME_UNIT_DISPLAY_NAMES:
            return ""
        total = float(n_timepoints) * interval
        return f"= {total:g} {_TIME_UNIT_DISPLAY_NAMES[unit]}"
    except Exception:
        return ""


def element_to_dict(element):
    """
    Convert an ElementTree Element object to a dictionary.
    """
    result = {}
    result.update(element.attrib)

    # Process the element's children recursively
    for child in element:
        # Recursively convert the child element to a dictionary
        child_dict = element_to_dict(child)
        # If the child tag already exists in the dictionary, convert it to a list
        if child.tag in result:
            if not isinstance(result[child.tag], list):
                result[child.tag] = [result[child.tag]]
            result[child.tag].append(child_dict)
        else:
            result[child.tag] = child_dict

    # If the element has text, add it to the dictionary after stripping whitespace
    if element.text:
        result[element.tag] = element.text.strip()

    def remove_namespace(tag):
        return tag.split('}', 1)[-1] if '}' in tag else tag
    
    def clean_dict(d):
        if isinstance(d, dict):
            return {remove_namespace(k): clean_dict(v) for k, v in d.items()}
        elif isinstance(d, list):
            return [clean_dict(i) for i in d]
        else:
            return d
    
    result = clean_dict(result)
    return result

def rel_elsize(elsize):
    """
    Calculate the relative element sizes for an array, smallest is set to 1
    """
    min_size=min(elsize)
    elsize_scaled= [round(x/min_size, 2) for x in elsize]
    return(elsize_scaled)

def expand_column_patterns(selected, available_columns):
    """
    Expand glob-style patterns (e.g. 'mean_intensity_*') against available column names.
    - Accepts 'selected' as either a single string or an iterable of strings.
    - Works when 'available_columns' is a pandas.Index, list, tuple, etc.
    - Deduplicates while preserving order.
    """
    # Normalize selected -> list[str]
    if selected is None:
        selected_list = []
    elif isinstance(selected, str):
        selected_list = [selected]
    else:
        selected_list = list(selected)

    # Normalize available_columns -> list[str]
    if available_columns is None:
        avail = []
    else:
        # pandas.Index has .tolist()
        if hasattr(available_columns, "tolist"):
            avail = list(available_columns.tolist())
        else:
            avail = list(available_columns)

    # If nothing to match against, return input (deduped)
    if len(avail) == 0:
        # de-dup preserving order
        seen, out = set(), []
        for name in selected_list:
            if name not in seen:
                out.append(name); seen.add(name)
        return out

    # Expand patterns
    out = []
    for name in selected_list:
        name = str(name)
        if any(ch in name for ch in "*?["):
            matches = [c for c in avail if fnmatch.fnmatchcase(str(c), name)]
            out.extend(matches if matches else [name])  # keep pattern if nothing matched
        else:
            out.append(name)

    # De-dup while preserving order
    seen, uniq = set(), []
    for f in out:
        if f not in seen:
            uniq.append(f); seen.add(f)
    return uniq
