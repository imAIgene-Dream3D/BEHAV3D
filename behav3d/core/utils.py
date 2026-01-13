from datetime import datetime
import fnmatch
import numpy as np

def format_time(
    start_time,
    end_time
):
    elapsed_time = end_time - start_time
    hours = int(elapsed_time // 3600)
    minutes = int((elapsed_time % 3600) // 60)
    seconds = int(elapsed_time % 60)
    return(hours, minutes, seconds)

def convert_time(time_interval, time_unit):
        time_conversions={
            "s": 3600,
            "m": 60,
            "h": 1
        }
        assert time_unit in time_conversions.keys(), f"time unit needs to be one of: {time_conversions.keys()}, is {time_unit}"
        time_interval = time_interval/time_conversions[time_unit]
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
