import pandas as pd
from pathlib import Path
from datetime import datetime
import fnmatch
import re

# def predict_classes(args):
#     clf_path, path, outpath, idx = args
#     # clf = joblib.load(clf_path)
#     return(1)

def detect_organoid_types_from_metadata(metadata):
    """
    Detect organoid types from metadata column names.
    Looks for columns matching pattern: or_{organoidtype}_line_condition, or_{organoidtype}_segments_image_path, etc.
    Returns list of organoid type names (e.g., ['organoid1', 'organoid2'])
    """
    if metadata is None:
        return []
    
    organoid_types = set()
    # Pattern: columns starting with 'or_' prefix
    for col in metadata.columns:
        if col.startswith('or_'):
            # Extract cell type name: or_organoid1_line_condition -> organoid1
            parts = col[3:].split('_', 1)  # Remove 'or_' prefix and split
            if parts:
                cell_type = parts[0]
                organoid_types.add(cell_type)
    
    return sorted(list(organoid_types))

def detect_immune_cell_types_from_metadata(metadata):
    """
    Detect immune cell types from metadata column names.
    Looks for columns matching pattern: im_{celltype}_line_condition, im_{celltype}_segments_image_path, etc.
    Returns list of immune cell type names (e.g., ['tcell', 'macro', 'nk'])
    """
    if metadata is None:
        return []
    
    immune_types = set()
    # Pattern: columns starting with 'im_' prefix
    for col in metadata.columns:
        if col.startswith('im_'):
            # Extract cell type name: im_tcell_line_condition -> tcell
            parts = col[3:].split('_', 1)  # Remove 'im_' prefix and split
            if parts:
                cell_type = parts[0]
                immune_types.add(cell_type)
    
    return sorted(list(immune_types))

def detect_other_cell_types_from_metadata(metadata):
    """
    Detect 'other' cell types from metadata column names.
    Looks for columns matching pattern: ot_{celltype}_line_condition, ot_{celltype}_segments_image_path, etc.
    Returns list of other cell type names (e.g., ['tumor1', 'fibroblast'])
    """
    if metadata is None:
        return []
    
    other_types = set()
    # Pattern: columns starting with 'ot_' prefix
    for col in metadata.columns:
        if col.startswith('ot_'):
            # Extract cell type name: ot_tumor1_line_condition -> tumor1
            parts = col[3:].split('_', 1)  # Remove 'ot_' prefix and split
            if parts:
                cell_type = parts[0]
                other_types.add(cell_type)
    
    return sorted(list(other_types))

def has_dead_channel(metadata):
    """
    Check if dead_channel column exists and has non-null values in metadata.
    Returns True if dead channel is present, False otherwise.
    """
    if metadata is None:
        return False
    
    if 'dead_channel' not in metadata.columns:
        return False
    
    # Check if any row has a non-null dead_channel value
    return metadata['dead_channel'].notna().any()

def load_behav3d_metadata(
    metadata_path
    ):
    """
    Load metadata CSV with dynamic cell type support.
    Uses prefixed columns: or_*, im_*, ot_* for organoid, immune, and other cell types.
    """
    # Basic dtype dict for common columns (cell-type specific columns are dynamic)
    dtype_dict = {
        "sample_name": str,
        "exp_nr": "Int64",
        "well": str,
        "dead_channel": "Int64",
        "pixel_distance_xy": float,
        "pixel_distance_z": float,
        "distance_unit": str,
        "time_interval": float,
        "time_unit": str,
        "dimension_order": str, 
        "raw_image_path": str,
        "signal_unmixing_image_path": str, #FUNC
        "original_raw_image_path": str, # FUNC
    }

    metadata = pd.read_csv(metadata_path)
    
    # Apply dtypes only to columns that exist
    for col, dtype in dtype_dict.items():
        if col in metadata.columns:
            try:
                metadata[col] = metadata[col].astype(dtype)
            except:
                pass  # Skip if conversion fails

    metadata = metadata.dropna(how="all").reset_index(drop=True)
    return metadata 

def check_behav3d_metadata(
    metadata,
    func=False
    ):
    """
    Validate metadata with dynamic cell type support.
    Uses prefixed columns: or_*, im_*, ot_* for organoid, immune, and other cell types.
    """
    # Detect cell types dynamically
    organoid_types = detect_organoid_types_from_metadata(metadata)
    immune_types = detect_immune_cell_types_from_metadata(metadata)
    other_types = detect_other_cell_types_from_metadata(metadata)
    
    # Basic required columns (always needed)
    required_columns = [
        "sample_name", 
        "exp_nr", 
        "well", 
        "pixel_distance_xy", 
        "pixel_distance_z", 
        "distance_unit", 
        "time_interval", 
        "time_unit",
        "dimension_order", 
        "raw_image_path",
    ]
    
    # Build dynamic required columns for each detected cell type
    path_columns = []  # Columns that can be empty (filled by pipeline)
    
    for org_type in organoid_types:
        required_columns.append(f"or_{org_type}_line_condition")
        path_columns.extend([
            f"or_{org_type}_segments_image_path",
            f"or_{org_type}_tracks_image_path",
            f"or_{org_type}_tracks_csv_path"
        ])
    
    for immune_type in immune_types:
        required_columns.append(f"im_{immune_type}_line_condition")
        path_columns.extend([
            f"im_{immune_type}_segments_image_path",
            f"im_{immune_type}_tracks_image_path",
            f"im_{immune_type}_tracks_csv_path"
        ])
    
    for other_type in other_types:
        required_columns.append(f"ot_{other_type}_line_condition")
        path_columns.extend([
            f"ot_{other_type}_segments_image_path",
            f"ot_{other_type}_tracks_image_path",
            f"ot_{other_type}_tracks_csv_path"
        ])
    
    # Optional func columns
    func_columns = []
    if func:
        if "signal_unmixing_image_path" in metadata.columns:
            func_columns.append("signal_unmixing_image_path")
        if "original_raw_image_path" in metadata.columns:
            func_columns.append("original_raw_image_path")
    
    # Check required columns exist
    missing_columns = set(required_columns) - set(metadata.columns)
    if missing_columns:
        missing_string = '\n'.join(missing_columns)
        raise AssertionError(f"Not all required columns are present in the metadata .csv file\n{missing_string}")
    
    # Check that non-path columns are filled
    columns_to_check = [col for col in required_columns if col not in path_columns]
    for col in columns_to_check:
        if metadata[col].isna().any():
            raise AssertionError(f"Column '{col}' has missing values. All required fields must be filled.")
    
    # Validate paths per row
    ok = True
    all_cell_types = [(f"or_{t}", t) for t in organoid_types] + \
                     [(f"im_{t}", t) for t in immune_types] + \
                     [(f"ot_{t}", t) for t in other_types]
    
    for rowidx, sample_metadata in metadata.iterrows():
        print(f"Row {rowidx+1}: {sample_metadata['sample_name']}")
        sample_name = sample_metadata['sample_name']
        
        # Check raw image exists
        assert Path(sample_metadata["raw_image_path"]).exists(), \
            f"The raw_image_path supplied for 'row {rowidx+1}: {sample_name}' does not exist"
        
        # Check each cell type's paths
        for prefix_type, display_type in all_cell_types:
            segments_col = f"{prefix_type}_segments_image_path"
            tracks_img_col = f"{prefix_type}_tracks_image_path"
            tracks_csv_col = f"{prefix_type}_tracks_csv_path"
            
            # Check segments
            if not pd.isna(sample_metadata[segments_col]) and sample_metadata[segments_col]:
                if not Path(sample_metadata[segments_col]).exists():
                    print(f"⚠️ {segments_col} path does not exist")
            else:
                print(f"!!! No segmented {display_type} image. Please run segmentation below.")
                ok = False
            
            # Check tracks image
            if not pd.isna(sample_metadata[tracks_img_col]) and sample_metadata[tracks_img_col]:
                if not Path(sample_metadata[tracks_img_col]).exists():
                    print(f"⚠️ {tracks_img_col} path does not exist")
            else:
                print(f"!!! No tracked {display_type} image. Please run tracking below.")
                ok = False
            
            # Check tracks CSV
            if not pd.isna(sample_metadata[tracks_csv_col]) and sample_metadata[tracks_csv_col]:
                if not Path(sample_metadata[tracks_csv_col]).exists():
                    print(f"⚠️ {tracks_csv_col} path does not exist")
            else:
                print(f"!!! No tracked {display_type} CSV. Please run tracking below.")
                ok = False
        
        # Check func columns if needed
        if func:
            if "signal_unmixing_image_path" in func_columns:
                if not pd.isna(sample_metadata["signal_unmixing_image_path"]):
                    assert Path(sample_metadata["signal_unmixing_image_path"]).exists(), \
                        f"signal_unmixing_image_path for Row {rowidx+1} '{sample_name}' does not exist"
                else:
                    print(f"!!! No signal unmixed image supplied.")
                    ok = False
        
        print("")
    
    if ok:
        print("✅ Metadata file is complete and correct")
    else:
        print("⚠️ Some processing steps need to be run (see warnings above)")

    
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