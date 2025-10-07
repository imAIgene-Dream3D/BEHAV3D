import pandas as pd
from pathlib import Path
from datetime import datetime
import fnmatch

# def predict_classes(args):
#     clf_path, path, outpath, idx = args
#     # clf = joblib.load(clf_path)
#     return(1)

def load_behav3d_metadata(
    metadata_path
    ):
    dtype_dict = {
        "sample_name": str,
        "organoid_line": str,
        "tcell_line": str,
        "exp_nr": "Int64",
        "well": str,
        "tcell_channel": "Int64",
        "live_channel": "Int64",
        "dead_channel": "Int64",
        "dead_dye_threshold": float,
        "pixel_distance_xy": float,
        "pixel_distance_z": float,
        "distance_unit": str,
        "time_interval": float,  # Assuming it could be a float
        "time_unit": str,
        "dimension_order": str, 
        "raw_image_path": str,  # Keeping as str for easy handling
        "tcell_segments_image_path": str,
        "tcell_tracks_image_path": str,
        "tcell_tracks_csv_path": str,
        "organoid_segments_image_path": str,
        "organoid_tracks_image_path": str,
        "organoid_tracks_csv_path": str,
        "signal_unmixing_image_path": str, #FUNC
        "organoid_2_segments_image_path": str, #FUNC
        "organoid_2_tracks_image_path": str, #FUNC
        "organoid_2_tracks_csv_path": str, #FUNC
        "original_raw_image_path": str, # FUNC
    }

    metadata = pd.read_csv(metadata_path, dtype=dtype_dict)

    # Apply dtypes only to columns that exist --> Avoid ValueError, later we check_behav3d_metadata
    # for col, dtype in dtype_dict.items():
    #     if col in metadata.columns:
    #         metadata[col] = metadata[col].astype(dtype)

    metadata = metadata.dropna(how="all").reset_index(drop=True)
    return metadata 

def check_behav3d_metadata(
    metadata,
    func=False
    ):
    
    required_columns = [
        "sample_name", 
        "organoid_line", 
        "tcell_line", 
        "exp_nr", 
        "well", 
        "tcell_channel", 
        "live_channel", 
        "dead_channel", 
        "pixel_distance_xy", 
        "pixel_distance_z", 
        "distance_unit", 
        "time_interval", 
        "time_unit",
        "dimension_order", 
        "raw_image_path",
        "tcell_segments_image_path", 
        "tcell_tracks_image_path", 
        "tcell_tracks_csv_path", 
        "organoid_segments_image_path",
        "organoid_tracks_image_path", 
        "organoid_tracks_csv_path"
    ]

    func_columns = [
            "signal_unmixing_image_path",
            "organoid_2_segments_image_path",
            "organoid_2_tracks_image_path",
            "organoid_2_tracks_csv_path",
            # "original_raw_image_path"  # Not included bacouse it is added later if signal unmixing 
        ]

    missing_columns = set(required_columns) - set(metadata.columns)
    missing_string = '\n'.join(missing_columns)

    columns = [
        "tcell_segments_image_path", 
        "tcell_tracks_image_path", 
        "tcell_tracks_csv_path", 
        "organoid_segments_image_path",
        "organoid_tracks_image_path", 
        "organoid_tracks_csv_path"]
    
    if func:
        assert all(col in metadata.columns for col in required_columns+func_columns), f"Not all required columns are present in the metadata .csv file\n{missing_string}"
        assert not any(metadata.drop(columns=columns+func_columns).isna().any()), "Some column values have not been supplied. Make sure you correctly supply values for all columns in the metadata .csv"
    else:
        assert all(col in metadata.columns for col in required_columns), f"Not all required columns are present in the metadata .csv file\n{missing_string}"
        assert not any(metadata.drop(columns=columns).isna().any()), "Some column values have not been supplied. Make sure you correctly supply values for all columns in the metadata .csv"
    
    ok = True
    for rowidx, sample_metadata in metadata.iterrows():
        print(f"Row {rowidx+1}: {sample_metadata['sample_name']}")
        sample_name = sample_metadata['sample_name']
        
        assert Path(sample_metadata["raw_image_path"]).exists(), f"The image_path supplied for 'row {rowidx+1}: {sample_name}' does not exist"
        
        # elsizes = load_elsizes(sample_metadata["raw_image_path"])
        # assert sample_metadata["pixel_distance_xy"] == elsizes["x"], f"Pixel distance xy supplied for 'row {rowidx+1}: {sample_name}' does not match the x pixel distance {elsizes['x']} retrieved from the raw image metadata"
        # assert sample_metadata["pixel_distance_xy"] == elsizes["x"], f"Pixel distance xy supplied for 'row {rowidx+1}: {sample_name}' does not match the y pixel distance {elsizes['y']} retrieved from the raw image metadata"
        # assert sample_metadata["pixel_distance_z"] == elsizes["z"], f"Pixel distance z supplied for 'row {rowidx+1}: {sample_name}' does not match the z pixel distance {elsizes['z']} in the image"
        
        ### T cell paths
        if not pd.isna(sample_metadata["tcell_segments_image_path"]):
            assert Path(sample_metadata["tcell_segments_image_path"]).exists(), f"The tcell_segments_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        elif pd.isna(sample_metadata["tcell_tracks_image_path"]):
            print(f"!!! No segmented or tracked tcell image is supplied. Please run segmentation and tracking below.")
            ok=False
        if not pd.isna(sample_metadata["tcell_tracks_image_path"]):
            assert Path(sample_metadata["tcell_tracks_image_path"]).exists(), f"The tcell_tracks_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        else:
            print(f"!!! No tracked tcell image is supplied. Please run tracking below.")
            ok=False
        if not pd.isna(sample_metadata["tcell_tracks_csv_path"]):
            assert Path(sample_metadata["tcell_tracks_csv_path"]).exists(), f"The tcell_tracks_csv_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        else:
            print(f"!!! No tracked tcell .csv is supplied. Please run tracking below.")
            ok=False
            
        ### Organoids paths
        if not pd.isna(sample_metadata["organoid_segments_image_path"]):
            assert Path(sample_metadata["organoid_segments_image_path"]).exists(), f"The organoid_segments_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        elif pd.isna(sample_metadata["organoid_tracks_image_path"]):
            print(f"!!! No segmented or tracked organoid image is supplied. Please run segmentation and tracking below.")
            ok=False
        if not pd.isna(sample_metadata["organoid_tracks_image_path"]):
            assert Path(sample_metadata["organoid_tracks_image_path"]).exists(), f"The organoid_tracks_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        else:
            print(f"!!! No tracked organoid image is supplied. Please run tracking below.")
            ok=False
        if not pd.isna(sample_metadata["organoid_tracks_csv_path"]):
            assert Path(sample_metadata["organoid_tracks_csv_path"]).exists(), f"The organoid_tracks_csv_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        else:
            print(f"!!! No tracked organoid .csv is supplied. Please run tracking below.")
            ok=False

        ### FUNC
        if func:        
            ### Organoids 2 paths
            if not pd.isna(sample_metadata["organoid_2_segments_image_path"]):
                assert Path(sample_metadata["organoid_2_segments_image_path"]).exists(), f"The organoid_2_segments_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
            elif pd.isna(sample_metadata["organoid_2_tracks_image_path"]):
                print(f"!!! No segmented or tracked organoid_2 image is supplied. If 2 organoid types, please run segmentation and tracking below.")
                ok=False
            if not pd.isna(sample_metadata["organoid_2_tracks_image_path"]):
                assert Path(sample_metadata["organoid_2_tracks_image_path"]).exists(), f"The organoid_2_tracks_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
            else:
                print(f"!!! No tracked organoid_2 image is supplied. If 2 organoid types, please run tracking below.")
                ok=False
            if not pd.isna(sample_metadata["organoid_2_tracks_csv_path"]):
                assert Path(sample_metadata["organoid_2_tracks_csv_path"]).exists(), f"The organoid_2_tracks_csv_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
            else:
                print(f"!!! No tracked organoid_2 .csv is supplied.If 2 organoid types, please run tracking below.")
                ok=False
            ### Signal unmixing
            if not pd.isna(sample_metadata["signal_unmixing_image_path"]):
                print(f"--------{sample_metadata["signal_unmixing_image_path"]}--------------")
                assert Path(sample_metadata["signal_unmixing_image_path"]).exists(), f"The signal_unmixing_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
            else:
                print(f"!!! No signal unmixed image is supplied. If including signal unmixing, please run it below.")
                ok=False

        print("")
    if ok:
        print("Metadata file is complete and correct")

    
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