import pandas as pd
from pathlib import Path

def load_behav3d_metadata(
    metadata_path
    ):
    dtype_dict = {
        "sample_name": str,
        "organoid_line": str,
        "tcell_line": str,
        "exp_nr": int,
        "well": str,
        "tcell_channel": int,
        "live_channel": int,
        "dead_channel": int,
        "dead_dye_threshold": float,
        "contact_threshold": float,
        "pixel_distance_xy": float,
        "pixel_distance_z": float,
        "distance_unit": str,
        "time_interval": float,  # Assuming it could be a float
        "time_unit": str,
        "raw_image_path": str,  # Keeping as str for easy handling
        "tcell_segments_image_path": str,
        "tcell_tracks_image_path": str,
        "organoid_tracks_image_path": str,
        "tcell_tracks_csv_path": str
    }
    metadata = pd.read_csv(metadata_path, dtype=dtype_dict)
    return metadata 

def check_behav3d_metadata(
    metadata
    ):
    assert not any(metadata.drop(columns=["tcell_tracks_image_path", "organoid_tracks_image_path", "tcell_tracks_csv_path"]).isna().any()), "Some column values have not been supplied. Make sure you correctly supply values for all columns in the metadata .csv"
    ok = True
    for rowidx, sample_metadata in metadata.iterrows():
        sample_name = sample_metadata['sample_name']
        
        assert Path(sample_metadata["raw_image_path"]).exists(), f"The image_path supplied for 'row {rowidx+1}: {sample_name}' does not exist"
        
        if not pd.isna(sample_metadata["tcell_segments_image_path"]):
            assert Path(sample_metadata["tcell_segments_image_path"]).exists(), f"The tcell_segments_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        elif not pd.isna(sample_metadata["tcell_tracks_image_path"]):
            print(f"!!! No segmented or tracked tcell image is supplied for 'row {rowidx+1}: {sample_name}'. Please run segmentation and tracking below.")
            ok=False
        if not pd.isna(sample_metadata["tcell_tracks_image_path"]):
            assert Path(sample_metadata["tcell_tracks_image_path"]).exists(), f"The tcell_tracks_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        else:
            print(f"!!! No tracked tcell image is supplied for 'row {rowidx+1}: {sample_name}'. Please run tracking below.")
            ok=False
            
        if not pd.isna(sample_metadata["organoid_tracks_image_path"]):
            assert Path(sample_metadata["organoid_tracks_image_path"]).exists(), f"The organoid_tracks_image_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        else:
            print(f"!!! No tracked organoid image is supplied for 'row {rowidx+1}: {sample_name}'. Please run tracking below.")
            ok=False
            
        if not pd.isna(sample_metadata["tcell_tracks_csv_path"]):
            assert Path(sample_metadata["tcell_tracks_csv_path"]).exists(), f"The tcell_tracks_csv_path supplied for Row {rowidx+1} '{sample_name}' does not exist"
        else:
            print(f"!!! No .csv of tcell tracks is supplied for 'row {rowidx+1}: {sample_name}'. Please run tracking below.")
            ok=False
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