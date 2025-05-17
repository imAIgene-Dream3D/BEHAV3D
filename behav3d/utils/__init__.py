import pandas as pd
from pathlib import Path
from datetime import datetime

def predict_classes(args):
    clf_path, path, outpath, idx = args
    # clf = joblib.load(clf_path)
    return(1)

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
        "tcell_tracks_csv_path": str,
        "organoid_segments_image_path": str,
        "organoid_tracks_image_path": str,
        "organoid_tracks_csv_path": str,
        
    }
    metadata = pd.read_csv(metadata_path, dtype=dtype_dict)
    return metadata 

def check_behav3d_metadata(
    metadata
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
        "contact_threshold", 
        "pixel_distance_xy", 
        "pixel_distance_z", 
        "distance_unit", 
        "time_interval", 
        "time_unit",
        "raw_image_path",
        "tcell_segments_image_path", 
        "tcell_tracks_image_path", 
        "tcell_tracks_csv_path", 
        "organoid_segments_image_path",
        "organoid_tracks_image_path", 
        "organoid_tracks_csv_path"
    ]
    missing_columns = set(required_columns) - set(metadata.columns)
    assert all(col in metadata.columns for col in required_columns), f"Not all required columns are present in the metadata .csv file\n{'\n'.join(missing_columns)}"
    assert not any(metadata.drop(columns=[
        "tcell_segments_image_path", 
        "tcell_tracks_image_path", 
        "tcell_tracks_csv_path", 
        "organoid_segments_image_path",
        "organoid_tracks_image_path", 
        "organoid_tracks_csv_path"]).isna().any()), "Some column values have not been supplied. Make sure you correctly supply values for all columns in the metadata .csv"
    
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