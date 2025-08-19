import time
from behav3d.utils.fileio import convert_file_to_zarr
from pathlib import Path

def convert_input_files_to_zarr(
    output_dir,
    metadata,
    ):
    
    for idx, sample in metadata.iterrows():
        print(f"Processing sample: {sample['sample_name']}") 
        start_time = time.time()
        
        sample_name = sample['sample_name']
        raw_image_path = Path(sample['raw_image_path'])
        raw_image_zarr =  Path(output_dir, "images", sample_name, f"{sample_name}.zarr")
        
        # organoid_segments_path = sample_metadata['organoid_tracks_image_path']
        # tcell_segments_path = sample_metadata['tcell_tracks_image_path']
        # organoid_tracks_path = sample_metadata['organoid_tracks_image_path']
        # tcell_tracks_path = sample_metadata['tcell_tracks_image_path']
        
        convert_file_to_zarr(
            path=raw_image_path, 
            outpath=raw_image_zarr, 
            overwrite=False
        )
                
        metadata.at[idx, "raw_image_path"] = str(raw_image_zarr)
    return(metadata)