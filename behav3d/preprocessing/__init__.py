import time
from behav3d.utils.fileio import convert_file_to_zarr, get_image_shape
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

        convert_file_to_zarr(
            path=raw_image_path, 
            outpath=raw_image_zarr, 
            overwrite=False
        )
                
        metadata.at[idx, "raw_image_path"] = str(raw_image_zarr)
    
    return(metadata)

def get_all_image_shapes(metadata):
    """
    Get the shapes of all images in the metadata.
    """
    shapes = []
    for _, row in metadata.iterrows():
        path = Path(row['raw_image_path'])
        if path.exists():
            shape = get_image_shape(path)
            shapes.append((row['sample_name'], shape))
        else:
            shapes.append((row['sample_name'], "File not found"))
    return shapes