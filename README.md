# BEHAV3D

## Installation

If required, install [miniforge](https://github.com/conda-forge/miniforge)

```
mamba env create --file env_*.yml
mamba activate behav3d
pip install -e <path/to/this/behav3/repository>
```


## How to run

BEHAV3D is run through a ipython notebook, in for example:
- Visual Code with an ipython extension
- Ipython notebook in a web browser

The script to run BEHAV3D: [run_behav3d.ipynb](./notebooks/run_behav3d.ipynb) (BEHAV3D/notebooks/run_behav3d.ipynb)

### Required files

[REQUIRED]
- A [metadata .csv](./configs/metadata.csv) that contains all samples that should be analyzed
- Per sample, a raw microscopy image (.czi, .tiff, .zarr)

[Optional]
- Already segmented data for each sample
- Already tracked data for each sample
  
## Metadata.csv structure

### Required
| Column Name                    | Explanation                                                                                         |
|-------------------------------|-----------------------------------------------------------------------------------------------------|
| sample_name                   | Name to assign to the sample; used for naming output files/folders.                                |
| organoid_line                 | Name of organoid cell line used in the experiment (e.g. 10T, 162M, etc.). Analysis visual results will be split on organoid lines                                     |
| tcell_line                    | Name or type of T cell used (e.g. TEG/CAR-T/etc.). Analysis visual results will be split on tcell lines                                 |
| exp_nr                        | Experiment number for tracking and reproducibility.                                                |
| well                          | Well identifier on the experimental plate (e.g. `well01`).                                          |
| tcell_channel                 | Index of the fluorescence channel representing T cells in the raw image. (Channel indexing starts at 0)                                      |
| live_channel                  | Index of the fluorescence channel representing alive cells in the raw image. (Channel indexing starts at 0)             |
| dead_channel                  | Index of the fluorescence channel representing dead signal in the raw image. (Channel indexing starts at 0)                                  |
| contact_threshold             | Distance threshold (in `distance_unit`) to define T cell-organoid or T cell - T cell contact.                        |
| pixel_distance_xy             | Physical pixel spacing in the XY plane (usually in micrometers per pixel).                         |
| pixel_distance_z              | Physical spacing between Z slices (z-step), in same unit as `pixel_distance_xy`.                   |
| distance_unit                 | Unit for all provided spatial distances (e.g. `μm` for micrometers).                                        |
| time_interval                 | Time between consecutive frames (usually in seconds).                                              |
| time_unit                     | Unit of the time interval (e.g. `s` for seconds, 'm' for minutes, 'h' for hours).                                                  |
| raw_image_path                | Full path to the original CZI or raw imaging file.                                                 |

### Optional
| Column Name                    | Explanation                                                                                         |
|-------------------------------|-----------------------------------------------------------------------------------------------------|
| tcell_segments_image_path     | Path to the segmented T cell image (e.g. .zarr file).                                              |
| tcell_tracks_image_path       | Path to the image showing T cell tracks over time.                                                 |
| tcell_tracks_csv_path         | Path to the CSV file containing T cell tracking data.                                              |
| organoid_segments_image_path  | Path to the segmented organoid image (e.g. .zarr file).                                            |
| organoid_tracks_image_path    | Path to the image showing organoid tracks over time.                                               |
| organoid_tracks_csv_path      | Path to the CSV file containing organoid tracking data.                                            |