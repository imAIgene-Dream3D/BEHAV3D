# BEHAV3D

## Installation

**Step 0** If required, install [miniforge](https://github.com/conda-forge/miniforge) or [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

https://github.com/conda-forge/miniforge/releases/latest/
https://conda-forge.org/download/



 - Mamba is a faster version of conda but both work the same
 - On Windows, isntallation with the mamba .exe might be blocked. If so, open a command terminal (cmd.exe) and run:
```
start /wait "" Miniforge3-Windows-x86_64.exe /InstallationType=JustMe /RegisterPython=0 /S /D=%UserProfile%\Miniforge3
```

**Step 1** Open a miniforge terminal or conda terminal by typing *miniforge* in your application search tab

**Step 2** Download the BEHAV3D repository from this github page, by clicking *Code* > *Download zip*. Unzip the file and place its content to a folder that from now on we call BEHAV3D reposotory.

**Step 3** Adapt location of this directory to your own one as indicated below and paste the following in the terminal to navigate to the BEHAV3D repository
```
cd <path/to/this/behav3/repository>
```
**Step 4**
This creates a python environment and installes required packages
Use the code below [env_mac.yml](./env_mac.yml) or [env_windows.yml](./env_windows.yml) dependent on what operating system you are on

*Mac install*:
In the terminal copy paste and enter the following:

```
mamba env create --file env_mac.yml
mamba activate behav3d
pip install -e .
python -m ipykernel install --user --name=behav3d --display-name "behav3d"
```

*Windows install*
In the terminal copy paste and enter the following:
```
mamba env create --file env_windows.yml
mamba activate behav3d
pip install -e .
python -m ipykernel install --user --name=behav3d --display-name "behav3d"
```

## How to run

BEHAV3D is run through a ipython notebook, in for example:

- Visual Code with an ipython extension (Recommended). Install VisualCode from here https://code.visualstudio.com/Download. Once installed open the program and go to extensions and install the ipython extension. Once this is done you can go to *select kernel* > *python* > *behav3d*     
  Open the following file that is in your BEHAV3D repository : [run_behav3d.ipynb](./notebooks/run_behav3d.ipynb) (BEHAV3D/notebooks/run_behav3d.ipynb)

- Ipython notebook in a web browser

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
