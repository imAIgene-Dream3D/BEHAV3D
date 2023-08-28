# BEHAV3D pipeline

## Important
***This code is currently under development***
_____

## Overview
BEHAV3D is dynamic immuno-organoid 3D imaging-transcriptomics platform to study tumor death dynamics; immune cell behavior and behavior-guided transcriptomics.
_____
## What type of data does BEHAV3D work with?
- Any type of multispectral time-lapse 3D (or 2D) imaging data, where objects such as tumor cells or tumor organoids are in co-culture with immune cells of interest.
_____
## What output can BEHAV3D provide?
- Behavioral clustering/classification of different types of cell dynamics
- Any type of change of cell state that can be detected by a change in fluorescent intensity e.g. cell death, reporter, Ca2+ signalling (In Progress)
- Tumor death dynamics quantification (In Progress)
- Backprojection of behavioral phenotype in Imaris 3D image visualization software (In Progress)
- Correlation between tumor death dynamics and behavioral phenotypes (In Progress)
_____
## How to cite this pipeline
Dekkers JF*, Alieva M*, Cleven A, Keramati F, Wezenaar AKL, van Vliet EJ, Puschhof J, Brazda P, Johanna I, Meringa AD, Rebel HG, Buchholz MB, Barrera Román M, Zeeman AL, de Blank S, Fasci D, Geurts MH, Cornel AM, Driehuis E, Millen R, Straetemans T, Nicolasen MJT, Aarts-Riemens T, Ariese HCR, Johnson HR, van Ineveld RL, Karaiskaki F, Kopper O, Bar-Ephraim YE, Kretzschmar K, Eggermont AMM, Nierkens S, Wehrens EJ, Stunnenberg HG, Clevers H, Kuball J, Sebestyen Z, Rios AC. **Uncovering the mode of action of engineered T cells in patient cancer organoids**. * *equal contibution* Nat Biotechnol. 2023 Jan https://doi.org/10.1038/s41587-022-01397-w
_____
## Software and Hardware requirements
BEHAV3D runs can be run from both command line or through the jupyter notebook and was tested on MacOS Big Sur with R version 4.1.1 and on WIndows 10 with R version 4.3.0 .

For segmentation, a computer with decent hardware is required dependent on the sizes of images that need to be processed:
- Ilastik segmentation of a timeseries with dimensions [4, 390, 36, 512, 512] (ctzyx) can be segmented in X minutes on a xxx workstation

For tracking, a computer with lower specs is often sufficiant:
- Tracking with TrackMate using Linear Assignment Problem tracking on a timeseries with dimensions [4, 390, 36, 512, 512] (ctzyx) and assigning TrackIDs back to the segmented image takes approx. 15 minutes

BEHAV3D runs in R studio, jupyter notebook or from command line and was tested on MacOS Big Sur with R version 4.1.1 and on Windows 10 with R version 4.3.0 .
_____
## Installation

```
git clone https://github.com/RiosGroup/BEHAV3D
cd ./BEHAV3D

install ilastik ()
conda create -n BEHAV3D python=3
conda activate BEHAV3D
pip install .
```

for open-source processing:
- Install Ilastik (https://www.ilastik.org/download.html)
- Install Fiji (https://imagej.net/software/fiji/downloads)

Set *ilastik_path* in your 'config.yml' file
``` 
For Mac:
.../ilastik-1.3.2-OSX.app/Contents/ilastik-release/run_ilastik.sh

For Linux:
.../ilastik-1.3.2-Linux/run_ilastik.sh

For Windows (untested):
.../Program Files/ilastik-1.3.2/ilastik.exe
```
_____

## Set-up
BEHAV3D contains several pipelines for the processing of time-lapse imaging. <br>

For all steps, BEHAV3D requires 2 files:<br>

### **config.yml** <br>
A yaml file that contains general (non-experimentally bound) parameters for BEHAV3D and paths to the used software.
An example version can be found in [...BEHAV3D/configs/config_template.yml](https://github.com/RiosGroup/BEHAV3D/blob/main/configs/config_template.yml) <br>
Explanation on what each variable changes is commented in that template

### **metadata.csv** <br>
To correctly import experimental data for analysis with BEHAV3D, it is required to fill in a .csv that contains information per experiment performed, requiring information on:<br>
For an example see: [...BEHAV3D/configs/metadata_template.tsv](https://github.com/RiosGroup/BEHAV3D/blob/main/configs/metadata_template.tsv)<br>

**Parameters that have to be specified:**
- **sample_name** <br>
The name of the sample which will be used for the naming of output files <br>
- **organoid_line** <br>
The name of the organoid line that has been used in the experiment <br>
- **tcell_line** <br>
The name of the T cell line that has been used in the experiment <br>
- **exp_nr** <br>
The number of the experiment <br>
- **well** <br>
The well in which the co-culture experiment has been performed
- **dead_dye_channel** <br>
The index of the channel that contains the dead dye for analysis of death mechanics<br>
- **dead_dye_threshold** <br>
The threshold of the mean dead dye channel inside a cell to register the cell as dead<br>
- **contact_threshold** <br>
The threshold distance to determine if cells are in contact with eachother, calculated from the borders of a cell <br>
- **pixel_distance_xy** <br>
The real-life distance of each pixel (often in µm) for the x and y dimension <br>
- **pixel_distance_z** <br>
The real-life distance of each pixel (often in µm) for the z dimension <br>
- **distance_unit** <br>
The unit of the distance parameters (e.g. µm) which will be used to noramlize each experiment to the same distance unit  <br>
- **time_interval** <br>
The real-life time inbetween each timepoint in the experiment  <br>
- **time_unit** <br>
The unit of the time parameters (e.g. m for minutes) which will be used to noramlize each experiment to the same time unit <br>
- **image_path** <br>
The path to the image containing the T cells, organoids and dead dyes, saved as a channel in .h5 format<br>
- **image_internal_path** <br>
The internal path of the .h5 file containing the image (e.g. "./image") <br>
- **tcell_track_csv** <br>
The path to the .csv containing the positions of each T cell track <br>

---
## Running
<br>

To run BEHAV3D, the following notebooks can be followed:
- **run_behav3d.ipynb** (Link to ipynb)<br>
    Run the full pipeline of BEHAV3D (including segmentation, tracking and feature extraction). Each block can be ran separetly to run the various parts of the BEHAV3D pipeline. It requires you to only change the config.yml path to the config you supply
- **run_behav3d_imaris.ipynb** (Link to ipynb)<br>
    Run the imaris pipeline of BEHAV3D (using preprocessed imaris statistics files (.csv) and converting it to track input that can be used for BEHAV3D feature extraction.
    It requires the multiple .csv files.
### Segmentation and tracking of input data
BEHAV3D provides an open-source approach to the segmentaion and tracking of time-lapse data.<br>
This step is ***optional***, as you can also process the data using other software (e.g. Imaris, Fiji)

**Imaris**
- Statistics files (.csv) containing features of tracks from objects (cells or organoids) from time-lapse statistics that are aquired by tracking these objects in a commercially available software (Imaris, Oxford Instruments).
However any type of time-lapse data can be processed with the pipeline, 
including measurements extracted from MTrackJ (Fiji) or others. 
Required is that the .csv contains a unique TrackID and the position of each segmnent in the track in t, z, y and x. 
- Aditional statistics describing the cell behavior such as speed, displacement are calculated by Imaris, however they can also be calculate by pre-processing algorithms from the cell coordinates. Statistics related to the expression of markers of interest (e.g live-dead cell dye) should be included to study the dynamic expression of these overtime. For statistics related to distance to organoids, use the *min_intensity in ch X* (corresponding to the channel number created by the Distance transformation Xtension. Rename it to be called *dist_org*.

The Imaris data then has to be converted to a .csv that is compatible with BEHAV3D processing, which can be easily done with a single function
### Feature extraction
This step extracts various track features from the tracked data. This includes movement features, dead dye features and contact features. Imaris processed statistics have to first be processed as previously described but can then be processed in almost the same manner as images processed with BEHAV3D skipping a few feature calculations that Imaris has already processed. <br>
The BEHAV3D pipeline then filters the tracks and additionally summarizes the features for analysis

### T cell track analysis
The T cell track analysis is currently still in R, although being converted into python.

For now, the filtered tracks can be used as input to the "predict_tcell_behavior.R" script (link to script)

### Organoid track analysis
The organoid track analysis is currently still in R, although being converted into python.

For now, you need Imaris processed files as input to the "organoid_death_dynamics.R" script (link to script)
_____

## Dataset example
In this repository we provide example datasets consisting of a multispectral time-lapse 3D imaging dataset originated from a co-culture of engeneered T cells and Tumor derived organoids from the BEHAV3D [original paper](https://www.nature.com/articles/s41587-022-01397-w). Multispectral imaging allows to identify: Live/dead T cells; Live/Dead organoids. For downstream analysis of organoids: Either individual tumor derived organoids are tracked overtime or the total organoid volume per well is tracked. For each generated object we acquire information on the dead cell dye intensity and position and volume of individual organoids. For downstream analysis of T cell: T cells are tracked overtime. For each Tracked T cell object we aquire, position per timepoint, speed, square displacement, distance to an organoid, dead dye intensity, major and minor axis length (used in some downstream analysis).
_____

## Repository
This repository contains a collection of scripts and example datasets enabling the following dowstream analysis. Follow the structure in the script folder for each module and each analysis type. Introduce the corresponding folder/ file direction on your own computer where required (note that to specify directory paths in R (/) forward slash is recommended):
_____

## Demo

You can run BEHAV3D on demo data to see examples of the results.\
\
There are 2 demos:
- tcell_demo    (For 'tcell_dynamics_classification' )
- organoid_demo (For 'organoid_death_dynamics')

**>Step 1** To set up the demo on you local PC, run [BEHAV3D/demos/set_up_demo.R](https://github.com/RiosGroup/BEHAV3D/blob/main/demos/set_up_demo.R)\
This sets up the paths in the BEHAV3D config file for the demo, then run the different modules on the demo (look below).

## Modules
### Organoids death dynamics module

This module examines the organoid death over time (individual organoids and per well statistics)

***To run from command line:***
```
Rscript ...BEHAV3D/scripts/organoid_death_dynamics/organoid_death_dynamics.R -c </Path/to/BEHAV3D/config>
```

***To run from Jupyter notebook:***<br>
***"Currently In Progress transfer to python"***<br>
