## BEHAV3D pipeline
### Overview
BEHAV3D is dynamic immuno-organoid 3D imaging-transcriptomics platform to study tumor death dynamics; immune cell behavior and behavior-guided transcriptomics. 

#### What type of data does BEHAV3D work with?
- Any type of multispectral time-lapse 3D (or 2D) imaging data, where objects such as cells or organoids are in co-culture or single culture. 
#### What output can BEHAV3D provide?
- Any type of change of cell state that can be detected by a change in fluorescent intensity e.g. cell death, reporter, Ca2+ signalling
- Classification of different types of cell dynamics
- Correlation between dynamics of different cell types
- Interaction between cell types
- Predicted behavior states infered to transcriptomic data

#### Software and Hardware requirements
BEHAV3D runs on R studio (version 4.0.2) and was tested on Windows 10. See package requirements for each module. 
#### Installation
Download the repository to your PC via direct dowload or git clone https://github.com/alievakrash/BEHAV3D.git in Git Bash.
#### Input data
The current version of the pipeline works with objects (cells or organoids) time-lapse statistics that are aquired by tracking these objects in a commercially available software (Imaris, Oxford Instruments). However any type of time-lapse data can be processed with the pipeline, including measruements extract from MTrackJ (Fiji) or others. Main feature that is needed are coordinates for the objects and a common ID for the same object that is tracked over time. Aditional statistics describing the cell behavior such as speed, displacement are calculated by Imaris, however they can also be calculate by pre-processing algorithms from the cell coordinates. Statistics related to the expression of markers of interest (e.g live-dead cell dye) should be included to study the dynamic expression of these overtime.
#### Dataset example
In this repository we provide an example dataset consisting of a multispectral time-lapse 3D imaging dataset originated from a co-culture of engeneered T cells and Tumor derived organoids. Multispectral imaging allows to identify: Live/dead T cells; Live/Dead organoids. For downstream analysis of organoids: either individual tumor derived organoids are tracked overtime or the total organoid volume per well is tracked. For each generated object we acquire information on the dead cell dye intensity and position and volume of individual organoids. For downstream analysis of T cell: T cells are tracked overtime. For each Tracked T cell object we aquire, position per timepoint, speed, square displacement, distance to an organoid, dead dye intensity, major and minor axis length (used in some downstream analysis).
## Repository
This repository contains a collection of scripts and example datasets enabling the following dowstream analysis. Follow the structure in the script folder for each module and each analysis:
### (1) Organoids death dynamics module
- Batch import tracked organoids data: 

-Run script [batch_import_organoids_data_for_each_n](https://github.com/alievakrash/BEHAV3D/blob/main/scripts/Organoids%20death%20dynamics/Batch%20import%20organoids%20data/batch_import_organoids_data_for_each_n.R)

-Insert the direction of the example dataset on your PC [here](https://github.com/alievakrash/BEHAV3D/blob/c0a285b105d6ee237dfc7f9b7bf912caa6f3e1cb/scripts/Organoids%20death%20dynamics/Batch%20import%20organoids%20data/batch_import_organoids_data_for_each_n.R#L9)

> Output files: Full_well_death_dynamics.rds and Individual_orgs_death_dynamics.rds

- Script for death dynamics per well

- Script for individual organoid death dynamics
### (2) T cell behavior classification module
- Script to import T cell data

- Script to Behavioral reference map generation (Dynamic time warping)

-Script for T cell behavior classification based on the Behavioral reference map (Random Forest classifier)

### (3) Behavior-guided transcriptomics module
- Script to create a behavioral probability map for scRNA seq data
