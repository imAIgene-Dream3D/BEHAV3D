# BEHAV3D pipeline
## Overview
BEHAV3D is dynamic immuno-organoid 3D imaging-transcriptomics platform to study tumor death dynamics; immune cell behavior and behavior-guided transcriptomics. 

## What type of data does BEHAV3D work with?
- Any type of multispectral time-lapse 3D (or 2D) imaging data, where objects such as cells or organoids are in co-culture or single culture. 
## What output can BEHAV3D provide?
- Any type of change of cell state that can be detected by a change in fluorescent intensity e.g. cell death, reporter, Ca2+ signalling
- Classification of different types of cell dynamics
- Correlation between dynamics of different cell types
- Interaction between cell types
- Predicted behavior states infered to transcriptomic data

## Software and Hardware requirements
BEHAV3D runs in R studio or from command line and was tested on MacOS Big Sur with R version 4.1.1. 

## Installation
Download the repository to your PC via direct dowload or git clone https://github.com/alievakrash/BEHAV3D.git in Git Bash.

For the required libraries, it is possible to run the "install_required_packages.R" script in the main BEHAV3D folder.

BEHAV3D uses the following libraries:
- abind
- dplyr
- dtwclust
- ggplot2
- gplots
- MESS
- optparse
- parallel
- pheatmap
- plyr
- randomForest
- readr
- reshape2
- scales
- Seurat
- spatstat
- sp
- stats
- tibble
- tidyr
- umap
- viridis
- xlsx
- yaml
- zoo

## Input data
The current version of the pipeline works with objects (cells or organoids) time-lapse statistics that are aquired by tracking these objects in a commercially available software (Imaris, Oxford Instruments). 
However any type of time-lapse data can be processed with the pipeline, including measruements extract from MTrackJ (Fiji) or others. Main feature that is needed are coordinates for the objects and a common ID for the same object that is tracked over time. Aditional statistics describing the cell behavior such as speed, displacement are calculated by Imaris, however they can also be calculate by pre-processing algorithms from the cell coordinates. Statistics related to the expression of markers of interest (e.g live-dead cell dye) should be included to study the dynamic expression of these overtime. For statistics related to distance to organoids, use the *min_intensity in ch X* (corresponding to the channel number created by the Distance transformation Xtension. Rename it to be called *dist_org*.

## Dataset example
In this repository we provide example datasets consisting of a multispectral time-lapse 3D imaging dataset originated from a co-culture of engeneered T cells and Tumor derived organoids. Multispectral imaging allows to identify: Live/dead T cells; Live/Dead organoids. For downstream analysis of organoids: Either individual tumor derived organoids are tracked overtime or the total organoid volume per well is tracked. For each generated object we acquire information on the dead cell dye intensity and position and volume of individual organoids. For downstream analysis of T cell: T cells are tracked overtime. For each Tracked T cell object we aquire, position per timepoint, speed, square displacement, distance to an organoid, dead dye intensity, major and minor axis length (used in some downstream analysis).

## Repository
This repository contains a collection of scripts and example datasets enabling the following dowstream analysis. Follow the structure in the script folder for each module and each analysis type. Introduce the corresponding folder/ file direction on your own computer where required (note that to specify directory paths in R (/) forward slash is recommended):

## Set-up

BEHAV3D uses 2 specific fiels to set up analysis:\

### **BEHAV3D config**
Contains all experiment-specific settings and paths to data for all modules in BEHAV3D\
An example version can be found in [...BEHAV3D/configs/config_template.yml](https://github.com/RiosGroup/BEHAV3D/blob/main/configs/config_template.yml)\
Explanation for the different variables are explained in that template

### **Experimental metadata template**
To correctly import data for BEHAV3D, it is required to fill in a .tsv that contains information per experiment performed, requiring information on:
- Experiment basename
- organoid_line
- tcell_line
- exp_nr
- well
- date
- dead_dye_channel (Channel that contains the channel with the dead dye)
- organoid_distance_channel (Channel that contains the distance to organoids)
- tcell_contact_threshold (threshold of distance to other tcells to be considered touching))
- tcell_dead_dye_threshold (threshold to consider an tcell "dead")
- tcell_stats_folder (path to folder with tcell track statistics)
- organoid_contact_threshold (threshold of distance to organoid to be considered touching)
- organoid_dead_dye_threshold (threshold to consider an organoid "dead")
- organoid_stats_folder (path to folder with organoid track statistics)


## Modules
### (1) Organoids death dynamics module

This module examines the organoid death over time (individual organoids and per well statistics)

***To run from command line:***
```
Rscript ...BEHAV3D/scripts/organoid_death_dynamics/organoid_death_dynamics.R -c </Path/to/BEHAV3D/config>
```

***To run from Rstudio:***

Change the config path on [line 18](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/Organoids%20death%20dynamics/organoid_death_dynamics.R#L13)

***Output_files***

All organoid death dynamics unfiltered
- Full_well_death_dynamics.pdf  (Mean of dead dye intensity over the full well)
- Full_well_death_dynamics.rds
- Full_individual_orgs_death_dynamics.pdf   (Mean of dead dye intensity over per organod track in each well)
- Full_individual_orgs_death_dynamics.rds
- Full_percentage_dead_org_over_time.pdf    (Percentage of overall organoid death occuring over time)
- Full_percentage_dead_org_over_time.rds

Filtered organoid death dynamics based on BEHAV3D config
- Exp_well_death_dynamics.pdf
- Exp_well_death_dynamics.rds
- Exp_individual_orgs_death_dynamics.pdf
- Exp_individual_orgs_death_dynamics.rds
- Exp_percentage_dead_org_over_time.pdf
- Exp_percentage_dead_org_over_time.rds

### (2) T cell behavior classification module

This module examines the tcell dynamics and either performs clustering (no randomforest supplied) or classifcation (randomforest supplied)

***To run from command line:***\
If running from command line, provide the config.yml file for your experiment, like so:\
```
Rscript ...BEHAV3D/scripts/tcell_dynamics_classification/predict_tcell_behavior.R -c </Path/to/BEHAV3D/config>
```
optionally, you can supply two additional parameters:
- -t/--tracks_rds : Path to RDS file containing processed T cell track data if not supplied in the BEHAV3D config
- -f/--force_redo : Supply if you want to re-import and reprocess the track data, default BEHAV3D checks if file already exists in the output folder structure and skips if available

```
Rscript ...BEHAV3D/scripts/tcell_dynamics_classification/predict_tcell_behavior.R -c </Path/to/BEHAV3D/config> -t </Path/to/TrackRDSfile> -f
```

***To run from Rstudio***\
Change the config path on [line 23](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/T%20cell%20dynamics%20classification/predict_tcell_behavior.R#L23)\
(Optional) Change the force_redo parameter on [line 16](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/T%20cell%20dynamics%20classification/predict_tcell_behavior.R#L16)\
(Optional) Change the for parameter on [line 17](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/T%20cell%20dynamics%20classification/predict_tcell_behavior.R#L17)

### ***(Optional) You can (re)train the randomforest with the following steps***
- Run ...BEHAV3D/scripts/tcell_dynamics_classification/predict_tcell_behavior.R 
- Run ...BEHAV3D/scripts/tcell_dynamics_classification/train_randomforest/train_random_forest_classifier.R from command line
```
Rscript ...BEHAV3D/scripts/tcell_dynamics_classification/train_randomforest/train_random_forest_classifier.R -i </Path/to/behavioral/reference/map> -o </Path/to/output/randomForest>
```
- or change the [input parameter](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/T%20cell%20dynamics%20classification/train_randomforest/train_random_forest_classifier.R#L14) and [output parameter](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/T%20cell%20dynamics%20classification/train_randomforest/train_random_forest_classifier.R#L15)


## (3) Behavior-guided transcriptomics module
This module integrates information from single cell sequencing and behavioral profiling, by predicting in a behavioral phenotype of single cells in scRNA seq data. For more information see Figure 4 in https://www.biorxiv.org/content/10.1101/2021.05.05.442764v2

Predict in silico the proportions of cells with different behavioral signatures in different experimental groups [non-engaged, non-engaged enriched, engaged, super-engaged]

***To run from command line:***\
```
Rscript ...BEHAV3D/scripts/behavior_guided_transcriptomics/1.in_silico_engager-superengager_selection.R -c </Path/to/BEHAV3D/config> [-t </Path/to/trackRDS>]
Rscript ...BEHAV3D/scripts/behavior_guided_transcriptomics/2.behavioral-guided_transcriptomics.R -c </Path/to/BEHAV3D/config> [-t </Path/to/trackRDS>]
```

***To run from Rstudio:***\

Change the config path in [1.in_silico_engager-superengager_selection.R](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/Behavior-guided%20transcriptomics/1.in_silico_engager-superengager_selection.R#L15) and [2.behavioral-guided_transcriptomics.R](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/Behavior-guided%20transcriptomics/2.behavioral-guided_transcriptomics.R#L16)

(Optional) Supply a tracks_rds in [1.in_silico_engager-superengager_selection.R](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/Behavior-guided%20transcriptomics/1.in_silico_engager-superengager_selection.R#L10) and [2.behavioral-guided_transcriptomics.R](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/Behavior-guided%20transcriptomics/2.behavioral-guided_transcriptomics.R#L11)

