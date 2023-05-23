# BEHAV3D pipeline
## Overview
BEHAV3D is dynamic immuno-organoid 3D imaging-transcriptomics platform to study tumor death dynamics; immune cell behavior and behavior-guided transcriptomics.

## What type of data does BEHAV3D work with?
- Any type of multispectral time-lapse 3D (or 2D) imaging data, where objects such as tumor cells or tumor organoids are in co-culture with immune cells of interest.
## What output can BEHAV3D provide?
- Any type of change of cell state that can be detected by a change in fluorescent intensity e.g. cell death, reporter, Ca2+ signalling
- Classification of different types of cell dynamics
- Tumor death dynamics quantification
- Backprojection of behavioral phenotype in Imaris 3D image visualization software 
- Correlation between tumor death dynamics and behavioral phenotypes
## How to cite this pipeline
Dekkers JF*, Alieva M*, Cleven A, Keramati F, Wezenaar AKL, van Vliet EJ, Puschhof J, Brazda P, Johanna I, Meringa AD, Rebel HG, Buchholz MB, Barrera Román M, Zeeman AL, de Blank S, Fasci D, Geurts MH, Cornel AM, Driehuis E, Millen R, Straetemans T, Nicolasen MJT, Aarts-Riemens T, Ariese HCR, Johnson HR, van Ineveld RL, Karaiskaki F, Kopper O, Bar-Ephraim YE, Kretzschmar K, Eggermont AMM, Nierkens S, Wehrens EJ, Stunnenberg HG, Clevers H, Kuball J, Sebestyen Z, Rios AC. **Uncovering the mode of action of engineered T cells in patient cancer organoids**. * *equal contibution* Nat Biotechnol. 2023 Jan https://doi.org/10.1038/s41587-022-01397-w
## Software and Hardware requirements
BEHAV3D runs in R studio or from command line and was tested on MacOS Big Sur with R version 4.1.1.

## Installation
Download the repository to your PC via direct dowload or git clone https://github.com/RiosGroup/BEHAV3D in Git Bash.

BEHAV3D uses the following R libraries:
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
In this repository we provide example datasets consisting of a multispectral time-lapse 3D imaging dataset originated from a co-culture of engeneered T cells and Tumor derived organoids from the BEHAV3D [original paper](https://www.nature.com/articles/s41587-022-01397-w). Multispectral imaging allows to identify: Live/dead T cells; Live/Dead organoids. For downstream analysis of organoids: Either individual tumor derived organoids are tracked overtime or the total organoid volume per well is tracked. For each generated object we acquire information on the dead cell dye intensity and position and volume of individual organoids. For downstream analysis of T cell: T cells are tracked overtime. For each Tracked T cell object we aquire, position per timepoint, speed, square displacement, distance to an organoid, dead dye intensity, major and minor axis length (used in some downstream analysis).

## Repository
This repository contains a collection of scripts and example datasets enabling the following dowstream analysis. Follow the structure in the script folder for each module and each analysis type. Introduce the corresponding folder/ file direction on your own computer where required (note that to specify directory paths in R (/) forward slash is recommended):

## Set-up

BEHAV3D uses 2 specific fiels to set up analysis:\

### **BEHAV3D config**
Contains all experiment-specific settings and paths to data for all modules in BEHAV3D\
An example version can be found in [...BEHAV3D/configs/config_template.yml](https://github.com/RiosGroup/BEHAV3D/blob/main/configs/config_template.yml)\
Explanation on what each variable changes is commented in that template

### **Experimental metadata template**
To correctly import data for BEHAV3D, it is required to fill in a .tsv that contains information per experiment performed, requiring information on:
- Experiment basename
- organoid_line
- tcell_line
- exp_nr
- well
- date
- dead_dye_channel (Channel # that contains the dead dye intensities)
- organoid_distance_channel (Channel # that contains the distance to organoids)
- tcell_contact_threshold (threshold of distance to other tcells to be considered touching, usually set to average cell diameter (10um), change if cells are of different size )
- tcell_dead_dye_threshold (threshold to consider an tcell "dead", can be visually estimates from the Imaris file based on mean intentisity of ch "dead cells")
- tcell_stats_folder (path to folder with tcell track statistics)
- organoid_contact_threshold (threshold of distance to organoid to be considered touching,estimate from the Imaris file based on min distance to ch "distance transformation" in T cell touching the organoid. Note that while this value is usually fixed in the same experiment, it might change if image resolution is different or between Imaris versions)
- organoid_dead_dye_threshold (threshold to consider an organoid "dead", same as for tcell_dead_dye_threshold)
- organoid_stats_folder (path to folder with organoid track statistics)

For an example see: [...BEHAV3D/configs/metadata_template.tsv](https://github.com/RiosGroup/BEHAV3D/blob/main/configs/metadata_template.tsv)\

## Demo

You can run BEHAV3D on demo data to see examples of the results.\
\
There are 2 demos:
- tcell_demo    (For 'tcell_dynamics_classification' )
- organoid_demo (For 'organoid_death_dynamics')

To set up the the demo, run [BEHAV3D/demos/set_up_demo.R](https://github.com/RiosGroup/BEHAV3D/blob/main/demos/set_up_demo.R)\
This sets up the paths in the BEHAV3D config file for the demo, then run the different modules on the demo (look below).

## Modules
### (1) Organoids death dynamics module

This module examines the organoid death over time (individual organoids and per well statistics)

***To run from command line:***
```
Rscript ...BEHAV3D/scripts/organoid_death_dynamics/organoid_death_dynamics.R -c </Path/to/BEHAV3D/config>
```

***To run from Rstudio:***

Change the config path on [line 18](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/organoid_death_dynamics/organoid_death_dynamics.R#L18)

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
Change the config path on [line 27](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/tcell_dynamics_classification/predict_tcell_behavior.R#L27)\
(Optional) Change the force_redo parameter on [line 16](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/tcell_dynamics_classification/predict_tcell_behavior.R#L16)\
(Optional) Change the for parameter on [line 17](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/tcell_dynamics_classification/predict_tcell_behavior.R#L17)

***Output_files***

output rds:
- raw_tcell_track_data.rds (Combined raw track data for all experiments)
- processed_tcell_track_data.rds (Combined processed track data for all experiments; Added contact, death etc.)
- behavioral_reference_map.rds   (Output of the t_cell behavior that can be used to train your own randomForest)
- classified_tcell_track_data.rds (Combined classified track data - Either randomForest or clustering - for all experiments)
- classified_tcell_track_data_summary.rds (Summary of classified track data - Either randomForest or clustering- for all experiments)
- cluster_perc_tcell_track_data.rds (Cluster percentages for the track data - used for creation of RF_ClassProp_WellvsCelltype.pdf)
- (If no randomForest) tcell_track_features.rds (Track features used to create Cluster_heatmap.pdf)

output pdf:
- (If randomForest supplied) RF_ClassProp_WellvsCelltype.pdf (Proportion of each randomForest predetermined cluster for all experiments)
- (If no randomForest) Umap_unclustered.pdf
- (If no randomForest) Umap_clustered.pdf
- (If no randomForest) Cluster_heatmap.pdf (Heatmap of track features for created clusters)
- (If no randomForest) umap_cluster_percentage_bars_separate.pdf (Proportion of each cluster for each separate experiment)
- (If no randomForest) umap_cluster_percentage_bars_combined.pdf (Proportion of each cluster for combiend experiments - based on organoid_lines and tcell_lines)

quality control:
- NrCellTracks_filtering_perExp.pdf (Shows the number of tracks remaining after each filtering step, shown per experiment)
- NrCellTracks_filtering_perFilt.pdf    (Shows the number of tracks remaining after each filtering step, shown per filtering step)
- TouchingvsNontouching_distribution.pdf    (Shows the number of Touching vs. Non-Touching T cells based on the defined "tcell_contact_threshold")
- DeadDye_distribution.pdf (Shows the distribution of mean dead dye intensity in Tcells per experiment, can be used to set correct "tcell_dead_dye_threshold" for dead cells)

### ***(Optional) You can (re)train the randomforest with the following steps***
- Run ...BEHAV3D/scripts/tcell_dynamics_classification/predict_tcell_behavior.R
- Run ...BEHAV3D/scripts/tcell_dynamics_classification/train_randomforest/train_random_forest_classifier.R from command line
```
Rscript ...BEHAV3D/scripts/tcell_dynamics_classification/train_randomforest/train_random_forest_classifier.R -i </Path/to/behavioral/reference/map> -o </Path/to/output/randomForest>
```
- or change the [input parameter](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/tcell_dynamics_classification/train_randomforest/train_random_forest_classifier.R#L14) and [output parameter](https://github.com/RiosGroup/BEHAV3D/blob/main/scripts/tcell_dynamics_classification/train_randomforest/train_random_forest_classifier.R#L15)


