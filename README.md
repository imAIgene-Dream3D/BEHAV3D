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
BEHAV3D runs on R studio (version 4.0.2) and was tested on Windows 10. 
#### Installation
Download the repository to your PC via direct dowload or git clone https://github.com/alievakrash/BEHAV3D.git in Git Bash.

See in each module and script the required libraries or install in R Studio libraries for all modules: devtools, dplyr, dtwclust, eulerr, ggplot2, ggpubr, ggrepel, gridExtra, hypergeo, kmodR, lme4, lmerTest, MESS, nlme, openxlsx, parallel, patternplot, pheatmap, plotly, plyr, png, purr, RColorBrewer, readxl, reshape2, rgeos, scales, Seurat, sp, spatstat, stats, tidyr, tidyverse, umap, VennDiagram, viridis, xlsx, zoo.
#### Input data
The current version of the pipeline works with objects (cells or organoids) time-lapse statistics that are aquired by tracking these objects in a commercially available software (Imaris, Oxford Instruments). However any type of time-lapse data can be processed with the pipeline, including measruements extract from MTrackJ (Fiji) or others. Main feature that is needed are coordinates for the objects and a common ID for the same object that is tracked over time. Aditional statistics describing the cell behavior such as speed, displacement are calculated by Imaris, however they can also be calculate by pre-processing algorithms from the cell coordinates. Statistics related to the expression of markers of interest (e.g live-dead cell dye) should be included to study the dynamic expression of these overtime. For statistics related to distance to organoids, use the *min_intensity in ch X* (corresponding to the channel number created by the Distance transformation Xtension. Rename it to be called *dist_org*.
#### Dataset example
In this repository we provide an example dataset consisting of a multispectral time-lapse 3D imaging dataset originated from a co-culture of engeneered T cells and Tumor derived organoids. Multispectral imaging allows to identify: Live/dead T cells; Live/Dead organoids. For downstream analysis of organoids: either individual tumor derived organoids are tracked overtime or the total organoid volume per well is tracked. For each generated object we acquire information on the dead cell dye intensity and position and volume of individual organoids. For downstream analysis of T cell: T cells are tracked overtime. For each Tracked T cell object we aquire, position per timepoint, speed, square displacement, distance to an organoid, dead dye intensity, major and minor axis length (used in some downstream analysis).
## Repository
This repository contains a collection of scripts and example datasets enabling the following dowstream analysis. Follow the structure in the script folder for each module and each analysis type. Introduce the corresponding folder/ file direction on your own computer where required (note that to specify directory paths in R (/) forward slash is recommended):
## (1) Organoids death dynamics module
- Batch import tracked organoids data: 

-Run script [batch_import_organoids_data_for_each_n](https://github.com/alievakrash/BEHAV3D/blob/main/scripts/Organoids%20death%20dynamics/Batch%20import%20organoids%20data/batch_import_organoids_data_for_each_n.R)

-Insert [here](https://github.com/alievakrash/BEHAV3D/blob/c0a285b105d6ee237dfc7f9b7bf912caa6f3e1cb/scripts/Organoids%20death%20dynamics/Batch%20import%20organoids%20data/batch_import_organoids_data_for_each_n.R#L9) the direction of the example dataset on your PC 

Output files: [Full_well_death_dynamics.rds](https://github.com/alievakrash/BEHAV3D/blob/553db58a0116559817b9f2109333cf4f7e58f4da/scripts/Organoids%20death%20dynamics/Test%20dataset/Full_well_death_dynamics) and [Individual_orgs_death_dynamics.rds](https://github.com/alievakrash/BEHAV3D/blob/553db58a0116559817b9f2109333cf4f7e58f4da/scripts/Organoids%20death%20dynamics/Test%20dataset/Individual_organoids_death_dynamics)

- Compute death dynamics per well:

-Run script [Calculate_mean_dead_dye_perwell.R](https://github.com/alievakrash/BEHAV3D/blob/5a2aed55ede54f2f14a987d7ab37480b8d15e038/scripts/Organoids%20death%20dynamics/Death%20dynamics%20per%20well/Calculate_mean_dead_dye_perwell.R)

-Insert [here](https://github.com/alievakrash/BEHAV3D/blob/fa49556265fc14c2da2355ef99561884ce65c807/scripts/Organoids%20death%20dynamics/Death%20dynamics%20per%20well/Calculate_mean_dead_dye_perwell.R#L5) the direction of the processed dataframe [Full_well_death_dynamics.rds](https://github.com/alievakrash/BEHAV3D/blob/553db58a0116559817b9f2109333cf4f7e58f4da/scripts/Organoids%20death%20dynamics/Test%20dataset/Full_well_death_dynamics) on your PC


- Compute individual organoid death dynamics:

-Run script [Individual organoids death dynamics.R](https://github.com/alievakrash/BEHAV3D/blob/5a2aed55ede54f2f14a987d7ab37480b8d15e038/scripts/Organoids%20death%20dynamics/Individual%20organoids%20death%20dynamics/Individual%20organoids%20death%20dynamics.R)

-Insert [here](https://github.com/alievakrash/BEHAV3D/blob/81ab7207a48fb60f8467737e5aa2e85f643f054d/scripts/Organoids%20death%20dynamics/Individual%20organoids%20death%20dynamics/Individual%20organoids%20death%20dynamics.R#L4) the direction of the processed dataframe [Individual_orgs_death_dynamics.rds](https://github.com/alievakrash/BEHAV3D/blob/553db58a0116559817b9f2109333cf4f7e58f4da/scripts/Organoids%20death%20dynamics/Test%20dataset/Individual_organoids_death_dynamics)

## (2) T cell behavior classification module
The T cell behavior classification module requires a configuration files that describes the BEHAV3D experiment:\
BEHAV3D/configs/config_tempalte.yml

- output_dir, where output files will be stored
- metadata_csv_path, a .txt file that contains a table of metadata associated with each experiment in the analysis pipeline (see below for additonal information)
- randomforest_path, the path to the trained randomforest to use for prediction
- reference_map_path, the path to the reference map that is used for randomforest training
- exp_duration, the maximum length of the experiment, all timepoints after will be truncated
- max_track_length, the maximum track length to use for analysis, longer tracks will be truncated
- tcell_contact_thr, nearest distance for T cells to be to each other that define as contact
- dead_tcell_thr, intensity cut-off to define dead cells

Optionally and recommended, for every BEHAV3D experiment a metadata .tsv that contains information on each sample used should be created.
This contains information such as the tcell line used, the organoid line used, the well that can be used to group samples. This metadata file should have the same filename as the rest of the statistics of that experiment, but ending in _metadata.csv. 

Example [20201027_(6)10T_b_metadata.csv](https://github.com/alievakrash/BEHAV3D/blob/dev_01/scripts/T%20cell%20dynamics%20classification/Example_dataset_T_cell_tracking/20201027_(6)10T_b_metadata.csv) 

A template can be found in:
- BEHAV3D/configs/metadata_template.tsv

Import and prediction can be run by running the following scripts:
- BEHAV3D/scripts/T cell dynamics classification/predict_tcell_behavior.R

If running from command line, provide the config.yml file for your experiment, like so:\
```
predict_tcell_behavior.R /a/b/config.yml
```
When running interactively, manually enter the path to the config at the start of the scripts
```
pars=<PATH/TO/CONFIG>
```

#### (Optional) You can (re)train the randomforest with the following steps:
Compute the Behavioral reference map generation (by Dynamic time warping):
- Run script [Create_Behavioral_Reference_map.R](https://github.com/alievakrash/BEHAV3D/blob/57c67317eea1af74d9aa82b33a9fab795d0a2dcc/scripts/T%20cell%20dynamics%20classification/Create_Behavioral_Reference_map.R)
- Insert [here](https://github.com/alievakrash/BEHAV3D/blob/57c67317eea1af74d9aa82b33a9fab795d0a2dcc/scripts/T%20cell%20dynamics%20classification/Create_Behavioral_Reference_map.R#L2) the direction of your dataset used for reference map. In this case we use the example dataset, consisting of only two different wells. For the creation of a new reference map use a compilation of datasets with different cell types of interest (likely to have different behavioral signatures). If you are creating a new reference map, insert [here](https://github.com/alievakrash/BEHAV3D/blob/c516cafc900cb71e8d33ba6d125b457923915bdb/scripts/T%20cell%20dynamics%20classification/Create_Behavioral_Reference_map.R#L92) the direction for your output dataframe.

Backproject behavioral signatures in the imaging dataset:
- Run script [Create_Behavioral_Reference_map.R](https://github.com/alievakrash/BEHAV3D/blob/57c67317eea1af74d9aa82b33a9fab795d0a2dcc/scripts/T%20cell%20dynamics%20classification/Create_Behavioral_Reference_map.R)
- Insert [here](https://github.com/alievakrash/BEHAV3D/blob/c516cafc900cb71e8d33ba6d125b457923915bdb/scripts/T%20cell%20dynamics%20classification/Create_Behavioral_Reference_map.R#L118) the direction of [master_example_data](https://github.com/alievakrash/BEHAV3D/blob/18f9332a54adf0b0d8e00d688802edc980aabdc9/scripts/T%20cell%20dynamics%20classification/example_dataset_T_cell_tracking/master_example_data) that is used to reconvert the unique TrackIDs that are created for processing back into the original TrackIDs, that will be used for backprojection.
- For each well of interest [adapt here](https://github.com/alievakrash/BEHAV3D/blob/c516cafc900cb71e8d33ba6d125b457923915bdb/scripts/T%20cell%20dynamics%20classification/Create_Behavioral_Reference_map.R#L127-L131) with the corresponding "ranks" and output direction.
- Predict T cell behavior classification for new datasets, based on the [Behavioral reference map](https://github.com/alievakrash/BEHAV3D/blob/57c67317eea1af74d9aa82b33a9fab795d0a2dcc/scripts/T%20cell%20dynamics%20classification/Behavioral%20reference%20map/Behavioral_Referance_map_git) :

Train the random forest, using the path to the behavioral reference map provided in the config.yml:
- BEHAV3D/scripts/T cell dynamics classification/train_randomforest/train_random_forest_classifier.R

## (3) Behavior-guided transcriptomics module
This module integrates information from single cell sequencing and behavioral profiling, by predicting in a behavioral phenotype of single cells in scRNA seq data. For more information see Figure 4 in https://www.biorxiv.org/content/10.1101/2021.05.05.442764v2

- Predict in silico the proportions of cells with different behavioral signatures in different experimental groups [non-engaged, non-engaged enriched, engaged, super-engaged]

-Run script [CD4_CD8_4_6_hours_beh_import](https://github.com/alievakrash/BEHAV3D/blob/dev_01/scripts/Behavior-guided%20transcriptomics/CD4_CD8_4_6_hours_beh_import.R) to select the tracks of a length that coincides with the length of the simulated experiment. If you are processing new experimenta conditions (e.g. T cells incubated with a new set of organoids), these first need to be imported according to [Import T cells data.R](https://github.com/alievakrash/BEHAV3D/blob/dev_01/scripts/T%20cell%20dynamics%20classification/Import%20T%20cells%20data.R) and the [master output](https://github.com/alievakrash/BEHAV3D/blob/555a3f59f6c8e5dca2f0ed877deea4055e947c9a/scripts/T%20cell%20dynamics%20classification/Import%20T%20cells%20data.R#L76) should be saved to be used as an input for [this script](https://github.com/alievakrash/BEHAV3D/blob/dev_01/scripts/Behavior-guided%20transcriptomics/CD4_CD8_4_6_hours_beh_import.R).

-Run script [CD4_CD8_4_6_hours_beh_predict](https://github.com/alievakrash/BEHAV3D/blob/dev_01/scripts/Behavior-guided%20transcriptomics/CD4_CD8_4_6_hours_beh_predict.R)

Store [output](https://github.com/alievakrash/BEHAV3D/blob/dev_01/scripts/Behavior-guided%20transcriptomics/example%20dataset/CD8_engagement_behavior_freq.csv) containing the proportion of behavioral signatures per experimental condition and per [cell type](https://github.com/alievakrash/BEHAV3D/blob/dev_01/scripts/Behavior-guided%20transcriptomics/example%20dataset/CD4_engagement_behavior_freq.csv) to be used below.

- Compute behavioral probability map for scRNA seq data:

-Run script [Behavioral-guided transcriptomics](https://github.com/alievakrash/BEHAV3D/blob/60a653b51d1417b0374225f720f5655c59f980ca/scripts/Behavior-guided%20transcriptomics/Behavioral-guided%20transcriptomics.R)

-Insert [here](https://github.com/alievakrash/BEHAV3D/blob/d3857d0ddebd6d1c3d88fe45c5a4aef8f648faf9/scripts/Behavior-guided%20transcriptomics/Behavioral-guided%20transcriptomics.R#L9) the directory containing the [scRNA seq dataset](https://github.com/alievakrash/BEHAV3D/blob/d3857d0ddebd6d1c3d88fe45c5a4aef8f648faf9/scripts/Behavior-guided%20transcriptomics/example%20dataset/scRNA_seq_dataset.rds) containing pseudotime trajectory on your PC.

-Insert [here](https://github.com/alievakrash/BEHAV3D/blob/beda6096a61b5795dde6316f0aa87dd3a1f1c1ae/scripts/Behavior-guided%20transcriptomics/Behavioral-guided%20transcriptomics.R#L18) the direction of [csv file](https://github.com/alievakrash/BEHAV3D/blob/96cf1b1e283354b41d578db8599a7be6f76f32b0/scripts/Behavior-guided%20transcriptomics/example%20dataset/CD8_behav_sig_per_exp_condition.csv) containing the proportion of behavioral signatures per experimental condition, calculated in silico based on imaging data.

-Insert [here](https://github.com/alievakrash/BEHAV3D/blob/96cf1b1e283354b41d578db8599a7be6f76f32b0/scripts/Behavior-guided%20transcriptomics/Behavioral-guided%20transcriptomics.R#L99) the direction of the [folder](https://github.com/alievakrash/BEHAV3D/tree/main/scripts/Behavior-guided%20transcriptomics/Probability_map) containing a set of probability maps quantified with different cluster resolutions.

-Check [here](https://github.com/alievakrash/BEHAV3D/blob/96cf1b1e283354b41d578db8599a7be6f76f32b0/scripts/Behavior-guided%20transcriptomics/Behavioral-guided%20transcriptomics.R#L161-L173) if any of the behavioral signatures has extreme outlier probability values due to skewed distributions and transformed to a maximal cutoff value as shown [here](https://github.com/alievakrash/BEHAV3D/blob/96cf1b1e283354b41d578db8599a7be6f76f32b0/scripts/Behavior-guided%20transcriptomics/Behavioral-guided%20transcriptomics.R#L170-L173) 

-Insert [here](https://github.com/alievakrash/BEHAV3D/blob/96cf1b1e283354b41d578db8599a7be6f76f32b0/scripts/Behavior-guided%20transcriptomics/Behavioral-guided%20transcriptomics.R#L236) the directory to store the output Seurat object. 
