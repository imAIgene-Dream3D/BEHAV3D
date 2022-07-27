## Import the processed tracks with a length that is coinciding the the length of the experiment in Fig 4a.

library(yaml)
library(dplyr)
library(stats)
library(tidyr)
library(ggplot2)
library(optparse)

force_redo=TRUE
tracks_provided=NULL

### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  ### !!!!!! Change the path to the BEHAV3D_config file here if running the code in RStudio !!!!!!
  ### Demo path
  BEHAV3D_dir = paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/")
  pars = yaml.load_file(paste0(BEHAV3D_dir, "/demos/behavioral_transcriptomics_demo/BEHAV3D_config.yml"))
  
  ### For your own file
  # pars = yaml.load_file("")
  
    } else {
  option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, 
                help="Path to the BEHAV3D config file", metavar="character"),
    make_option(c("-f", "--force_redo"), action="store_true", default=FALSE, 
                help="Force the pipeline to re-import data even if files exists"),
    make_option(c("-t", "--tracks_rds"), type="character", default=NULL, 
                help="(Optional) Path to RDS file containing clustered T cell track data", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  if (is.null(opt$config)){
    print_help(opt_parser)
    stop("Config file -c|--config, must be supplied", call.=FALSE)
  }
  pars = yaml.load_file(opt$config)
  tracks_provided=opt$tracks_rds
  force_redo=opt$force_redo
  BEHAV3D_dir = paste0(dirname(dirname(dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))]))))), "/")
}
  
main_output_dir=pars$output_dir
output_dir=paste0(main_output_dir,"/transcriptomics/results/")
dir.create(output_dir, recursive=TRUE)

pars$tcell_exp_duration=pars$second_2
  pars$output_dir=output_dir
write_yaml(pars, file=paste0(output_dir,"in_silico_BEHAV3D_config.yml"))

if (is.null(tracks_provided)){
  tcelL_class_script=paste0(BEHAV3D_dir, "/scripts/tcell_dynamics_classification/predict_tcell_behavior.R")
  if (force_redo==TRUE){
    system(paste0("Rscript '", tcelL_class_script, "' -c ", paste0("'",output_dir,"in_silico_BEHAV3D_config.yml' -f")))
  } else{
    system(paste0("Rscript '", tcelL_class_script, "' -c ", paste0("'",output_dir,"in_silico_BEHAV3D_config.yml'")))
  }
  master_clust_Live6 = readRDS(file=paste0(output_dir,"/tcell_behavior/results/","classified_tcell_track_data.rds"))
} else {
  master_clust_Live6 = readRDS(file=opt$tracks_rds)
}

#######################################################################################################################################################
######### Reproduce in silico the experiment performed in Fig 4a. Separate cells at each washing step and compute the behaviors of the cells that were in contact with an organoid vs the ones that weren't
#################################################################################################################################################
#### Define the timepoints where the washing steps occured. Instead of selecting just one timepoint, that gives a very small amount of tracks, we define the time as a time-range around 3 hours or 5 hours of imaging. these correspond to 4 and 6 hours of culture respectively
around_5h<-c(pars$second_1:pars$second_2)
around_3h<-c(pars$first_1:pars$first_2)

### calculate the frequency of behaviors per conditions for CD8 only:

### Simulate the selection at 5 hours of T cells that were or were not in contact with organoids at that timepoint
master_clust_Live6_CD4_CD8<-filter(master_clust_Live6, grepl("CD4|CD8", tcell_line))

### Select cells at timepoint "around 5 hours" and remove cluster 1 since in the real experiment dead cells are excluded by FACS
master_clust_Live6_h_CD8<-subset(master_clust_Live6_CD4_CD8, Time%in%around_5h & grepl("CD8", tcell_line) & !cluster=="1") 
### For each trackID select the last timepoint in the 5 hours range
master_clust_Live6_h_CD8<-master_clust_Live6_h_CD8%>%group_by(TrackID)%>%filter(Time==max(Time)) 

## Plot the engager and the non engager populations and their cluster behavior distribution
Plot_Eng_CD8<-ggplot(master_clust_Live6_h_CD8, aes(x=factor(contact), fill=as.factor(cluster))) + 
  geom_bar(position="fill")+
  # scale_x_discrete(breaks=c(0,1), minor_breaks = NULL, n.breaks=2, labels = c("True", "False"))+
  ylab("Percentage")+
  xlab("Contact")+
  scale_x_discrete(labels=c("No Contact", "Contact"))+
  scale_fill_manual(
    values=c(
    "gold3",
    "darkolivegreen3",
    "seagreen3",
    "forestgreen",
    "dodgerblue",
    "cyan1",
    "indianred",
    "firebrick",
    "brown1"),drop = FALSE)+
  ggtitle(paste0("Engager & non-engager CD8 between \ntimepoint ", pars$second_1, " and ", pars$second_2))+
  theme_void()+
  theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
  coord_flip() +
  # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
  facet_grid(organoid_line  ~ tcell_line)
Plot_Eng_CD8

### Now plot the super engagers and the never engagers. These cells under go two washing and separation steps:
### 1) at 3 hours they are separated into  contacting organoids and non-contacting organoids and then put back in culture. For the non-contacting organoids condition new organoids are added.
### 2) at 5 hours they are separated again into contact organoids and non=contacting organoids T cells.

### Orverall we recover two populations: Super-engagers> T cells that were in contact with organoids at 3 and 5 hours
### Never engagers> T cells that were not in contact with organoids at 3 hours and also not at 5 hours.

### Simulate the selection at 3 hours of T cells that were in contact with an orgnaoid
master_clust_Live4_h_EN_CD8<-subset(master_clust_Live6_CD4_CD8, Time%in%around_3h& grepl("CD8", tcell_line)& contact==1)
master_clust_Live4_h_EN_CD8<-master_clust_Live4_h_EN_CD8%>%group_by(TrackID)%>%filter(Time==min(Time))
### Simulate the selection at 3 hours of T cells that were not in contact with an orgnaoid
master_clust_Live4_h_NEN_CD8<-subset(master_clust_Live6_CD4_CD8, Time%in%around_3h & grepl("CD8", tcell_line) &contact==0)
master_clust_Live4_h_NEN_CD8<-master_clust_Live4_h_NEN_CD8%>%group_by(TrackID)%>%filter(Time==min(Time))

###Simulate the separation at 5 hours of T cells in contact with organoid or not, based on the previously selected cells
master_clust_Live6_h_SEN_CD8<-subset(master_clust_Live6_h_CD8, contact==1 &TrackID%in%master_clust_Live4_h_EN_CD8$TrackID)
master_clust_Live6_h_NEN_CD8<-subset(master_clust_Live6_h_CD8, contact==0 &TrackID%in%master_clust_Live4_h_NEN_CD8$TrackID)
master_clust_Live6_h_NEN_SEN_CD8<-rbind(master_clust_Live6_h_NEN_CD8,master_clust_Live6_h_SEN_CD8)
master_clust_Live6_h_NEN_SEN_CD8<-subset(master_clust_Live6_h_NEN_SEN_CD8, !cluster=="1") ## remove dying cell cluster, since will be deleted by FACS

Plot_SEng_CD8<-ggplot(master_clust_Live6_h_NEN_SEN_CD8, aes(fill=as.factor(cluster), x=as.factor(contact))) + 
  geom_bar( position="fill")+ 
  ylab("Percentage")+
  xlab("Contact")+
  scale_x_discrete(labels=c("No Contact", "Contact"))+
  scale_fill_manual(
    values=c(
      "gold3",
      "darkolivegreen3",
      "seagreen3",
      "forestgreen",
      "dodgerblue",
      "cyan1",
      "indianred",
      "firebrick",
      "brown1"),drop = FALSE)+
  ggtitle(paste0("Super-engager & Never-engager CD8"))+
  theme_void()+
  theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
  coord_flip() +
  # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
  facet_grid(organoid_line  ~ tcell_line)
Plot_SEng_CD8

#### Repeat the same for the CD4 cells
### calculate the frequency of behaviors per conditions for CD4 only:
### Simulate the selection at 5 hours of T cells that were or were not in contact with organoids at that timepoint
master_clust_Live6_h_CD4<-subset(master_clust_Live6_CD4_CD8, Time%in%around_5h & grepl("CD4", tcell_line) & !cluster=="1")
master_clust_Live6_h_CD4<-master_clust_Live6_h_CD4%>%group_by(TrackID)%>%filter(Time==max(Time))

Plot_Eng_CD4<-ggplot(master_clust_Live6_h_CD4, aes(fill=as.factor(cluster), x=as.factor(contact))) + 
  geom_bar( position="fill")+ 
  geom_bar(position="fill")+
  # scale_x_discrete(breaks=c(0,1), minor_breaks = NULL, n.breaks=2, labels = c("True", "False"))+
  ylab("Percentage")+
  xlab("Contact")+
  scale_x_discrete(labels=c("No Contact", "Contact"))+
  scale_fill_manual(
    values=c(
      "gold3",
      "darkolivegreen3",
      "seagreen3",
      "forestgreen",
      "dodgerblue",
      "cyan1",
      "indianred",
      "firebrick",
      "brown1"),drop = FALSE)+
  ggtitle(paste0("Engager & non-engager CD4 between \ntimepoint ", pars$second_1, " and ", pars$second_2))+
  theme_void()+
  theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
  coord_flip() +
  # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
  facet_grid(organoid_line  ~ tcell_line)
Plot_Eng_CD4

### Never engagers> T cells that were not in contact with organoids at 3 hours and also not at 5 hours.
### Simulate the selection at 3 hours of T cells that were in contact with an orgnaoid
master_clust_Live4_h_EN_CD4<-subset(master_clust_Live6_CD4_CD8, Time%in%around_3h& grepl("CD4", tcell_line)& contact==1)
master_clust_Live4_h_EN_CD4<-master_clust_Live4_h_EN_CD4%>%group_by(TrackID)%>%filter(Time==min(Time))
### Simulate the selection at 3 hours of T cells that were not in contact with an orgnaoid
master_clust_Live4_h_NEN_CD4<-subset(master_clust_Live6_CD4_CD8, Time%in%around_3h & grepl("CD4", tcell_line) &contact==0)
master_clust_Live4_h_NEN_CD4<-master_clust_Live4_h_NEN_CD4%>%group_by(TrackID)%>%filter(Time==min(Time))

###Simulate the separation at 5 hours of T cells in contact with organoid or not, based on the previously selected cells
master_clust_Live6_h_SEN_CD4<-subset(master_clust_Live6_h_CD4, contact==1 &TrackID%in%master_clust_Live4_h_EN_CD4$TrackID)
master_clust_Live6_h_NEN_CD4<-subset(master_clust_Live6_h_CD4, contact==0 &TrackID%in%master_clust_Live4_h_NEN_CD4$TrackID)
master_clust_Live6_h_NEN_SEN_CD4<-rbind(master_clust_Live6_h_NEN_CD4,master_clust_Live6_h_SEN_CD4)
master_clust_Live6_h_NEN_SEN_CD4<-subset(master_clust_Live6_h_NEN_SEN_CD4, !cluster=="1") ## remove dying cell cluster, since will be deleted by FACS

Plot_SEng_CD4<-ggplot(master_clust_Live6_h_NEN_SEN_CD4, aes(fill=as.factor(cluster), x=as.factor(contact))) + 
  geom_bar( position="fill")+ 
  ylab("Percentage")+
  xlab("Contact")+
  scale_x_discrete(labels=c("No Contact", "Contact"))+
  scale_fill_manual(
    values=c(
      "gold3",
      "darkolivegreen3",
      "seagreen3",
      "forestgreen",
      "dodgerblue",
      "cyan1",
      "indianred",
      "firebrick",
      "brown1"),drop = FALSE)+
  ggtitle(paste0("Super-engager & Never-engager CD4"))+
  theme_void()+
  theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
  coord_flip() +
  # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
  facet_grid(organoid_line  ~ tcell_line)
Plot_SEng_CD4

########################### CD8 ########################################################
### calculare the proportions for combination with scRNA seq:
### create a condition for each experimental condition
master_clust_Live6_h_CD8$engagement<-ifelse(master_clust_Live6_h_CD8$contact==1, "engager","nonengager")
master_clust_Live6_h_NEN_SEN_CD8$engagement<-ifelse(master_clust_Live6_h_NEN_SEN_CD8$contact==1,"super-engaged","never-engaged")
## join both datasets:
CD8_engagement<-rbind(master_clust_Live6_h_CD8,master_clust_Live6_h_NEN_SEN_CD8)
library(dplyr)
CD8_behav<-CD8_engagement %>% group_by(cluster,engagement) %>%summarize(frec = n())
engagement_n<-CD8_engagement %>% group_by(engagement) %>%summarize(total_n = n())
CD8_behav<-left_join(CD8_behav,engagement_n)
CD8_behav$cluster_prop<-CD8_behav$frec/CD8_behav$total_n
colnames(CD8_behav)=c("behavioral_cluster", "exp_condition", "count", "total_n", "cluster_proportion")

CD8_engagement$engagement <- factor(CD8_engagement$engagement, levels = c("never-engaged", "nonengager", "engager", "super-engaged"))
Plot_x<-ggplot(CD8_engagement, aes(fill=as.factor(cluster), x=engagement)) + 
  geom_bar( position="fill")+ 
  ylab("Percentage")+
  xlab("Contact")+
  # scale_x_discrete(breaks=labels=c("No Contact", "Contact"))+
  scale_fill_manual(
    values=c(
      "gold3",
      "darkolivegreen3",
      "seagreen3",
      "forestgreen",
      "dodgerblue",
      "cyan1",
      "indianred",
      "firebrick",
      "brown1"),drop = FALSE)+
  ggtitle(paste0("Proportions of the four different types of engagers CD8"))+
  theme_void()+
  theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
  coord_flip() +
  # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
  facet_grid(organoid_line  ~ tcell_line)
Plot_x
### save data to csv same folder where the pseudotime clustering from Farid
write.csv(CD8_behav,paste0(output_dir, "CD8_engagement_behavior_freq.csv"), row.names = FALSE)

#### Calculate the mean contact time per condition for information (used in the paper)
## SE mean contact time:
master_clust_Live6_h_CD8_se_pop<-subset(CD8_engagement, engagement=="super-engaged")
###calculate the mean engagement time for each signature:
mean_cont_SE<- master_clust_Live6_CD4_CD8%>%filter(TrackID%in%master_clust_Live6_h_CD8_se_pop$TrackID)%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)
mean_cont_SE$exp<-"super-engaged"


## E mean contact time:
master_clust_Live6_h_CD8_e_pop<-subset(CD8_engagement, engagement=="engaged")
###calculate the mean engagement time for each signature:
mean_cont_E<- master_clust_Live6_CD4_CD8%>%filter(TrackID%in%master_clust_Live6_h_CD8_e_pop$TrackID)%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)
mean_cont_E$exp<-"engaged"


## NonE mean contact time:
master_clust_Live6_h_CD8_noe_pop<-subset(CD8_engagement, engagement=="non-engaged")
###calculate the mean engagement time for each signature:
mean_cont_NoE<- master_clust_Live6_CD4_CD8%>%filter(TrackID%in%master_clust_Live6_h_CD8_noe_pop$TrackID)%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)
mean_cont_NoE$exp<-"non-engaged"


## NonE mean contact time:
master_clust_Live6_h_CD8_NE_pop<-subset(CD8_engagement, engagement=="never-engaged")
###calculate the mean engagement time for each signature:
mean_cont_NE<- master_clust_Live6_CD4_CD8%>%filter(TrackID%in%master_clust_Live6_h_CD8_NE_pop$TrackID)%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)

mean_cont_NE$exp<-"never-engaged"

mean_cont<-rbind(mean_cont_NE,mean_cont_NoE,mean_cont_E,mean_cont_SE)

saveRDS(mean_cont, paste0(output_dir,"min_contact_per_hour_13T_per_exp"))



########################### CD4 ########################################################
### calculare the proportions for combination with scRNA seq:
### create a condition for each experimental condition
master_clust_Live6_h_CD4$engagement<-ifelse(master_clust_Live6_h_CD4$contact==1, "engager","nonengager")
master_clust_Live6_h_NEN_SEN_CD4$engagement<-ifelse(master_clust_Live6_h_NEN_SEN_CD4$contact==1,"super-engaged","never-engaged")
## join both datasets:
CD4_engagement<-rbind(master_clust_Live6_h_CD4,master_clust_Live6_h_NEN_SEN_CD4)
library(dplyr)
CD4_behav<-CD4_engagement %>% group_by(cluster,engagement) %>%summarize(frec = n())
engagement_n<-CD4_engagement %>% group_by(engagement) %>%summarize(total_n = n())
CD4_behav<-left_join(CD4_behav,engagement_n)
CD4_behav$cluster_prop<-CD4_behav$frec/CD4_behav$total_n
colnames(CD4_behav)=c("behavioral_signature", "exp_condition", "count", "total_n", "cluster_proportion")

CD4_engagement$engagement <- factor(CD4_engagement$engagement, levels = c("never-engaged", "nonengager", "engager", "super-engaged"))
Plot_y<-ggplot(CD4_engagement, aes(fill=as.factor(cluster), x=engagement)) + 
  geom_bar( position="fill")+ 
  ylab("Percentage")+
  xlab("Contact")+
  # scale_x_discrete(breaks=labels=c("No Contact", "Contact"))+
  scale_fill_manual(
    values=c(
      "gold3",
      "darkolivegreen3",
      "seagreen3",
      "forestgreen",
      "dodgerblue",
      "cyan1",
      "indianred",
      "firebrick",
      "brown1"),drop = FALSE)+
  ggtitle(paste0("Proportions of the four different types of engagers CD4"))+
  theme_void()+
  theme(aspect.ratio = 0.2, axis.text.y = element_text())+  
  coord_flip() +
  # facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line) #For plotting all wells separately
  facet_grid(organoid_line  ~ tcell_line)
Plot_y

### save data to csv same folder where the pseudotime clustering from Farid
write.csv(CD4_behav,paste0(output_dir,"CD4_engagement_behavior_freq.csv"), row.names = FALSE)

## Save output for figure Fig 4a
pdf(paste0(output_dir,"Engager_Super_ENG_proportions_CD4_CD8.pdf"))
Plot_x
Plot_y
dev.off()

