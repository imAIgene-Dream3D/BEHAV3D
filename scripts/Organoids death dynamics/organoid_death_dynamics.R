set.seed(1234)
getwd()
library(plyr)
library(readr)
library(dplyr)
library(yaml)

### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  pars = yaml.load_file("/Users/samdeblank/Documents/1.projects/tcell_paper/20220330_fixes/10_per_cells/BEHAV3D_configOrg.yml")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  pars <- yaml.load_file(args[1])
}

## directory where the files are located
# working_directory <- pars$data_dir
output_dir=paste0(pars$output_dir,"/")

## function to import organoid data from csv files extracted by Imaris
read_imaris_csv <- function(flnm) {
  read_csv(flnm, skip = 3, col_types = cols()) %>% 
    mutate(filename = flnm)
}

# Import file-specific metadata for all images used in this analysis.
pat = pars$metadata_csv
metadata=read.csv(pars$metadata_csv, sep="\t")

if ( any(is.na(metadata$stat_folder)) ){
  metadata$stat_folder=apply(metadata, 1, function(x) paste0(pars$data_dir, x["basename"], "_Statistics"))
}

# Function to import organoid data specifically from Imaris generated csv files
read_ims_csv <- function(stat_folder, pattern) {
  read_plus <- function(flnm) {
    read_csv(flnm, skip = 3, col_types = cols()) %>% 
      mutate(filename = flnm)
  }
  pattern_file <- list.files(path = stat_folder, pattern = pattern, full.names=TRUE)
  print(pattern_file)
  ims_csv <- read_plus(pattern_file)
  return(ims_csv)
}

stat_folders <- metadata$stat_folder
pat = "*Volume"
volume_csv <- ldply(stat_folders, read_ims_csv, pattern=pat)
# import sum_red
pat = paste0("*Intensity_Mean_Ch=", pars$dead_dye_channel, "_Img=1")
sum_red_csv <- ldply(stat_folders, read_ims_csv, pattern=pat)
# import area
pat = "*Area"
area_csv <- ldply(stat_folders, read_ims_csv, pattern=pat)
# import position
pat = "*Position"
pos_csv <- ldply(stat_folders, read_ims_csv, pattern=pat)


live_deadROI <- cbind(volume_csv[,c("Volume","Time", "TrackID", "ID")], 
                      sum_red_csv[,c("Intensity Mean")], 
                      area_csv[,c("Area")],
                      pos_csv[,c("Position X","Position Y","Position Z","filename")])
colnames(live_deadROI) <- c("Volume","Time","TrackID","ID","red_sum","area", "pos_x","pos_y","pos_z", "filename")


### Convert the filename to the same format in both datasets (master and metadata)
live_deadROI$basename <- gsub("_org_Position.csv", "",live_deadROI$filename, perl=TRUE)
live_deadROI$basename=basename(live_deadROI$basename)
##Join the information of metadata to master:
live_deadROI<-left_join(live_deadROI, metadata)

## make TrackID unique:
category <- as.factor(live_deadROI$basename)
ranks <- rank(-table(category), ties.method="first")
ranks <- as.data.frame(ranks)
ranks$basename <- row.names(ranks)
live_deadROI <- left_join(live_deadROI, ranks)  ## plot with all the tracks together 

live_deadROI$Track2 <- with(live_deadROI, interaction(TrackID, ranks))
live_deadROI$Track2 <- gsub(".", '', live_deadROI$Track2, fixed = T)
live_deadROI$Track2 <- as.numeric(as.character(live_deadROI$Track2))


detach(package:plyr)
library(dplyr)
live_deadROI$red_sum <- live_deadROI$red_sum*live_deadROI$Volume
live_deadROI1 <-live_deadROI %>% 
  group_by(Track2, Time, organoid_line, tcell_line,exp_nr, well, date) %>% 
  summarise(Volume = sum(Volume), sum_red= sum(red_sum), area=sum(area), pos_x=mean(pos_x), pos_y=mean(pos_y), pos_z=mean(pos_z))

live_deadROI3 <- live_deadROI1 ## plot with all the tracks together 
live_deadROI3$red <- live_deadROI3$sum_red/live_deadROI3$Volume
live_deadROI3 <- live_deadROI3[complete.cases(live_deadROI3), ]
## set time in hours
live_deadROI3$Time2<-(live_deadROI3$Time-1)/2


### Quantify the dead cell dye intensity per well
live_deadROI32<-live_deadROI3%>% 
  group_by(Track2) %>%arrange(Track2)%>% filter(n() > 40) ## filter cells with at least 40 timepoints
### pulled all the dead cell dye signal from the well
live_deadROI7 <-live_deadROI32 %>% 
  group_by(Time2, organoid_line, tcell_line,exp_nr, well, date) %>% 
  summarise(Volume = sum(Volume), red= sum(sum_red))
live_deadROI7$red<-live_deadROI7$red/live_deadROI7$Volume

# setwd(paste0( file_location,"analysis/"))  ## set here your working directory that contains the example dataset

### SAVE dataframe with all the well values for processing in a different script
saveRDS(live_deadROI7, file = paste0(output_dir,"Full_well_death_dynamics"))

### Process the death dynamics per individual organoid
library(scales)
live_deadROI3<-live_deadROI3%>% ungroup()%>%
  mutate(red_rs=rescale(red, to=c(0,100)))
live_deadROI3<-live_deadROI3%>% 
  group_by(Track2) %>%arrange(Track2)%>% filter(n() > 40)
temp1 <- aggregate(red_rs ~ Track2+organoid_line+tcell_line+exp_nr+well+date, data = live_deadROI3, max)  ##calculate the max red of each track
colnames(temp1) [length(names(temp1))] <- "max_red"
live_deadROI4 <- merge(temp1, live_deadROI3)
live_deadROI6 <- subset(live_deadROI4 , Volume>1000)
live_deadROI6deadT0 <-live_deadROI6%>%group_by(Track2)%>%filter(Time==min(Time) & red_rs<10 )
live_deadROI6 <-live_deadROI6%>%filter(Track2 %in% live_deadROI6deadT0$Track2)
## Filter For the organoids that increase in red dead cell dye, subsitute by the max:
library(tibble)  # for `rownames_to_column` and `column_to_rownames`
temp1 <- live_deadROI6%>%rownames_to_column('row') %>%arrange(Time) %>% group_by(Track2) %>% filter(row_number() <= which.max(red))%>%column_to_rownames('row')
# Filter the rows that are missing (after reaching the max):
temp2<-subset(live_deadROI6, !row.names(live_deadROI6)%in%row.names(temp1))
temp2$red_rs<-temp2$max_red
live_deadROI6<-rbind(temp1, temp2)
### plot to check outcome
library(ggplot2)
Plot <- ggplot(live_deadROI6, aes(Time2,red_rs, color = Track2, group = Track2)) + 
  geom_smooth(method="loess", size = 1, se=F, span=1) +
  theme_bw() + 
  ylab("dead dye intensity") + 
  xlab("Time (hours)") +
  theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), legend.text=element_text(size= 10))+
  labs(color = "Organoid")+
  facet_grid(organoid_line~well, scales = "free")
ggtitle("Individual org increase in dead dye intensity TEG")

Plot
saveRDS(live_deadROI6, file = paste0(output_dir,"Individual_orgs_death_dynamics"))  ### save here a dataframe with all the organoids values
ggsave(paste0(output_dir,"Individual_orgs_death_dynamics.png"), device="png")


library(dplyr)
library(scales)
library(MESS)
##Import dataframe live_deadROI7 for each experiment. Combine several if necessary
Combi<-live_deadROI7

## Rescale the mean dead cell dye per experiment
Combi<-Combi%>% 
  group_by(exp_nr) %>% 
  mutate(red_sc=rescale(red, to=c(0,100)))%>%ungroup()
## for each well normalize the dead cell dye mean intensity to timepoint 1
Combi<-Combi%>%group_by( organoid_line,tcell_line , exp_nr , well, date)%>%arrange(Time2)%>%mutate(red_norm=red_sc - first(red_sc))%>%filter(Time2<23)

## calculate area under the curve for each donor/T cell co-culture
Area_under_Curve<-Combi%>%
  group_by( organoid_line,tcell_line , exp_nr , well, date)%>%summarise(auc=auc(Time2,red_norm))%>%ungroup()
## calculate delta increase for each donor/T cell co-culture
delta <-Combi%>%group_by( organoid_line,tcell_line , exp_nr , well, date)%>%arrange(Time2)%>%summarize(delta= last(red_norm) - first(red_norm))%>%ungroup()
## combine both delta and AUC in one dataframe
AUC_delta<-left_join(Area_under_Curve,delta,by=c("organoid_line" ,"tcell_line" ,"exp_nr"  , "well" , "date"  ))
library(xlsx)
## export calculation for plot 
write.xlsx(AUC_delta, file = paste0(output_dir,"AUC_delta.xlsx"))
