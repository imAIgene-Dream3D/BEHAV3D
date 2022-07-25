set.seed(1234)
getwd()
library(plyr)
library(readr)
library(dplyr)
library(yaml)
library(optparse)
library(ggplot2)

### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  ### !!!!!! Change the path to the BEHAV3D_config file here if running the code in RStudio !!!!!!
  pars = yaml.load_file("/Users/samdeblank/OneDrive - Prinses Maxima Centrum/github/BEHAV3D-2.0/demos/combined_demo_data/BEHAV3D_config.yml")
} else {
  option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, 
                help="Path to the BEHAV3D config file", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  if (is.null(opt$config)){
    print_help(opt_parser)
    stop("Config file -c|--config, must be supplied", call.=FALSE)
  }
  pars = yaml.load_file(opt$config)
}

## Create output directory
pars$data_dir = paste0(pars$data_dir,"/")
output_dir=paste0(pars$output_dir,"/organoid_dynamics/results/")
dir.create(output_dir, recursive=TRUE)

### Import file-specific metadata for all images used in this analysis.
pat = pars$metadata_csv
metadata=read.csv(pars$metadata_csv, sep="\t")

if ( any(is.na(metadata$organoid_stats_folder)) ){
  metadata$organoid_stats_folder=apply(metadata, 1, function(x) paste0(pars$data_dir, x["basename"], "_Statistics"))
}

### Function to import organoid data specifically from Imaris generated csv files
read_ims_csv <- function(metadata_row, pattern) {
  read_plus <- function(flnm, stat_folder) {
    read_csv(flnm, skip = 3, col_types = cols(TrackID= col_character())) %>% 
      mutate(filename = flnm, stat_folder=stat_folder) 
  }
  basename=gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", metadata_row[['basename']])
  pattern <- paste0(basename,".*", pattern, ".*")
  print(pattern)
  pattern_file <- list.files(path = metadata_row[['stats_folder']], pattern = pattern, full.names=TRUE)
  if (identical(pattern_file, character(0))){
    print(paste("No file with pattern '", pattern, "' found for", metadata_row['stats_folder']))
  } 
  ims_csv <- read_plus(pattern_file, metadata_row['stats_folder'])
  return(ims_csv)
}

stat_folders <- metadata[c("basename", "organoid_stats_folder")]
colnames(stat_folders) <- c("basename", "stats_folder")

pat = "Volume"
volume_csv <- do.call("rbind", apply(stat_folders, 1, read_ims_csv, pattern=pat))
# import sum_red
datalist = list()
for (i in 1:length(stat_folders$stats_folder)){
  pat=paste0("Intensity_Mean_Ch=", metadata$dead_dye_channel[i], "_Img=1")
  img_csv = read_ims_csv(stat_folders[i,], pattern=pat)
  if (!identical(img_csv, character(0))){
    datalist[[i]]=img_csv
  }
}
sum_red_csv=do.call(rbind, datalist)

# import area
pat = "Area"
area_csv <- do.call("rbind", apply(stat_folders, 1, read_ims_csv, pattern=pat))
# import position
pat = "Position"
pos_csv <- do.call("rbind", apply(stat_folders, 1, read_ims_csv, pattern=pat))

live_deadROI <- cbind(volume_csv[,c("Volume","Time", "TrackID", "ID")], 
                      sum_red_csv[,c("Intensity Mean")], 
                      area_csv[,c("Area")],
                      pos_csv[,c("Position X","Position Y","Position Z","filename", "stat_folder")])
colnames(live_deadROI) <- c("Volume","Time","TrackID","ID","dead_dye_mean","area", "pos_x","pos_y","pos_z", "filename", "organoid_stats_folder")

### Join the information of metadata to master:
live_deadROI<-left_join(live_deadROI, metadata)

### Make TrackID unique for each file:
category <- as.factor(live_deadROI$basename)
ranks <- rank(-table(category), ties.method="first")
ranks <- as.data.frame(ranks)
ranks$basename <- row.names(ranks)
live_deadROI <- left_join(live_deadROI, ranks)  ## plot with all the tracks together 
live_deadROI$TrackID2 <- factor(paste(live_deadROI$ranks, live_deadROI$TrackID, sep="_"))

live_deadROI$Original_TrackID <- live_deadROI$TrackID
live_deadROI$TrackID<-live_deadROI$TrackID2
live_deadROI$TrackID2<-NULL

detach(package:plyr)
library(dplyr)

### Some organoids fall apart into separated segments, combine them and calculate total dead dye signal
live_deadROI$dead_dye_sum <- live_deadROI$dead_dye_mean*live_deadROI$Volume
live_deadROI1 <-live_deadROI %>% 
  group_by(TrackID, Time, organoid_line, tcell_line,exp_nr, well, date) %>% 
  summarise(Volume = sum(Volume), dead_dye_sum= sum(dead_dye_sum), area=sum(area), pos_x=mean(pos_x), pos_y=mean(pos_y), pos_z=mean(pos_z))

live_deadROI3 <- live_deadROI1 ## plot with all the tracks together 
live_deadROI3$dead_dye_mean <- live_deadROI3$dead_dye_sum/live_deadROI3$Volume
live_deadROI3 <- live_deadROI3[complete.cases(live_deadROI3), ]

### Set time in hours
live_deadROI3$Time2<-(live_deadROI3$Time-1)/2

### Quantify the dead cell dye intensity per well
### Filter organoids with at least pars$organoid_min_track_length timepoints
live_deadROI32<-live_deadROI3%>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(n() > pars$organoid_min_track_length) 

### Combine all the dead cell dye signal from a well
live_deadROI7 <-live_deadROI32 %>% 
  group_by(Time, Time2, organoid_line, tcell_line,exp_nr, well, date) %>% 
  summarise(Volume = sum(Volume), dead_dye_sum= sum(dead_dye_sum))
live_deadROI7$dead_dye_mean<-live_deadROI7$dead_dye_sum/live_deadROI7$Volume

### Plot the dead dye intensity per well over time

ggplot(live_deadROI7, aes(Time,dead_dye_sum)) + 
  geom_smooth(method="loess", size = 1, se=F, span=1) +
  geom_point(size=0.5) +
  theme_bw() + 
  ylab("Sum of dead dye intensity") + 
  xlab("Timepoints") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), legend.text=element_text(size= 10))+
  labs(color = "Organoid")+
  facet_grid(organoid_line~well, scales = "free")+
  ggtitle("Total sum of dead dye intensity in organoids per well")

ggplot(live_deadROI7, aes(Time,dead_dye_mean)) + 
  geom_smooth(method="loess", size = 1, se=F, span=1) +
  geom_point(size=0.5) +
  theme_bw() + 
  ylab("Mean of dead dye intensity") + 
  xlab("Timepoints") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), legend.text=element_text(size= 10))+
  labs(color = "Organoid")+
  facet_grid(organoid_line~well, scales = "free")+
  ggtitle("Mean of dead dye intensity in organoids per well")

ggsave(paste0(output_dir,"Full_well_death_dynamics.pdf"), device="pdf")
### SAVE dataframe with all the well values for processing in a different script
saveRDS(live_deadROI7, file = paste0(output_dir,"Full_well_death_dynamics"))

live_deadROI7_exp_dur = live_deadROI7[live_deadROI7$Time <= pars$organoid_exp_duration,]
ggplot(live_deadROI7_exp_dur, aes(Time,dead_dye_mean)) + 
  geom_smooth(method="loess", size = 1, se=F, span=1) +
  geom_point(size=0.5) +
  theme_bw() + 
  ylab("Mean of dead dye intensity") + 
  xlab("Timepoints") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10), axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), legend.text=element_text(size= 10))+
  labs(color = "Organoid")+
  facet_grid(organoid_line~well, scales = "free")+
  ggtitle("Mean of dead dye intensity in organoids per well")

saveRDS(live_deadROI7_exp_dur, file = paste0(output_dir,"Exp_well_death_dynamics"))
ggsave(paste0(output_dir,"Exp_well_death_dynamics.pdf"), device="pdf")

### Process the death dynamics per individual organoid
library(scales)
live_deadROI3<-live_deadROI3%>% ungroup()%>%
  mutate(dead_dye_mean_rescaled=rescale(dead_dye_mean, to=c(0,100)))
live_deadROI3<-live_deadROI3%>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(n() > pars$organoid_min_track_length)
temp1 <- aggregate(dead_dye_mean_rescaled ~ TrackID+organoid_line+tcell_line+exp_nr+well+date, data = live_deadROI3, max)  ##calculate the max red of each track
temp2 <- aggregate(dead_dye_mean ~ TrackID+organoid_line+tcell_line+exp_nr+well+date, data = live_deadROI3, max)
colnames(temp1) [length(names(temp1))] <- "max_dead_dye_mean_rescaled"
colnames(temp2) [length(names(temp1))] <- "max_dead_dye_mean"

temp_merge <-merge(temp1, temp2)
### Filter out organoids that have dead dye signal above 'organoid_dead_dye_threshold' (e.g. dead) and are bigger than 'organoid_min_volume'
live_deadROI4 <- merge(temp_merge, live_deadROI3)

live_deadROI6 <- left_join(live_deadROI4, metadata)
live_deadROI6aliveT0 <-live_deadROI6%>%group_by(TrackID)%>%filter(Time==min(Time) & dead_dye_mean<organoid_dead_dye_threshold & Volume>pars$organoid_min_volume)
live_deadROI6 <-live_deadROI6%>%filter(TrackID %in% live_deadROI6aliveT0$TrackID)

## Filter for the organoids that increase in red dead cell dye, substitute by the max:
library(tibble)  # for `rownames_to_column` and `column_to_rownames`
temp1 <- live_deadROI6%>%rownames_to_column('row') %>%arrange(Time) %>% group_by(TrackID) %>% filter(row_number() <= which.max(dead_dye_mean))%>%column_to_rownames('row')
# Filter the rows that are missing (after reaching the max):
temp2<-subset(live_deadROI6, !row.names(live_deadROI6)%in%row.names(temp1))
temp2$dead_dye_mean_rescaled<-temp2$max_dead_dye_mean_rescaled
live_deadROI6<-rbind(temp1, temp2)

### plot to check outcome
ggplot(live_deadROI6, aes(Time,dead_dye_mean_rescaled, color = TrackID, group = TrackID)) + 
  geom_smooth(method="loess", size = 1, se=F, span=1) +
  theme_bw() + 
  ylab("dead dye intensity") + 
  xlab("Time (hours)") +
  theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), legend.text=element_text(size= 10))+
  labs(color = "Organoid")+
  facet_grid(interaction(exp_nr,well,organoid_line)  ~ tcell_line, scales = "free")+
  ggtitle("Individual org increase in dead dye intensity TEG")

saveRDS(live_deadROI6, file = paste0(output_dir,"Full_individual_orgs_death_dynamics"))  ### save here a dataframe with all the organoids values
ggsave(paste0(output_dir,"Full_individual_orgs_death_dynamics.pdf"), device="pdf")

live_deadROI6_exp_dur = live_deadROI6[live_deadROI6$Time<=pars$organoid_exp_duration,]
### plot to check outcome
ggplot(live_deadROI6_exp_dur, aes(Time,dead_dye_mean_rescaled, color = TrackID, group = TrackID)) + 
  geom_smooth(method="loess", size = 1, se=F, span=1) +
  theme_bw() + 
  ylab("dead dye intensity") + 
  xlab("Time (hours)") +
  theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), legend.text=element_text(size= 10))+
  labs(color = "Organoid")+
  facet_grid(interaction(exp_nr,well,organoid_line)  ~ tcell_line, scales = "free")+
  ggtitle("Individual org increase in dead dye intensity TEG")

saveRDS(live_deadROI6_exp_dur, file = paste0(output_dir,"Exp_individual_orgs_death_dynamics"))  ### save here a dataframe with all the organoids values
ggsave(paste0(output_dir,"Exp_individual_orgs_death_dynamics.pdf"), device="pdf")

library(dplyr)
library(scales)
library(MESS)

##Import dataframe live_deadROI7 for each experiment. Combine several if necessary
Combi<-live_deadROI7

## Rescale the mean dead cell dye per experiment
Combi<-Combi%>% 
  group_by(exp_nr) %>% 
  mutate(dead_dye_mean_rescaled=rescale(dead_dye_mean, to=c(0,100)))%>%ungroup()
## for each well normalize the dead cell dye mean intensity to timepoint 1
Combi<-Combi%>%group_by( organoid_line,tcell_line , exp_nr , well, date)%>%arrange(Time2)%>%mutate(dead_dye_mean_norm=dead_dye_mean_rescaled - first(dead_dye_mean_rescaled))%>%filter(Time<pars$organoid_exp_duration)

## calculate area under the curve for each donor/T cell co-culture
Area_under_Curve<-Combi%>%
  group_by( organoid_line,tcell_line , exp_nr , well, date)%>%summarise(auc=auc(Time2,dead_dye_mean_norm))%>%ungroup()
## calculate delta increase for each donor/T cell co-culture
delta <-Combi%>%group_by( organoid_line,tcell_line , exp_nr , well, date)%>%arrange(Time2)%>%summarize(delta= last(dead_dye_mean_norm) - first(dead_dye_mean_norm))%>%ungroup()
## combine both delta and AUC in one dataframe
AUC_delta<-left_join(Area_under_Curve,delta,by=c("organoid_line" ,"tcell_line" ,"exp_nr"  , "well" , "date"  ))
library(xlsx)
## export calculation for plot 
write.xlsx(AUC_delta, file = paste0(output_dir,"AUC_delta.xlsx"))

live_deadROI6$dead<-ifelse(live_deadROI6$dead_dye_mean<live_deadROI6$organoid_dead_dye_threshold,0,1)
library(tidyr)
### fill missing values with NA
live_deadROI6<-live_deadROI6 %>%group_by(organoid_line, tcell_line,exp_nr, well, date)%>%arrange(Time, .by_group = TRUE)%>% complete(Time, nesting(TrackID))
## refill NA with previous value
live_deadROI6<-live_deadROI6 %>%group_by(organoid_line, tcell_line,exp_nr, well, date)%>%arrange(Time, .by_group = TRUE)%>%fill(names(live_deadROI6))  
## select the dead cells
live_deadROI4_dead <-live_deadROI6%>%group_by(organoid_line, tcell_line,exp_nr, well, date, TrackID)%>%filter(max_dead_dye_mean>organoid_dead_dye_threshold )
# live_deadROI4_dead$track_exp <- with(live_deadROI4_dead , interaction(TrackID,Time))
live_deadROI4_live <-live_deadROI6%>%group_by(TrackID, organoid_line, tcell_line,exp_nr, well, date)%>%filter(max_dead_dye_mean<=organoid_dead_dye_threshold ) ## select the live cells
temp2 <-live_deadROI4_dead%>%group_by(TrackID, organoid_line, tcell_line,exp_nr, well, date)## group by track ID
temp2 <- temp2%>%arrange(Time, .by_group = TRUE) 
### After timepoint of dying set further timepoints to dead=1 as well
live_deadROI5_dead <- temp2 %>% group_by(TrackID, organoid_line, tcell_line,exp_nr, well, date) %>% filter(row_number() >= which.max(dead)) ## stop at max dead
live_deadROI6 <- temp2 %>% group_by(TrackID, organoid_line, tcell_line,exp_nr, well, date) %>% mutate(dead=replace(dead, row_number() > which.max(dead), 1))
live_deadROI6 <- live_deadROI6 %>% group_by(TrackID, organoid_line, tcell_line,exp_nr, well, date) %>% arrange(Time, .by_group = TRUE) %>% mutate(dead=replace(dead, row_number() > which.max(dead), 1))

# test1 <- live_deadROI4_dead%>%filter(!track_exp%in% live_deadROI5_dead$track_exp)  ###find all the rows that are beyond max red
# test1$dead<-1 ## set as dead beyond max red
# live_deadROI6_2<-rbind(live_deadROI4_live,live_deadROI5_dead[,c(-16)],test1[,c(-16)])  ## bind dataframe again
## calculate the number of organoids at T0 to calculate the percentage
live_deadROI6_n_t0<-live_deadROI6%>%filter(Time==1)%>%group_by(organoid_line, tcell_line,exp_nr, well, date)%>%summarise(starting_nr_organoids= n())
live_deadROI6_2<-left_join(live_deadROI6,live_deadROI6_n_t0)
live_deadROI6_2$starting_nr_organoids<-as.numeric(live_deadROI6_2$starting_nr_organoids)
# for each dead calculate the percentage of dead cells, if no, then 0 
live_deadROI6_per<-live_deadROI6_2%>%group_by(Time, organoid_line, tcell_line,exp_nr, well, date, .drop=FALSE)%>%summarize(n=sum(dead==1), starting_nr_organoids=mean(starting_nr_organoids))%>%mutate(perc = n*100 / starting_nr_organoids)

library(ggplot2)
# live_deadROI6_per$Time<-live_deadROI6_per$Time.x/2
### normalize to the initial number of dead organoids
perc_dead<- live_deadROI6_per%>%
  group_by(organoid_line, tcell_line,exp_nr, well, date)%>%
  arrange(Time)%>%
  mutate(perc_dead_norm=perc - first(perc))

p1 <- ggplot(perc_dead, aes(as.numeric(Time), perc_dead_norm)) + 
  geom_line()+
  # geom_smooth(method = "loess",size = 0.5, se = T, alpha=0.3, span=1)+
  theme_bw() + 
  ylab(paste0("Percentage of dead organoids (n=", live_deadROI6_n_t0$starting_nr_organoids, ")")) + 
  # coord_cartesian(ylim=c(0,100), xlim=c(0,pars$organoid_exp_duration))+
  scale_y_continuous(limit=c(0,100),oob=squish)+
  scale_x_continuous(limit=c(0,max(perc_dead$Time)))+
  xlab("Time") +
  facet_grid(interaction(exp_nr,organoid_line)  ~ tcell_line, scales = "free")+
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_text(size = 10), 
        legend.text=element_text(size= 10))+
  ggtitle("Percentage of dead organoids over time")+theme(aspect.ratio=1)

p1

saveRDS(perc_dead, file = paste0(output_dir,"Full_percentage_dead_org_over_time.rds"))  ### save here a dataframe with all the organoids values
ggsave(paste0(output_dir,"Full_percentage_dead_org_over_time.pdf"), device="pdf")

perc_dead_exp_dur <- perc_dead %>%
  filter(Time<=pars$organoid_exp_duration)

p1 <- ggplot(perc_dead_exp_dur, aes(as.numeric(Time), perc_dead_norm)) + 
  geom_line()+
  # geom_smooth(method = "loess",size = 0.5, se = T, alpha=0.3, span=1)+
  theme_bw() + 
  ylab(paste0("Percentage of dead organoids (n=", live_deadROI6_n_t0$starting_nr_organoids, ")")) + 
  # coord_cartesian(ylim=c(0,100), xlim=c(0,pars$organoid_exp_duration))+
  scale_y_continuous(limit=c(0,100),oob=squish)+
  scale_x_continuous(limit=c(0,min(max(perc_dead_exp_dur$Time),pars$organoid_exp_duration)))+
  xlab("Time") +
  facet_grid(interaction(exp_nr,organoid_line)  ~ tcell_line, scales = "free")+
  theme(axis.text.x = element_text(size=10), 
        axis.text.y = element_text(size=10), 
        axis.title.y = element_text(size = 10), 
        axis.title.x = element_text(size = 10), 
        legend.text=element_text(size= 10))+
  ggtitle("Percentage of dead organoids over time")+theme(aspect.ratio=1)

p1

saveRDS(perc_dead_exp_dur, file = paste0(output_dir,"Exp_percentage_dead_org_over_time.rds"))  ### save here a dataframe with all the organoids values
ggsave(paste0(output_dir,"Exp_percentage_dead_org_over_time.pdf"), device="pdf")
