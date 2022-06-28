set.seed(12)
getwd()
library(plyr)
library(readr)
library(dplyr)
library(reshape2)
library(zoo)
library(spatstat)
library(sp)
library(stats)
library(yaml)

### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  ### !!!!!! Change the path to the BEHAV3D_config file here if running the code in RStudio !!!!!!
  pars = yaml.load_file("/Users/samdeblank/surfdrive/Shared/Dream3DLab (Groupfolder)/1.Projects/AIM_ALLImmune/3.Analysis/BEHAV3D_analysis/AIM_MB2_Exp21_Tcellstats/BEHAV3D_config_combined.yml")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  pars <- yaml.load_file(args[1])
}

### Function to count the number of tracks in the dataset
count_tracks = function (track_table){
  nr_tracks_unfilt=track_table
  nr_tracks_unfilt$name = paste(nr_tracks_unfilt$organoid_line, nr_tracks_unfilt$exp_nr, nr_tracks_unfilt$well)
  nr_tracks_unfilt = nr_tracks_unfilt %>% 
    group_by(name, organoid_line) %>% 
    dplyr::summarize(
      nr_tracks=length(unique(TrackID))
    ) %>% 
    ungroup()
}

pars$data_dir = paste0(pars$data_dir,"/")
output_dir=paste0(pars$output_dir,"/")
model_path <- pars$randomforest

###############################
######### Data import #########
###############################

### Import file-specific metadata for all images used in this analysis.
pat = pars$metadata_csv
metadata=read.csv(pars$metadata_csv, sep="\t", check.names=FALSE)

track_counts=metadata
track_counts$name = paste(metadata$organoid_line, metadata$exp_nr, metadata$well)
track_counts=track_counts[,c("name", "organoid_line")]

### Check if folder with statistics is named in the metadata table, if not, try default naming of data basename + "_Statistics"
if ( any(is.na(metadata$stat_folder)) ){
  metadata$stat_folder=apply(metadata, 1, function(x) paste0(pars$data_dir, x["basename"], "_Statistics"))
}

### Function to import organoid data specifically from Imaris generated csv files
read_ims_csv <- function(stat_folder, pattern) {
  read_plus <- function(flnm) {
    read_csv(flnm, skip = 3, col_types = cols()) %>% 
      mutate(filename = flnm) 
  }
  pattern_file <- list.files(path = stat_folder, pattern = pattern, full.names=TRUE)
  if (identical(pattern_file, character(0))){
    print(paste("No file with pattern '", pattern, "' found for", stat_folder))
  } 
  ims_csv <- read_plus(pattern_file)
  return(ims_csv)
}

stat_folders <- metadata$stat_folder

# import Displacement^2
pat = "*Displacement\\^2"
displacement=ldply(stat_folders, read_ims_csv, pattern=pat)

# import Speed
pat = "*Speed"
speed <- ldply(stat_folders, read_ims_csv, pattern=pat)

# import mean dead dye intensity values
datalist = list()
for (i in 1:length(stat_folders)){
  pat=paste0("*Intensity_Mean_Ch=", metadata$dead_dye_channel[i], "_Img=1")
  print(pat)
  img_csv = read_ims_csv(stat_folders[i], pat)
  if (!identical(img_csv, character(0))){
    datalist[[i]]=img_csv
  }
}
red_lym=do.call(rbind, datalist)

# import Minimal distance to organoids
datalist = list()
for (i in 1:length(stat_folders)){
  pat=paste0("*Intensity_Min_Ch=", metadata$organoid_distance_channel[i], "_Img=1")
  # print(pat)
  datalist[[i]]=read_ims_csv(stat_folders[i], pat)
}
dist_org=do.call(rbind, datalist)

# import Position
pat = "*Position.csv"
pos <- ldply(stat_folders, read_ims_csv, pattern=pat)

### Join all Imaris information
master <- cbind(
  displacement[,c("Displacement^2","Time","TrackID" ,"ID")], 
  speed[,c("Speed" )], 
  dist_org[,c("Intensity Min")], 
  red_lym[,c("Intensity Mean")], 
  pos[,c("Position X" ,"Position Y" ,"Position Z","filename")]
)

### Get the basename from the filename for combination with metadata
master$basename <- gsub("_Position.csv", "", master$filename, perl=TRUE)
master$basename=basename(master$basename)
colnames(master) <- c("displacement","Time","TrackID","ID","speed","dist_org","red_lym","X-pos","Y-pos","Z-pos", "filename", "basename")

### Join the information of metadata to master:
master<-left_join(master, metadata)

### Create a unique TRACKID. 
### Each file processes with Imaris has separate TRACKIDs and these must be made unique before merging
category <- as.factor(master$filename)
ranks <- rank(-table(category), ties.method="first")
ranks <- as.data.frame(ranks)
ranks$filename <- row.names(ranks)
master <- left_join(master, ranks) 
master$TrackID2 <- with(master, interaction(ranks, TrackID, sep="_"))

### Remove the variable TrackID and only use unique TrackID2 (unique identifier instead)
master$TrackID<-master$TrackID2
master$TrackID2<-NULL
### Remove filename
master$filename<-NULL

### save RDS for later use (e.g. Backprojection of classified TrackIDs)
saveRDS(master, paste0(output_dir,"raw_tcell_track_data.rds"))

track_counts=left_join(track_counts, count_tracks(master))
colnames(track_counts)[colnames(track_counts)=="nr_tracks"]="unfiltered"

###############################
####### Data processing #######
###############################

detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)
pars$exp_duration

master <- master[which(master$Time<=pars$exp_duration), ] ##Make sure that all the time-series have the same length, in this case 10hours

track_counts=left_join(track_counts, count_tracks(master))
colnames(track_counts)[colnames(track_counts)=="nr_tracks"]="filt_exp_duration"

### Perform check for duplicates, should be empty
data_dup <- master%>%group_by(Time)%>%
  dplyr::count(TrackID) %>% 
  filter(n > 1) %>% 
  select(-n)

if (dim(data_dup)[1]!=0){
  stop("There are duplicates in the data, which should not be the case. Stopping execution...")
}

### For each well separately:
### Estimate if a T cells interacts with another T cell 
### This is done by calculating the minimal distance to the nearest neighboring T cell. 

List = list() ## create list for each timepoint
List2 = list() ## create list for each experiment

### For loop to look for the distance to the nearest neighbor at each timepoint and each experiment
for ( m in unique(master$well)){
  distance_1<-master[which(master$well==m),]
  List = list()
  for (i in unique(distance_1$Time)){
    distanceDFi <- distance_1[distance_1$Time==i,]
    distanceDF2<-as.data.frame(distanceDFi[,c("Time","X-pos","Y-pos","Z-pos")])
    coordin<-ppx(distanceDF2, coord.type=c("t","s","s", "s")) ## create a multidimensional space-time pattern
    dist<- nndist(coordin)
    distanceDFi_dist<-cbind(distanceDFi, dist)
    List[[length(List)+1]] <-distanceDFi_dist   ## store to list object
  }
  master_distance <- data.frame(do.call(rbind, List)) ## convert List to dataframe
  List2[[length(List2)+1]] <-master_distance 
}

master_dist<-do.call(rbind, List2)
colnames(master_dist)[which(names(master_dist) == "dist")] <- "nearest_Tcell"
master_dist$contact_lym <- ifelse(master_dist$nearest_Tcell<master_dist$tcell_contact_thr,1,0)

library(reshape2)
library(zoo)

### Since not all the tracks are tracked at all timepoints,
### Interpolate missing values and fill the NA (for each variable)

### Select the variables for which we need to interpolate NAs (numeric)
column_names<-names(master_dist)
column_names <- c("displacement", "dist_org", "red_lym", "contact_lym")

### Create a first dataset with refilled values for speed:
time_series<-acast(master_dist, Time ~ TrackID, value.var='speed',fun.aggregate = mean)

### rownames timepoints:
row.names(time_series)<-unique(master_dist$Time)

### Get rid of NA by interpolation
time_series_zoo<-zoo(time_series, row.names(time_series))
time_series_zoo<-na.approx(time_series_zoo) ## replace by interpolated value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series)
data<-time_series2[complete.cases(time_series2), ] 
colnames(data)<-c("Time", "TrackID", "speed")

### Store this data for calculating lagged speed later:
time_series2_speed<-data

for (i in column_names){
  time_series<-acast(master_dist, Time ~ TrackID, value.var=i,fun.aggregate = mean)
  row.names(time_series)<-unique(master_dist$Time)
  ### get rid of NA
  time_series_zoo<-zoo(time_series,row.names(time_series))
  time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
  time_series<-as.matrix(time_series_zoo)
  time_series2<-melt(time_series)
  new<-time_series2[complete.cases(time_series2), ] 
  data[ , ncol(data) + 1] <- new[3]                  # Append new column
  colnames(data)[ncol(data)] <- paste0(i)
}

library(dplyr)
### For cell interaction we need to consider the following:
### When two cells interact it is often the one cell moves and interacts with another one that is static
### In this case one might consider that only the motile cell is actively interacting and the static cell is just passively interacting
### To determine when a cell is actively interacting we measure for each cell what was its mean speed over the last 10 timepoints (20 mins)
time_series2_meanspeed <-time_series2_speed %>% 
  group_by(TrackID) %>%mutate(meanspeed=rollapply(speed,10,mean,align='right',fill=NA))

### Refill all missing values with the last value
time_series<-acast(time_series2_meanspeed, Time ~ TrackID, value.var='meanspeed',fun.aggregate = mean)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.locf(time_series_zoo, fromLast=T) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2_meanspeed<-melt(time_series)
colnames(time_series2_meanspeed)<-c("Time", "TrackID", "meanspeed")

### Remove last NAs
time_series2_meanspeed<-na.omit(time_series2_meanspeed)

### Create a dataframe with all the variables with corrected missing values
master_corrected <- data

### Join the information on the cell line, experiment number and well number
master_temp<- master[c("TrackID", colnames(metadata))]
#master_temp<- master[c("TrackID2", colnames(metadata)[-1])]
master_temp<-master_temp[!duplicated(master_temp$TrackID),]
master_corrected<- left_join(master_corrected, master_temp, by=c("TrackID"))

### Merge the information for the mean speed over the last 20 mins
master_corrected1<- merge(master_corrected, time_series2_meanspeed, by = c("Time","TrackID"))

### Update the binary variable for contact with organoids
### It can vary between experiments depending on the intensity of the T cells or organoids. 
### Check the threshold of contact in the imaging data and update in the metadata csv
master_corrected1$contact <- ifelse(master_corrected1$dist_org>master_corrected1$organoid_contact_threshold, 0,1)

### Plot the number of touching vs. non-touching T cells
library(ggplot2)
ggplot(master_corrected1, aes(x=contact, color=as.factor(exp_nr))) +
  geom_histogram(fill="white", alpha=0.5, position="identity")+facet_grid(organoid_line~well, scales = "free")

ggsave(
  paste0(output_dir,"TouchingvsNontouching_distribution.png"), 
  device="png", height=210, width=297, units="mm"
)

### Remove contact threshold variable
master_corrected1$organoid_contact_threshold<-NULL

### For clustering it is necessary to compare T cell tracks that have a similar length. 
### For that we select cell track that have at least 100 timepoints. 
### Detach package  'plyr' as it can interfere with 'dplyr'
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

master_corrected2<-master_corrected1 %>% 
  group_by(TrackID) %>% arrange(TrackID)%>% filter(Time>00&Time<pars$exp_duration)%>% filter(n() >= pars$min_track_length)

track_counts=left_join(track_counts, count_tracks(master_corrected2))
colnames(track_counts)[colnames(track_counts)=="nr_tracks"]="filt_minLength"

### Create a variable for the relative Time
master_corrected2<-master_corrected2 %>% 
  group_by(TrackID) %>%arrange(Time)%>%mutate(Time2 = Time - first(Time))
### For the Tracks that have more then 100 timepoints filter only the first 100.
master_corrected2<-master_corrected2 %>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(Time2<pars$max_track_length)

### To exclude noise due to dead cells remove the dead t cells from the beginning
# master_corrected3 <- master_corrected2
red_lym_over_time=master_corrected2
red_lym_over_time$name=paste(master_corrected2$organoid_line, master_corrected2$exp_nr, master_corrected2$well)
ggplot(red_lym_over_time[red_lym_over_time$Time==1,], aes(x=Time, y=red_lym))+
  geom_violin(aes(fill=name))+
  geom_jitter()+
  facet_grid(~name)+
  theme_bw()

ggsave(
  paste0(output_dir,"RedLym_distribution.png"), 
  device="png", height=210, width=297, units="mm"
)

### Filter out T cells that are dead at the start of the experiment
master_corrected3deadT0 <-master_corrected2%>%group_by(TrackID)%>%filter((Time2==0) & red_lym<dead_dye_threshold )

master_corrected3 <-master_corrected2%>%filter(TrackID %in% master_corrected3deadT0$TrackID )

track_counts=left_join(track_counts, count_tracks(master_corrected3))
colnames(track_counts)[colnames(track_counts)=="nr_tracks"]="filt_DeadCellStart"

### Create a binary variable for live or dead cells:
master_corrected3$death<- ifelse(master_corrected3$red_lym<master_corrected3$dead_dye_threshold,0,1)

### Create a variable for cumulative interaction with organoids
master_corrected3<-master_corrected3 %>% 
  group_by(TrackID) %>%mutate(contact2=(ave(contact, cumsum(!contact), FUN = cumsum)))
### Create a variable for T cells interact with other T cells while in the environment
master_corrected3$contact_lym<- ifelse(master_corrected3$contact==1,0,master_corrected3$contact_lym)
### For T cells inteacting in the environment keep as "interacting" only cells that had a mean speed in the last 20 minutes that is in the upper quantile.
master_corrected3<-master_corrected3%>%group_by(exp_nr)%>%mutate(contact_lym=ifelse(meanspeed<quantile(meanspeed,p=0.75),0,contact_lym))

### Save processed data on the tcells for possible further analysis
saveRDS(master_corrected3, file = paste0(output_dir,"processed_tcell_track_data.rds"))

write.table(track_counts, file=paste0(output_dir, "nr_of_tracks_after_filtering.tsv"), sep="\t", row.names=FALSE)

library(reshape2)
melted_track_counts = melt(track_counts)
detach("package:reshape2", unload=TRUE)
ggplot(melted_track_counts, aes(x=name, y=value, fill=variable)) + 
  geom_col(width=0.75, position="dodge") + 
  theme_bw() + 
  # facet_wrap(~name, ncol=4)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), aspect.ratio=0.75) +
  ylab("# Tracks") +
  xlab("Experiment")

ggsave(
  paste0(output_dir,"NrCellTracks_filtering_perExp.png"), 
  device="png", height=210, width=297, units="mm"
)

ggplot(melted_track_counts, aes(x=name, y=value, fill=organoid_line)) + 
  facet_grid(~variable)+
  geom_bar(stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(paste0(output_dir,"NrCellTracks_filtering_perFilt.png"), device="png", height=210, width=297, units="mm")


####### IF MODEL_PATH IS DEFINED PERFORM CLASSIFICATION BASED ON EXISTING DATA
####### OTHERWISE, PERFORM UNSUPERVISED UMAP CLUSTERING
if (model_path != ""){
  print("#### Model_path defined, performing random forest classification...")
  ###############################
  ##### Behavior prediction #####
  ###############################
  
  library(scales)
  library(randomForest)
  
  # Load in the previously trained Random Forest model
  tryCatch(
    {
      load(model_path)
    },
    error=function(e) {
      message("Provided model path does not exist") 
      message("Please set the path to the provided Random Forest model or \ntrain it on your own reference set (see: 'train_randomforest')")
      message(paste("Provided model path:", model_path))
      message("\nHere's the original error message:")
      message(e)
      stop()
      # Choose a return value in case of error
    }
  )
  
  #### Import new dataset to predict behavior (e.g. import dataset called "master_corrected3_example")
  master_test<-master_corrected3

  ### Normalize the data
  master_test2<-master_test%>% ungroup()%>%
    # group_by(tcell_line, organoid_line, exp_nr, well) %>%
    group_by(exp_nr) %>%
    mutate(z.disp = (displacement-mean(displacement))/sd(displacement),z.speed = (speed-mean(speed))/sd(speed), z.red = (red_lym-mean(red_lym))/sd(red_lym))%>%
    mutate(q.disp=ifelse(z.disp>(quantile(z.disp, p=0.75)),z.disp,min(z.disp)), q.speed=ifelse(z.speed>(quantile(z.speed, p=0.75)),z.speed,min(z.speed)),q.red=ifelse(z.red>(quantile(z.red, p=0.75)),z.red,min(z.red)))%>%
    mutate(q.disp=rescale(q.disp, to=c(0,100)),q.speed=rescale(q.speed, to=c(0,100)),q.red=rescale(q.red, to=c(0,100)),s.contact=rescale(contact, to=c(0,1)),s.contact_lym=rescale(contact_lym, to=c(0,1)))%>%
    mutate(q.disp=q.disp/mean(quantile(q.disp, p=0.9999)),q.speed=q.speed/mean(quantile(q.speed, p=0.9999)),q.red=q.red/mean(quantile(q.red, p=0.9999)))%>%
    ungroup()
  
  ### Calculate time-series descriptive statistics
  test_dataset <- master_test2%>%group_by(TrackID)%>% arrange(Time)%>%
    summarise(mean_speed= mean(q.speed),median_speed= median(q.speed),speed_sd= sd(q.speed),q3_speed= quantile(q.speed,0.90),
              mean_displacement = mean(q.disp),median_displacement = median(q.disp),
              displacement_sd=sd(q.disp),q3_disp= quantile(q.disp,0.90),
              mean_red_lym = mean(q.red),red_lym_sd=sd(q.red),q3_red= quantile(q.red,0.90),
              contact=mean(s.contact),mean_contact2=mean(contact2),contact2=max(contact2))
  
  ### create a  matrix from the predictors
  test1<-as.matrix(test_dataset[,2:15]) 
  model_predict_test<-predict(model,test1,type="response")
  test_dataset_predicted<-cbind(model_predict_test,test_dataset)
  test_dataset_predicted$cluster<-test_dataset_predicted$model_predict
  
  ### Plot the proportion of cells for each behavioral signature
  ## Join information of cell ID
  classification<-test_dataset_predicted[,c(2,17)]
  master_classified<-left_join(master_test,classification, by="TrackID")
  cell_ID<-master_classified[!duplicated(master_classified$TrackID),c("TrackID","organoid_line","tcell_line","exp_nr", "well")] 
  
  
  classified_tracks<-test_dataset_predicted
  classified_tracks$cluster2<-classified_tracks$cluster
  classified_tracks<-left_join(classified_tracks,cell_ID)
  classified_tracks<-classified_tracks%>%arrange(cluster2)
  
  saveRDS(classified_tracks, file = paste0(output_dir,"classified_tcell_track_data.rds"))
  
  ### Quantify the number of cells per well
  Number_cell_exp<-classified_tracks%>%group_by(well, exp_nr, tcell_line, organoid_line)%>%
    summarise(total_cell = n())
  Percentage_clus<-left_join(Number_cell_exp,classified_tracks)
  Percentage_clus <- Percentage_clus%>%group_by(cluster2,tcell_line, well,exp_nr, organoid_line)%>%
    summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()
  
  ### Add clusters that don't exist for experiments to the percentage table
  for (cluster in unique(Percentage_clus$cluster2)){
    for (row in 1:nrow(Percentage_clus)) {
      tcell_line = as.character(Percentage_clus[row, "tcell_line"])
      organoid_line = as.character(Percentage_clus[row, "organoid_line"])
      well = as.character(Percentage_clus[row, "well"])
      exp = as.integer(Percentage_clus[row, "exp_nr"])
      total_cell = as.integer(Number_cell_exp[
        Number_cell_exp$tcell_line==tcell_line &
          Number_cell_exp$well==well &
          Number_cell_exp$exp_nr==exp &
          Number_cell_exp$organoid_line==organoid_line, "total_cell"
      ])
      if (
        dim(
          Percentage_clus[Percentage_clus$cluster2==cluster &
                          Percentage_clus$tcell_line==tcell_line &
                          Percentage_clus$well==well &
                          Percentage_clus$exp_nr==exp &
                          Percentage_clus$organoid_line==organoid_line,]
        )[1]==0
      ){
        extra_row=data.frame(
          cluster2=cluster,
          tcell_line=tcell_line,
          well=well,
          exp_nr=exp,
          organoid_line=organoid_line,
          total_cell=total_cell,
          num_cluster=0,
          percentage=0
        )
        Percentage_clus<-rbind(Percentage_clus, extra_row)
      }
    }
  }
  
  saveRDS(Percentage_clus, file = paste0(output_dir,"cluster_perc_tcell_track_data.rds"))
  
  ### Plot proportion per well and per cell type
  Percentage_clus$tcell_line = as.factor(Percentage_clus$tcell_line)
  
  Per<-ggplot(Percentage_clus, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
    geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
  Per <- Per + facet_grid(interaction(exp_nr, well, organoid_line)  ~ tcell_line)

  Per<-Per+theme_void() + scale_fill_manual(values=c("gold3",
                                                     "darkolivegreen3",
                                                     "seagreen3",
                                                     "forestgreen",
                                                     "dodgerblue",
                                                     "cyan1",
                                                     "indianred",
                                                     "firebrick",
                                                     "brown1"))+
    theme(aspect.ratio = 0.2,strip.text.x = element_text(angle = 90))
  
  Per
  
  ggsave(paste0(output_dir,"RF_ClassProp_WellvsCelltype.png"), device="png")
  ggsave(paste0(output_dir,"RF_ClassProp_WellvsCelltype.pdf"), device="pdf")
} else {
  print("#### Model_path NOT defined, performing UMAP clustering...")
  #######################################################
  ###########     UNSUPERVISED CLUSTERING     ###########
  #######################################################
  
  library(dplyr)
  library(dtwclust)
  library(stats)
  library(scales)
  
  master_processed<-master_corrected3%>% 
    group_by(exp_nr) %>% 
    mutate(z.disp = (displacement-mean(displacement))/sd(displacement),z.speed = (speed-mean(speed))/sd(speed), z.red = (red_lym-mean(red_lym))/sd(red_lym))%>%
    mutate(q.disp=ifelse(z.disp>(quantile(z.disp, p=0.75)),z.disp,min(z.disp)), q.speed=ifelse(z.speed>(quantile(z.speed, p=0.75)),z.speed,min(z.speed)),q.red=ifelse(z.red>(quantile(z.red, p=0.75)),z.red,min(z.red)))%>%
    mutate(q.disp=rescale(q.disp, to=c(0,1)),q.speed=rescale(q.speed, to=c(0,1)),q.red=rescale(q.red, to=c(0,1)),s.contact=rescale(contact, to=c(0,1)),s.contact_lym=rescale(contact_lym, to=c(0,1))) %>%
    mutate(q.disp=q.disp/mean(quantile(q.disp, p=0.9999999)),q.speed=q.speed/mean(quantile(q.speed, p=0.9999999)),q.red=q.red/mean(quantile(q.red, p=0.9999999)))%>%ungroup()
  
  ### Arrange data by time:
  master_processed<-master_processed%>%group_by(TrackID)%>%arrange(Time)
  master_processed$TrackID = as.character(master_processed$TrackID)
  # master_processed_ref<-master_processed_ref%>%group_by(TrackID)%>%arrange(Time)
  
  ###Split the data in a list of TrackIDs with multivariate data for each Track overtime
  list_multivariate <- split(master_processed[,c("q.disp", "q.speed", "q.red","s.contact", "s.contact_lym")],master_processed$TrackID) 
  # list_multivariate_ref <- split(master_processed_ref[,c("q.disp", "q.speed", "q.red","s.contact", "s.contact_lym")],master_processed_ref$TrackID) 
  
  ##Set up parallel working for big datasets
  # load parallel
  library(parallel)
  # create multi-process workers
  workers <- makeCluster(detectCores()-2)
  # load dtwclust in each one, and make them use 1 thread per worker
  invisible(clusterEvalQ(workers, {
    library(dtwclust)
    RcppParallel::setThreadOptions(1L)
  }))
  # register your workers, e.g. with doParallel
  require(doParallel)
  registerDoParallel(workers)
  
  ###MULTIVARIATE cross-distance matrix calculate for different tracks
  distmat <- proxy::dist(list_multivariate, method = "dtw")
  matrix_distmat<-as.matrix(distmat)
  ## Store TrackID names
  TrackID<-as.character(names(list_multivariate))
  
  library(umap)
  ## Project cross-distance matrix in a UMAP
  umap_dist<- umap(matrix_distmat,n_components=2,input="dist",init = "random", 
                   n_neighbors=pars$umap_n_neighbors, min_dist=pars$umap_minimal_distance, spread=1)  ### adjust parameters
  #Visualize plot
  pdf(file=paste0(output_dir,"Umap_unclustered.pdf"))
  
  plot(umap_dist$`layout`, 
       col=rgb(0,0,0,alpha=0.1), 
       pch=19,
       asp=0.4)
  
  dev.off()
  
  umap_1 <- as.data.frame(umap_dist$`layout`) 
  
  Track2_umap<-cbind(TrackID,umap_1)
  temp_df<-master_corrected3[,c("TrackID", "well","exp_nr","organoid_line")]
  temp_df<- temp_df[!duplicated(temp_df$TrackID),]
  umap_2 <- left_join(Track2_umap ,temp_df)
  
  ## Perform clustering. Select clusterig type that suits more your dataset
  #library(mclust)
  library(ggplot2)
  #mc.norm <- Mclust(umap_dist$`layout`, 9, start =123) ## slower but better clustering
  #umap_3 <- cbind(mc.norm$classification, umap_2)
  
  # Try dbscan
  library(dbscan)
  cl <- dbscan(umap_dist$`layout`, eps= 0.8, minPts = 10)  ## minPts is the min amount of values per cluster
  umap_3 <- cbind(cl$cluster, umap_2)
  
  # Try kmeans
  #km.norm <- kmeans(umap_dist$`layout`,10, nstart = 100)
  #umap_3 <- cbind(km.norm$cluster, umap_2)
  #outliers<- rownames(km.mod[["L"]]) ##find outliers
  #umap_3 <-umap_3%>%filter(!TrackID %in% outliers) ##remove outliers
  colnames(umap_3)[1]<- "cluster2"
  ## remove the outlier cluster if necessary
  umap_3<-subset(umap_3, cluster2!=0)
  
  ##plot
  # umap_3$cluster2 = factor(umap_3$cluster2, levels=c("4","7","5","8","3", "2","9","6","1"))
  # levels(umap_3$cluster2) <- c("1","2","3","4","5","6", "7","8","9")  ##reorder clusters
  ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(cluster2))) +  
    geom_point(size=2, alpha=0.6) + labs(color="cluster")+
    xlab("") + ylab("")+
    
    # +scale_color_manual(values = c("gold3",
    #                                                      "darkolivegreen3",
    #                                                      "seagreen3",
    #                                                      "blue3",
    #                                                      "dodgerblue",
    #                                                      "cyan1",
    #                                                      "indianred",
    #                                                      "firebrick",
    #                                                      "brown1" ))+
    ggtitle("umap Cluster ") +
    theme_light(base_size=20) +theme_bw()+
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(), aspect.ratio=1)+coord_fixed()
  
  ggsave(
    paste0(output_dir,"Umap_clustered.pdf"), 
    device="pdf", height=210, width=297, units="mm"
  )
  
  
  #### To the original dataset add information on cluster type
  master_clustered <- merge(master_processed ,umap_3[c("TrackID","cluster2")], by.x = "TrackID", by.y = "TrackID")
  # master_clustered$cluster = master_clustered$cluster2
  ### Save Reference UMAP for training of random forest classifier
  # saveRDS(master_clustered, file = "New_Behavioral_Referance_map") ### store here your reference map that can be used to predict behaviors in new experiments
  
  ## Plot a heatmap to show the relative values of each behavior parameter
  ## Create a dataframe the summarizes the mean values for each parameter
  sum_all <- master_clustered%>% select( speed, displacement, red_lym, contact2, contact_lym, cluster2, contact)%>% group_by(cluster2)%>%
    summarise(contact_len=mean(contact2),n_contact_org= mean(contact),displacement2 = median(displacement), speed = median(speed), interaction_T_cells= mean(contact_lym),death = median(red_lym))
  ## Rescale the values from each parameter
  sum_all <- sum_all%>%mutate(contact_len= rescale(contact_len, to=c(0,100)) ,n_contact_org= rescale(n_contact_org, to=c(0,100)),displacement2 = rescale(displacement2, to=c(0,100)), speed = rescale(speed, to=c(0,100)),interaction_T_cells= rescale(interaction_T_cells, to=c(0,100)), death =rescale(death, to=c(0,100)))
  
  ## Reshape dataframe
  library(reshape2)
  library(viridis)
  sum_all<-melt(sum_all,id.vars = "cluster2")
  sum_all$cluster2<-as.factor(sum_all$cluster2)
  ## Plot heatmap
  gg <- ggplot(data = sum_all, aes(x = variable, cluster2, fill = value))
  gg <- gg + geom_tile()
  gg <- gg + scale_fill_viridis(option="C", name="AU")
  gg <- gg + labs(x=NULL, y="Cluster", title="cluster represention")+theme(aspect.ratio=1.7,axis.text.x = element_text(angle = 45, hjust = 1))+ylim(rev(levels(sum_all$cluster2)))
  gg
  
  ggsave(
    paste0(output_dir,"Cluster_heatmap.pdf"), 
    device="pdf", height=210, width=297, units="mm"
  )
  
  master_clustered_info<-master_clustered[!duplicated(master_clustered$TrackID),c("TrackID","organoid_line","tcell_line","exp_nr", "well", "cluster2")] 
  
  master_clustered2 =master_clustered%>%group_by(TrackID)%>% arrange(Time)%>%
                      summarise(mean_speed= mean(q.speed),median_speed= median(q.speed),speed_sd= sd(q.speed),q3_speed= quantile(q.speed,0.90),
                                mean_displacement = mean(q.disp),median_displacement = median(q.disp),
                                displacement_sd=sd(q.disp),q3_disp= quantile(q.disp,0.90),
                                mean_red_lym = mean(q.red),red_lym_sd=sd(q.red),q3_red= quantile(q.red,0.90),
                                contact=mean(s.contact),mean_contact2=mean(contact2),contact2=max(contact2))
  
  master_clustered2<-left_join(master_clustered2,master_clustered_info)
  
  
  ### Quantify the number of cells per well
  Number_cell_exp<-master_clustered2%>%group_by(well, exp_nr, tcell_line, organoid_line)%>%
    summarise(total_cell = n())
  Percentage_clus<-left_join(Number_cell_exp, master_clustered2)
  Percentage_clus <- Percentage_clus%>%group_by(cluster2,tcell_line, well,exp_nr, organoid_line)%>%
    summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()
  
  Per<-ggplot(Percentage_clus, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
    geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
  Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ tcell_line)
  # Per <- Per + facet_grid(exp_nr + well + organoid_line  ~ tcell_line)
  
  # Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ interaction(organoid_line,tcell_line))
  Per<-Per+theme_void() +
    theme(aspect.ratio = 0.2,strip.text.x = element_text(angle = 90))
  
  Per
  
  ggsave(
    paste0(output_dir,"umap_cluster_percentage_bars_separate.pdf"), 
    device="pdf", height=210, width=297, units="mm"
  )
  
  ############### Combine all wells
  Number_cell_exp_combined<-master_clustered2%>%group_by(exp_nr, tcell_line, organoid_line)%>%
    summarise(total_cell = n())
  Percentage_clus_combined<-left_join(Number_cell_exp_combined,master_clustered2)
  Percentage_clus_combined <- Percentage_clus_combined%>%group_by(cluster2,exp_nr, tcell_line, organoid_line)%>%
    summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()
  
  Per<-ggplot(Percentage_clus_combined, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
    geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
  Per <- Per + facet_grid(interaction(exp_nr,organoid_line)  ~ tcell_line)
  # Per <- Per + facet_grid(exp_nr + well + organoid_line  ~ tcell_line)
  
  # Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ interaction(organoid_line,tcell_line))
  Per<-Per+theme_void() +
    theme(aspect.ratio = 0.2,strip.text.x = element_text(angle = 90))
  
  Per
  
  ggsave(
    paste0(output_dir,"umap_cluster_percentage_bars_combined.pdf"), 
    device="pdf", height=210, width=297, units="mm"
  )
}
