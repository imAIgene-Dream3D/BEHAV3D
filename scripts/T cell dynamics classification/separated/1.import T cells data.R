set.seed(1234)
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

pars = yaml.load_file("/Users/samdeblank/surfdrive/Shared/T cell paper/Stats reports/t cells/2021-08-11_ror1_CART_n3/BEHAV3D_config.yml")

## directory where the files are located
working_directory <- pars$data_dir
setwd(working_directory)
output_dir=paste0(pars$output_dir,"/")

## function to import organoid data from csv files extracted by Imaris
read_imaris_csv <- function(flnm) {
  read_csv(flnm, skip = 3, col_types = cols()) %>% 
    mutate(filename = flnm)
}

# Import file-specific metadata for all images used in this analysis.
pat = pars$metadata_csv
metadata=read.csv(pars$metadata_csv, sep="\t")

read_ims_csv <- function(pattern, recursive=TRUE) {
  read_plus <- function(flnm) {
    read_csv(flnm, skip = 3, col_types = cols()) %>% 
      mutate(filename = flnm)
  }
  files <- list.files(path = working_directory, pattern = pattern,  recursive = recursive)
  ims_csv <- ldply(files, read_plus)
  return(ims_csv)
}

# import Displacement^2
pat = "*Displacement"
displacement <- read_ims_csv(pattern=pat)

# import Speed
pat = "*Speed"
speed <- read_ims_csv(pattern=pat)

# import mean dead dye intensity values
pat = "*Intensity_Mean_Ch=3_Img=1"
red_lym <- read_ims_csv(pattern=pat)

# import Minimal distance to organoids
pat = "*dist_org"
dist_org <- read_ims_csv(pattern=pat)

# import Position
pat = "*Position"
pos <- read_ims_csv(pattern=pat)

### join all
master <- cbind(
          displacement[,c("Displacement^2","Time","TrackID" ,"ID")], 
          speed[,c("Speed" )], 
          dist_org[,c("Intensity Min")], 
          red_lym[,c("Intensity Mean")], 
          pos[,c("Position X" ,"Position Y" ,"Position Z","filename")]
          )

### Convert the filename to the same format in both datasets (master and metadata)
master$basename <- gsub("_Position.csv", "", master$filename, perl=TRUE)
master$basename=basename(master$basename)

colnames(master) <- c("displacement","Time","TrackID","ID","speed","dist_org","red_lym","X-pos","Y-pos","Z-pos", "filename", "basename")

##Join the information of metadata to master:
master<-left_join(master, metadata)

### create a unique TRACKID. Each file processes with Imaris must have a unique track ID.
category <- as.factor(master$filename)
ranks <- rank(-table(category), ties.method="first")
ranks <- as.data.frame(ranks)
ranks$filename <- row.names(ranks)
master <- left_join(master, ranks) 
# master$TrackID2 <- with(master, interaction(TrackID, ranks))
master$TrackID2 <- with(master, interaction(ranks, TrackID, sep="_"))

# master$TrackID2 <- gsub(".", '', master$TrackID2, fixed = T)
# master$TrackID2 <- as.numeric(as.character(master$TrackID2))

## Remove filename
master$filename<-NULL

## COLnames

### save RDS for later use (e.g. backprojection of classified TrackIDs)
saveRDS(master, paste0(output_dir,"raw_tcell_track_data.rds"))
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

### VARIABLE
master <- master[which(master$Time<=300), ] ##Make sure that all the time-series have the same length, in this case 10hours

## check for duplicates: gives 0 duplicates
data_dup <- master%>%group_by(Time)%>%
  count(TrackID2) %>% 
  filter(n > 1) %>% 
  select(-n)

if (dim(data_dup)[1]!=0){
  stop("There are duplicates in the data, which should not be the case. Stopping execution...")
}

## For each well separately:
## Estimate if a T cells interacts with another T cell interact. This is done by calculating the minmal distance between to the nearest neighbor. 

List = list() ## create list for each timepoint
List2 = list() ## create list for each experiment
## For loop to look for the distance to the nearest neighbor at each timepoint and each experiment
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

## create a binary variable that defines if a well is interacting with another wells nased on its distance
# VARIABLE
master_dist$contact_lym<- ifelse(master_dist$nearest_Tcell<10,1,0)

###remove the varible TrackID and use only TrackID2 (unique identifier instead)
master_dist$TrackID<-master_dist$TrackID2
master_dist$TrackID2<-NULL

library(reshape2)
library(zoo)
## Since not all the tracks are tracked at all timepoints interpolate missing values and fill the NA (for each variable)
## select the variables for which we need to interpolate NAs (numeric)
column_names<-names(master_dist)
column_names <- c("displacement", "dist_org", "red_lym", "contact_lym")

## create a first dataset with refilled values for speed:
time_series<-acast(master_dist, Time ~ TrackID, value.var='speed',fun.aggregate = mean)
## rownames timepoints:
row.names(time_series)<-unique(master_dist$Time)
## get rid of NA
time_series_zoo<-zoo(time_series, row.names(time_series))
time_series_zoo<-na.approx(time_series_zoo) ## replace by interpolated value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series) # ERROR HERE
data<-time_series2[complete.cases(time_series2), ] 
colnames(data)<-c("Time", "TrackID", "speed")
## store this data for calculating lagged speed later:
time_series2_speed<-data

### ----------
for (i in column_names){
  time_series<-acast(master_dist, Time ~ TrackID, value.var=i,fun.aggregate = mean)
  row.names(time_series)<-unique(master_dist$Time)
  ## get rid of NA
  time_series_zoo<-zoo(time_series,row.names(time_series))
  time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
  time_series<-as.matrix(time_series_zoo)
  time_series2<-melt(time_series)
  new<-time_series2[complete.cases(time_series2), ] 
  data[ , ncol(data) + 1] <- new[3]                  # Append new column
  colnames(data)[ncol(data)] <- paste0(i)
}

library(dplyr)
##For cell interaction we need to consider the following:
## When two cells interact it is often the one cell moves and interacts with another one that is static
## In this case one might consider that only the motile cell is actively interacting and the static cell is just passively interacting
## To determine when a cell is actively interacting we measure for each cell what was its mean speed over the last 10 timepoints (20 mins)
time_series2_meanspeed <-time_series2_speed %>% 
  group_by(TrackID) %>%mutate(meanspeed=rollapply(speed,10,mean,align='right',fill=NA))
## refill all missing values with the last value
time_series<-acast(time_series2_meanspeed, Time ~ TrackID, value.var='meanspeed',fun.aggregate = mean)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.locf(time_series_zoo, fromLast=T) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2_meanspeed<-melt(time_series)
colnames(time_series2_meanspeed)<-c("Time", "TrackID", "meanspeed")

## remove last NAs
time_series2_meanspeed<-na.omit(time_series2_meanspeed)
### Create a dataframe with all the variables with corrected missing values
master_corrected <- data
## join the information on the cell type/ experiment number and group number
master_temp<- master[c("TrackID2", colnames(metadata)[-1])]
master_temp<-master_temp[!duplicated(master_temp$TrackID2),]
master_corrected<- left_join(master_corrected ,master_temp, by=c("TrackID"="TrackID2"))
#Merge the information for the mean speed over the last 20 mins
master_corrected1<- merge(master_corrected, time_series2_meanspeed, by = c("Time","TrackID"))

## Update the binary variable for contact with organoids (it can vary between experiments depending on the intensity of the T cells or organoids. Check the threshold of contact in the imaging data and update in the metadata csv):

master_corrected1$contact <- ifelse(master_corrected1$dist_org>master_corrected1$contact_threshold, 0,1)

library(ggplot2)
ggplot(master_corrected1, aes(x=contact, color=as.factor(exp_nr))) +
  geom_histogram(fill="white", alpha=0.5, position="identity")+facet_grid(organoid_line~well, scales = "free")

## Remove contact threhold viarable:
master_corrected1$contact_threshold<-NULL

## For clustering it is necessary to compare T cell tracks that have a similar length. 
## For that we select cell track that have at least 100 timepoints. 
## Detach package  'plyr' as it can interfere with 'dplyr'
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)

master_corrected2<-master_corrected1 %>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(Time>00&Time<300)%>% filter(n() > 99)
## Create a variable for the relative Time
master_corrected2<-master_corrected2 %>% 
  group_by(TrackID) %>%arrange(Time)%>%mutate(Time2 = Time - first(Time))
## For the Tracks that have more then 100 timepoints filter only the first 100.
master_corrected2<-master_corrected2 %>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(Time2<100)
## To exclude noise due to dead cells remove the dead t cells from the beginning
master_corrected3deadT0 <-master_corrected2%>%group_by(TrackID)%>%filter((Time2==0) & red_lym<10 )
master_corrected3 <-master_corrected2%>%filter(TrackID %in% master_corrected3deadT0$TrackID )
## Create a binary variable for live or dead cells:
master_corrected3$death<- ifelse(master_corrected3$red_lym<10,0,1)
## Create a variable for cumulative interaction with organoids
master_corrected3<-master_corrected3 %>% 
  group_by(TrackID) %>%mutate(contact2=(ave(contact, cumsum(!contact), FUN = cumsum)))
## Create a variable for T cells interact with other T cells while in the environment
master_corrected3$contact_lym<- ifelse(master_corrected3$contact==1,0,master_corrected3$contact_lym)
## For T cells inteacting in the environment keep as "interacting" only cells that had a mean speed in the last 20 minutes that is in the upper quantile.
master_corrected3<-master_corrected3%>%group_by(exp_nr)%>%mutate(contact_lym=ifelse(meanspeed<quantile(meanspeed,p=0.75),0,contact_lym))


saveRDS(master_corrected3, file = paste0(output_dir,"processed_tcell_track_data.rds"))

