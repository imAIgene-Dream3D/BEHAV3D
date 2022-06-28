## Import T cell statistics for CD4 and CD8 cells only

### To Simulate the distribution of the behaviors in silico in our experiment (Fig4 BioXve Dekkers &Alieva, 2021) where we enrich for different behaviors by separating and washing organoids
### We need to reprocess the tracks, this time selecting tracks of a longuer length (the time that the experiment lasts). We will only select the tracks that fall in the range of the experiment that will
### then be processed and classified with the script CD4_CD8_4_6_hours_beh_predict

library(dplyr)
library(stats)
library(tidyr)
library(scales)
library(reshape2)
library(zoo)
library(spatstat)
library(sp)
library(stats)

getwd()

### Import processed tracks from the Reference map tracks
master<-readRDS("master_reference_map1")

colnames(master)[c(14)]<-"contact"

#### cluster incidence per timepoint
master$cell_type<-master$cell
master$cell_type <- gsub(".*CD8.*", "CD8", master$cell_type)
master$cell_type <- gsub(".*CD4.*", "CD4", master$cell_type)

## subset CD4 CD8 cells:
master<-subset(master, cell_type%in%c("CD8","CD4"))
## For each well separately:
## Estimate if a T cells interacts with another T cell interact. This is done by calculating the minmal distance between to the nearest neighbor. 
List = list() ## create list for each timepoint
List2 = list() ## create list for each experiment
## For loop to look for the distance to the nearest neighbor at each timepoint and each experiment
for ( m in unique(master$cell)){
  distance_1<-master[which(master$cell==m),]
  distanceDF<-as.data.frame(distance_1[,c(2,8,9,10)])
  List = list()
  for (i in unique(distanceDF$Time)){
    distanceDFi <- distanceDF[distanceDF$Time==i,]
    coordin<-ppx(distanceDFi, coord.type=c("t","s","s", "s")) ## create a multidimensional space-time pattern
    dist<- nndist(coordin)
    List[[length(List)+1]] <-dist   ## store to list object
  }
  nearest_n <- data.frame(matrix(unlist(List))) ## convert List to dataframe
  master_distance <- data.frame(distance_1,nearest_n)
  List2[[length(List2)+1]] <-master_distance 
}

master_dist<-do.call(rbind, List2)
colnames(master_dist)[c(16)] <- c("nearest_Tcell")
## create a binary variable that defines if a cell is interacting with another cells nased on its distance
master_dist$contact_lym<- ifelse(master_dist$nearest_Tcell<10,1,0)


## fill the NA missing values (for each variable)
##for speed
time_series<-acast(master_dist, Time ~ TrackID, value.var='speed',fun.aggregate = mean)
## get rid of NA
library(zoo)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series)
time_series2_speed<-time_series2[complete.cases(time_series2), ] 
colnames(time_series2_speed)<-c("Time", "TrackID", "speed")

##for displacement
time_series<-acast(master_dist, Time ~ TrackID, value.var='displacement',fun.aggregate = mean)
## get rid of NA
library(zoo)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series)
time_series2_disp<-time_series2[complete.cases(time_series2), ] 
colnames(time_series2_disp)<-c("Time", "TrackID", "displacement")

##for distorg
time_series<-acast(master_dist, Time ~ TrackID, value.var='dist_org',fun.aggregate = mean)
## get rid of NA
library(zoo)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series)
time_series2_dist_org<-time_series2[complete.cases(time_series2), ] 
colnames(time_series2_dist_org)<-c("Time", "TrackID", "dist_org")

##for red_lym
time_series<-acast(master_dist, Time ~ TrackID, value.var='red_lym',fun.aggregate = mean)
## get rid of NA
library(zoo)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series)
time_series2_red_lym<-time_series2[complete.cases(time_series2), ] 
colnames(time_series2_red_lym)<-c("Time", "TrackID", "red_lym")

##for contact
time_series<-acast(master_dist, Time ~ TrackID, value.var='contact',fun.aggregate = mean)
## get rid of NA
library(zoo)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series)
time_series2_contact<-time_series2[complete.cases(time_series2), ] 
colnames(time_series2_contact)<-c("Time", "TrackID", "contact")

##for contact lym
time_series<-acast(master_dist, Time ~ TrackID, value.var='contact_lym',fun.aggregate = mean)
## get rid of NA
library(zoo)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.approx(time_series_zoo) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2<-melt(time_series)
time_series2_contact_lym<-time_series2[complete.cases(time_series2), ] 
colnames(time_series2_contact_lym)<-c("Time", "TrackID", "contact_lym")

library(dplyr)
##for meanspeed
time_series2_meanspeed <-time_series2_speed %>% 
  group_by(TrackID) %>%mutate(meanspeed=rollapply(speed,10,mean,align='right',fill=NA))
## refill all missing values with the last value
time_series<-acast(time_series2_meanspeed, Time ~ TrackID, value.var='meanspeed',fun.aggregate = mean)
time_series_zoo<-zoo(time_series)
time_series_zoo<-na.locf(time_series_zoo, fromLast=T) ## replace by last value
time_series<-as.matrix(time_series_zoo)
time_series2_meanspeed<-melt(time_series)
colnames(time_series2_meanspeed)<-c("Time", "TrackID", "meanspeed")

### rebind all the files with NAs 
master_clust1 <- cbind(time_series2_speed,time_series2_disp,time_series2_dist_org,time_series2_red_lym, time_series2_contact, time_series2_contact_lym)
master_clust1 <- master_clust1[,!duplicated(colnames(master_clust1))]
master_cell<- master[c(3,11)]
master_cell<-master_cell[!duplicated(master_cell$TrackID),]
master_exp<- master[c(3,12)]
master_exp<-master_exp[!duplicated(master_exp$TrackID),]
master_exp2<- master[c(3,13)]
master_exp2<-master_exp2[!duplicated(master_exp2$TrackID),]
master_cell_type<- master[c(3,15)]
master_cell_type<-master_cell_type[!duplicated(master_cell_type$TrackID),]
library(dplyr)
master_clust1<- left_join(master_clust1 ,master_cell)
master_clust1<- left_join(master_clust1 ,master_exp)
master_clust1<- left_join(master_clust1 ,master_exp2)
master_clust1<- left_join(master_clust1 ,master_cell_type)
#add the meanspeed values
master_clust1<- merge(master_clust1, time_series2_meanspeed, by = c("Time","TrackID"))
## count the time of contact:

master_clust1$contact<-ifelse(master_clust1$dist_org>0.0015,0,1)


##### select the cells that are interacting at 5 hours (6 hours of incubation) 150timepoints
# select tracks with more then 2 hours
master_clust1<-master_clust1 %>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(Time<155)%>%filter(n() > 59)

### take the last 100 timepoint sor each
master_clust<-master_clust1 %>% 
  group_by(TrackID)%>%arrange(Time)%>%slice_tail(n = 60)


master_clust<-master_clust %>% 
  group_by(TrackID) %>%arrange(Time)%>%mutate(Time2 = Time - first(Time))

master_clust<-master_clust %>% 
  group_by(TrackID) %>%arrange(TrackID)%>% filter(Time2<100)
## remove the dead t cells from the beginning
master_clust3deadT0 <-master_clust%>%group_by(TrackID)%>%filter((Time2==0) & red_lym<10 )
master_clust3 <-master_clust%>%filter(TrackID %in% master_clust3deadT0$TrackID )
##make death binary:
master_clust3$death<- ifelse(master_clust3$red_lym<10,0,1)

## Create a variable for cumulative interaction with organoids
master_corrected3<-master_clust3 %>% 
  group_by(TrackID) %>%mutate(contact2=(ave(contact, cumsum(!contact), FUN = cumsum)))
## Create a variable for T cells interact with other T cells while in the environment
master_corrected3$contact_lym<- ifelse(master_corrected3$contact==1,0,master_corrected3$contact_lym)
## For T cells inteacting in the environment keep as "interacting" only cells that had a mean speed in the last 20 minutes that is in the upper quantile.
master_corrected3<-master_corrected3%>%group_by(exp)%>%mutate(contact_lym=ifelse(meanspeed<quantile(meanspeed,p=0.75),0,contact_lym))


##check how many cells do we have at timepoint 1:
View(subset(master_corrected3, Time%in%c(150:150)))

saveRDS(master_corrected3,"master_CD3_scRNA_seq_inference_4")
###################################Predict behavior up to 4 hours #############################################################
