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

library(dplyr)
## function to import organoid data from csv files extracted by Imaris
read_plus <- function(flnm) {
  read_csv(flnm, skip = 3) %>% 
    mutate(filename = flnm)
}
## directory where the files are located
working_directory <- "D:/R/scripts/T_cell paper/FINAL SCRIPTS_20210408/Fig2/b/github/example_t cell_data"
setwd(working_directory)
# import Displacement^2
pat = "*Displacement"
files <- list.files(path = working_directory, pattern = pat,  recursive = TRUE)
displacement <- ldply(files, read_plus)
# import Speed
pat = "*Speed"
files <- list.files(path = working_directory, pattern = pat,  recursive = TRUE)
speed <- ldply(files, read_plus)
# import mean dead dye intensity values
pat = "*Intensity_Mean_Ch=3_Img=1"
files <- list.files(path = working_directory, pattern = pat, recursive = TRUE)
red_lym <- ldply(files, read_plus)
# import Minimal distance to organoids
pat = "*dist_org"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
dist_org <- ldply(files, read_plus)
# import Position
pat = "*Position"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
pos <- ldply(files, read_plus)

### join all
master <- cbind(displacement[,c(1,4,5,6)], speed[,c(1)], dist_org[,c(1)], red_lym[,c(1)], pos[,c(1,2,3,11)])

### remove unnecesary signs from filename:
master$filename<- gsub("/", "", master$filename)
master$filename <- gsub("\\(", "", master$filename)
master$filename <- gsub("\\)", "", master$filename)
## rename well // date // cell type according to filename
master$well<-master$filename  ### well name
master$well <- gsub(".*4_13T_CD4teg_CD8TEG.*", "well4_13T_exp1", master$well, perl=TRUE)
master$well <- gsub(".*20201027_610T.*", "well6_10T_exp2", master$well, perl=TRUE)

master$exp<-master$filename  ### well name
master$exp <- gsub(".*2020-07-28.*", "20200728", master$exp, perl=TRUE)
master$exp <- gsub(".*20201027.*", "20201027", master$exp, perl=TRUE)


master$cell_type<-master$filename  ### well name
master$cell_type <- gsub(".*610T_g_.*", "CD4_TEG_10T", master$cell_type, perl=TRUE)
master$cell_type <- gsub(".*610T_b_.*", "CD8_TEG_10T", master$cell_type, perl=TRUE)
master$cell_type <- gsub(".*13T_CD4teg_CD8TEG_green.*", "CD4_TEG_13T", master$cell_type, perl=TRUE)
master$cell_type <- gsub(".*13T_CD4teg_CD8TEG_blue.*", "CD8_TEG_13T", master$cell_type, perl=TRUE)


### create a unique TRACKID. Each file processes with Imaris must have a unique track ID.
category <- as.factor(master$filename)
ranks <- rank(-table(category), ties.method="first")
ranks <- as.data.frame(ranks)
ranks$filename <- row.names(ranks)
master <- left_join(master, ranks) 
master$TrackID2 <- with(master, interaction(TrackID, ranks))
master$TrackID2 <- gsub(".", '', master$TrackID2, fixed = T)
master$TrackID2 <- as.numeric(as.character(master$TrackID2))

## Remove filename
master$filename<-NULL

## COLnames
colnames(master) <- c("displacement","Time","TrackID","ID","speed","dist_org","red_lym","X-pos","Y-pos","Z-pos","well","exp","cell_type","ranks","TrackID2")

### save RDS for later use (e.g. backprojection of classified TrackIDs)
saveRDS(master, "master_example_data")
detach("package:reshape2", unload=TRUE)
detach("package:plyr", unload=TRUE)
master <- master[which(master$Time<=300), ] ##Make sure that all the time-series have the same length, inthis case 10hours

## check for duplicates: gives 0 duplicates
data_dup <- master%>%group_by(Time)%>%
  count(TrackID2) %>% 
  filter(n > 1) %>% 
  select(-n)
data_dup ###must be empty


## For each well separately:
## Estimate if a T cells interacts with another T cell interact. This is done by calculating the minmal distance between to the nearest neighbor. 
List = list() ## create list for each timepoint
List2 = list() ## create list for each experiment
## For loop to look for the distance to the nearest neighbor at each timepoint and each experiment
for ( m in unique(master$well)){
  distance_1<-master[which(master$well==m),]
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
## create a binary variable that defines if a well is interacting with another wells nased on its distance
master_dist$contact_lym<- ifelse(master_dist$nearest_Tcell<10,1,0)

###remove the varible TrackID and use only TrackID2 (unique identifier instead)
master_dist$TrackID<-master_dist$TrackID2
master_dist$TrackID2<-NULL

library(reshape2)
## Since not all the tracks are tracked at all timepoints interpolate missing values and fill the NA (for each variable)
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

### Create a dataframe with all the variables with corrected missing values
master_corrected <- cbind(time_series2_speed,time_series2_disp,time_series2_dist_org,time_series2_red_lym, time_series2_contact_lym)
## remove duplicated column names
master_corrected <- master_corrected[,!duplicated(colnames(master_corrected))]
## join the information on the cell type/ experiment number and group number
master_temp<- master[c("TrackID2", "well","exp", "cell_type")]
master_temp<-master_temp[!duplicated(master_temp$TrackID2),]
master_corrected<- left_join(master_corrected ,master_temp, by=c("TrackID"="TrackID2"))
#Merge the information for the mean speed over the last 20 mins
master_corrected1<- merge(master_corrected, time_series2_meanspeed, by = c("Time","TrackID"))

## Update the binary variable for contact with organoids (it can vary between experiments depending on the intensity of the T cells or organoids. Check the threshold of contact in the imaging data and update):
master_corrected1_1<-subset(master_corrected1, exp%in%c("20200728"))
master_corrected1_2<-subset(master_corrected1, exp%in%c("20201027"))

master_corrected1_1$contact <- ifelse(master_corrected1_1$dist_org>0.0012, 0,1)
master_corrected1_2$contact <- ifelse(master_corrected1_2$dist_org>0.0001, 0,1)


master_corrected1<-rbind(master_corrected1_1,master_corrected1_2)

library(ggplot2)
ggplot(master_corrected1, aes(x=contact, color=exp)) +
  geom_histogram(fill="white", alpha=0.5, position="identity")+facet_grid(.~well, scales = "free")

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
master_corrected3<-master_corrected3%>%group_by(exp)%>%mutate(contact_lym=ifelse(meanspeed<quantile(meanspeed,p=0.75),0,contact_lym))


saveRDS(master_corrected3, file = "master_corrected3_example")
