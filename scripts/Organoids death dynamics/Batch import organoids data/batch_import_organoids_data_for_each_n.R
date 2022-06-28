library(plyr)
library(readr)
library(dplyr)
library(yaml)

read_plus <- function(flnm) {
  read_csv(flnm, skip=3) %>% 
    mutate(filename = flnm)
}

pars = yaml.load_file("/Users/samdeblank/OneDrive - Prinses Maxima Centrum/github/BEHAV3D/configs/config_template.yml")
setwd(pars$output_dir)  ## set here your working directory that contains the example dataset
working_directory <- pars$data_dir

# import volumes
pat = "*Volume"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
volume_csv <- ldply(files, read_plus)
# import sum_red
pat = "*Intensity_Mean_Ch=3_Img=1"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
sum_red_csv <- ldply(files, read_plus)
# import area
pat = "*Area"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
area_csv <- ldply(files, read_plus)
# import position
pat = "*Position"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
pos_csv <- ldply(files, read_plus)


live_deadROI <- cbind(volume_csv[,c(1,4,5,6)], 
                      sum_red_csv[,c(1)], 
                      area_csv[,c(1)],
                      pos_csv[,c(1,2,3,11)])
colnames(live_deadROI ) <- c("Volume","Time","TrackID","ID","red_sum","area", "pos_x","pos_y","pos_z", "Org")

## make TrackID unique:
category <- as.factor(live_deadROI$Org)
ranks <- rank(-table(category), ties.method="first")
ranks <- as.data.frame(ranks)
ranks$Org <- row.names(ranks)
live_deadROI <- left_join(live_deadROI, ranks)  ## plot with all the tracks together 

live_deadROI$Track2 <- with(live_deadROI, interaction(TrackID, ranks))
live_deadROI$Track2 <- gsub(".", '', live_deadROI$Track2, fixed = T)
live_deadROI$Track2 <- as.numeric(as.character(live_deadROI$Track2))


detach(package:plyr)
library(dplyr)
live_deadROI$red_sum <- live_deadROI$red_sum*live_deadROI$Volume
live_deadROI1 <-live_deadROI %>% 
  select(Track2, Time, Volume, red_sum, area, pos_x, pos_y, pos_z, Org) %>% 
  group_by(Track2, Time) %>% 
  summarise(Volume = sum(Volume), red= sum(red_sum), area=sum(area), pos_x=mean(pos_x), pos_y=mean(pos_y), pos_z=mean(pos_z))

live_deadROI2<- live_deadROI %>% 
  group_by(Track2, Time) %>% 
  distinct(Org)

live_deadROI3 <- left_join(live_deadROI1, live_deadROI2)  ## plot with all the tracks together 
live_deadROI3$sumred <- live_deadROI3$red
live_deadROI3$red <- live_deadROI3$red/live_deadROI3$Volume
live_deadROI3 <- live_deadROI3[complete.cases(live_deadROI3), ] 


## change name of the cell type
live_deadROI3$Org <- gsub("/", "", live_deadROI3$Org)
live_deadROI3$Org <- gsub("\\(", "", live_deadROI3$Org)
live_deadROI3$Org <- gsub("\\)", "", live_deadROI3$Org)
live_deadROI3$Org <- gsub("\\[", "", live_deadROI3$Org)
live_deadROI3$Org <- gsub("\\]", "", live_deadROI3$Org)

live_deadROI3$exp <- live_deadROI3$Org  ## for rthe experiemtnal batch
live_deadROI3$type <- live_deadROI3$exp  ## for the Orgtype LM1 or TEG

names<-live_deadROI3[!duplicated(live_deadROI3$Org),]
names$Org
names
all_donors <- c("20t", "20T","10t","10T", "36t", "36T", "100t", "100T", "38t", "38T", "13t", "13T", "14t", "14T", "62t", "62T", "25t", "25T", "169m", "169M", "209m", "209M", "27t", "27T", "34t", "34T", "1837m", "1837M", "100tnew", "100Tnew")
#
tcell_names <- list("001","LM1")
### loop to identify tumor donor, TEG type and exp date
library(stringr)
length_Names= length(names$Org)
length_Donors= length(all_donors)
for (i in 1:length_Names){
  #this part will search for n= and the number 
  #that comes afterwards and replaces exp type bz exp and the number
  strings <- c(names$Org[i])
  if (stringr::str_detect(strings, '(n\\=)')){
    exp_number <- stringr::str_match(strings, '(n\\=)(\\d)')[,3]
    exp_name <- paste0("exp", exp_number)
    live_deadROI3$exp <- gsub(strings,exp_name,live_deadROI3$exp)}
  # else {
  #exp_number <- stringr::str_match(strings, '(_n)(\\d)')[,3]
  #exp_name <- paste0("exp", exp_number)
  #live_deadROI3$exp <- gsub(strings,exp_name,live_deadROI3$exp)}
  else if (stringr::str_detect(strings, '(_n)')){
    exp_number <- stringr::str_match(strings, '(_n)(\\d)')[,3]
    exp_name <- paste0("exp", exp_number)
    live_deadROI3$exp <- gsub(strings,exp_name,live_deadROI3$exp)}
  else{live_deadROI3$exp <- gsub(strings,"weird naming",live_deadROI3$exp)}
  
  # This part searches for 001 and LM1 in the file name and replaces it in the table
  if (stringr::str_detect(strings, "001")){
    live_deadROI3$type <- gsub(strings,"TEG",live_deadROI3$type)}
  else if (stringr::str_detect(strings, "LM1|lm1")){
    live_deadROI3$type <- gsub(strings,"LM1",live_deadROI3$type)}
  else {live_deadROI3$type <- gsub(strings,"Non",live_deadROI3$type)}
  
  #this part searches for each item of the list "donor_names" within all files and replaces Org by the donor name
  for (j in 1:length_Donors){
    temp = all_donors[j]
    if(stringr::str_detect(strings, temp)){
      temp_upper = str_to_upper(temp)
      live_deadROI3$Org <- gsub(strings,temp_upper,live_deadROI3$Org)
    }}
}
names$Org[i] 
exp_number
exp_name
temp
## set time in hours
live_deadROI3$Time2<-(live_deadROI3$Time-1)/2
library(ggplot2)
### Quantify the dead cell dye intensity per well
live_deadROI32<-live_deadROI3%>% 
  group_by(Track2) %>%arrange(Track2)%>% filter(n() > 40) ## filter cells with at least 40 timepoints
### pulled all the dead cell dye signal from the well
live_deadROI7 <-live_deadROI32 %>% 
  select(Track2, Time2, Volume, sumred, Org, type) %>% 
  group_by(Org, Time2, type) %>% 
  summarise(Volume = sum(Volume), red= sum(sumred))
live_deadROI7$red<-live_deadROI7$red/live_deadROI7$Volume

### set an experiment ID if processing several experimental replicates
live_deadROI7$exp<-"1"
### SAVE dataframe for processing in a different script
saveRDS(live_deadROI7, file = "Full_well_death_dynamics")   ### save here a dataframe with all the well values



### Process the death dynamics per individual organoid
library(scales)
live_deadROI3<-live_deadROI3%>% ungroup()%>%
  mutate(red=rescale(red, to=c(0,100)))
live_deadROI3<-live_deadROI3%>% 
  group_by(Track2) %>%arrange(Track2)%>% filter(n() > 40)
temp1 <- aggregate(red ~ Track2+Org+type, data = live_deadROI3, max)  ##calculate the max red of each track
colnames(temp1) [4] <- "max_red"
live_deadROI4 <- merge(temp1, live_deadROI3)
live_deadROI6 <- subset(live_deadROI4 , Volume>1000)
live_deadROI6deadT0 <-live_deadROI6%>%group_by(Track2)%>%filter((Time==1) & red<10 )
live_deadROI6 <-live_deadROI6%>%filter(Track2 %in% live_deadROI6deadT0$Track2)


### plot to check outcome
library(ggplot2)
Plot <- ggplot(live_deadROI6, aes(Time2,red, color = Track2, group = Track2)) + 
  geom_smooth(method="loess", size = 1, se=F, span=1) +
  theme_bw() + 
  ylab("dead dye intensity") + 
  xlab("Time (hours)") +
  theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 20), legend.text=element_text(size= 10))+
  labs(color = "Organoid")
ggtitle("Individual org increase in dead dye intensity TEG")

Plot
### SAVE Dataframe for later processing
saveRDS(live_deadROI6, file = "Individual_orgs_death_dynamics")  ### save here a dataframe with all the organoids values
