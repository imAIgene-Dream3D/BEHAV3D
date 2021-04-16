library(plyr)
library(readr)
library(dplyr)
## function to import organoid data from csv files extracted by Imaris
read_plus <- function(flnm) {
  read_csv(flnm, skip=3) %>% 
    mutate(filename = flnm)
}

## set directory where csv files are located
setwd('E:/Processed/Floppy extracts/test screen Mak/Organoid stats/test')
working_directory <- "E:/Processed/Floppy extracts/test screen Mak/Organoid stats/test"

# import volume per organoid object
pat = "*Volume"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
volume_csv <- ldply(files, read_plus)
# import mean dead dye intensity statistics per organoid object
pat = "*Intensity_Mean_Ch=3_Img=1"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
sum_red_csv <- ldply(files, read_plus)
# import area per organoid object
pat = "*Area"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
area_csv <- ldply(files, read_plus)
# import position per organoid object
pat = "*Position"
files <- list.files(path = working_directory, pattern = pat, full.names = T, recursive = TRUE)
pos_csv <- ldply(files, read_plus)

##create dataframe containing the statistics of interest
master <- cbind(volume_csv[,c(1,4,5,6)], 
                      sum_red_csv[,c(1)], 
                      area_csv[,c(1)],
                      pos_csv[,c(1,2,3,11)])
colnames(master ) <- c("Volume","Time","TrackID","ID","red_sum","area", "pos_x","pos_y","pos_z", "Org")

## make each TrackID unique (TrackIDs imported from different imaging files can be repeated):
category <- as.factor(master$Org)
ranks <- rank(-table(category), ties.method="first")
ranks <- as.data.frame(ranks)
ranks$Org <- row.names(ranks)
master <- left_join(master, ranks) 

master$Track2 <- with(master, interaction(TrackID, ranks))
master$Track2 <- gsub(".", '', master$Track2, fixed = T)
master$Track2 <- as.numeric(as.character(master$Track2))


detach(package:plyr) ## detach plyr as 
library(dplyr)
## "one organoid" can be composed by several detached particles, specially when it start to die.
## Unify the statistics for the full organoid per timepoint
master$red_sum <- master$red_sum*master$Volume
master1 <-master %>% 
  select(Track2, Time, Volume, red_sum, area, pos_x, pos_y, pos_z, Org) %>% 
  group_by(Track2, Time) %>% 
  summarise(Volume = sum(Volume), red= sum(red_sum), area=sum(area), pos_x=mean(pos_x), pos_y=mean(pos_y), pos_z=mean(pos_z))

master2<- master %>% 
  group_by(Track2, Time) %>% 
  distinct(Org)

master3 <- left_join(master1, master2)  
master3$sumred <- master3$red
master3$red <- master3$red/master3$Volume
master3 <- master3[complete.cases(master3), ] 


## read file name to extract information about the type of organoid, type of  T cell and experiment per well analysed
master3$Org <- gsub("/", "", master3$Org)
master3$Org <- gsub("\\(", "", master3$Org)
master3$Org <- gsub("\\)", "", master3$Org)

master3$exp <- master3$Org  ## create a variable for the experiment number
master3$type <- master3$exp  ## create a variable for the cell type
## make a vector with all the filenames from which information needs to be extracted
names<-master3[!duplicated(master3$Org),]
## make a vector with all the organoid names that have to be matched in the filename
all_donors <- c("20t", "20T","10t","10T", "36t", "36T", "100t", "100T", "38t", "38T", "13t", "13T", "14t", "14T", "62t", "62T", "25t", "25T", "169m", "169M", "209m", "209M", "27t", "27T", "34t", "34T", "1837m", "1837M", "100tnew", "100Tnew")
## make a vector with all the T cell types that have to be matched in the filename
tcell_names <- list("001","LM1")
## for loop to detect specific strings and substite exp, Org and type by the appropriate name
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
    master3$exp <- gsub(strings,exp_name,master3$exp)}
  # else {
  #exp_number <- stringr::str_match(strings, '(_n)(\\d)')[,3]
  #exp_name <- paste0("exp", exp_number)
  #master3$exp <- gsub(strings,exp_name,master3$exp)}
  else if (stringr::str_detect(strings, '(_n)')){
    exp_number <- stringr::str_match(strings, '(_n)(\\d)')[,3]
    exp_name <- paste0("exp", exp_number)
    master3$exp <- gsub(strings,exp_name,master3$exp)}
  else{master3$exp <- gsub(strings,"weird naming",master3$exp)}
  
  # This part searches for 001 and LM1 in the file name and replaces it in the table
  if (stringr::str_detect(strings, "001")){
    master3$type <- gsub(strings,"TEG",master3$type)}
  else if (stringr::str_detect(strings, "LM1|lm1")){
    master3$type <- gsub(strings,"LM1",master3$type)}
  else {master3$type <- gsub(strings,"Non",master3$type)}
  
  #this part searches for each item of the list "donor_names" within all files and replaces Org by the donor name
  for (j in 1:length_Donors){
    temp = all_donors[j]
    if(stringr::str_detect(strings, temp)){
      temp_upper = str_to_upper(temp)
      master3$Org <- gsub(strings,temp_upper,master3$Org)
    }}
}

### save here a dataframe with all the organoids statistics
saveRDS(master3, file = "Organoid_data") 
detach("package:plyr", unload=TRUE)