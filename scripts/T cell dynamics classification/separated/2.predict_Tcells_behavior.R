library(dplyr)
library(stats)
library(tidyr)
library(scales)
library(randomForest)
library(ggplot2)
library(yaml)
## Import reference map that will be used as the ground truth dataset to train and test a random forest classifier

pars = yaml.load_file("/Users/samdeblank/surfdrive/Shared/T cell paper/Stats reports/t cells/2021-08-11_ror1_CART_n3/BEHAV3D_config.yml")

reference_map <- pars$reference_map
output_dir <- pars$output_dir
model_path <- pars$randomforest

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
master_test<-readRDS(paste0(output_dir, "processed_tcell_track_data.rds"))
### Normalize the data
master_test2<-master_test%>% ungroup()%>%
  # group_by(exp_nr) %>% 
  mutate(z.disp = (displacement-mean(displacement))/sd(displacement),z.speed = (speed-mean(speed))/sd(speed), z.red = (red_lym-mean(red_lym))/sd(red_lym))%>%
  mutate(q.disp=ifelse(z.disp>(quantile(z.disp, p=0.75)),z.disp,min(z.disp)), q.speed=ifelse(z.speed>(quantile(z.speed, p=0.75)),z.speed,min(z.speed)),q.red=ifelse(z.red>(quantile(z.red, p=0.75)),z.red,min(z.red)))%>%
  mutate(q.disp=rescale(q.disp, to=c(0,100)),q.speed=rescale(q.speed, to=c(0,100)),q.red=rescale(q.red, to=c(0,100)),s.contact=rescale(contact, to=c(0,1)),s.contact_lym=rescale(contact_lym, to=c(0,1)))%>%
  mutate(q.disp=q.disp/mean(quantile(q.disp, p=0.9999)),q.speed=q.speed/mean(quantile(q.speed, p=0.9999)),q.red=q.red/mean(quantile(q.red, p=0.9999)))%>%ungroup()

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

saveRDS(Percentage_clus, file = paste0(output_dir,"cluster_perc_tcell_track_data.rds"))
### Plot proportion per well and per cell type
Percentage_clus$tcell_line = as.factor(Percentage_clus$tcell_line)

Per<-ggplot(Percentage_clus, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
Per <- Per + facet_grid(interaction(tcell_line, exp_nr,well,organoid_line)  ~ tcell_line)
# Per <- Per + facet_grid(exp_nr + well + organoid_line  ~ tcell_line)

# Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ interaction(organoid_line,tcell_line))
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

### Run a chisq.test to see if the distibution is different between conditions
library(reshape2)
table1<-acast(Percentage_clus, cluster2~tcell_line, value.var = "percentage")
table1[is.na(table1)]<-0  ## If any cluster are missing convert NA to 0

## Run Pearson Chi sq test to see if there is a different distribution of clusters between cell types
chisq.test(table1)
### Normalize the data
master_test2<-master_test%>% ungroup()%>%
  group_by(exp_nr) %>% 
  mutate(z.disp = (displacement-mean(displacement))/sd(displacement),z.speed = (speed-mean(speed))/sd(speed), z.red = (red_lym-mean(red_lym))/sd(red_lym))%>%
  mutate(q.disp=ifelse(z.disp>(quantile(z.disp, p=0.75)),z.disp,min(z.disp)), q.speed=ifelse(z.speed>(quantile(z.speed, p=0.75)),z.speed,min(z.speed)),q.red=ifelse(z.red>(quantile(z.red, p=0.75)),z.red,min(z.red)))%>%
  mutate(q.disp=rescale(q.disp, to=c(0,100)),q.speed=rescale(q.speed, to=c(0,100)),q.red=rescale(q.red, to=c(0,100)),s.contact=rescale(contact, to=c(0,1)),s.contact_lym=rescale(contact_lym, to=c(0,1)))%>%
  mutate(q.disp=q.disp/mean(quantile(q.disp, p=0.9999)),q.speed=q.speed/mean(quantile(q.speed, p=0.9999)),q.red=q.red/mean(quantile(q.red, p=0.9999)))%>%ungroup()

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

### Quantify the number of cells per well
Number_cell_exp<-classified_tracks%>%group_by(well)%>%
  summarise(total_cell = n())
Percentage_clus<-left_join(Number_cell_exp,classified_tracks)
Percentage_clus <- Percentage_clus%>%group_by(cluster2,tcell_line, well,exp_nr, organoid_line)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()

### Plot proportion per well and per cell type
Per<-ggplot(Percentage_clus, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ interaction(organoid_line,tcell_line))
Per<-Per+theme_void()+ scale_fill_manual(values=c("gold3",
                                                  "darkolivegreen3",
                                                  "seagreen3",
                                                  "forestgreen",
                                                  "dodgerblue",
                                                  "cyan1",
                                                  "indianred",
                                                  "firebrick",
                                                  "brown1"))+theme(aspect.ratio = 1,strip.text.x = element_text(angle = 90))
Per

ggsave(paste0(output_dir,"RF_ClassProp_WellvsCelltype.png"), device="png")


### Run a chisq.test to see if the distibution is different between conditions
library(reshape2)
table1<-acast(Percentage_clus, cluster2~tcell_line, value.var = "percentage")
table1[is.na(table1)]<-0  ## If any cluster are missing convert NA to 0

## Run Pearson Chi sq test to see if there is a different distribution of clusters between cell types
chisq.test(table1)

