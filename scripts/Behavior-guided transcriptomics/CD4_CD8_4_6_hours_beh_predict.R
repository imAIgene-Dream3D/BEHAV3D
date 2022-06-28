## Impor the processed tracks with a length that is coinciding the the length of the experiment in Fig 4a.
setwd()
library(dplyr)
master_CD3<-readRDS("master_CD3_scRNA_seq_inference_4")
## Import the training dataset
library(stats)
library(tidyr)
train_dataset<-readRDS("master_clust_Live7_REF_1")
set.seed(123)
train_dataset$cluster2<-as.numeric(train_dataset$cluster2)
train_dataset_grouped <- train_dataset %>% 
  group_by(TrackID) %>% 
  nest()%>%ungroup
train_dataset_1<- train_dataset_grouped%>%
  sample_frac(0.95)
train_dataset_grouped<-as.data.frame(train_dataset_grouped)%>% 
  unnest(cols = c(data))
train_dataset_1<-as.data.frame(train_dataset_1)%>% 
  unnest(cols = c(data))

train_dataset_2 <-train_dataset_grouped%>%filter(!TrackID%in% train_dataset_1$TrackID)

set.seed(321)
rsq <- function(x, y) summary(lm(y~x))$r.squared

train_dataset1 <- train_dataset_1%>% group_by(TrackID)%>%
  summarise(mean_speed= mean(q.speed),median_speed= median(q.speed),speed_sd= sd(q.speed),q3_speed= quantile(q.speed,0.90),
            mean_displacement = mean(q.disp),median_displacement = median(q.disp),q3_disp= quantile(q.disp,0.90),
            displacement_sd=sd(q.disp),
            mean_red_lym = mean(q.red),red_lym_sd=sd(q.red),q3_red= quantile(q.red,0.90),
            touch=mean(s.touch),mean_touch2=mean(touch2),touch2=max(touch2), cluster=mean(cluster2), V1=mean(V1), V2=mean(V2))

train_dataset1_2 <- train_dataset_2%>% group_by(TrackID)%>%
  summarise(mean_speed= mean(q.speed),median_speed= median(q.speed),speed_sd= sd(q.speed),q3_speed= quantile(q.speed,0.90),
            mean_displacement = mean(q.disp),median_displacement = median(q.disp),q3_disp= quantile(q.disp,0.90),
            displacement_sd=sd(q.disp),
            mean_red_lym = mean(q.red),red_lym_sd=sd(q.red),q3_red= quantile(q.red,0.90),
            touch=mean(s.touch),mean_touch2=mean(touch2),touch2=max(touch2), cluster=mean(cluster2), V1=mean(V1), V2=mean(V2))


### rename the dataset to predict:
library(scales)
master_CD32<-master_CD3%>% ungroup()%>%
  mutate(z.disp = (displacement-mean(displacement))/sd(displacement),z.speed = (speed-mean(speed))/sd(speed), z.red = (red_lym-mean(red_lym))/sd(red_lym))%>%
  mutate(q.disp=ifelse(z.disp>(quantile(z.disp, p=0.75)),z.disp,min(z.disp)), q.speed=ifelse(z.speed>(quantile(z.speed, p=0.75)),z.speed,min(z.speed)),q.red=ifelse(z.red>(quantile(z.red, p=0.75)),z.red,min(z.red)))%>%
  mutate(q.disp=rescale(q.disp, to=c(0,100)),q.speed=rescale(q.speed, to=c(0,100)),q.red=rescale(q.red, to=c(0,100)),s.contact=rescale(contact, to=c(0,1)),s.contact_lym=rescale(contact_lym, to=c(0,1)))%>%
  mutate(q.disp=q.disp/mean(quantile(q.disp, p=0.9999)),q.speed=q.speed/mean(quantile(q.speed, p=0.9999)),q.red=q.red/mean(quantile(q.red, p=0.9999)))%>%ungroup()

## calculate stats of the test dataset
test_dataset <- master_CD32%>%group_by(TrackID)%>% arrange(Time)%>%
  summarise(mean_speed= mean(q.speed),median_speed= median(q.speed),speed_sd= sd(q.speed),q3_speed= quantile(q.speed,0.90),
            mean_displacement = mean(q.disp),median_displacement = median(q.disp),
            displacement_sd=sd(q.disp),q3_disp= quantile(q.disp,0.90),
            mean_red_lym = mean(q.red),red_lym_sd=sd(q.red),q3_red= quantile(q.red,0.90),
            contact=mean(s.contact),mean_contact2=mean(contact2),contact2=max(contact2))
library(randomForest)
test1<-as.matrix(test_dataset[,2:15]) ### sataset to assign to a class
train_dataset1<-as.data.frame(train_dataset1)  ### training dataset
##set names of columns to the same:
colnames(train_dataset1)[2:15]<-colnames(test_dataset)[2:15]
y <-as.factor(train_dataset1[,16]) ##strack ID
x <- as.matrix(train_dataset1[,2:15]) ##data for training

model <- randomForest(y=y,x=x,ntree=100, importance=TRUE)
predictors_importance<-importance(model)
predictors_importance
model$confusion[, 'class.error']
model_predict<-predict(model,test1,type="response")


## check the classification of the train dataset, where do error occur
library(ggplot2)
library(scales)
test_dataset_predicted<-cbind(model_predict,test_dataset)
test_dataset_predicted$cluster<-test_dataset_predicted$model_predict
sum_all <- test_dataset_predicted%>% select( mean_speed, mean_displacement, mean_red_lym, cluster, contact)%>% group_by(cluster)%>%
  summarise(contact= mean(contact),displacement2 = median(mean_displacement), speed = median(mean_speed), death = median(mean_red_lym))

sum_all <- sum_all%>%mutate(contact= rescale(contact, to=c(0,100)),displacement2 = rescale(displacement2, to=c(0,100)), speed = rescale(speed, to=c(0,100)), death =rescale(death, to=c(0,100)))
library(reshape2)
sum_all<-melt(sum_all,id.vars = "cluster")
library(viridis)
library(plotly)
library(ggplot2)
sum_all$cluster<-as.factor(sum_all$cluster)
sum_all$cluster <- factor(sum_all$cluster,levels(sum_all$cluster)[9:1])
gg <- ggplot(data = sum_all, aes(x = variable, as.factor(cluster), fill = value))
gg <- gg + geom_tile()
gg <- gg + scale_fill_viridis(option="C", name="AU")
gg <- gg + labs(x=NULL, y="Cluster", title="cluster represenation")+theme(aspect.ratio=1.7,axis.text.x = element_text(angle = 45, hjust = 1))
gg




master_CD3_class<-test_dataset_predicted[,c(2,17)]

master_clust_Live6<-left_join(master_CD3,master_CD3_class, by="TrackID")


###calculate the mean engagement time for each signature:
mean_cont<- master_clust_Live6%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)


###calculate the mean engagement time for each signature:
mean_cont<- master_clust_Live6%>%group_by(cluster,TrackID)%>%
  summarise(contact_per_h= sum(contact))


min_max_cont<- mean_cont%>%group_by(cluster)%>%
  summarise(min= quantile(contact_per_h,0.05),max= quantile(contact_per_h,0.95))

saveRDS(min_max_cont, "min_max_cont_13T")

mean_cont_length<- master_clust_Live6%>%group_by(cluster, TrackID)%>%
  summarise(max_length= max(contact2))%>%group_by(cluster)%>%summarise(contact_length_per_h= mean(max_length)/2)

mean_cont<-left_join(mean_cont,mean_cont_length)

saveRDS(mean_cont, "min_contact_per_hour_13T")


which( colnames(master_clust_Live6)=="cell_type" ) 
cell_cell<-master_clust_Live6[!duplicated(master_clust_Live6$TrackID),c(2,12)] 
d_tsne_3<-test_dataset_predicted
d_tsne_3$cluster2<-d_tsne_3$cluster
d_tsne_3<-left_join(d_tsne_3,cell_cell)
d_tsne_3<-d_tsne_3%>%arrange(cluster)
#########################################################################################################
#Plot pie chart
library(dplyr)
Number_cell_exp<-d_tsne_3%>%group_by(cell_type)%>%
  summarise(total_cell = n())
Percentage_clus<-left_join(Number_cell_exp,d_tsne_3)
Percentage_clus <- Percentage_clus%>%group_by(cluster,cell_type)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()

### Plot circular map
Per<-ggplot(Percentage_clus, aes(fill=as.factor(cluster), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
Per <- Per + facet_grid(cell_type ~.)+theme_void()+theme(strip.text.x = element_text(angle = 90),strip.text.y = element_text(angle = 0))
Per<-Per+ scale_fill_manual(values=c("gold3",
                                     "darkolivegreen3",
                                     "seagreen3",
                                     "forestgreen",
                                     "dodgerblue",
                                     "cyan1",
                                     "indianred",
                                     "firebrick",
                                     "brown1"))
Per




#######################################################################################################################################################
######### Reproduce in silico the experiment performed in Fig 4a. Separate cells at each washing step and compute the behavirs of the cells that were in contact with an organoid vs the ones that weren't
#################################################################################################################################################
#### Define the timepoints where the washing steps occured. Instead of selecting just one timepoint, that gives a very small amount of tracks, we define the time as a time-range around 3 hours or 5 hours of imaging. these correspond to 4 and 6 hours of culture respectively
around_5h<-c(145:154)
around_3h<-c(86:95)

### calculate the frequency of behaviors per conditions for CD8 only:

### Simulate the selection at 5 hours of T cells that were or were not in contact with organoids at that timepoint
master_clust_Live6_h_CD8<-subset(master_clust_Live6, Time%in%around_5h& cell_type=="CD8" &!cluster=="1") ### Select cells at timepoint "around 5 hours" and remove cluster 1 since in the real experiment dead cells are excluded by FACS
master_clust_Live6_h_CD8<-master_clust_Live6_h_CD8%>%group_by(TrackID)%>%filter(Time==max(Time)) ## For each trackID select teh last timepoint in the 5 hours range

## Plor the engager and the non engager populations and their cluster behavior distribution
Plot_Eng_CD8<-ggplot(master_clust_Live6_h_CD8, aes(fill=as.factor(cluster), x=contact)) + 
  geom_bar( position="fill")+ coord_flip()+ scale_fill_manual(values=c(
    "darkolivegreen3",
    "seagreen3",
    "forestgreen",
    "dodgerblue",
    "cyan1",
    "indianred",
    "firebrick",
    "brown1"))+ggtitle("Engager & non-engager CD8")
Plot_Eng_CD8

### Now plot the super engagers and the never engagers. These cells under go two washing and separation steps:
### 1) at 3 hours they are separated into  contacting organoids and non-contacting organoids and then put back in culture. For the non-contacting organoids condition new organoids are added.
### 2) at 5 hours they are separated again into contact organoids and non=contacting organoids T cells.

### Orverall we recover two populations: Super-engagers> T cells that were in contact with organoids at 3 and 5 hours
### Never engagers> T cells that were not in contact with organoids at 3 hours and also not at 5 hours.

### Simulate the selection at 3 hours of T cells that were in contact with an orgnaoid
master_clust_Live4_h_EN_CD8<-subset(master_clust_Live6, Time%in%around_3h& cell_type=="CD8"& contact==1)
master_clust_Live4_h_EN_CD8<-master_clust_Live4_h_EN_CD8%>%group_by(TrackID)%>%filter(Time==min(Time))
### Simulate the selection at 3 hours of T cells that were not in contact with an orgnaoid
master_clust_Live4_h_NEN_CD8<-subset(master_clust_Live6, Time%in%around_3h & cell_type=="CD8" &contact==0)
master_clust_Live4_h_NEN_CD8<-master_clust_Live4_h_NEN_CD8%>%group_by(TrackID)%>%filter(Time==min(Time))

###Simulate the separation at 5 hours of T cells in contact with organoid or not, based on the previously selected cells
master_clust_Live6_h_SEN_CD8<-subset(master_clust_Live6_h_CD8, contact==1 &TrackID%in%master_clust_Live4_h_EN_CD8$TrackID)
master_clust_Live6_h_NEN_CD8<-subset(master_clust_Live6_h_CD8, contact==0 &TrackID%in%master_clust_Live4_h_NEN_CD8$TrackID)
master_clust_Live6_h_NEN_SEN_CD8<-rbind(master_clust_Live6_h_NEN_CD8,master_clust_Live6_h_SEN_CD8)
master_clust_Live6_h_NEN_SEN_CD8<-subset(master_clust_Live6_h_NEN_SEN_CD8, !cluster=="1") ## remove dying cell cluster, since will be delted by FACS

Plot_SEng_CD8<-ggplot(master_clust_Live6_h_NEN_SEN_CD8, aes(fill=as.factor(cluster), x=contact)) + 
  geom_bar( position="fill")+ coord_flip()+ scale_fill_manual(values=c(
    "darkolivegreen3",
    "seagreen3",
    "forestgreen",
    "dodgerblue",
    "cyan1",
    "indianred",
    "firebrick",
    "brown1"))+ggtitle("Super-engager & Never-engager CD8")
Plot_SEng_CD8




#### Repeat the same for the CD4 cells
### calculate the frequency of behaviors per conditions for CD4 only:
### Simulate the selection at 5 hours of T cells that were or were not in contact with organoids at that timepoint
master_clust_Live6_h_CD4<-subset(master_clust_Live6, Time%in%around_5h & cell_type=="CD4"&!cluster=="1")
master_clust_Live6_h_CD4<-master_clust_Live6_h_CD4%>%group_by(TrackID)%>%filter(Time==max(Time))

Plot_Eng_CD4<-ggplot(master_clust_Live6_h_CD4, aes(fill=as.factor(cluster), x=contact)) + 
  geom_bar( position="fill")+ coord_flip()+ scale_fill_manual(values=c(
    "darkolivegreen3",
    "seagreen3",
    "forestgreen",
    "dodgerblue",
    "cyan1",
    "indianred",
    "firebrick",
    "brown1"))+ggtitle("Engager & non-engager CD4")
Plot_Eng_CD4

### Simulate the selection at 3 hours of T cells that were in contact with an orgnaoid
master_clust_Live4_h_EN_CD4<-subset(master_clust_Live6, Time%in%around_3h & cell_type=="CD4" &contact==1)
master_clust_Live4_h_EN_CD8<-master_clust_Live4_h_EN_CD8%>%group_by(TrackID)%>%filter(Time==min(Time))
### Simulate the selection at 3 hours of T cells that were not in contact with an orgnaoid
master_clust_Live4_h_NEN_CD4<-subset(master_clust_Live6, Time%in%around_3h & cell_type=="CD4"&contact==0)
master_clust_Live4_h_NEN_CD4<-master_clust_Live4_h_NEN_CD4%>%group_by(TrackID)%>%filter(Time==min(Time))
###Simulate the separation at 5 hours of T cells in contact with organoid or not
master_clust_Live6_h_SEN_CD4<-subset(master_clust_Live6_h, contact==1 &TrackID%in%master_clust_Live4_h_EN_CD4$TrackID)
master_clust_Live6_h_NEN_CD4<-subset(master_clust_Live6_h, contact==0 &TrackID%in%master_clust_Live4_h_NEN_CD4$TrackID)
master_clust_Live6_h_NEN_SEN_CD4<-rbind(master_clust_Live6_h_NEN_CD4,master_clust_Live6_h_SEN_CD4)
master_clust_Live6_h_NEN_SEN_CD4<-subset(master_clust_Live6_h_NEN_SEN_CD4, !cluster=="1")

Plot_SEng_CD4<-ggplot(master_clust_Live6_h_NEN_SEN_CD4, aes(fill=as.factor(cluster), x=contact)) + 
  geom_bar( position="fill")+ coord_flip()+ scale_fill_manual(values=c(
    "darkolivegreen3",
    "seagreen3",
    "forestgreen",
    "dodgerblue",
    "cyan1",
    "indianred",
    "firebrick",
    "brown1"))+ggtitle("Super-engager & Never-engager CD4")
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
CD8_behav$cluster_prop<-100*CD8_behav$frec/CD8_behav$total_n

Plot_x<-ggplot(CD8_engagement, aes(fill=as.factor(cluster), x=engagement)) + 
  geom_bar( position="fill")+ coord_flip()+ scale_fill_manual(values=c(
    "darkolivegreen3",
    "seagreen3",
    "forestgreen",
    "dodgerblue",
    "cyan1",
    "indianred",
    "firebrick",
    "brown1"))+ggtitle("Super-engager & Never-engager CD8")
Plot_x
### save data to csv same folder where the pseudotime clustering from Farid
write.csv(CD8_behav,"CD8_engagement_behavior_freq.csv")


#### Calculate the mean contact time per condition for information (used in the paper)
## SE mean contact time:
master_clust_Live6_h_CD8_se_pop<-subset(CD8_engagement, engagement=="super-engaged")
###calculate the mean engagement time for each signature:
mean_cont_SE<- master_clust_Live6%>%filter(TrackID%in%master_clust_Live6_h_CD8_se_pop$TrackID)%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)
mean_cont_SE$exp<-"super-engaged"


## E mean contact time:
master_clust_Live6_h_CD8_e_pop<-subset(CD8_engagement, engagement=="engager")
###calculate the mean engagement time for each signature:
mean_cont_E<- master_clust_Live6%>%filter(TrackID%in%master_clust_Live6_h_CD8_e_pop$TrackID)%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)
mean_cont_E$exp<-"engager"


## NonE mean contact time:
master_clust_Live6_h_CD8_noe_pop<-subset(CD8_engagement, engagement=="nonengager")
###calculate the mean engagement time for each signature:
mean_cont_NoE<- master_clust_Live6%>%filter(TrackID%in%master_clust_Live6_h_CD8_noe_pop$TrackID)%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)
mean_cont_NoE$exp<-"nonengager"


## NonE mean contact time:
master_clust_Live6_h_CD8_NE_pop<-subset(CD8_engagement, engagement=="never-engaged")
###calculate the mean engagement time for each signature:
mean_cont_NE<- master_clust_Live6%>%filter(TrackID%in%master_clust_Live6_h_CD8_NE_pop$TrackID)%>%group_by(cluster)%>%
  summarise(contact_per_h= mean(contact)*60)

mean_cont_NE$exp<-"never-engaged"

mean_cont<-rbind(mean_cont_NE,mean_cont_NoE,mean_cont_E,mean_cont_SE)

saveRDS(mean_cont, "min_contact_per_hour_13T_per_exp")



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
CD4_behav$cluster_prop<-100*CD4_behav$frec/CD4_behav$total_n


Plot_y<-ggplot(CD4_engagement, aes(fill=as.factor(cluster), x=engagement)) + 
  geom_bar( position="fill")+ coord_flip()+ scale_fill_manual(values=c(
    "darkolivegreen3",
    "seagreen3",
    "forestgreen",
    "dodgerblue",
    "cyan1",
    "indianred",
    "firebrick",
    "brown1"))+ggtitle("Super-engager & Never-engager CD4")
Plot_y



### save data to csv same folder where the pseudotime clustering from Farid
write.csv(CD4_behav,"CD4_engagement_behavior_freq.csv")



## Save output for figure Fig 4a
pdf("Engager_Super_ENG_proportions_CD4_CD8.pdf")
Plot_x
Plot_y
dev.off()
