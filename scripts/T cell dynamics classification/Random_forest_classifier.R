library(dplyr)
library(stats)
library(tidyr)
library(scales)
library(randomForest)
library(ggplot2)
## Import reference map that will be used as the ground truth dataset to train and test a random forest classifier

train_dataset<-readRDS("Behavioral_Referance_map_git")
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
            contact=mean(s.contact),mean_contact2=mean(contact2),contact2=max(contact2), cluster=mean(cluster2))

validation_dataset <- train_dataset_2%>% group_by(TrackID)%>%
  summarise(mean_speed= mean(q.speed),median_speed= median(q.speed),speed_sd= sd(q.speed),q3_speed= quantile(q.speed,0.90),
            mean_displacement = mean(q.disp),median_displacement = median(q.disp),q3_disp= quantile(q.disp,0.90),
            displacement_sd=sd(q.disp),
            mean_red_lym = mean(q.red),red_lym_sd=sd(q.red),q3_red= quantile(q.red,0.90),
            contact=mean(s.contact),mean_contact2=mean(contact2),contact2=max(contact2), cluster=mean(cluster2))


train_dataset1<-as.data.frame(train_dataset1)  ### training dataset
##set names of columns to the same:
y <-as.factor(train_dataset1[,16]) ##strack ID
x <- as.matrix(train_dataset1[,2:15]) ##data for training

model <- randomForest(y=y,x=x,ntree=100, importance=TRUE)
predictors_importance<-importance(model)
predictors_importance
model$confusion[, 'class.error']

###error of the test dataset
test_error2 <- randomForest(y=y,x=x,ntree=100, importance=TRUE, 
                            xtest=validation_dataset[,2:15], ytest=as.factor(validation_dataset$cluster))



###Test for error per cluster using the validation dataset
## predict the cluster types in the validation dataset using the trained classifier
model_predict_train<-predict(model,validation_dataset[,2:15],type="response")
## Merge predicted and ground truth values
ground_truth<-validation_dataset[!duplicated(validation_dataset$TrackID),]
ground_truth_classified<-cbind(ground_truth,model_predict_train)
## Plot the ground truth vs predicted values
p2 <-ggplot(ground_truth_classified, aes(x=as.factor(model_predict_train),y=as.factor(cluster), color = as.factor(cluster))) + 
  geom_jitter(width = 0.2, height = 0.2)+ scale_colour_manual(values = c("gold3",
                                                                         "darkolivegreen3",
                                                                         "seagreen3",
                                                                         "blue3",
                                                                         "dodgerblue",
                                                                         "cyan1",
                                                                         "indianred",
                                                                         "firebrick",
                                                                         "brown1"))+
  theme_bw() + theme(aspect.ratio = 1)+
  xlab("model predicted") +ylab("ground_truth") 

p2

####Plot tree error evolution per cluster
model_df<-model[["err.rate"]]
model_df<-as.data.frame(model_df)
model_df$tree<-rownames(model_df)
model_df$tree<-as.numeric(model_df$tree)
model_df <- data.frame(tree= model_df$tree,stack(model_df,select=-tree))
model_df$ind = factor(model_df$ind, levels=c("OOB","1","2","3","4","5","6","7","8","9"))
p2 <-ggplot(model_df , aes(x=tree,y=values, group = ind, color = as.factor(ind))) + 
  geom_line(size=1)+ scale_colour_manual(values = c("grey","gold3",
                                              "darkolivegreen3",
                                              "seagreen3",
                                              "blue3",
                                              "dodgerblue",
                                              "cyan1",
                                              "indianred",
                                              "firebrick",
                                              "brown1"))+
  theme_classic() + theme(aspect.ratio = 0.5)+
  xlab("tree") +ylab("error") 

p2




###################################################################################################################################################################
#### Import new dataset to predict behavior (e.g. import dataset called "master_corrected3_example")
master_test<-readRDS("D:/R/scripts/T_cell paper/FINAL SCRIPTS_20210408/Fig2/b/github/example_t cell_data/master_corrected3_example")
### Normalize the data
master_test2<-master_test%>% ungroup()%>%
  group_by(exp) %>% 
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
cell_ID<-master_classified[!duplicated(master_classified$TrackID),c(2,8,9,10)] 


classified_tracks<-test_dataset_predicted
classified_tracks$cluster2<-classified_tracks$cluster
classified_tracks<-left_join(classified_tracks,cell_ID)
classified_tracks<-classified_tracks%>%arrange(cluster2)

### Quantify the number of cells per well
Number_cell_exp<-classified_tracks%>%group_by(well)%>%
  summarise(total_cell = n())
Percentage_clus<-left_join(Number_cell_exp,classified_tracks)
Percentage_clus <- Percentage_clus%>%group_by(cluster2,cell_type, well,exp)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()

### Plot proportion per well and per cell type
Per<-ggplot(Percentage_clus, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
Per <- Per + facet_grid(interaction(well,exp) ~ cell_type)
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
