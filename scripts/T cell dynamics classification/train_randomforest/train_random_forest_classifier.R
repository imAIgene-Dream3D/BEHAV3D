library(dplyr)
library(stats)
library(tidyr)
library(scales)
library(randomForest)
library(ggplot2)
library(yaml)
## Import reference map that will be used as the ground truth dataset to train and test a random forest classifier

pars = yaml.load_file("/Users/samdeblank/surfdrive/Shared/T cell paper/Revision_nature_biotech/Sam_analysis/WT1_pooled/BEHAV3D_config.yml")

reference_map <- pars$reference_map
model_path <- pars$randomforest_path

# reference_map <- "/Users/samdeblank/OneDrive - Prinses Maxima Centrum/github/BEHAV3D/scripts/T cell dynamics classification/Behavioral reference map/"
# output_dir="/Users/samdeblank/Documents/tcell_paper/unprocessed/2021-06-10_WT1_n1/1.split_data/2021-06-10_WT1_n1(2)_36T_[ims1_2021-09-09T14-14-13.594]_dt_Statistics/0.Analysis/"

  
train_dataset<-readRDS(reference_map)
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
ggsave(paste0(output_dir,"RF_GTvsPred.png"), device="png")

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
ggsave(paste0(output_dir,"RF_TreeErrorEvolution.png"), device="png")

save(model, file = model_path)
