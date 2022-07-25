set.seed(123)
library(dplyr)
library(stats)
library(tidyr)
library(scales)
library(randomForest)
library(ggplot2)
library(yaml)

## Import reference map that will be used as the ground truth dataset to train and test a random forest classifier
### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  ### !!!!!! Change the path to the BEHAV3D_config file here if running the code in RStudio !!!!!!
  reference_map <- "/Users/samdeblank/surfdrive/Shared/Dream3DLab (Groupfolder)/1.Projects/AIM_ALLImmune/3.Analysis/BEHAV3D_analysis/AIM_MB2_Exp21_Tcellstats/results/tcell_behavior/results/behavioral_reference_map.rds"
  model_output_path <- "/Users/samdeblank/surfdrive/Shared/Dream3DLab (Groupfolder)/1.Projects/AIM_ALLImmune/3.Analysis/BEHAV3D_analysis/AIM_MB2_Exp21_Tcellstats/results/train_randomforest/TrainedRandomForest.Rdata"
} else {
  option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL, 
                help="Path to behavioral_reference_map", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL, 
                help="Output path for the RandomForest model", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  if (is.null(opt$input)){
    print_help(opt_parser)
    stop("-i|--input, must be supplied", call.=FALSE)
  }
  if (is.null(opt$output)){
    print_help(opt_parser)
    stop("-o|--output, must be supplied", call.=FALSE)
  }
  reference_map <- readRDS(file = opt$input)
  model_output_path <- opt$output
}

output_dir = paste0(dirname(model_output_path), "/")
train_dataset<-readRDS(reference_map)

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

y <-as.factor(train_dataset1[,"cluster"]) ##strack ID
x <- as.matrix(train_dataset1[,-which(names(train_dataset1) == "TrackID" | names(train_dataset1) == "cluster")]) ##data for training

model <- randomForest(y=y,x=x,ntree=100, importance=TRUE)
predictors_importance<-importance(model)
predictors_importance
model$confusion[, 'class.error']

###error of the test dataset
test_error2 <- randomForest(y=y,x=x,ntree=100, importance=TRUE, 
                            xtest=validation_dataset[,2:15], ytest=as.factor(validation_dataset$cluster))



###Test for error per cluster using the validation dataset
## predict the cluster types in the validation dataset using the trained classifier
model_predict_train<-predict(model,validation_dataset[,-which(names(train_dataset1) == "TrackID" | names(train_dataset1) == "cluster")],type="response")
## Merge predicted and ground truth values
ground_truth<-validation_dataset[!duplicated(validation_dataset$TrackID),]
ground_truth_classified<-cbind(ground_truth,model_predict_train)
## Plot the ground truth vs predicted values
p2 <-ggplot(ground_truth_classified, aes(x=as.factor(model_predict_train),y=as.factor(cluster), color = as.factor(cluster))) + 
  geom_jitter(width = 0.2, height = 0.2)+ 
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
model_df$ind = as.factor(model_df$ind)
p2 <-ggplot(model_df , aes(x=tree,y=values, group = ind, color = as.factor(ind))) + 
  geom_line(size=1)+ 
  theme_classic() + theme(aspect.ratio = 0.5)+
  xlab("tree") +ylab("error") 

p2
ggsave(paste0(output_dir,"RF_TreeErrorEvolution.png"), device="png")

save(model, file = model_output_path)
