set.seed(123)
master_corrected3<-readRDS(file = "master_corrected3_example")
##import T cell data:
library(dplyr)
library(dtwclust)
library(stats)
library(scales)

### Normalize per experiment
## Normalize and rescale all variables for multivariate clustering. For continuos values such as speed, displacement and mean dead dye intensity. Only values of the upper quantile are taken. This ensures a clearer separation of static and dead cells.

master_processed<-master_corrected3%>% 
  group_by(exp) %>% 
  mutate(z.disp = (displacement-mean(displacement))/sd(displacement),z.speed = (speed-mean(speed))/sd(speed), z.red = (red_lym-mean(red_lym))/sd(red_lym))%>%
  mutate(q.disp=ifelse(z.disp>(quantile(z.disp, p=0.75)),z.disp,min(z.disp)), q.speed=ifelse(z.speed>(quantile(z.speed, p=0.75)),z.speed,min(z.speed)),q.red=ifelse(z.red>(quantile(z.red, p=0.75)),z.red,min(z.red)))%>%
  mutate(q.disp=rescale(q.disp, to=c(0,1)),q.speed=rescale(q.speed, to=c(0,1)),q.red=rescale(q.red, to=c(0,1)),s.contact=rescale(contact, to=c(0,1)),s.contact_lym=rescale(contact_lym, to=c(0,1))) %>%
  mutate(q.disp=q.disp/mean(quantile(q.disp, p=0.9999999)),q.speed=q.speed/mean(quantile(q.speed, p=0.9999999)),q.red=q.red/mean(quantile(q.red, p=0.9999999)))%>%ungroup()

### Arrange data by time:
master_processed<-master_processed%>%group_by(TrackID)%>%arrange(Time)
###Split the data in a list of TrackIDs with multivariate data for each Track overtime
list_multivariate <- split(master_processed[,c("q.disp", "q.speed", "q.red","s.contact", "s.contact_lym")],master_processed$TrackID) 

##Set up parallel working for big datasets
# load parallel
library(parallel)
# create multi-process workers
workers <- makeCluster(detectCores()-2)
# load dtwclust in each one, and make them use 1 thread per worker
invisible(clusterEvalQ(workers, {
  library(dtwclust)
  RcppParallel::setThreadOptions(1L)
}))
# register your workers, e.g. with doParallel
require(doParallel)
registerDoParallel(workers)

###MULTIVARIATE cross-distance matrix calculate for different tracks
distmat <- proxy::dist(list_multivariate, method = "dtw")
matrix_distmat<-as.matrix(distmat)
## Store TrackID names
TrackID<-as.numeric(names(list_multivariate))
library(umap)
set.seed(123)
## Project cross-distance matrix in a UMAP
umap_dist<- umap(matrix_distmat,n_components=2,input="dist",init = "random", 
                 n_neighbors=3, min_dist=0.1, spread=1)  ### adjust parameters
#Visualize plot
plot(umap_dist$`layout`, 
     col=rgb(0,0,0,alpha=0.1), 
     pch=19,
     asp=0.4)
umap_1 <- as.data.frame(umap_dist$`layout`) 
Track2_umap<-cbind(TrackID,umap_1)
temp_df<-master_corrected3[,c("TrackID", "well","exp","cell_type")]
temp_df<- temp_df[!duplicated(temp_df$TrackID),]
umap_2 <- left_join(Track2_umap ,temp_df)

## Perform clustering. Select clusterig type that suits more your dataset
library(kmodR)
library(ggplot2)
set.seed(12)
km.mod <- kmod(umap_dist$`layout`,k=9, l=80) ## adjust parameters
umap_3 <- cbind(km.mod$XC_dist_sqr_assign, umap_2)
outliers<- rownames(km.mod[["L"]]) ##find outliers
umap_3 <-umap_3%>%filter(!TrackID %in% outliers) ##remove outliers
colnames(umap_3)[2]<- "cluster2"
##plot
umap_3$cluster2 = factor(umap_3$cluster2, levels=c("4","7","5","8","3", "2","9","6","1"))
levels(umap_3$cluster2) <- c("1","2","3","4","5","6", "7","8","9")  ##reorder clusters
ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(cluster2))) +  
  geom_point(size=2, alpha=0.6) + labs(color="cluster")+
  xlab("") + ylab("") +scale_color_manual(values = c("gold3",
                                                     "darkolivegreen3",
                                                     "seagreen3",
                                                     "blue3",
                                                     "dodgerblue",
                                                     "cyan1",
                                                     "indianred",
                                                     "firebrick",
                                                     "brown1" ))+
  ggtitle("umap Cluster ") +
  theme_light(base_size=20) +theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(), aspect.ratio=1)+coord_fixed()



#### To the original dataset add information on cluster type
master_clustered <- merge(master_processed ,umap_3[c("TrackID","cluster2")], by.x = "TrackID", by.y = "TrackID")
### Save Reference UMAP for training of random forest classifier
saveRDS(master_clustered, file = "New_Behavioral_Referance_map") ### store here your reference map that can be used to predict behaviors in new experiments

## Plot a heatmap to show the relative values of each behavior parameter
## Create a dataframe the summarizes the mean values for each parameter
sum_all <- master_clustered%>% select( speed, displacement, red_lym, contact2, contact_lym, cluster2, contact)%>% group_by(cluster2)%>%
  summarise(n_contact_org= mean(contact),displacement2 = median(displacement), speed = median(speed), interaction_T_cells= mean(contact_lym),death = median(red_lym))
## Rescale the values from each parameter
sum_all <- sum_all%>%mutate(n_contact_org= rescale(n_contact_org, to=c(0,100)),displacement2 = rescale(displacement2, to=c(0,100)), speed = rescale(speed, to=c(0,100)),interaction_T_cells= rescale(interaction_T_cells, to=c(0,100)), death =rescale(death, to=c(0,100)))

## Reshape dataframe
library(reshape2)
library(viridis)
sum_all<-melt(sum_all,id.vars = "cluster2")
sum_all$cluster2<-as.factor(sum_all$cluster2)
## Plot heatmap
gg <- ggplot(data = sum_all, aes(x = variable, cluster2, fill = value))
gg <- gg + geom_tile()
gg <- gg + scale_fill_viridis(option="C", name="AU")
gg <- gg + labs(x=NULL, y="Cluster", title="cluster represention")+theme(aspect.ratio=1.7,axis.text.x = element_text(angle = 45, hjust = 1))+ylim(rev(levels(sum_all$cluster2)))
gg




### Backproject the clustered data to the imaging dataset. Each TrackID (was coded to be made unique) has to be converted back to its original TrackID
## read master dataset that contains both TrackID and TrackID2
master<-readRDS("master_example_data")
### keep only the TrackID2 that were classified
master2 <-master%>%filter(TrackID2 %in% master_clustered$TrackID )
clustertype<-master_clustered[,c("TrackID", "cluster2")]
clustertype<- clustertype[!duplicated(clustertype$TrackID),]
master3 <- left_join(master2 ,clustertype, by=c("TrackID2"="TrackID"))

### Now for each well of interest export the corresponding TrackID assigned to a cluster. Ranks are specific for each well analyzed
##E.g
Ranks_2<-subset(master3,ranks==2) 
Ranks_2<-Ranks_2[!duplicated(Ranks_2$TrackID),c("TrackID","cluster2")] 
Ranks_2_list<-split(Ranks_2,Ranks_2$cluster2)
### Save this list that allows to identify in the imaging dataset to which cluster does each cell belong to.
write(paste(as.character(Ranks_2_list), sep="' '", collapse=", "), "Backproject.txt")

