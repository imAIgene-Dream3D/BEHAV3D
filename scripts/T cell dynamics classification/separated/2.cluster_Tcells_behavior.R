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
  tracks_provided=""
} else {
  option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, 
                help="Path to the BEHAV3D config file", metavar="character"),
    make_option(c("-t", "--tracks_rds"), type="character", default=NULL, 
                help="(Optional) Path to RDS file containing processed T cell track data", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  if (is.null(opt$config)){
    print_help(opt_parser)
    stop("Config file -c|--config, must be supplied", call.=FALSE)
  }
  tracks_provided=opt$tracks_rds
  
}

master_processed<-readRDS(tracks_provided)

#######################################################
###########     UNSUPERVISED CLUSTERING     ###########
#######################################################

library(dplyr)
library(dtwclust)
library(stats)
library(scales)

master_processed<-master_corrected3%>% 
  group_by(exp_nr) %>% 
  mutate(z.disp = (displacement-mean(displacement))/sd(displacement),z.speed = (speed-mean(speed))/sd(speed), z.red = (red_lym-mean(red_lym))/sd(red_lym))%>%
  mutate(q.disp=ifelse(z.disp>(quantile(z.disp, p=0.75)),z.disp,min(z.disp)), q.speed=ifelse(z.speed>(quantile(z.speed, p=0.75)),z.speed,min(z.speed)),q.red=ifelse(z.red>(quantile(z.red, p=0.75)),z.red,min(z.red)))%>%
  mutate(q.disp=rescale(q.disp, to=c(0,1)),q.speed=rescale(q.speed, to=c(0,1)),q.red=rescale(q.red, to=c(0,1)),s.contact=rescale(contact, to=c(0,1)),s.contact_lym=rescale(contact_lym, to=c(0,1))) %>%
  mutate(q.disp=q.disp/mean(quantile(q.disp, p=0.9999999)),q.speed=q.speed/mean(quantile(q.speed, p=0.9999999)),q.red=q.red/mean(quantile(q.red, p=0.9999999)))%>%ungroup()

### Arrange data by time:
master_processed<-master_processed%>%group_by(TrackID)%>%arrange(Time)
master_processed$TrackID = as.character(master_processed$TrackID)
# master_processed_ref<-master_processed_ref%>%group_by(TrackID)%>%arrange(Time)

saveRDS(master_processed, file = paste0(output_dir,"tcell_track_features.rds"))

###Split the data in a list of TrackIDs with multivariate data for each Track overtime
list_multivariate <- split(master_processed[,c("q.disp", "q.speed", "q.red","s.contact", "s.contact_lym")],master_processed$TrackID) 
# list_multivariate_ref <- split(master_processed_ref[,c("q.disp", "q.speed", "q.red","s.contact", "s.contact_lym")],master_processed_ref$TrackID) 

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
TrackID<-as.character(names(list_multivariate))

library(umap)
## Project cross-distance matrix in a UMAP
umap_dist<- umap(matrix_distmat,n_components=2,input="dist",init = "random", 
                 n_neighbors=pars$umap_n_neighbors, min_dist=pars$umap_minimal_distance, spread=1)  ### adjust parameters
#Visualize plot
pdf(file=paste0(output_dir,"Umap_unclustered.pdf"))

plot(umap_dist$`layout`, 
     col=rgb(0,0,0,alpha=0.1), 
     pch=19,
     asp=0.4)

dev.off()

umap_1 <- as.data.frame(umap_dist$`layout`) 

Track2_umap<-cbind(TrackID,umap_1)
temp_df<-master_corrected3[,c("TrackID", "well","exp_nr","organoid_line")]
temp_df<- temp_df[!duplicated(temp_df$TrackID),]
umap_2 <- left_join(Track2_umap ,temp_df)

## Perform clustering. Select clusterig type that suits more your dataset
library(ggplot2)

km.norm <- kmeans(umap_dist$`layout`,pars$nr_of_clusters, nstart = 100)
umap_3 <- cbind(km.norm$cluster, umap_2)
colnames(umap_3)[1]<- "cluster2"
## remove the outlier cluster if necessary
umap_3<-subset(umap_3, cluster2!=0)


ggplot(umap_3, aes(x=V1, y=V2, color=as.factor(cluster2))) +  
  geom_point(size=2, alpha=0.6) + labs(color="cluster")+
  xlab("") + ylab("")+
  ggtitle("umap Cluster ") +
  theme_light(base_size=20) +theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(), aspect.ratio=1)+coord_fixed()

ggsave(
  paste0(output_dir,"Umap_clustered.pdf"), 
  device="pdf", height=210, width=297, units="mm"
)

#### To the original dataset add information on cluster type
master_clustered <- merge(master_processed ,umap_3[c("TrackID","cluster2")], by.x = "TrackID", by.y = "TrackID")
# master_clustered$cluster = master_clustered$cluster2
### Save Reference UMAP for training of random forest classifier
saveRDS(master_clustered, file = paste0(output_dir,"behavioral_reference_map.rds")) ### store here your reference map that can be used to predict behaviors in new experiments

## Plot a heatmap to show the relative values of each behavior parameter
## Create a dataframe the summarizes the mean values for each parameter
sum_all <- master_clustered%>% select( speed, displacement, red_lym, contact2, contact_lym, cluster2, contact)%>% group_by(cluster2)%>%
  summarise(contact_len=mean(contact2),n_contact_org= mean(contact),displacement2 = median(displacement), speed = median(speed), interaction_T_cells= mean(contact_lym),death = median(red_lym))
## Rescale the values from each parameter
sum_all <- sum_all%>%mutate(contact_len= rescale(contact_len, to=c(0,100)) ,n_contact_org= rescale(n_contact_org, to=c(0,100)),displacement2 = rescale(displacement2, to=c(0,100)), speed = rescale(speed, to=c(0,100)),interaction_T_cells= rescale(interaction_T_cells, to=c(0,100)), death =rescale(death, to=c(0,100)))

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

ggsave(
  paste0(output_dir,"Cluster_heatmap.pdf"), 
  device="pdf", height=210, width=297, units="mm"
)

master_clustered_info<-master_clustered[!duplicated(master_clustered$TrackID),c("TrackID","organoid_line","tcell_line","exp_nr", "well", "cluster2")] 

master_clustered2 =master_clustered%>%group_by(TrackID)%>% arrange(Time)%>%
  summarise(mean_speed= mean(q.speed),median_speed= median(q.speed),speed_sd= sd(q.speed),q3_speed= quantile(q.speed,0.90),
            mean_displacement = mean(q.disp),median_displacement = median(q.disp),
            displacement_sd=sd(q.disp),q3_disp= quantile(q.disp,0.90),
            mean_red_lym = mean(q.red),red_lym_sd=sd(q.red),q3_red= quantile(q.red,0.90),
            contact=mean(s.contact),mean_contact2=mean(contact2),contact2=max(contact2))

master_clustered2<-left_join(master_clustered2,master_clustered_info)


### Quantify the number of cells per well
Number_cell_exp<-master_clustered2%>%group_by(well, exp_nr, tcell_line, organoid_line)%>%
  summarise(total_cell = n())
Percentage_clus<-left_join(Number_cell_exp, master_clustered2)
Percentage_clus <- Percentage_clus%>%group_by(cluster2,tcell_line, well,exp_nr, organoid_line)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()

Per<-ggplot(Percentage_clus, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ tcell_line)
# Per <- Per + facet_grid(exp_nr + well + organoid_line  ~ tcell_line)

# Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ interaction(organoid_line,tcell_line))
Per<-Per+theme_void() +
  theme(aspect.ratio = 0.2,strip.text.x = element_text(angle = 90))

Per

ggsave(
  paste0(output_dir,"umap_cluster_percentage_bars_separate.pdf"), 
  device="pdf", height=210, width=297, units="mm"
)

############### Combine all wells
Number_cell_exp_combined<-master_clustered2%>%group_by(exp_nr, tcell_line, organoid_line)%>%
  summarise(total_cell = n())
Percentage_clus_combined<-left_join(Number_cell_exp_combined,master_clustered2)
Percentage_clus_combined <- Percentage_clus_combined%>%group_by(cluster2,exp_nr, tcell_line, organoid_line)%>%
  summarise(total_cell = mean(total_cell), num_cluster=n())%>%mutate(percentage=num_cluster*100/total_cell)%>%ungroup()

Per<-ggplot(Percentage_clus_combined, aes(fill=as.factor(cluster2), y=percentage, x="")) + 
  geom_bar( stat="identity", position="fill")+ coord_flip()+ scale_y_reverse()
Per <- Per + facet_grid(interaction(exp_nr,organoid_line)  ~ tcell_line)
# Per <- Per + facet_grid(exp_nr + well + organoid_line  ~ tcell_line)

# Per <- Per + facet_grid(interaction(exp_nr,well,organoid_line)  ~ interaction(organoid_line,tcell_line))
Per<-Per+theme_void() +
  theme(aspect.ratio = 0.2,strip.text.x = element_text(angle = 90))

Per

ggsave(
  paste0(output_dir,"umap_cluster_percentage_bars_combined.pdf"), 
  device="pdf", height=210, width=297, units="mm"
)