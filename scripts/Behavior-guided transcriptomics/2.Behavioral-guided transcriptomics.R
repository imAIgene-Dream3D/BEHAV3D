library(scales)
library(ggplot2)
library(dplyr)
library(viridis)
library(Seurat)
library(reshape2)
library(pheatmap)
library(gplots)
library(yaml)
## Import scRNA seq dataset
if (interactive()) {
  ### !!!!!! Change the path to the BEHAV3D_config file here if running the code in RStudio !!!!!!
  pars = yaml.load_file("/Users/samdeblank/OneDrive - Prinses Maxima Centrum/github/BEHAV3D-2.0/demos/behavioral_transcriptomics_demo/BEHAV3D_config.yml")
} else {
  option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, 
                help="Path to the BEHAV3D config file", metavar="character"),
    make_option(c("-t", "--tracks_rds"), type="character", default=NULL, 
                help="(Optional) Path to RDS file containing clustered T cell track data", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  if (is.null(opt$config)){
    print_help(opt_parser)
    stop("Config file -c|--config, must be supplied", call.=FALSE)
  }
  pars = yaml.load_file(opt$config)
  tracks_provided=opt$tracks_rds
}

scRNA_seq_dataset<-readRDS(pars$scRNAseq_dataset)
main_output_dir=pars$output_dir
output_dir=paste0(main_output_dir,"/transcriptomics/results/")
dir.create(output_dir, recursive=TRUE)

p1<-DimPlot(scRNA_seq_dataset, reduction = "umap", pt.size=1, label.size = 8)+theme(aspect.ratio = 1)
p1


###Import the proportion of behavioral signatures per exp condition calculated in silico based on imaging data:
CD8_behav<-read.csv(paste0(output_dir, "CD8_engagement_behavior_freq.csv"))
CD8_behav$behavioral_cluster=factor(CD8_behav$behavioral_cluster, levels=c(1,2,3,4,5,6,7,8,9))

prob_output_dir=paste0(output_dir,"/probability_map/")
dir.create(prob_output_dir)
###### For behavioral inference we will build a probability map for cells are similar to each other. For this we cluster the scRNA seq data with different resolutions (creating a range of bigger to smaller clusters)
###### This allows to create an accurate probability map, assuming that cells that are similar to each other also have a similar probabilty map.
## For loop for different resoltuon of clusters.
for (r in seq_along(1:7)){
  scRNA_seq_dataset_1 <- FindClusters(scRNA_seq_dataset,method=RunUMAP, resolution = r)
  
  p2<-DimPlot(scRNA_seq_dataset_1, reduction = "umap", pt.size=1, label.size = 8)+theme(aspect.ratio = 1)
  p2
  
  ### for each cell identify which UMAP_cl and which exp_condition does it belong to
  CD8_metadata<-scRNA_seq_dataset_1@meta.data
  UMAP_cl<-scRNA_seq_dataset_1@active.ident
  CD8_metadata<-cbind(CD8_metadata,UMAP_cl)
  CD8_metadata$cell_ID<-row.names(CD8_metadata)
  
  ## calculate the proportion of each exp_condition pero UMAP_cluster
  
  CD8_freq<-CD8_metadata %>% group_by(UMAP_cl,exp_condition) %>%summarize(frec = n())
  engagement_n<-CD8_metadata %>% group_by(UMAP_cl) %>%summarize(total_n = n())
  CD8_freq<-left_join(CD8_freq,engagement_n)
  CD8_freq$perc_eng_UMAP<-CD8_freq$frec/CD8_freq$total_n
  
  
  ### Join the Information on proportions of UMAP_cl // behavioral signature // exp_condition
  CD8_behav_joined<-merge(CD8_behav[c("behavioral_cluster","exp_condition","cluster_percentage")],CD8_freq[c("UMAP_cl","exp_condition","perc_eng_UMAP")])
  
  
  ### For each UMAP_cl calculate the proportion of each Behavioral signature:
  
  List3 = list()
  for (m in unique(CD8_behav_joined$UMAP_cl)){
    CD8_behav_joined_m<-subset(CD8_behav_joined, UMAP_cl==m)  ### for each UMAP cluster calculate what is the chance of each behav_signature
    List = list()
    List2 = list()
    for (i in unique(CD8_behav_joined_m$behavioral_cluster)){
      behavioral_cluster<-subset(CD8_behav_joined_m,behavioral_cluster==i) ### check for each beahv_signature what is the chance of it being present in a certain engagement state
      behavioral_cluster<-as.data.frame(sum(behavioral_cluster$cluster_percentage*behavioral_cluster$perc_eng_UMAP))
      behav_sign_n<-i
      List[[length(List)+1]] <-behavioral_cluster 
      List2[[length(List2)+1]] <-behav_sign_n
      ## store to list object
    }
    CD8_behav_joined_m_pred <- data.frame(matrix(unlist(List))) ## convert List to dta
    behavioral_cluster <- data.frame(matrix(unlist(List2))) ## convert List to dta
    rm(List)
    rm(List2)
    CD8_behav_joined_m_pred<-cbind(CD8_behav_joined_m_pred,behavioral_cluster)
    colnames(CD8_behav_joined_m_pred)<-c("composition","behavioral_cluster")
    CD8_behav_joined_m_pred$UMAP_cl<-m
    List3[[length(List3)+1]] <-CD8_behav_joined_m_pred
  }
  CD8_behav_joined_pred <- do.call(rbind, List3)
  CD8_behav_joined_pred$UMAP_cl=factor(CD8_behav_joined_pred$UMAP_cl, levels=levels(scRNA_seq_dataset_1@active.ident))
  CD8_behav_joined_pred$behavioral_cluster= factor(CD8_behav_joined_pred$behavioral_cluster,levels=c(1,2,3,4,5,6,7,8,9))
  g5<-ggplot(data=CD8_behav_joined_pred, aes(x=UMAP_cl, y=composition, fill=as.factor(behavioral_cluster))) +
    geom_bar(stat="identity")+ 
    scale_fill_manual(
      values=c(
        "gold3",
        "darkolivegreen3",
        "seagreen3",
        "forestgreen",
        "dodgerblue",
        "cyan1",
        "indianred",
        "firebrick",
        "brown1"),drop = FALSE)+
    ggtitle("Behavior enrichment")+
    scale_x_discrete(breaks=levels(scRNA_seq_dataset_1@active.ident),drop=FALSE)+
    scale_y_continuous(breaks=c(0,25,50,75,100), labels=c(0,25,50,75,100), limits=c(0,100))
  plot(g5)
  
  ggsave(paste0(prob_output_dir, "behavior_enrichment_v_transcriptomics_umap_reslevel", r, ".pdf"))
  CD8_beh_run_MED1<-as.data.frame(acast(CD8_behav_joined_pred, UMAP_cl~behavioral_cluster,value.var="composition"))
  
  CD8_beh_run_MED1$UMAP_cl<-row.names(CD8_beh_run_MED1)
  ##Join with scRNA dataset metadata to make a and save dataframe with cell ID assigned to their probability of falling into each behavioral_sign:
  CD8_metadata2_run_MED1<-left_join(CD8_metadata[c(11,12)],CD8_beh_run_MED1)
  
  write.csv(CD8_metadata2_run_MED1[-c(1:2)],paste0(prob_output_dir, "Probability_map_UMAP_cl_", r,".csv"), row.names = F)
}

######################################################################################################################################
###Calculate the average probability map per cell based on probability obtained with difference clustering resolution##################################

## Import all the probability maps obtained with different cluster resolutions
temp = list.files(path=prob_output_dir, pattern="*.csv", full.names=T)
named.list <- lapply(temp, read.csv) 
library(abind)
arr <-abind(named.list, along = 3) ### create a 3dimensional array with each cellID as a row
##Where probability is set as NA change it to 0:
arr[is.na(arr)] <- 0
### Calculate the average probabability for all the 7 runs, for each individual cells and each Behavioral signature
mean_prob<-rowMeans(arr, dims = 2)
mean_prob<-as.data.frame(mean_prob)  ### Convert to a dataframe
for (i in levels(CD8_behav$behavioral_cluster)){
  if (! paste0("X",i) %in% colnames(mean_prob)){
    mean_prob[paste0("X",i)] <- 0
  }
}


row.names(mean_prob)<-CD8_metadata2_run_MED1$cell_ID  ## set cellID as row.names
mean_prob$cell_ID<-row.names(mean_prob)
###Incorporate pseudotime information:
mean_prob<-left_join(mean_prob,CD8_metadata[c(7,12)])
### rescale the pseudotime
mean_prob$Pseudotime<-rescale(mean_prob$Pseudotime, to=c(0,1))

### Multiply the probability of engagement along the Pseudotime. So the further along the trajectory more strength for the engager populations.

mean_prob$X9<-mean_prob$X9*mean_prob$Pseudotime

mean_prob$X8<-mean_prob$X8*mean_prob$Pseudotime
mean_prob$X7<-mean_prob$X7*mean_prob$Pseudotime
### Add the probability information to the seurat object, creating new metadata for each behavioral signature
### First normalize the probability
mean_prob2<-as.data.frame(scale(mean_prob[c(1:9)]))

scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X0, col.name = "Medium_exposed")
p3<-FeaturePlot(scRNA_seq_dataset, 'Medium_exposed')+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)


scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X2, col.name = "Static")
p4<-FeaturePlot(scRNA_seq_dataset, 'Static')+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)

scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X3, col.name = "Lazy")
p5<-FeaturePlot(scRNA_seq_dataset, 'Lazy')+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)


scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X4, col.name = "Slow_scanner")
p6<-FeaturePlot(scRNA_seq_dataset, 'Slow_scanner')+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)


scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X5, col.name = "Medium_scanner")
p7<-FeaturePlot(scRNA_seq_dataset, 'Medium_scanner')+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)


scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X6, col.name = "Super_scanner")
p8<-FeaturePlot(scRNA_seq_dataset, 'Super_scanner')+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)


scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X7, col.name = "Tickler")
p9<-FeaturePlot(scRNA_seq_dataset, 'Tickler',max.cutoff = 2.8)+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)

scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X8, col.name = "Engager")
p10<-FeaturePlot(scRNA_seq_dataset, 'Engager')+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)

scRNA_seq_dataset<-AddMetaData(scRNA_seq_dataset, mean_prob2$X9, col.name = "Super_Engager")
p11<-FeaturePlot(scRNA_seq_dataset, 'Super_Engager')+ scale_colour_viridis(option = "B")+theme(aspect.ratio = 1)

###Update the metadata object:
scRNA_seq_dataset_meta<-scRNA_seq_dataset@meta.data
### Create a Histogram of probability of belonging to a certain population:
scRNA_seq_dataset_meta$Medium_exposed<-rescale(scRNA_seq_dataset_meta$Medium_exposed, to=c(0,1))
scRNA_seq_dataset_meta$Static<-rescale(scRNA_seq_dataset_meta$Static, to=c(0,1))
scRNA_seq_dataset_meta$Lazy<-rescale(scRNA_seq_dataset_meta$Lazy, to=c(0,1))
scRNA_seq_dataset_meta$Slow_scanner<-rescale(scRNA_seq_dataset_meta$Slow_scanner, to=c(0,1))
scRNA_seq_dataset_meta$Medium_scanner<-rescale(scRNA_seq_dataset_meta$Medium_scanner, to=c(0,1))
scRNA_seq_dataset_meta$Super_scanner<-rescale(scRNA_seq_dataset_meta$Super_scanner, to=c(0,1))
scRNA_seq_dataset_meta$Engager<-rescale(scRNA_seq_dataset_meta$Engager, to=c(0,1))
scRNA_seq_dataset_meta$Super_Engager<-rescale(scRNA_seq_dataset_meta$Super_Engager, to=c(0,1))
### ticklers have several outlier values above: 2.4. So they are converted to the maximal value before rescaling
hist(scRNA_seq_dataset_meta$Tickler)
sort(boxplot.stats(scRNA_seq_dataset_meta$Tickler)$out)
scRNA_seq_dataset_meta$Tickler<-ifelse(scRNA_seq_dataset_meta$Tickler>quantile(scRNA_seq_dataset_meta$Tickler,0.98), quantile(scRNA_seq_dataset_meta$Tickler,0.98), scRNA_seq_dataset_meta$Tickler)
scRNA_seq_dataset_meta$Tickler<-rescale(scRNA_seq_dataset_meta$Tickler, to=c(0,1))



###heatmap per individual behavior and with pheatmap

### Include rescaled pseudotime in the heatmap in order to organize the cells by pseudotime
mat_exp = as.data.frame(cbind(rescale(scRNA_seq_dataset_meta$Pseudotime, to=c(0,1)),scRNA_seq_dataset_meta$Medium_exposed,rowMeans(as.matrix(scRNA_seq_dataset_meta[c("Static","Lazy","Slow_scanner","Medium_scanner","Super_scanner")])),scRNA_seq_dataset_meta$Tickler,scRNA_seq_dataset_meta$Engager, scRNA_seq_dataset_meta$Super_Engager))
mat_exp2 =mat_exp[order(mat_exp$V1),]
mat<-cbind(mat_exp2$V1,mat_exp2$V2, mat_exp2$V3,mat_exp2$V4,mat_exp2$V5,mat_exp2$V6)
### Plot heatmap
heatmap.2(t(mat),scale = c("none"), trace="n", col=viridis,Colv = NA, 
          dendrogram = "none", labCol = "", labRow = c("Pseudotime","No target control","Static to Scanner","Tickler","Engager","Super-engager"),Rowv=FALSE, cexRow = 0.5)


### Plot the pseudotime line on top of heatmap. in order to see at which pseudotime the behavior changes
df_hm = as.data.frame(cbind(scRNA_seq_dataset_meta$Pseudotime,scRNA_seq_dataset_meta$Medium_exposed,rowMeans(as.matrix(scRNA_seq_dataset_meta[c("Static","Lazy","Slow_scanner","Medium_scanner","Super_scanner")])),scRNA_seq_dataset_meta$Tickler,scRNA_seq_dataset_meta$Engager, scRNA_seq_dataset_meta$Super_Engager))

df_hm<- df_hm[with(df_hm, order(V1)), ]
df_hm$x<-seq(1:1814)


###smooth the line:
plot_df_hm_smooth <- as.data.frame(predict(smooth.spline(df_hm$x, 
                                                         df_hm$V1, 
                                                         cv=TRUE,lambda=.0015), x = df_hm$x),
                                   col.names = c('cells', 'Pseudotime'))



df_hm<-cbind(df_hm,plot_df_hm_smooth[c("Pseudotime")])

p16<-ggplot(df_hm) + 
  geom_line( mapping=aes(x = x, y = Pseudotime),color="chartreuse3", size=2) +theme_classic()+theme(aspect.ratio = 0.5)+
  scale_y_continuous(name = "Pseudotime") 


p16


### Save the heatmap and pseudotime line that goes on top
pdf("CD8_probability.pdf")
heatmap.2(t(mat),scale = c("none"), trace="n", col=viridis,Colv = NA, 
          dendrogram = "none", labCol = "", labRow = c("Pseudotime","No target control","Static to Scanner","Tickler","Engager","Super-engager"),Rowv=FALSE, cexRow = 0.5)
p16

dev.off()



### separate the Pseudotime into groups based on the change of behavioral probability:
scRNA_seq_dataset$PT_cluster<- ifelse(scRNA_seq_dataset$Pseudotime<14, "TEGs_alone", "tumor_exposed")

scRNA_seq_dataset$PT_cluster<- ifelse(scRNA_seq_dataset$Pseudotime>21, "engagement", scRNA_seq_dataset$PT_cluster)

scRNA_seq_dataset$PT_cluster<- ifelse(scRNA_seq_dataset$Pseudotime>41, "prolongued_engagement", scRNA_seq_dataset$PT_cluster)

### Add PT_cluster as a new metadata object
scRNA_seq_dataset <- SetIdent(object = scRNA_seq_dataset, value = "PT_cluster")

levels(x = scRNA_seq_dataset) <- c( "TEGs_alone" ,"tumor_exposed","engagement",  "prolongued_engagement" )

### Save the processed Seurat object
setwd("DIRECTORY WHERE TO SAVE THE PSEUDOTIME")
saveRDS(scRNA_seq_dataset,"scRNA_seq_dataset_processed.rds")
