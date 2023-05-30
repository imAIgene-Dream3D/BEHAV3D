### Backproject the clustered data to the imaging dataset. Each TrackID (was coded to be made unique) has to be converted back to its original TrackID

library(yaml)
library(dplyr)
### Checks if being run in GUI (e.g. Rstudio) or command line
if (interactive()) {
  ### !!!!!! Change the path to the BEHAV3D_config file here if running the code in RStudio !!!!!!
  ### Demo path
  BEHAV3D_dir = paste0(dirname(dirname(dirname(rstudioapi::getSourceEditorContext()$path))),"/")
  pars = yaml.load_file(paste0(BEHAV3D_dir, "/demos/tcell_demo/BEHAV3D_config.yml"))
  
  ### For your own file, uncomment following line and add own path to the BEHAV3D_config.yml
  # pars = yaml.load_file("")
  
} else {
  option_list = list(
    make_option(c("-c", "--config"), type="character", default=NULL, 
                help="Path to the BEHAV3D config file", metavar="character"),
    make_option(c("-f", "--force_redo"), action="store_true", default=FALSE, 
                help="Force the pipeline to re-import data even if files exists"),
    make_option(c("-t", "--tracks_rds"), type="character", default=NULL, 
                help="(Optional) Path to RDS file containing processed T cell track data", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  if (is.null(opt$config)){
    print_help(opt_parser)
    stop("Config file -c|--config, must be supplied", call.=FALSE)
  }
  pars = yaml.load_file(opt$config)
  force_redo=opt$force_redo
  tracks_provided=opt$tracks_rds
}

output_dir=paste0(pars$output_dir,"/tcell_behavior/results/")

#import metadata
pat = pars$metadata_csv
metadata=read.csv(pars$metadata_csv, sep="\t", check.names=FALSE)
## read master dataset that contains both TrackID and TrackID2
master<-readRDS(paste0(output_dir,"raw_tcell_track_data.rds"))
classified_tracks<-readRDS(paste0(output_dir,"classified_tcell_track_data_summary.rds"))
### keep only the TrackID2 that were classified
master2 <-master%>%filter(TrackID %in% classified_tracks$TrackID )
clustertype<-classified_tracks[,c("TrackID", "cluster2")]
clustertype<- clustertype[!duplicated(clustertype$TrackID),]
master3 <- left_join(master2 ,clustertype, by=c("TrackID"))

### Now for each well of interest export the corresponding TrackID assigned to a cluster. 
for(to_backproject in unique(master3$basename)){
  To_export<-subset(master3,basename==to_backproject)
  To_export<-To_export[!duplicated(To_export$Original_TrackID),c("Original_TrackID","cluster2")]
  To_export_list<-split(To_export,To_export$cluster2)
  ### Save this list that allows to identify in the imaging dataset to which cluster does each cell belong to.
  backproject_dir=paste0(output_dir,"backprojection/",to_backproject)
  dir.create(backproject_dir, recursive=TRUE)
  write(paste(as.character(To_export_list), sep="' '", collapse=", "), paste0(backproject_dir,"/Backproject.txt"))
}
