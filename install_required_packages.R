packages <- c(
  "abind",
  "dplyr",
  "dtwclust",
  "ggplot2",
  "gplots",
  "MESS",
  "optparse",
  "parallel",
  "pheatmap",
  "plyr",
  "randomForest",
  "readr",
  "reshape2",
  "scales",
  "Seurat",
  "spatstat",
  "sp",
  "stats",
  "tibble",
  "tidyr",
  "umap",
  "viridis",
  "xlsx",
  "yaml",
  "zoo"
  )

packagecheck <- match( packages, utils::installed.packages()[,1] )
packagestoinstall <- packages[ is.na( packagecheck ) ]
  
if(length( packagestoinstall ) > 0) {
  utils::install.packages(packagestoinstall, dependencies = TRUE, repos = "https://cloud.r-project.org/")
} else {
  print( "All requested packages already installed" )
}