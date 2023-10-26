library(yaml)


# Function to set up demo configuration
setup_demo_config <- function(demo_name) {
  demo_config_path <- paste0(BEHAV3D_dir, "/demos/", demo_name, "_demo/BEHAV3D_config.yml")
  demo_config <- yaml.load_file(demo_config_path)
  
  demo_config$data_dir <- paste0(BEHAV3D_dir, "/demos/", demo_name, "_demo/data")
  demo_config$output_dir <- paste0(BEHAV3D_dir, "/demos/", demo_name, "_demo/example_output")
  demo_config$metadata_csv_path <- paste0(BEHAV3D_dir, "/demos/", demo_name, "_demo/BEHAV3D_metadata.tsv")
  demo_config$randomforest_path <- paste0(BEHAV3D_dir, "/references/TrainedRandomForest.Rdata")
  
  write_yaml(demo_config, file = demo_config_path)
}


if (interactive()) {
  BEHAV3D_dir = paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),"/")
} else {
  BEHAV3D_dir = paste0(dirname(dirname(normalizePath(sub("--file=",
                                                         "",
                                                        commandArgs(trailingOnly = FALSE)[grep("--file=",
                                                        commandArgs(trailingOnly = FALSE))])))),
                                                        "/")
}

# Set up Tcell demo
setup_demo_config("tcell")

# Set up Organoid demo
setup_demo_config("organoid")

