library(yaml)
if (interactive()) {
  BEHAV3D_dir = paste0(dirname(dirname(rstudioapi::getSourceEditorContext()$path)),"/")
} else {
  BEHAV3D_dir = paste0(dirname(dirname(normalizePath(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])))), "/")
}

### Set up Tcell demo
tcell_demo_config_path=paste0(BEHAV3D_dir, "/demos/tcell_demo/BEHAV3D_config.yml")
tcell_demo_config = yaml.load_file(tcell_demo_config_path)

tcell_demo_config$data_dir = paste0(BEHAV3D_dir, "/demos/tcell_demo/data")
tcell_demo_config$output_dir = paste0(BEHAV3D_dir, "/demos/tcell_demo/example_output")
tcell_demo_config$metadata_csv_path = paste0(BEHAV3D_dir, "/demos/tcell_demo/BEHAV3D_metadata.tsv")
tcell_demo_config$randomforest_path = paste0(BEHAV3D_dir, "/references/TrainedRandomForest.Rdata")

write_yaml(tcell_demo_config, file=tcell_demo_config_path)

### Set up Organoid demo
organoid_demo_config_path=paste0(BEHAV3D_dir, "/demos/organoid_demo/BEHAV3D_config.yml")
organoid_demo_config = yaml.load_file(organoid_demo_config_path)

organoid_demo_config$data_dir = paste0(BEHAV3D_dir, "/demos/organoid_demo/data")
organoid_demo_config$output_dir = paste0(BEHAV3D_dir, "/demos/organoid_demo/example_output")
organoid_demo_config$metadata_csv_path = paste0(BEHAV3D_dir, "/demos/organoid_demo/BEHAV3D_metadata.tsv")
organoid_demo_config$randomforest_path = paste0(BEHAV3D_dir, "/references/TrainedRandomForest.Rdata")

write_yaml(organoid_demo_config, file=organoid_demo_config_path)
