# Step 01: download TCGA data and save a SummarizedExperiment.
source("R/utils.R", local = TRUE)
args <- parse_args()
cfg <- load_cfg(args$config)

source("R/step_tcga_download.R", local = TRUE)
log_info("Running step 01 (config=%s)", args$config)
tcga_download(cfg)
