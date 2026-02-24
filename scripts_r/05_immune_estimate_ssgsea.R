# Step 05: immune scoring (ESTIMATE + optional ssGSEA).
source("R/utils.R", local = TRUE)
args <- parse_args()
cfg <- load_cfg(args$config)

source("R/step_immune.R", local = TRUE)
log_info("Running step 05 (config=%s)", args$config)
immune_scoring(cfg)
