# Step 02: differential expression using edgeR (safe default design).
source("R/utils.R", local = TRUE)
args <- parse_args()
cfg <- load_cfg(args$config)

source("R/step_deg_edger.R", local = TRUE)
log_info("Running step 02 (config=%s)", args$config)
deg_edger(cfg)
