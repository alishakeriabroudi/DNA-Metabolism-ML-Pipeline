# Step 04: consensus clustering (ConsensusClusterPlus).
source("R/utils.R", local = TRUE)
args <- parse_args()
cfg <- load_cfg(args$config)

source("R/step_clustering.R", local = TRUE)
log_info("Running step 04 (config=%s)", args$config)
consensus_clustering(cfg)
