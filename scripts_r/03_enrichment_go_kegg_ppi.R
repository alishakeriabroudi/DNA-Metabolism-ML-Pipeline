# Step 03: GO/KEGG enrichment (+ offline PPI note).
source("R/utils.R", local = TRUE)
args <- parse_args()
cfg <- load_cfg(args$config)

source("R/step_enrichment.R", local = TRUE)
log_info("Running step 03 (config=%s)", args$config)
enrichment_go_kegg_ppi(cfg)
