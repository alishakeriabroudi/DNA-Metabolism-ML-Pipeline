# Simple smoke test: verifies key outputs exist after running the core pipeline.
stopifnot(file.exists("data/tcga/counts.rds"))
stopifnot(file.exists("results/deg/edger_all.csv"))
stopifnot(file.exists("results/enrichment/go_bp.csv"))
stopifnot(file.exists("results/clustering/subtype_C1C2.csv"))
stopifnot(file.exists("results/immune/estimate_scores.csv"))
message("Smoke test passed.")
