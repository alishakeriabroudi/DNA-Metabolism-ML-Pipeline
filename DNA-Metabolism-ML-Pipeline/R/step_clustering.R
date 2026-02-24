source("R/utils.R", local = TRUE)

consensus_clustering <- function(cfg) {
  require_pkg("ConsensusClusterPlus")
  require_pkg("SummarizedExperiment")

  se_path <- file.path(cfg$outputs$tcga_dir, "tcga_se.rds")
  if (!file.exists(se_path)) stopf("Missing input: %s (run 01_tcga_download first)", se_path)

  ensure_dir(cfg$outputs$clustering_dir)

  se <- readRDS(se_path)
  mat <- SummarizedExperiment::assay(se)

  # Optional: restrict to gene set
  genes <- read_gene_list(cfg$inputs$dna_metabolism_genes)
  if (!is.null(genes)) {
    keep <- intersect(rownames(mat), genes)
    if (length(keep) >= 10) {
      mat <- mat[keep, , drop = FALSE]
      log_info("Clustering restricted to gene set (n=%d).", nrow(mat))
    } else {
      log_info("Gene set provided but too small after intersection; using all genes.")
    }
  }

  # variance filter to keep top features
  vars <- apply(mat, 1, stats::var)
  top_n <- min(2000, length(vars))
  mat <- mat[order(vars, decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]

  max_k <- cfg$parameters$max_k %||% 6
  reps <- cfg$parameters$reps %||% 100
  p_item <- cfg$parameters$p_item %||% 0.8
  p_feature <- cfg$parameters$p_feature %||% 0.8
  alg <- cfg$parameters$cluster_alg %||% "hc"
  dist <- cfg$parameters$distance %||% "pearson"

  log_info("Consensus clustering: maxK=%d reps=%d pItem=%.2f pFeature=%.2f alg=%s dist=%s",
           max_k, reps, p_item, p_feature, alg, dist)

  res <- ConsensusClusterPlus::ConsensusClusterPlus(
    as.matrix(mat),
    maxK = max_k,
    reps = reps,
    pItem = p_item,
    pFeature = p_feature,
    clusterAlg = alg,
    distance = dist,
    seed = 123,
    plot = "pdf",
    title = cfg$outputs$clustering_dir
  )

  # Save a simple label file using bestK = 2 (safe default)
  best_k <- min(2, max_k)
  cl <- res[[best_k]]$consensusClass
  labels <- data.frame(sample = names(cl), subtype = paste0("C", as.integer(cl)))
  out_path <- file.path(cfg$outputs$clustering_dir, sprintf("subtypes_k%d.csv", best_k))
  safe_write_csv(labels, out_path)
  log_info("Saved subtype labels: %s", out_path)
}
