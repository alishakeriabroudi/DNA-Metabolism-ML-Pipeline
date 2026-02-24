source("R/utils.R", local = TRUE)

deg_edger <- function(cfg) {
  require_pkg("edgeR")
  require_pkg("SummarizedExperiment")

  se_path <- file.path(cfg$outputs$tcga_dir, "tcga_se.rds")
  if (!file.exists(se_path)) stopf("Missing input: %s (run 01_tcga_download first)", se_path)

  ensure_dir(cfg$outputs$deg_dir)
  se <- readRDS(se_path)

  counts <- SummarizedExperiment::assay(se)
  meta <- as.data.frame(SummarizedExperiment::colData(se))

  # Minimal design: if no group column exists, create a single group (safe default)
  group <- meta$shortLetterCode %||% rep("TUMOR", ncol(counts))
  group <- as.factor(group)

  y <- edgeR::DGEList(counts = counts, group = group)
  keep <- edgeR::filterByExpr(y)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y)

  design <- stats::model.matrix(~ group)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit, coef = 2)

  tab <- edgeR::topTags(qlf, n = Inf)$table
  tab$gene <- rownames(tab)
  tab <- tab[, c("gene", setdiff(colnames(tab), "gene"))]

  out_path <- file.path(cfg$outputs$deg_dir, "deg_edger_all.csv")
  safe_write_csv(tab, out_path)
  log_info("Saved DEG table: %s", out_path)

  # filtered
  fdr <- cfg$parameters$fdr_cutoff %||% 0.05
  lfc <- cfg$parameters$logfc_cutoff %||% 1.0
  tab_sig <- subset(tab, FDR <= fdr & abs(logFC) >= lfc)
  out_sig <- file.path(cfg$outputs$deg_dir, "deg_edger_sig.csv")
  safe_write_csv(tab_sig, out_sig)
  log_info("Saved filtered DEG table: %s (n=%d)", out_sig, nrow(tab_sig))
}
