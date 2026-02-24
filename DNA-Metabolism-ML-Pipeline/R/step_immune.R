source("R/utils.R", local = TRUE)

immune_scoring <- function(cfg) {
  require_pkg("estimate")
  require_pkg("SummarizedExperiment")

  se_path <- file.path(cfg$outputs$tcga_dir, "tcga_se.rds")
  if (!file.exists(se_path)) stopf("Missing input: %s (run 01_tcga_download first)", se_path)

  ensure_dir(cfg$outputs$immune_dir)
  se <- readRDS(se_path)
  expr <- SummarizedExperiment::assay(se)

  # ESTIMATE expects genes as rownames and samples as colnames; write to file
  gct_path <- file.path(cfg$outputs$immune_dir, "expr.gct")
  estimate::write.gct(expr, gct_path)
  filt_path <- file.path(cfg$outputs$immune_dir, "expr_filtered.gct")
  estimate::filterCommonGenes(input.f = gct_path, output.f = filt_path, id = "GeneSymbol")

  scores_path <- file.path(cfg$outputs$immune_dir, "estimate_scores.gct")
  estimate::estimateScore(filt_path, scores_path, platform = "illumina")

  # Convert to CSV for convenience
  scores <- estimate::read.gct(scores_path)
  out_csv <- file.path(cfg$outputs$immune_dir, "estimate_scores.csv")
  safe_write_csv(scores, out_csv)
  log_info("Saved ESTIMATE scores: %s", out_csv)

  # Optional ssGSEA if a GMT is provided
  if (isTRUE(cfg$parameters$ssgsea) && file.exists(cfg$inputs$immune_gmt)) {
    require_pkg("GSVA")
    require_pkg("GSEABase")
    gmt <- GSEABase::getGmt(cfg$inputs$immune_gmt)
    ssg <- GSVA::gsva(as.matrix(expr), gmt, method = "ssgsea", kcdf = "Poisson", ssgsea.norm = TRUE)
    out_ssg <- file.path(cfg$outputs$immune_dir, "ssgsea_scores.csv")
    safe_write_csv(as.data.frame(t(ssg)), out_ssg)
    log_info("Saved ssGSEA scores: %s", out_ssg)
  } else {
    log_info("ssGSEA skipped (set parameters.ssgsea=true and provide inputs.immune_gmt).")
  }
}
