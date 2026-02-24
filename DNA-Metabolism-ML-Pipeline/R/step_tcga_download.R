source("R/utils.R", local = TRUE)

tcga_download <- function(cfg) {
  require_pkg("TCGAbiolinks")

  ensure_dir(cfg$outputs$tcga_dir)

  log_info("TCGA download: project=%s, workflow=%s, sample_type=%s",
           cfg$project$cancer_project, cfg$project$rna_seq_workflow, cfg$project$sample_type)

  query <- TCGAbiolinks::GDCquery(
    project = cfg$project$cancer_project,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = cfg$project$rna_seq_workflow,
    sample.type = cfg$project$sample_type
  )

  TCGAbiolinks::GDCdownload(query, directory = cfg$outputs$tcga_dir)
  se <- TCGAbiolinks::GDCprepare(query, directory = cfg$outputs$tcga_dir)

  saveRDS(se, file = file.path(cfg$outputs$tcga_dir, "tcga_se.rds"))
  log_info("Saved: %s", file.path(cfg$outputs$tcga_dir, "tcga_se.rds"))
}
