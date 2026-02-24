suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(yaml)
})

cfg <- yaml::read_yaml("config/config.yaml")
proj <- cfg$tcga$project
workflow <- cfg$tcga$workflow

dir.create("data/tcga", recursive = TRUE, showWarnings = FALSE)

query <- GDCquery(project = proj,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = workflow)

message("Downloading TCGA data...")
GDCdownload(query)
se <- GDCprepare(query)

saveRDS(assay(se), "data/tcga/counts.rds")
saveRDS(as.data.frame(colData(se)), "data/tcga/meta.rds")
saveRDS(GDCquery_clinic(project = proj, type = "clinical"), "data/tcga/clinical.rds")

message("Saved TCGA counts/meta/clinical in data/tcga/.")
