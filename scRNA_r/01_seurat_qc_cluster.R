suppressPackageStartupMessages({ library(Seurat) })
dir.create("data/scrna", recursive=TRUE, showWarnings=FALSE)
dir.create("results/scrna", recursive=TRUE, showWarnings=FALSE)
stop("Edit this script to load your scRNA matrix (e.g., Read10X), then save a Seurat object to results/scrna/seurat_obj.rds")
