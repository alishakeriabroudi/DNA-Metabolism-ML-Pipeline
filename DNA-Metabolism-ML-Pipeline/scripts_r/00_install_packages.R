# One-time installer for required R packages.
# Keeps the core pipeline reproducible without forcing a heavy environment manager.

pkgs_cran <- c("yaml", "BiocManager", "data.table")
pkgs_bioc <- c(
  "TCGAbiolinks", "SummarizedExperiment", "edgeR",
  "clusterProfiler", "org.Hs.eg.db", "ConsensusClusterPlus",
  "estimate"
)
pkgs_optional <- c("GSVA", "GSEABase")

install_if_missing <- function(pkgs) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("Installing: ", p)
      install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing(pkgs_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

for (p in c(pkgs_bioc, pkgs_optional)) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Bioconductor install: ", p)
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}

message("Done. You can now run: bash run_core_pipeline.sh")
