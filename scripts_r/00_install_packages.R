if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor
bioc_pkgs <- c(
  "TCGAbiolinks","SummarizedExperiment","edgeR","limma",
  "clusterProfiler","org.Hs.eg.db",
  "GSVA","estimate","STRINGdb",
  "impute"   # required by some ML packages
)

# CRAN survival-ML ecosystem
cran_pkgs <- c(
  "dplyr","data.table","ggplot2","yaml",
  "survival","survminer","MASS",
  "glmnet","CoxBoost","randomForestSRC","gbm",
  "plsRcox","survivalsvm","superpc",
  "timeROC","Rtsne"
)

for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
for (p in bioc_pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)

message("Done.")
