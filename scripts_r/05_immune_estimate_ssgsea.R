suppressPackageStartupMessages({
  library(edgeR); library(GSVA); library(estimate)
})

dir.create("results/immune", recursive = TRUE, showWarnings = FALSE)
counts <- readRDS("data/tcga/counts.rds")
logcpm <- cpm(counts, log=TRUE, prior.count=1)

tmp_gct <- "results/immune/expr.gct"
tmp_common <- "results/immune/expr_common.gct"
tmp_scores <- "results/immune/estimate_scores.gct"

expr <- as.data.frame(logcpm)
expr <- cbind(NAME=rownames(expr), Description=rownames(expr), expr)
write.table(expr, file=tmp_gct, sep="\t", quote=FALSE, row.names=FALSE)

filterCommonGenes(input.f=tmp_gct, output.f=tmp_common, id="GeneSymbol")
estimateScore(input.ds=tmp_common, output.ds=tmp_scores, platform="illumina")

scores <- read.table(tmp_scores, skip=2, header=TRUE, sep="\t", check.names=FALSE)
write.csv(scores, "results/immune/estimate_scores.csv", row.names=FALSE)

gmt_path <- "data/gene_sets/immune_signatures.gmt"
if (file.exists(gmt_path)) {
  gene_sets <- getGmt(gmt_path)
  ssgsea <- gsva(as.matrix(logcpm), gene_sets, method="ssgsea", kcdf="Gaussian", abs.ranking=TRUE)
  write.csv(t(ssgsea), "results/immune/ssgsea_scores.csv", row.names=TRUE)
}

message("Immune analysis complete.")
