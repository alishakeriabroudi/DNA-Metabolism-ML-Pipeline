suppressPackageStartupMessages({
  library(edgeR); library(dplyr); library(yaml)
})

cfg <- yaml::read_yaml("config/config.yaml")
fdr_thr <- cfg$deg$fdr
fc_thr  <- cfg$deg$fc

counts <- readRDS("data/tcga/counts.rds")
meta   <- readRDS("data/tcga/meta.rds")

dir.create("results/deg", recursive = TRUE, showWarnings = FALSE)

# Attempt to form Tumor/Normal labels from TCGA metadata (PRAD may have few normals).
if (!("shortLetterCode" %in% colnames(meta)) && !("sample_type" %in% colnames(meta))) {
  stop("No sample type columns found. Adapt this script for your data.")
}

if ("shortLetterCode" %in% colnames(meta)) {
  group <- ifelse(meta$shortLetterCode == "TP", "Tumor",
           ifelse(meta$shortLetterCode == "NT", "Normal", NA))
} else {
  group <- ifelse(grepl("Tumor", meta$sample_type, ignore.case=TRUE), "Tumor",
           ifelse(grepl("Normal", meta$sample_type, ignore.case=TRUE), "Normal", NA))
}

keep <- !is.na(group)
counts <- counts[, keep, drop=FALSE]
group <- factor(group[keep])

if (length(unique(group)) < 2) stop("Only one group present. Use GEO tumor-vs-normal DEG or provide labels.")

dge <- DGEList(counts=counts, group=group)
dge <- dge[filterByExpr(dge), , keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)

design <- model.matrix(~group)
dge <- estimateDisp(dge, design)
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef=2)

tab <- topTags(qlf, n=Inf)$table
tab$FDR <- p.adjust(tab$PValue, method="BH")

deg <- tab %>%
  mutate(gene=rownames(tab)) %>%
  filter(FDR < fdr_thr & abs(logFC) >= log2(fc_thr)) %>%
  arrange(FDR)

write.csv(tab, "results/deg/edger_all.csv", row.names=TRUE)
write.csv(deg, "results/deg/edger_deg.csv", row.names=FALSE)

message(sprintf("Saved DEG: n=%d", nrow(deg)))
