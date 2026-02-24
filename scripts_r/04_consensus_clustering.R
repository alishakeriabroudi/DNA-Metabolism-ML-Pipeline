suppressPackageStartupMessages({
  library(ConsensusClusterPlus); library(edgeR); library(yaml)
})

cfg <- yaml::read_yaml("config/config.yaml")
counts <- readRDS("data/tcga/counts.rds")
logcpm <- cpm(counts, log=TRUE, prior.count=1)

dir.create("results/clustering", recursive = TRUE, showWarnings = FALSE)

gene_list_path <- "data/gene_sets/dna_metabolism_genes.txt"
if (file.exists(gene_list_path)) {
  genes <- readLines(gene_list_path)
  logcpm <- logcpm[rownames(logcpm) %in% genes, , drop=FALSE]
}

res <- ConsensusClusterPlus(as.matrix(logcpm),
                            maxK=cfg$clustering$max_k,
                            reps=cfg$clustering$reps,
                            pItem=cfg$clustering$p_item,
                            pFeature=cfg$clustering$p_feature,
                            clusterAlg="hc",
                            distance="pearson",
                            seed=cfg$project$seed,
                            plot="pdf",
                            title="results/clustering/consensus")

k <- cfg$clustering$k_selected
subtype <- res[[k]]$consensusClass
saveRDS(subtype, "results/clustering/subtype_C1C2.rds")
write.csv(data.frame(sample=names(subtype), subtype=subtype),
          "results/clustering/subtype_C1C2.csv", row.names=FALSE)

message(sprintf("Saved subtype labels at K=%d", k))
