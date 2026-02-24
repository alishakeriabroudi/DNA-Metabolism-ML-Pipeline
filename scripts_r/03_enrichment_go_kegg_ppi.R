suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(STRINGdb)
})

dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)

deg <- read.csv("results/deg/edger_deg.csv")
genes <- unique(deg$gene)

conv <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez <- unique(conv$ENTREZID)

ego <- enrichGO(gene=entrez, OrgDb=org.Hs.eg.db, keyType="ENTREZID",
                ont="BP", pAdjustMethod="BH", qvalueCutoff=0.05, readable=TRUE)
ekk <- enrichKEGG(gene=entrez, organism="hsa", pAdjustMethod="BH", qvalueCutoff=0.05)

write.csv(as.data.frame(ego), "results/enrichment/go_bp.csv", row.names=FALSE)
write.csv(as.data.frame(ekk), "results/enrichment/kegg.csv", row.names=FALSE)

string_db <- STRINGdb$new(version="11.5", species=9606, score_threshold=400, input_directory="")
mapped <- string_db$map(data.frame(gene=genes), "gene", removeUnmappedRows=TRUE)
ppi <- string_db$get_interactions(mapped$STRING_id)
write.csv(ppi, "results/enrichment/string_ppi_edges.csv", row.names=FALSE)

message("Saved enrichment results and STRING PPI edges.")
