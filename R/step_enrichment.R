source("R/utils.R", local = TRUE)

enrichment_go_kegg_ppi <- function(cfg) {
  require_pkg("clusterProfiler")
  require_pkg("org.Hs.eg.db")

  deg_sig <- file.path(cfg$outputs$deg_dir, "deg_edger_sig.csv")
  if (!file.exists(deg_sig)) stopf("Missing input: %s (run 02_deg_edger first)", deg_sig)

  ensure_dir(cfg$outputs$enrichment_dir)
  tab <- utils::read.csv(deg_sig, stringsAsFactors = FALSE)
  genes <- unique(tab$gene)
  if (length(genes) < 5) stopf("Not enough DE genes for enrichment (n=%d).", length(genes))

  eg <- clusterProfiler::bitr(
    genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db
  )
  if (nrow(eg) == 0) stopf("No genes mapped to ENTREZID.")

  ego <- clusterProfiler::enrichGO(
    gene = unique(eg$ENTREZID),
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  ekegg <- clusterProfiler::enrichKEGG(
    gene = unique(eg$ENTREZID),
    organism = "hsa",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2
  )

  go_path <- file.path(cfg$outputs$enrichment_dir, "enrich_go_bp.csv")
  kegg_path <- file.path(cfg$outputs$enrichment_dir, "enrich_kegg.csv")
  safe_write_csv(as.data.frame(ego), go_path)
  safe_write_csv(as.data.frame(ekegg), kegg_path)

  log_info("Saved GO: %s", go_path)
  log_info("Saved KEGG: %s", kegg_path)

  # PPI export note: actual STRING API calls are intentionally avoided for reproducibility/offline runs.
  note_path <- file.path(cfg$outputs$enrichment_dir, "ppi_string_NOTE.txt")
  writeLines(c(
    "STRING PPI export:",
    "- This pipeline writes a placeholder note instead of calling STRING over the internet.",
    "- To create a Cytoscape edge list: upload DE genes to STRING (string-db.org) and export edges.",
    "- Or implement STRINGdb package usage here if you accept online dependencies."
  ), con = note_path, useBytes = TRUE)
  log_info("Wrote PPI note: %s", note_path)
}
