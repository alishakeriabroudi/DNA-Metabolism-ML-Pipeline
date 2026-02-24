#!/usr/bin/env bash
set -euo pipefail

echo "[1/5] TCGA download"
Rscript scripts_r/01_tcga_download.R

echo "[2/5] DEG (edgeR)"
Rscript scripts_r/02_deg_edger.R

echo "[3/5] Enrichment + PPI"
Rscript scripts_r/03_enrichment_go_kegg_ppi.R

echo "[4/5] Consensus clustering"
Rscript scripts_r/04_consensus_clustering.R

echo "[5/5] Immune scoring"
Rscript scripts_r/05_immune_estimate_ssgsea.R

echo "Core pipeline completed."
