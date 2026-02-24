#!/usr/bin/env bash
set -euo pipefail

CONFIG_PATH="${1:-config/config.yaml}"

echo "[pipeline] Using config: ${CONFIG_PATH}"
Rscript scripts_r/01_tcga_download.R --config "${CONFIG_PATH}"
Rscript scripts_r/02_deg_edger.R --config "${CONFIG_PATH}"
Rscript scripts_r/03_enrichment_go_kegg_ppi.R --config "${CONFIG_PATH}"
Rscript scripts_r/04_consensus_clustering.R --config "${CONFIG_PATH}"
Rscript scripts_r/05_immune_estimate_ssgsea.R --config "${CONFIG_PATH}"

echo "[pipeline] Done."
