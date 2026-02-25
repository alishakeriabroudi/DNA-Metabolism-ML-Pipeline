# DNA-Metabolism-ML-Pipeline

A reproducible **TCGA prostate cancer** transcriptomics pipeline (bulk RNA‑seq) with:
- TCGA download (TCGAbiolinks)
- Differential expression (edgeR)
- GO/KEGG enrichment (clusterProfiler)
- STRING PPI export (for Cytoscape)
- Consensus clustering (ConsensusClusterPlus)
- Immune scoring (ESTIMATE; optional ssGSEA if you provide a GMT)

This repo is organized to be **reviewable**, **config‑driven**, and easy to run on a clean machine.

## Quickstart

### 1) R packages (one time)
Open R and run:

```r
source("scripts_r/00_install_packages.R")
```

### 2) Run the core pipeline (config‑driven)
```bash
Rscript scripts_r/01_tcga_download.R --config config/config.yaml
Rscript scripts_r/02_deg_edger.R --config config/config.yaml
Rscript scripts_r/03_enrichment_go_kegg_ppi.R --config config/config.yaml
Rscript scripts_r/04_consensus_clustering.R --config config/config.yaml
Rscript scripts_r/05_immune_estimate_ssgsea.R --config config/config.yaml
```

Or run everything:
```bash
bash run_core_pipeline.sh
```

## Inputs (optional)
- `data/gene_sets/dna_metabolism_genes.txt` (one gene symbol per line) — restrict clustering to a gene set
- `data/gene_sets/immune_signatures.gmt` — enable ssGSEA immune scoring

If missing, scripts fall back to safe defaults.

## Outputs
All outputs go to `results/` (see subfolders inside `config/config.yaml`).

## Notes
- Scripts are **thin wrappers** that call functions in `R/` (shared utilities + pipeline steps).
- The pipeline is designed to be robust: clear errors, deterministic folder structure, and minimal hard‑coding.

## Associated publication
- **Using machine learning to discover DNA metabolism biomarkers that direct prostate cancer treatment**  
- Ali Shakeri Abroudi, Melika Djamali, Hossein Azizi — *Scientific Reports* 15, 26117 (2025)  
- DOI: https://doi.org/10.1038/s41598-025-11457-1

## License
MIT (see `LICENSE`).
