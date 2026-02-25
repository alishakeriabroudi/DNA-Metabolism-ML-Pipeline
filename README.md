# DNA-Metabolism-ML-Pipeline

A reproducible **TCGA prostate cancer** transcriptomics pipeline (bulk RNA-seq) with:

- TCGA download (TCGAbiolinks)
- Differential expression (edgeR)
- GO/KEGG enrichment (clusterProfiler)
- STRING PPI export (for Cytoscape)
- Consensus clustering (ConsensusClusterPlus)
- Immune scoring (ESTIMATE; optional ssGSEA if you provide a GMT)

This repo is organized to be **reviewable**, **config-driven**, and easy to run on a clean machine.

---

## Associated publication

**[Using machine learning to discover DNA metabolism biomarkers that direct prostate cancer treatment](https://doi.org/10.1038/s41598-025-11457-1)**  
Ali Shakeri Abroudi, Melika Djamali, Hossein Azizi — *Scientific Reports* 15, 26117 (2025)  
DOI: https://doi.org/10.1038/s41598-025-11457-1

## Software DOI (Zenodo)

- **Cite all versions (Concept DOI):** https://doi.org/10.5281/zenodo.18772013  
- **This release (v1.0.1, Version DOI):** https://doi.org/10.5281/zenodo.18772014
  
---

## Requirements

- **R (>= 4.2)** recommended
- *(Optional but recommended)* `bash` + a Unix-like shell (macOS/Linux/WSL)
- Internet access for TCGA download

---

## Quickstart

### 1) Install R packages (one time)

Open R and run:

```r
source("scripts_r/00_install_packages.R")
```

### 2) Run the core pipeline (config-driven)

Run step-by-step:

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

---

## Configuration

Most paths and parameters live in:

- `config/config.yaml`

All outputs are written under `results/` using the folder structure defined in the config.

---

## Inputs (optional)

- `data/gene_sets/dna_metabolism_genes.txt` (one gene symbol per line)  
  Restricts clustering / downstream steps to a gene set (if enabled by config).

- `data/gene_sets/immune_signatures.gmt`  
  Enables ssGSEA immune scoring (if enabled by config).

If these files are missing, scripts fall back to safe defaults.

---

## Outputs

All outputs go to `results/` (see subfolders inside `config/config.yaml`). Typical outputs include:

- Differential expression tables (edgeR)
- GO/KEGG enrichment results
- STRING PPI exports for Cytoscape
- Consensus clustering assignments and plots
- Immune scoring outputs (ESTIMATE; optional ssGSEA)

---

## Notes

- Scripts in `scripts_r/` are **thin wrappers** that call functions in `R/` (shared utilities + pipeline steps).
- The pipeline aims to be robust: clear errors, deterministic folder structure, and minimal hard-coding.

---

## Citation

If you use this repository, please cite the associated paper above.  
A `CITATION.cff` file is recommended/provided for GitHub citation support.

---

## License

MIT (see `LICENSE`).
