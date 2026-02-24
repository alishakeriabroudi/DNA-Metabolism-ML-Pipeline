# Companion Analysis Repository (Clean Core + Experimental)

This repository is a **clean, reproducible companion codebase** for a prostate cancer transcriptomics workflow:
- TCGA bulk RNA-seq download (TCGAbiolinks)
- DEG discovery (edgeR)
- GO/KEGG enrichment (clusterProfiler)
- STRING PPI edge export (for Cytoscape)
- Consensus clustering (ConsensusClusterPlus)
- Immune scoring (ESTIMATE; optional ssGSEA if you provide a GMT)

The **core pipeline is intentionally kept simple and robust** so you can upload it to GitHub without fragile dependencies.

Optional, more complex modeling ideas (e.g., large survival-model ensembles) are placed in **`experimental/`**.

---

## Quickstart

### 1) Create Python environment (only used for optional tools)
```bash
conda env create -f environment.yml
conda activate ml-prostate-pipeline
```

### 2) Install R packages (one time)
Open R and run:
```r
source("scripts_r/00_install_packages.R")
```

### 3) Run the core pipeline
```bash
Rscript scripts_r/01_tcga_download.R
Rscript scripts_r/02_deg_edger.R
Rscript scripts_r/03_enrichment_go_kegg_ppi.R
Rscript scripts_r/04_consensus_clustering.R
Rscript scripts_r/05_immune_estimate_ssgsea.R
```

A convenience runner is also provided:
```bash
bash run_core_pipeline.sh
```

---

## Inputs you may optionally provide

Place optional gene sets in:
- `data/gene_sets/dna_metabolism_genes.txt` (one gene symbol per line) to restrict clustering
- `data/gene_sets/immune_signatures.gmt` for ssGSEA immune scoring

If these files are missing, the scripts fall back to safe defaults.

---

## Outputs

- `results/deg/` — DEG tables
- `results/enrichment/` — GO/KEGG + STRING PPI edge list
- `results/clustering/` — consensus clustering outputs + subtype labels
- `results/immune/` — ESTIMATE scores (+ ssGSEA scores if GMT provided)

---

## Experimental

See `experimental/README_EXPERIMENTAL.md`.

---

## License

MIT (see `LICENSE`).
