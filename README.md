# Pompe Disease Microarray Transcriptomics Analysis

A comprehensive R pipeline for normalization, differential expression analysis, annotation, functional enrichment (GO/KEGG), miRNA target analysis and downstream protein–protein interaction (PPI) exploration in Pompe disease microarray data. The pipeline also performs a MELAS sensitivity analysis by repeating the DEG workflow with the flagged MELAS sample removed.

## Key Features

- Robust pre-processing and RMA normalization
- Limma-based differential expression (DEG) analysis
- Functional enrichment via *clusterProfiler*
- miRNA target investigation with *multiMiR*
- PPI network construction using *STRINGdb*
- Drug–gene interaction mining through *rDGIdb*
- Publication-ready visualisations

## Usage

1. Clone the repository and install [Git LFS](https://git-lfs.com/).
2. Run `git lfs install` and `git lfs pull` to download large metadata files such as `pheno.csv`.
3. Ensure `pheno.csv` is located in your `~/Documents` directory or in `metadata/` at the repo root. The analysis script will copy it into place automatically.
4. Place your CEL files in the `metadata/` directory alongside `pheno.csv`.
5. Run `microarray_analysis_pompe.R` in R (you may need to set the environment variable `POMPE_ROOT` to the repo root).
6. Results, plots and sensitivity analyses will be written to the `results/` folder.


## Getting Started

1. Install [Git LFS](https://git-lfs.com/) and track the metadata file:
   ```bash
   git lfs install
   git lfs track "metadata/pheno.csv"
   git add .gitattributes metadata/pheno.csv
   git commit -m "Add pheno metadata"
   ```
2. The analysis script expects `pheno.csv` in your `Documents` folder. If it is
   absent, the script will download it using the URL defined inside
   `microarray_analysis_pompe.R`.
3. Place raw CEL files in `metadata/` (or update `cel_dir` in the script).
4. Run `microarray_analysis_pompe.R` in R with the repository root as the working directory:
   ```r
   source("microarray_analysis_pompe.R")
   ```
5. Results, figures, PPI networks and drug‑discovery tables will be written to
   the `results/` directory tree.

## Notes

- The script automatically detects the metadata file in `~/Documents` and performs a
  sensitivity analysis excluding the single MELAS-labelled sample.
- Ensure Git LFS is configured before running, otherwise the metadata file will be missing.

## Requirements

- R 4.5.1
- Bioconductor
- Packages: affy, limma, Biobase, hgu133plus2.db, clusterProfiler, org.Hs.eg.db, enrichplot, DOSE, GOSemSim, pheatmap, multiMiR, ggplot2, STRINGdb, ggraph, rDGIdb, tidygraph, igraph

## Citation

If you use this pipeline, please cite:

> Tosun F. & Bayrak H. "Transcriptomic and Functional Analysis of Pompe Disease" (2025). [GitHub link]

## License

MIT
