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

1. **Clone the repository**
   ```bash
   git clone <repo-url>
   cd microarray-analysis-Pompe
   ```
2. **Install Git LFS and fetch metadata**
   ```bash
   git lfs install
   git lfs pull --include="pheno.csv"
   ```
   The large `pheno.csv` will be downloaded to your `~/Documents` folder.
3. **Prepare CEL files**
   Place all CEL files referenced in `pheno.csv` in the same `~/Documents` directory.
4. **Install R dependencies** (R ≥4.5.1)
   ```r
   packages <- c("affy","limma","Biobase","hgu133plus2.db","clusterProfiler",
                 "org.Hs.eg.db","enrichplot","DOSE","GOSemSim","pheatmap",
                 "multiMiR","STRINGdb","rDGIdb","ggplot2")
   BiocManager::install(packages)
   ```
5. **Run the analysis**
   ```bash
   Rscript microarray_analysis_pompe.R
   ```
   Results and figures will be written to the `results/` folder.

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
