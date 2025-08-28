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
This repository hosts an end-to-end R pipeline for analysing Pompe disease microarray data.
The workflow covers normalisation, differential expression, functional enrichment, miRNA
integration, **protein–protein interaction (PPI) networks**, and **drug discovery based on
miRNA targets**.


## Repository layout
- `microarray_analysis_pompe.R` – main R script containing the full analysis workflow.
- `metadata/` – large phenotype and CEL metadata tracked via Git LFS.
- `results/` – directory created at run time to store QC reports, DEG lists, plots and
  downstream analyses.

## Getting the data
Large metadata files are versioned with [Git Large File Storage](https://git-lfs.com/).
After cloning the repo run:

```bash
git lfs install
git lfs pull
```

The analysis script expects a phenotype file `pheno.csv` under `metadata/`.
If the file is not present locally the script will automatically download it to your
`~/Documents` directory (see code for the URL placeholder).

## Running the pipeline
1. Install R (≥4.3) and the required packages listed inside the script
   (`affy`, `limma`, `STRINGdb`, `rDGIdb`, etc.).
2. Place CEL files and `pheno.csv` in `metadata/`.
3. From R/RStudio run:

```r
source("microarray_analysis_pompe.R")
```

The script performs:
- RMA normalisation and quality control plots
- Differential expression analysis with `limma`
- Functional enrichment (GO/KEGG)
- **Sensitivity analysis excluding the MELAS sample**
- miRNA target retrieval via `multiMiR`
- **PPI network construction with STRINGdb**
- **Drug–gene interaction search using rDGIdb for miRNA targets**

Figures resembling those in the referenced manuscript are generated throughout
(`results/plots/`, `results/ppi/`, `results/miRNA_analysis/`).

## Notes
- Only one individual is flagged as MELAS in `pheno.csv`; sensitivity analyses automatically
  exclude this sample.
- Ensure Git LFS is set up before committing or pulling large metadata files.
  
## Citation
If you use this pipeline, please cite:

> Tosun F. & Bayrak H. "Transcriptomic and Functional Analysis of Pompe Disease" (2025).

## License
MIT
