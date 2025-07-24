# microarray-analysis

# Pompe Disease Microarray Transcriptomics Analysis

A comprehensive R pipeline for normalization, differential expression analysis, annotation, functional enrichment (GO/KEGG), and miRNA target analysis in Pompe disease microarray data.

## Features

- Robust pre-processing and normalization
- Limma-based DEG analysis
- Annotation via biomaRt or annotation packages
- GO/KEGG pathway enrichment with clusterProfiler
- miRNA target analysis with multiMiR
- Publication-ready figures

## Usage

1. Place your CEL files in the `dataset` folder.
2. Run `main_analysis.R` in R (you may need to set working directory to repo root).
3. Results and figures will be saved in the `output/` folder.

## Requirements

- R 4.5.1
- Bioconductor
- Packages: affy, limma, Biobase, hgu133plus2.db, biomaRt, clusterProfiler, org.Hs.eg.db, enrichplot, DOSE, GOSemSim, pheatmap, WebGestaltR, multiMiR, ggplot2

## Citation

If you use this pipeline, please cite:

> Tosun F. & Bayrak H. "Transcriptomic and Functional Analysis of Pompe Disease" (2025). [GitHub link]

## License

MIT
