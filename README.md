# Pompe Transcriptomics – GSE38680

This repository hosts a reproducible R workflow for analysing the
Infantile-Onset Pompe Disease (IOPD) dataset **GSE38680** generated on the
Affymetrix U133 Plus 2.0 platform. The analysis covers normalisation,
differential expression, functional enrichment, miRNA integration and PPI
network exploration.

## Quick start

```bash
# 1) set up environment
conda env create -f environment/environment.yml
conda activate pompe

# 2) obtain data
# CEL files and `pheno.csv` are tracked via Git LFS under `metadata/`
git lfs install
git lfs pull

# 3) run analysis
Rscript microarray_analysis_pompe.R
```

Results and plots are written to the `results/` directory.

## Project structure
- `microarray_analysis_pompe.R` – main R script
- `config/` – `samplesheet.csv` and `params.yaml`
- `metadata/` – CEL files and phenotype table (LFS)
- `environment/` – Conda environment specification
- `docs/` – additional documentation such as methods

## Data and ethics
All input data originate from the public GEO accession **GSE38680**.
No personal health information is included.

  ## Data via Git LFS
After cloning:
  git lfs install
  git lfs pull
  
## Citation
Please cite:

>Tosun F., Bayrak H. 2025. *Transcriptomic and Functional Analysis of Pompe Disease*.

See `CITATION.cff` for machine‑readable metadata.
