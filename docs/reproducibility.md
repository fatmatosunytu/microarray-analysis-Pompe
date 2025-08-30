# Reproducibility

1. Install the software environment:
   ```bash
   conda env create -f environment/environment.yml
   conda activate pompe
   ```
2. Pull required metadata tracked with Git LFS:
   ```bash
   git lfs install
   git lfs pull
   ```
3. Ensure `config/samplesheet.csv` and `metadata/pheno.csv` reflect your setup.
4. Run the analysis script:
   ```bash
   Rscript microarray_analysis_pompe.R
   ```
Results are written to the `results/` directory.
