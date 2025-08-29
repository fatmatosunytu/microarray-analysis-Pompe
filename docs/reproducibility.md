 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a//dev/null b/docs/reproducibility.md
index 0000000000000000000000000000000000000000..947af31e6c26a1af3e06a65df68e8c77f9739253 100644
--- a//dev/null
+++ b/docs/reproducibility.md
@@ -0,0 +1,18 @@
+# Reproducibility
+
+1. Install the software environment:
+   ```bash
+   conda env create -f environment/environment.yml
+   conda activate pompe
+   ```
+2. Pull required metadata tracked with Git LFS:
+   ```bash
+   git lfs install
+   git lfs pull
+   ```
+3. Ensure `config/samplesheet.csv` and `metadata/pheno.csv` reflect your setup.
+4. Run the analysis script:
+   ```bash
+   Rscript microarray_analysis_pompe.R
+   ```
+Results are written to the `results/` directory.
 
EOF
)
