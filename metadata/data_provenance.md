 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a//dev/null b/metadata/data_provenance.md
index 0000000000000000000000000000000000000000..7368409a08251df5db8db50452a95644b33ca06d 100644
--- a//dev/null
+++ b/metadata/data_provenance.md
@@ -0,0 +1,6 @@
+# Data Provenance
+
+- **Dataset:** GSE38680 (Infantile-Onset Pompe Disease)
+- **Source:** [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38680)
+- **Accessed:** 2025-01-01
+- **Notes:** Publicly available expression data; no personal health information.
+- **Storage**: Raw CEL files and `pheno.csv` tracked via Git LFS.
 
EOF
)
