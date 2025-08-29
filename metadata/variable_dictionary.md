 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a//dev/null b/metadata/variable_dictionary.md
index 0000000000000000000000000000000000000000..ef49ed5c47858de2354f0032fcd13184b35daeab 100644
--- a//dev/null
+++ b/metadata/variable_dictionary.md
@@ -0,0 +1,10 @@
+# Variable Dictionary
+
+| column    | description                                                    |
+|-----------|----------------------------------------------------------------|
+| sample_id | GEO sample accession (e.g., GSM947461)                         |
+| filename  | CEL file name associated with the sample                       |
+| group     | Phenotype group: IOPD or Control                               |
+| age       | Age in months                                                  |
+| sex       | Biological sex (M/F)                                           |
+| sample   | Unique sample identifier (e.g., GEO accession)                  |
+| is_melas | `TRUE` if sample is flagged for MELAS sensitivity analysis      |
 
EOF
)
