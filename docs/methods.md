(cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a//dev/null b/docs/methods.md
index 0000000000000000000000000000000000000000..b7c8bef6c3c849629934d8ea1326007709b36316 100644
--- a//dev/null
+++ b/docs/methods.md
@@ -0,0 +1,8 @@
+# Methods
+
+This analysis reproduces the IOPD microarray workflow described in the manuscript.
+Data were retrieved from GEO accession [GSE38680](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38680).
+Affymetrix U133 Plus 2.0 CEL files are RMA normalised using the **affy** package.
+Differential expression is assessed with **limma** and Benjamini–Hochberg correction.
+Enrichment analyses rely on **clusterProfiler** with GO and KEGG databases.
+Functional enrichment for GO and KEGG terms is performed with `clusterProfiler` and `org.Hs.eg.db`.
+miRNA targets are obtained via `multiMiR`; PPI networks are constructed with `STRINGdb`; drug–gene interactions are queried with `rDGIdb`.
+No personal health information is included in this repository.
 
EOF
)
