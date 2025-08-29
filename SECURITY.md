 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a//dev/null b/SECURITY.md
index 0000000000000000000000000000000000000000..396d3dcc920c99c5ce291a86a6f4a8e5bf4d3a4e 100644
--- a//dev/null
+++ b/SECURITY.md
@@ -0,0 +1,4 @@
+# Security Policy
+
+- The dataset used in this project (GSE38680) is publicly available and contains no personal health information.
+- If you discover a security vulnerability or data privacy issue, please open an issue or contact the maintainers. We will address it promptly.
 
If you discover a security issue, please contact us at
[fatmatosun987@gmail.com](mailto:fatmatosun987@gmail.com).

We do not currently participate in a bug bounty program.

EOF
)
