 (cd "$(git rev-parse --show-toplevel)" && git apply --3way <<'EOF' 
diff --git a//dev/null b/CONTRIBUTING.md
index 0000000000000000000000000000000000000000..bcbde924d409cc504445ed2c97b4e702a3fc28ad 100644
--- a//dev/null
+++ b/CONTRIBUTING.md
@@ -0,0 +1,22 @@
+# Contributing
+
+Thanks for your interest in improving the Pompe microarray analysis pipeline!
+
+## Getting Started
+- Fork the repository and create your branch from `main`.
+- Ensure Git LFS is installed (`git lfs install`).
+- Run the analysis locally to confirm changes: `Rscript microarray_analysis_pompe.R`.
+
+## Coding Standards
+- Use R style compatible with `styler` and `lintr`.
+- Place data under `metadata/` and results under `results/`.
+- Update documentation and examples when changing interfaces.
+
+## Pull Requests
+- Reference related issues in commit messages.
+- Include a summary of changes and testing steps.
+- Ensure CI workflows pass.
+
+## Reporting Issues
+Please open an issue using the templates in `.github/ISSUE_TEMPLATE` and provide
+as much detail as possible.
 
EOF
)
