# =============================================================================
# Pompe Disease Microarray Transcriptomics Analysis Pipeline
# Author: Fatma Tosun- Harun Bayrak
# Date: 202-07-24
# Affiliation: [Your Affiliation]
# Description: Full pipeline for normalization, DEG analysis, annotation, functional enrichment (GO/KEGG), miRNA analysis, visualization
# License: MIT
# =============================================================================

# ---------------------------LOAD AND INSTALL REQUIRED PACKAGES-----------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("affy", "limma", "Biobase", "hgu133plus2.db"), ask = FALSE)

library(affy)
library(limma)
library(Biobase)
library(hgu133plus2.db)

# 1. Set the working directory and read the CEL files
setwd("C:/Users/User/Desktop/dataset")
gset <- ReadAffy()
gset <- rma(gset)  # Normalize et

# 2. Extract the normalized expression matrix
expr_matrix <- exprs(gset)

# RMA outputs are generally not needed, but if you want to add them:
if(any(is.na(expr_matrix))) {
  expr_matrix[is.na(expr_matrix)] <- rowMeans(expr_matrix, na.rm = TRUE)[row(expr_matrix)[is.na(expr_matrix)]]
}
if(any(is.infinite(expr_matrix))) {
  expr_matrix[is.infinite(expr_matrix)] <- quantile(expr_matrix, 0.99, na.rm = TRUE)
}

# Get the sample names (column names)
colnames_ex <- colnames(expr_matrix)

# Extract GSM IDs (example: GSM947461_H0013501 → GSM947461)
colnames_ex_gsm <- sapply(strsplit(colnames_ex, "_"), `[`, 1)

# Group definitions
control_gsm <- c("GSM947461", "GSM947462", "GSM947463", "GSM947464", "GSM947465", 
                 "GSM947466", "GSM947467", "GSM947468", "GSM947469", "GSM947470")
pompe_gsm <- c("GSM947471", "GSM947472", "GSM947473", "GSM947474", "GSM947475", 
               "GSM947476", "GSM947477", "GSM947478", "GSM947479")

# Create the group vector
groups <- ifelse(colnames_ex_gsm %in% control_gsm, "Control",
                 ifelse(colnames_ex_gsm %in% pompe_gsm, "Pompe", NA))

# Check for mismatches
if (any(is.na(groups))) {
  stop("Eşleşemeyen örnek(ler): ", paste(colnames_ex[is.na(groups)], collapse = ", "))
}

# Define as a factor
groups <- factor(groups, levels = c("Control", "Pompe"))

# Check the distribution
table(groups)

#----------------------------- ANALYSIS PREPARATION ----------------------------

# 1. Get and clean sample names
sample_names <- colnames(expr_matrix)
sample_names <- gsub("\\.CEL$", "", sample_names)

# 2. Define group information (required for analysis)
control_samples <- c("GSM947461_H0013501", "GSM947462_H0013502", "GSM947463_H0013504", 
                     "GSM947464_H0013505", "GSM947465_H0013506", "GSM947466_H0013507", 
                     "GSM947467_H0013508", "GSM947468_H0013509", "GSM947469_H0013510", 
                     "GSM947470_H0013503")
pompe_samples <- c("GSM947471_H0013492", "GSM947472_H0013493", "GSM947473_H0013494", 
                   "GSM947474_H0013495", "GSM947475_H0013496", "GSM947476_H0013497", 
                   "GSM947477_H0013498", "GSM947478_H0013499", "GSM947479_H0013500")

groups <- ifelse(sample_names %in% control_samples, "Control",
                 ifelse(sample_names %in% pompe_samples, "Pompe", NA))

# 3. Check for missing matches
if(any(is.na(groups))) {
  stop("Eşleşmeyen örnekler bulundu: ", paste(sample_names[is.na(groups)], collapse = ", "))
}

groups <- factor(groups, levels = c("Control", "Pompe"))

# 4. Design matrix and model
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

fit <- lmFit(expr_matrix, design)
contrast.matrix <- makeContrasts(Pompe - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 5. Get DEG results
deg_results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH", p.value = 0.05, lfc = 1)

# -----------------------------------DEG RESULTS--------------------------------

# 1. Install relevant packages
BiocManager::install(c("hgu133plus2.db", "biomaRt", "clusterProfiler", "org.Hs.eg.db", "enrichplot"))
install.packages("pheatmap")
library(hgu133plus2.db)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(pheatmap)


# 2. Annotation via biomaRt
View(deg_results)
probe_ids <- rownames(deg_results)

tryCatch({
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_symbols <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol"),
                        filters = "affy_hg_u133_plus_2",
                        values = probe_ids,
                        mart = mart)
  
  deg_results$Gene.Symbol <- gene_symbols$hgnc_symbol[match(probe_ids, gene_symbols$affy_hg_u133_plus_2)]
  
  write.csv(deg_results, "Pompe_vs_Control_Anlamli_Genler_Sembollerle.csv")
  cat("Gen sembolleri başarıyla eklendi.\n")
}, error = function(e) {
  cat("BiomaRt bağlantı hatası:", e$message, "\n")
})

# 3. Save and display results
write.csv(deg_results, "Pompe_vs_Control_limma_DEGs_with_genes.csv", row.names = TRUE)
View(deg_results)

# 4. Statistical summary
de_summary <- summary(decideTests(fit2, adjust.method = "BH", p.value = 0.05, lfc = 1))
print(de_summary)

# -----------------------------------VOLCANO PLOT-------------------------------

plotVolcano <- function(fit2, coef = 1, p.value = 0.05, lfc = 1) {
  full_results <- topTable(fit2, coef = coef, number = Inf)
  full_results$negLogP <- -log10(full_results$adj.P.Val)
  full_results$color <- "grey"
  full_results$color[full_results$adj.P.Val < p.value & full_results$logFC >  lfc] <- "red"
  full_results$color[full_results$adj.P.Val < p.value & full_results$logFC < -lfc] <- "blue"
  
  plot(full_results$logFC,
       full_results$negLogP,
       col = full_results$color,
       pch = 20, cex = 0.6,
       xlab = "log2 Fold Change",
       ylab = "-log10(adjusted P.Value)",
       main = "Volcano Plot - Pompe vs Control")
  abline(v = c(-lfc, lfc), col = "black", lty = 2)
  abline(h = -log10(p.value), col = "black", lty = 2)
  legend("topright", legend = c("Down", "Not Sig", "Up"), col = c("blue", "grey", "red"), pch = 20)
}
plotVolcano(fit2)

# -----------------------------------HEATMAP------------------------------------

# Top 50 genes for heatmap
top_genes <- rownames(deg_results)[1:50]
expr_top <- expr_matrix[top_genes, ]
annotation_col <- data.frame(Group = groups)
rownames(annotation_col) <- colnames(expr_top)

# If not loaded:
if (!require("pheatmap")) install.packages("pheatmap")
library(pheatmap)

pheatmap(expr_top,
         scale = "row",
         show_rownames = FALSE,
         annotation_col = annotation_col,
         main = "Top 50 DEG Heatmap")

# -----------------------------------GO AND KEGG ANALYSIS-----------------------

# -----------------------------------------------------
# 1. Install relevant packages
# -----------------------------------------------------
# You need the following Bioconductor and CRAN packages for this analysis:
# clusterProfiler, org.Hs.eg.db, enrichplot, ggplot2, DOSE, GOSemSim

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(GOSemSim)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

necessary_pkgs <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "ggplot2", "DOSE", "GOSemSim")
for (pkg in necessary_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# -----------------------------------------------------
# 2. Prepare the Significant Gene List
# -----------------------------------------------------
#We draw the symbols of genes found to be significant from differential expression analysis.
significant_genes <- deg_results$Gene.Symbol[!is.na(deg_results$Gene.Symbol)]
significant_genes <- unique(significant_genes)

# -----------------------------------------------------
# 3. GO Biological Process (BP) Enrichment Analysis
# -----------------------------------------------------
go_bp <- enrichGO(
  gene = significant_genes,
  OrgDb = org.Hs.eg.db,     # Human Database
  keyType = "SYMBOL",       # Gene Key: Symbol
  ont = "BP",               # BP: Biological Process ontology
  pvalueCutoff = 0.05,      # p-value threshold
  qvalueCutoff = 0.05,      # FDR corrected threshold
  readable = TRUE           # Convert gene IDs to readable names
)

#Save results as CSV (full table)
write.csv(as.data.frame(go_bp), "GO_Biological_Process.csv", row.names = FALSE)

# -----------------------------------------------------
# 4. KEGG Pathway Enrichment Analysis
# -----------------------------------------------------
# In KEGG analysis, genes must be in Entrez ID format!
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = significant_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)
entrez_ids <- na.omit(entrez_ids)

kegg <- enrichKEGG(
  gene = entrez_ids,
  organism = 'hsa',          # Homo sapiens
  pvalueCutoff = 0.05
)

write.csv(as.data.frame(kegg), "KEGG_Enrichment.csv", row.names = FALSE)
barplot(kegg, showCategory = 15, title = "KEGG Pathways")
browseKEGG(kegg, 'hsa04142')  # example: Lysosome pathway can be added to other pathwaysr


# -----------------------------------------------------
# 5. Visualization of GO Results
# -----------------------------------------------------

library(GOSemSim)

if (nrow(go_bp) > 0) {
  # Dotplot: The 10 most enriched GO terms
  p1 <- dotplot(go_bp, showCategory = 10, title = "GO Biological Process")
  ggsave("GO_BP_dotplot.png", plot = p1, width = 10, height = 8)
  
  # Barplot: Alternative visualization
  p2 <- barplot(go_bp, showCategory = 10, title = "GO Biological Process")
  ggsave("GO_BP_barplot.png", plot = p2, width = 10, height = 8)
  
  # Semantic Similarity (for emapplot): Semantic similarity matrix between terms
  hsGO <- godata('org.Hs.eg.db', ont = "BP")
  go_bp_sim <- pairwise_termsim(go_bp, semData = hsGO)
  
  # emapplot: Network of relationships between GO terms
  tryCatch({
    p3 <- emapplot(go_bp_sim, showCategory = 10)
    ggsave("GO_BP_emapplot.png", plot = p3, width = 10, height = 8)
  }, error = function(e) cat("Could not create emapplot:", e$message, "\n"))
  
  # cnetplot: Network showing Gen-GO term relationships
  tryCatch({
    p4 <- cnetplot(go_bp, categorySize = "pvalue", foldChange = NULL)
    ggsave("GO_BP_cnetplot.png", plot = p4, width = 12, height = 10)
  }, error = function(e) cat("Could not create cnetplot:", e$message, "\n"))
  
  # Custom Dotplot (with ggplot2)
  go_df <- as.data.frame(go_bp)
  top_go <- head(go_df[order(go_df$p.adjust), ], 15)
  
  p5 <- ggplot(top_go, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = "GO Biological Process (Custom)",
         x = "-log10(Adjusted P-value)", y = "GO Term") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10))
  ggsave("GO_BP_custom_dotplot.png", plot = p5, width = 12, height = 8)
  
  # Save the 10 most significant terms as a separate CSV
  write.csv(top_go[1:10, ], "Top_GO_Terms.csv", row.names = FALSE)
  
} else {
  cat("No GO enrichment results found.\n")
}

# -----------------------------------------------------
# 6. Visualizing KEGG Results
# -----------------------------------------------------
if (nrow(as.data.frame(kegg)) > 0) {
  p_kegg <- barplot(kegg, showCategory = 15, title = "KEGG Pathways")
  ggsave("KEGG_barplot.png", plot = p_kegg, width = 10, height = 8)
} else {
  cat("No KEGG enrichment results found.\n")
}

# -----------------------------------------------------
# 7. Examining the Results Table for
# -----------------------------------------------------
# If you want to open GO enrichment results as a table:
 go_df <- as.data.frame(go_bp)
 View(go_df)


 # -----------------------------------miRNA ANALYSIS----------------------------
 
 # — Completed (normalization, limma, annotation, heatmap, GO/KEGG…)
 
 # 1. NA / Infinite value check
 if (any(is.na(expr_matrix))) {
   expr_matrix[is.na(expr_matrix)] <- rowMeans(expr_matrix, na.rm = TRUE)[row(expr_matrix)[is.na(expr_matrix)]]
 }
 
 if (any(is.infinite(expr_matrix))) {
   expr_matrix[is.infinite(expr_matrix)] <- quantile(expr_matrix, 0.99, na.rm = TRUE)
 }
 
 cat("Expression matrix size:", dim(expr_matrix), "\n")
 cat("Design matrix size:", dim(design), "\n")
 
 # 2. Analysis with alternative simple model LM
 
 group_list <- factor(c(rep("Control", 10), rep("Pompe", 9)), levels = c("Control", "Pompe"))
 
 gs <- group_list
 design_simple <- model.matrix(~ gs)
 fit_simple <- lmFit(as.matrix(expr_matrix), design_simple)
 fit_simple <- eBayes(fit_simple)
 
 deg_simple <- topTable(fit_simple, coef = "gsPompe", number = Inf,
                        adjust.method = "BH", p.value = 0.05, lfc = 1)
 
 write.csv(deg_simple, "Pompe_vs_Control_simple_DEGs.csv", row.names = TRUE)
 deg_table <- read.csv("Pompe_vs_Control_simple_DEGs.csv", stringsAsFactors = FALSE)
 View(deg_table)  # Shows the table in a new tab in RStudio
 
 
 
 # Install the required package (hgu133plus2.db or Adding symbols with biomaRt annotation)
 
 if (!requireNamespace("hgu133plus2.db")) BiocManager::install("hgu133plus2.db")
 library(hgu133plus2.db)
 library(biomaRt)
 
 # … (your annotation code will be here)
 # 3. Retrieve the probe IDs in order.
 probe_ids <- rownames(deg_table)
 
 # 4. Match the symbol with select() (ENTREZID, GENENAME can also be added if desired)
 anno_df <- select(hgu133plus2.db,
                   keys = probe_ids,
                   columns = c("SYMBOL", "GENENAME"),
                   keytype = "PROBEID")
 
 # 5. Merge the resulting annotation table with the original DEG table.
 # (note the order! - probeID must be unique)
 deg_table$SYMBOL <- anno_df$SYMBOL[match(probe_ids, anno_df$PROBEID)]
 deg_table$GENENAME <- anno_df$GENENAME[match(probe_ids, anno_df$PROBEID)]
 write.csv(deg_table, "DEG_Table_Annotated.csv")  # Save the results
 
 # 6. WebGestaltR + miRNA analysis
 
 if (!require("WebGestaltR", quietly=TRUE)) {
   BiocManager::install("WebGestaltR")
   library(WebGestaltR)
 }
 significant_genes <- unique(deg_results$Gene.Symbol[!is.na(deg_results$Gene.Symbol)])
 available_dbs <- WebGestaltR::listGeneSet(organism="hsapiens")
 miRNA_dbs <- available_dbs[grep("miRNA", available_dbs$name, ignore.case=TRUE), ]
 print(miRNA_dbs)
 
 # 7. PATHWAY analysis via WebGestaltR
 pathway_result <- WebGestaltR(
   enrichMethod="ORA", organism="hsapiens",
   enrichDatabase="pathway_KEGG",
   interestGene=significant_genes,
   interestGeneType="genesymbol", referenceSet="genome",
   minNum=5, fdrThr=0.05,
   isOutput=TRUE, outputDirectory=getwd(),
   projectName="Pompe_Pathway_Analysis"
 )
 View(pathway_result)  # View in RStudio
 top_pw <- head(pathway_result[order(pathway_result$FDR), ], 10)
 barplot(-log10(top_pw$FDR), names.arg=top_pw$description, las=2, cex.names=0.7,
         main="Top KEGG pathways", ylab="-log10(FDR)", col="skyblue")
 
 
 # 8. validated (simplified) with multiMiR
 if (!require("multiMiR", quietly=TRUE)) BiocManager::install("multiMiR")
 library(multiMiR)
 
 
 # 8a. Query only reliable and working databases (except targetscan)
 mm <- tryCatch({
   get_multimir(
     org = "hsa",
     target = significant_genes,
     table = "validated", 
     summary = TRUE
   )
 }, error = function(e) {
   message("Error in multiMiR query:", e$message)
   return(NULL)
 })
 
 if (!is.null(mm)) {
   print(summary(mm))
   
   mmr_sum <- with(mm@data, tapply(type, mature_mirna_id, length))
   top_mmrs <- head(sort(mmr_sum, decreasing = TRUE), 10)
   barplot(top_mmrs, las = 2, col = "black", main = "miRNA Target Count")
 } else {
   cat("Could not retrieve multiMiR data, analysis skipped.\n")
 }
  
 # 8b. UMAP and quality control charts
 library(umap); library(maptools)
 ex2 <- na.omit(ex[!duplicated(rownames(ex)), ])
 ump <- umap(t(ex2), n_neighbors=8, random_state=123)
 plot(ump$layout, col=gs, pch=20, main="UMAP plot (n=8)")
 legend("topright", legend=levels(gs), col=1:length(levels(gs)), pch=20)
 
 
 
 
 
 # ---- multiMiR sorgusu sonrası ----
 # mm nesnesi (ya da mm@data) elinde olmalı!
 
 if (!is.null(mm)) {
   # Hedef miRNA sayısı tablosu
   miRNA_target_counts <- with(mm@data, tapply(type, mature_mirna_id, length))
   top_miRNAs <- sort(miRNA_target_counts, decreasing = TRUE)[1:10]
   df <- data.frame(
     miRNA = names(top_miRNAs),
     Target_Count = as.numeric(top_miRNAs)
   )
   
   # ggplot2 ile barplot (yatay)
   library(ggplot2)
   p <- ggplot(df, aes(x = reorder(miRNA, Target_Count), y = Target_Count)) +
     geom_bar(stat = "identity", fill = "deepskyblue4") +
     coord_flip() +
     labs(title = "Top Targeting miRNAs",
          x = "miRNA",
          y = "Number of Validated Targets") +
     theme_minimal(base_size = 15)
   
   print(p)
   
   # PNG olarak kaydet
   ggsave("top_targeting_miRNAs.png", plot = p, width = 8, height = 6, dpi = 300)
   # PDF olarak kaydet
   ggsave("top_targeting_miRNAs.pdf", plot = p, width = 8, height = 6)
 }
 