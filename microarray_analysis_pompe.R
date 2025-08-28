# =============================================================================
# Pompe Disease Microarray Transcriptomics Analysis Pipeline
# Author: Fatma Tosun- Harun Bayrak
# Date: 202-07-24
# Affiliation: [Your Affiliation]
# Description: Full pipeline for normalization, DEG analysis, annotation, functional enrichment (GO/KEGG), miRNA analysis, visualization
# License: MIT
# =============================================================================

#========================== Install and load required packages for analysis ===================================

packages <- c(␊
  "pkgbuild", "AnnotationDbi", "Biobase", "DOSE", "GEOquery", "GOSemSim", "R.utils",␊
  "affy", "annotate", "clusterProfiler", "data.table", "dplyr", "enrichplot",␊
  "genefilter", "ggplot2", "ggrepel", "grid", "gridExtra", "hgu133plus2.db",␊
  "jsonlite", "limma", "msigdbr", "multiMiR", "org.Hs.eg.db", "pheatmap",␊
  "sva", "umap", "STRINGdb", "rDGIdb", "igraph", "ggraph", "tidygraph"
)␊

#------------------------   Load required packages; assumes they are already installed   ----------------------

invisible(lapply(packages, library, character.only = TRUE))

  align_samples <- function(cel_files, pheno) {
   sample_names <- basename(cel_files)
   pheno_df <- as.data.frame(pheno, stringsAsFactors = FALSE)
   pheno_keys <- toupper(basename(trimws(pheno_df$filename)))
   sample_keys <- toupper(sample_names)
   idx <- match(sample_keys, pheno_keys)
   if (any(is.na(idx))) {
     stop("pheno$filename içinde bulunamayan örnek(ler): ",
         paste(sample_names[is.na(idx)], collapse = ", "))
   }
   aligned <- pheno_df[idx, , drop = FALSE]
   rownames(aligned) <- sample_names
   stopifnot(nrow(aligned) == length(sample_names))
   stopifnot(identical(rownames(aligned), sample_names))
   aligned
 }

pkgbuild::has_build_tools(debug = TRUE)

args <- commandArgs(trailingOnly = TRUE)
FDR_THRESH   <- if (length(args) >= 1) as.numeric(args[1]) else 0.1
LOGFC_THRESH <- if (length(args) >= 2) as.numeric(args[2]) else 0.5

#-----------------------------   Load CEL files and normalize with RMA    ---------------------------------------------

  root_dir <- normalizePath(Sys.getenv("POMPE_ROOT", getwd()), winslash = "/")
  dir.create(file.path(root_dir, "metadata"), showWarnings = FALSE, recursive = TRUE)␊
  meta_dir <- file.path(root_dir, "metadata")␊
␊
  stopifnot(dir.exists(meta_dir))␊

  dir.create(file.path(root_dir, "results", "qc"),    recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root_dir, "results", "deg"),   recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root_dir, "results", "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root_dir, "results", "miRNA_analysis"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root_dir, "logs"),             recursive = TRUE, showWarnings = FALSE)

  writeLines(c(capture.output(sessionInfo())), file.path(root_dir, "logs", "session_info.txt"))

  set.seed(20240724)

  doc_file <- file.path(Sys.getenv("HOME"), "Documents", "pheno.csv")
  if (!file.exists(doc_file)) {
    metadata_url <- "https://example.com/pheno.csv"  # replace with actual download link
    dir.create(dirname(doc_file), recursive = TRUE, showWarnings = FALSE)
    utils::download.file(metadata_url, doc_file, mode = "wb")
  } else {
    stop("pheno.csv not found. Ensure metadata is downloaded via Git LFS to ~/Documents.")
  }
  file.copy(doc_file, file.path(meta_dir, "pheno.csv"), overwrite = TRUE)
  pheno <- data.table::fread(file.path(meta_dir, "pheno.csv"),
                             na.strings = c("", "NA", "NaN"))

pheno$filename <- basename(trimws(pheno$filename))
cel_dir   <- meta_dir
cel_files <- file.path(cel_dir, pheno$filename)
missing   <- pheno$filename[!file.exists(cel_files)]
if (length(missing)) stop("Missing CEL files in metadata/: ", paste(missing, collapse = ", "))

stopifnot(all(c("sample","group","filename") %in% names(pheno)))
g <- tolower(trimws(as.character(pheno$group)))
g[g %in% c("control","healthy","normal")] <- "control"
g[grepl("pompe", g)] <- "pompe"
pheno$group <- factor(ifelse(g %in% c("control","pompe"),
                             ifelse(g == "control","Control","Pompe"), NA_character_),
                      levels = c("Control","Pompe"))
stopifnot(!any(is.na(pheno$group)))

if (!"is_melas" %in% names(pheno)) pheno$is_melas <- FALSE
pheno$is_melas <- tolower(trimws(as.character(pheno$is_melas))) %in%
  c("1","true","t","yes","y")

melas_idx <- which(pheno$is_melas)
if (length(melas_idx) > 1) {
  warning("More than one MELAS sample detected; using the first flagged sample.")
  melas_idx <- melas_idx[1]
}

if (!"batch" %in% names(pheno)) pheno$batch <- NA

# ----------------   ExpressionSet ID eşitleme (assayData / phenoData / protocolData)   -----------------------

pheno_aligned <- align_samples(cel_files, pheno)
if (inherits(pheno_aligned, "data.table")) data.table::setDF(pheno_aligned)

if (anyDuplicated(rownames(pheno_aligned))) {
  rownames(pheno_aligned) <- make.unique(rownames(pheno_aligned))
}

if (!exists("eset")) {
  stopifnot(exists("cel_files"))
  raw  <- ReadAffy(filenames = cel_files)
  eset <- rma(raw)
  expr <- exprs(eset)
}

Biobase::sampleNames(eset) <- rownames(pheno_aligned)

expr <- Biobase::exprs(eset)                 
colnames(expr) <- Biobase::sampleNames(eset) 

ord <- match(Biobase::sampleNames(eset), rownames(pheno_aligned))
stopifnot(!any(is.na(ord)))                
pheno_aligned <- pheno_aligned[ord, , drop = FALSE]

Biobase::pData(eset) <- pheno_aligned
Biobase::protocolData(eset) <- Biobase::AnnotatedDataFrame(
  data.frame(row.names = Biobase::sampleNames(eset))
)

methods::validObject(eset)  # TRUE dönmeli
expr_matrix <- Biobase::exprs(eset)
colnames(expr_matrix) <- Biobase::sampleNames(eset)

#========================== Conduct quality control and filter low-expression probes ===================================

  iqr <- apply(expr, 1, IQR)
thr <- quantile(iqr, 0.20, na.rm = TRUE)   # alt %20 IQR dışarı
keep <- iqr > thr
expr_f <- expr[keep, ]
message("Low-expression removed: ", sum(!keep), " probes (", round(mean(!keep)*100,2), "%)")
write.csv(data.frame(probe=rownames(expr)[!keep], IQR=iqr[!keep]),
          file.path(root_dir,"results","qc","low_expression_removed.csv"), row.names = FALSE)

pca <- prcomp(t(expr_f), scale. = TRUE)
pca_df <- data.frame(pca$x[,1:2], group = pheno_aligned$group,
                     batch = pheno_aligned$batch, sample = pheno_aligned$sample)
gg <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = factor(batch), label = sample)) +
  geom_point(size=3) + ggrepel::geom_text_repel(size=2.8) +
  theme_bw(14) + labs(title="PCA – RMA + IQR filtre", shape = "batch")
ggplot2::ggsave(file.path(root_dir,"results","qc","pca_groups_batches.png"), gg, width=9, height=7)

#========================== Identify differentially expressed genes with limma ===================================

expr_cb <- if (exists("expr_f")) expr_f else expr
stopifnot(is.matrix(expr_cb))
stopifnot(exists("pheno_aligned"))
groups <- factor(pheno_aligned$group, levels = c("Control","Pompe"))
design <- model.matrix(~0 + groups)
colnames(design) <- c("Control","Pompe")
stopifnot(ncol(expr_cb) == nrow(design))

print(dim(expr_cb)); print(dim(design)); head(design)

age_m <- sex_f <- rin_n <- NULL

if ("age_at_baseline_mounth" %in% names(pheno_aligned) ||
    "age_at_baseline_month"  %in% names(pheno_aligned)) {
  age_col <- intersect(c("age_at_baseline_mounth","age_at_baseline_month"),
                       names(pheno_aligned))[1]
  age_m <- suppressWarnings(as.numeric(gsub(",", ".", trimws(pheno_aligned[[age_col]]))))
}

if ("sex" %in% names(pheno_aligned)) {
  sex_f <- factor(toupper(trimws(as.character(pheno_aligned$sex))))
}

if ("rin" %in% names(pheno_aligned)) {
  rin_raw <- trimws(as.character(pheno_aligned$rin))
  rin_raw[rin_raw %in% c("TRUE","T","Yes","Y")] <- NA  # sayı olmayanları NA yap
  rin_n <- suppressWarnings(as.numeric(gsub(",", ".", rin_raw)))
}

design_df <- data.frame(groups = factor(pheno_aligned$group, levels=c("Control","Pompe")),
                        stringsAsFactors = FALSE)

if (!is.null(age_m)) {
  sda <- stats::sd(age_m, na.rm = TRUE)
  if (!is.na(sda) && sda > 0) design_df$age_m <- age_m
}

if (!is.null(sex_f)) {
  sex_f <- droplevels(sex_f)
  if (nlevels(sex_f) >= 2) design_df$sex_f <- sex_f
}

if (!is.null(rin_n)) {
  sdr <- stats::sd(rin_n, na.rm = TRUE)
  if (!is.na(sdr) && sdr > 0) design_df$rin_n <- rin_n
}

design_cov <- stats::model.matrix(~ 0 + ., data = design_df)
colnames(design_cov) <- sub("^groupsControl$", "Control",
                            sub("^groupsPompe$",   "Pompe", colnames(design_cov)))
rownames(design_cov) <- rownames(pheno_aligned)

  fit  <- limma::lmFit(expr_cb, design)
cont <- limma::makeContrasts(PompeVsControl = Pompe - Control, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

tt_all <- limma::topTable(fit2, coef = "PompeVsControl", number = Inf, adjust.method = "BH")

probes <- rownames(expr_cb)
ann <- AnnotationDbi::select(hgu133plus2.db, keys = probes,
                             columns = c("SYMBOL","GENENAME"),
                             keytype = "PROBEID")

#sym    <- ann$SYMBOL[match(probes, ann$PROBEID)]
#genenm <- ann$GENENAME[match(probes, ann$PROBEID)]

#tt_all$SYMBOL   <- sym[rownames(tt_all)]
#tt_all$GENENAME <- genenm[rownames(tt_all)]


ann <- AnnotationDbi::select(hgu133plus2.db,
                             keys = rownames(tt_all),
                             columns = c("SYMBOL","GENENAME"),
                             keytype = "PROBEID")

lkp_sym  <- setNames(ann$SYMBOL,  ann$PROBEID)
lkp_name <- setNames(ann$GENENAME, ann$PROBEID)

tt_all$SYMBOL   <- unname(lkp_sym[rownames(tt_all)])
tt_all$GENENAME <- unname(lkp_name[rownames(tt_all)])


tt_all$abs_logFC <- abs(tt_all$logFC)
tt_gene <- tt_all[!is.na(tt_all$SYMBOL) & nzchar(tt_all$SYMBOL), ]
tt_gene <- tt_gene[order(tt_gene$SYMBOL, -tt_gene$abs_logFC), ]
tt_gene <- tt_gene[!duplicated(tt_gene$SYMBOL), ]

stopifnot(exists("tt_all"))

table(is.na(tt_all$SYMBOL))   # artık çoğu FALSE olmalı
nrow(tt_gene)                 # >0 olmalı

deg <- subset(tt_gene, adj.P.Val < FDR_THRESH & abs(logFC) >= LOGFC_THRESH)
if (nrow(deg) == 0) stop("No significant DEGs; lower thresholds or check design.")

up_n   <- sum(deg$logFC > 0, na.rm = TRUE)
down_n <- sum(deg$logFC < 0, na.rm = TRUE)

write.csv(tt_all,  file.path(root_dir,"results","deg","all_probes_BH.csv"))
write.csv(tt_gene, file.path(root_dir,"results","deg","collapsed_by_gene.csv"), row.names = TRUE)

# ===================== PPI analysis via STRINGdb =====================
ppi_dir <- file.path(root_dir, "results", "ppi")
dir.create(ppi_dir, recursive = TRUE, showWarnings = FALSE)
if (nrow(deg) > 1) {
  try({
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
    deg_mapped <- string_db$map(deg, "SYMBOL", removeUnmappedRows = TRUE)
    if (nrow(deg_mapped) > 1) {
      ppi_edges <- string_db$get_interactions(deg_mapped$STRING_id)
      write.csv(ppi_edges, file.path(ppi_dir, "string_interactions.csv"), row.names = FALSE)
      png(file.path(ppi_dir, "string_network.png"), width = 800, height = 600)
      string_db$plot_network(deg_mapped$STRING_id)
      dev.off()
    }
  }, silent = TRUE)
}

deg_for_ppi <- if (exists("melas_res") && !is.null(melas_res$deg_noM_sig)) melas_res$deg_noM_sig else deg

if (nrow(deg_for_ppi) > 0) {
  try({
    string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
    deg_mapped <- string_db$map(deg_for_ppi, "ENTREZID", removeUnmappedRows = TRUE)
    inter <- string_db$get_interactions(deg_mapped$STRING_id)
    if (nrow(inter) > 0) {
      g_ppi <- igraph::graph_from_data_frame(inter, directed = FALSE)
      plt_ppi <- ggraph::ggraph(g_ppi, layout = "fr") +
        ggraph::geom_edge_link(alpha = 0.3) +
        ggraph::geom_node_point(color = "steelblue", size = 3) +
        ggraph::theme_void()
      ggsave(file.path(root_dir, "results", "plots", "ppi_network.png"), plt_ppi,
             width = 8, height = 6, dpi = 300)
    }
  }, silent = TRUE)

  try({
    mir_dir <- file.path(root_dir, "results", "miRNA_analysis")
    mir_hits <- multiMiR::get_multimir(target = rownames(deg_for_ppi), table = "validated")
    genes <- unique(mir_hits@data$target_symbol)
    drug_res <- rDGIdb::queryDGIdb(genes)
    write.csv(drug_res$matchedTerms,
              file.path(mir_dir, "drug_candidates.csv"), row.names = FALSE)
  }, silent = TRUE)
}

#======================= Protein-protein interaction network (STRINGdb) =======================
if (nrow(deg) > 0) {
  string_db <- STRINGdb::STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
  mapped_deg <- string_db$map(deg, "SYMBOL", removeUnmappedRows = TRUE)
  hits <- mapped_deg$STRING_id
  if (length(hits) > 0) {
    png(file.path(root_dir, "results", "deg", "ppi_network.png"), width = 10, height = 8, units = "in", res = 300)
    string_db$plot_network(hits)
    dev.off()
    ppi_edges <- string_db$get_interactions(hits)
    write.csv(ppi_edges, file.path(root_dir, "results", "deg", "ppi_edges.csv"), row.names = FALSE)
  }
}
nrow(tt_all)
table("SYMBOL_is_NA" = is.na(tt_all$SYMBOL))
mean(!is.na(tt_all$SYMBOL))  # kapsama oranı
probes <- rownames(expr_cb)
all(rownames(tt_all) %in% probes)

head(AnnotationDbi::select(hgu133plus2.db,
                           keys=head(probes, 5),
                           columns=c("SYMBOL","GENENAME"),
                           keytype="PROBEID"))


cat("\n=== ÖZET ===\n")
cat("Design boyutu: ", nrow(design), " örnek x ", ncol(design), " değişken\n", sep = "")
cat("Toplam gen (probe-level): ", nrow(tt_all), "\n", sep = "")
cat("Tekilleştirilmiş gen sayısı: ", nrow(tt_gene), "\n", sep = "")
write.csv(deg,     file.path(root_dir,"results","deg", sprintf("DEG_FDR_lt_%s.csv", FDR_THRESH)), row.names = TRUE)
cat("İlk 5 DEG:\n"); print(utils::head(deg[, c("SYMBOL","logFC","adj.P.Val","GENENAME")], 5))

de_sig <- deg             # FDR < FDR_THRESH gene-level
expr_use <- expr_cb       # batch/SVA sonrası matris (grafikler için)

if (!exists("deg")) {
  p_thresh   <- get0("p_thresh",   ifnotfound = get0("FDR_THRESH",   ifnotfound = 0.1))
  lfc_thresh <- get0("lfc_thresh", ifnotfound = get0("LOGFC_THRESH", ifnotfound = 0.5))
  stopifnot(exists("tt_gene"))
  deg <- subset(tt_gene, adj.P.Val < p_thresh & abs(logFC) >= lfc_thresh)
}
total_deg <- nrow(deg)
up_deg    <- sum(deg$logFC > 0, na.rm = TRUE)
down_deg  <- sum(deg$logFC < 0, na.rm = TRUE)
cat("Toplam DEG:", total_deg, "| Up:", up_deg, "| Down:", down_deg, "\n")

to_scalar_logical <- function(x) {
  if (is.logical(x)) return(isTRUE(any(x, na.rm = TRUE)))
  if (is.numeric(x)) return(isTRUE(any(x != 0, na.rm = TRUE)))
  if (is.character(x)) {
    y <- tolower(trimws(x))
    return(isTRUE(any(y %in% c("1","true","t","yes","y"), na.rm = TRUE)))
  }
  FALSE
}

prepare_melas <- function(pheno, expr, results_dir = file.path(getwd(), "results")) {
  ids <- colnames(expr)
  if (!"filename" %in% names(pheno))
    stop("pheno_aligned içinde 'filename' kolonu yok (GSM...CEL).")
  
  idx <- match(ids, pheno$filename)
  stopifnot(!any(is.na(idx)))
  pheno <- pheno[idx, , drop = FALSE]
  rownames(pheno) <- ids
  
  run_res <- run_deg_analysis(expr, pheno, run_melas = TRUE)
  deg_cov    <- run_res$deg_cov
  deg_cov_sig <- run_res$deg_cov_sig
  deg_noM    <- run_res$deg_noM
  deg_noM_sig <- run_res$deg_noM_sig
  
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(deg_cov,  file.path(results_dir, "DE_all_covariate_MELAS_included.csv"))
  write.csv(deg_noM,  file.path(results_dir, "DE_all_noMELAS.csv"))
  write.csv(deg_cov_sig, file.path(results_dir,
                                   sprintf("DE_sig_covariate_MELAS_included_FDR%s_LFC%s.csv",
                                           FDR_THRESH, LOGFC_THRESH)))
  write.csv(deg_noM_sig, file.path(results_dir,
                                   sprintf("DE_sig_noMELAS_FDR%s_LFC%s.csv",
                                           FDR_THRESH, LOGFC_THRESH)))
  
  list(pheno_aligned = pheno,
       deg_cov = deg_cov, deg_cov_sig = deg_cov_sig,
       deg_noM = deg_noM, deg_noM_sig = deg_noM_sig)
}

melas_flag <-
  ( "is_melas" %in% names(pheno_aligned) && to_scalar_logical(pheno_aligned$is_melas) ) ||
  ( "specialcase" %in% names(pheno_aligned) && to_scalar_logical(pheno_aligned$specialcase) )


melas_flag <-
  ( "is_melas" %in% names(pheno_aligned) &&
      isTRUE(any(tolower(trimws(as.character(pheno_aligned$is_melas))) %in% c("1","true","t","yes","y"), na.rm = TRUE)) ) ||
  ( "specialcase" %in% names(pheno_aligned) &&
      isTRUE(any(as.logical(pheno_aligned$specialcase), na.rm = TRUE)) )

pheno_aligned$is_melas <- factor(ifelse(melas_flag, "Yes","No"), levels = c("No","Yes"))

pheno_aligned$group <- factor(pheno_aligned$group, levels = c("Control","Pompe"))
if ("age_at_baseline_mounth" %in% names(pheno_aligned))
  pheno_aligned$age_at_baseline_mounth <- suppressWarnings(as.numeric(pheno_aligned$age_at_baseline_mounth))
if ("sex" %in% names(pheno_aligned))
  pheno_aligned$sex <- factor(pheno_aligned$sex)  # "M"/"F"
if ("rin" %in% names(pheno_aligned))
  pheno_aligned$rin <- suppressWarnings(as.numeric(pheno_aligned$rin))

stopifnot(nrow(pheno_aligned) == ncol(expr_f))
stopifnot(identical(rownames(pheno_aligned), colnames(expr_f)))

run_deg_analysis <- function(expr, pheno, run_melas = FALSE,
                             lfc_threshold = LOGFC_THRESH, fdr_threshold = FDR_THRESH) {
  design_cov <- model.matrix(~0 + group + is_melas + age_at_baseline_mounth + sex + rin,
                             data = pheno)
  colnames(design_cov) <- make.names(colnames(design_cov))
  fit_cov  <- limma::lmFit(expr, design_cov)
  cont_cov <- limma::makeContrasts(Pompe - Control, levels = design_cov)
  fit2_cov <- limma::eBayes(limma::contrasts.fit(fit_cov, cont_cov))
  cov_res  <- limma::topTable(fit2_cov, coef = 1, number = Inf, adjust.method = "BH")
  out <- list(deg_cov = cov_res,␊
              deg_cov_sig = subset(cov_res, adj.P.Val < fdr_threshold & abs(logFC) >= lfc_threshold))␊
  if (run_melas) {␊
    keep_idx <- which(
    keep_idx <- which(pheno$is_melas == "No")
    expr_noM <- expr[, keep_idx, drop = FALSE]
    ph_noM   <- droplevels(pheno[keep_idx, ])
    ph_noM$group <- factor(ph_noM$group, levels = c("Control","Pompe"))
    design_noM <- model.matrix(~0 + group + age_at_baseline_mounth + sex + rin, data = ph_noM)
    colnames(design_noM) <- make.names(colnames(design_noM))
    fit_noM  <- limma::lmFit(expr_noM, design_noM)
    cont_noM <- limma::makeContrasts(Pompe - Control, levels = design_noM)
    fit2_noM <- limma::eBayes(limma::contrasts.fit(fit_noM, cont_noM))
    deg_noM    <- limma::topTable(fit2_noM, coef = 1, number = Inf, adjust.method = "BH")
    out$deg_noM <- deg_noM
    out$deg_noM_sig <- subset(deg_noM, adj.P.Val < fdr_threshold & abs(logFC) >= lfc_threshold)
  }
  out
}

melas_res <- prepare_melas(pheno_aligned, expr_cb, file.path(root_dir, "results", "deg"))

plot_dir <- file.path(root_dir, "results", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

expr_plot <- if (exists("expr_use")) expr_use else if (exists("expr_cb")) expr_cb else expr_f
stopifnot(ncol(expr_plot) == nrow(pheno_aligned))

stopifnot(exists("tt_all"), exists("tt_gene"), exists("deg"))

volc_main <- within(tt_all, {
  negLog10FDR <- -log10(adj.P.Val)
  sig <- adj.P.Val < FDR_THRESH
  dir <- ifelse(sig & logFC > 0, "Up",
                ifelse(sig & logFC < 0, "Down", "NS"))
})

lab_main <- head(volc_main[order(volc_main$adj.P.Val, -abs(volc_main$logFC)), ], 15)
lab_main$SYMBOL[is.na(lab_main$SYMBOL) | lab_main$SYMBOL == ""] <- rownames(lab_main)[is.na(lab_main$SYMBOL) | lab_main$SYMBOL == ""]

volc <- ggplot(volc_main, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(shape = dir), alpha = 0.6) +
  geom_hline(yintercept = -log10(FDR_THRESH), linetype = 2) +
  geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = 3) +
  ggrepel::geom_text_repel(
    data = lab_main,
    aes(label = SYMBOL),
    size = 3, max.overlaps = 100
  ) +
  theme_bw(14) +
  labs(title = "Volcano (Pompe vs Control) — Main",
       x = "log2 Fold Change", y = "-log10(FDR)")
ggsave(file.path(plot_dir, "volcano_main.png"), volc, width = 8, height = 6, dpi = 300)

if (exists("tt_all_s") && exists("deg_sens")) {
  volc_s <- within(tt_all_s, {
    negLog10FDR <- -log10(adj.P.Val)
    sig <- adj.P.Val < FDR_THRESH
    dir <- ifelse(sig & logFC > 0, "Up",
                  ifelse(sig & logFC < 0, "Down", "NS"))
  })
  lab_s <- head(volc_s[order(volc_s$adj.P.Val, -abs(volc_s$logFC)), ], 15)
  lab_s$SYMBOL[is.na(lab_s$SYMBOL) | lab_s$SYMBOL == ""] <- rownames(lab_s)[is.na(lab_s$SYMBOL) | lab_s$SYMBOL == ""]
  
  p_s <- ggplot(volc_s, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(shape = dir), alpha = 0.6) +
    geom_hline(yintercept = -log10(FDR_THRESH), linetype = 2) +
    geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), linetype = 3) +
    ggrepel::geom_text_repel(
      data = lab_s,
      aes(label = SYMBOL),
      size = 3, max.overlaps = 100
    ) +
    theme_bw(14) +
    labs(title = "Volcano (Pompe vs Control) — MELAS excluded",
         x = "log2 Fold Change", y = "-log10(FDR)")
  ggsave(file.path(plot_dir, "volcano_noMELAS_main.png"), p_s, width = 8, height = 6, dpi = 300)
}


expr_plot <- if (exists("expr_use")) expr_use else if (exists("expr_cb")) expr_cb else expr_f  # batch/SVA sonrası

stopifnot(exists("tt_gene"), exists("deg"))
stopifnot(nrow(deg) > 0)

expr_plot <- if (exists("expr_use")) expr_use else if (exists("expr_cb")) expr_cb else expr_f

TOP_N <- min(50, nrow(deg))
stopifnot(TOP_N >= 1)

tt_gene$SYMBOL <- as.character(tt_gene$SYMBOL)
deg$SYMBOL     <- as.character(deg$SYMBOL)

top_syms <- head(deg$SYMBOL[order(deg$adj.P.Val)], TOP_N)

reps <- rownames(tt_gene)[ match(top_syms, tt_gene$SYMBOL) ]
reps <- unique(reps[!is.na(reps)])

sub_expr <- expr_plot[ intersect(rownames(expr_plot), reps), , drop = FALSE ]
stopifnot(nrow(sub_expr) > 0)   # mantıksal test

rownames(sub_expr) <- tt_gene$SYMBOL[ match(rownames(sub_expr), rownames(tt_gene)) ]
  
ann_col <- data.frame(Group = pheno_aligned$group,
                      row.names = rownames(pheno_aligned))
if ("batch" %in% names(pheno_aligned))
  ann_col$Batch <- pheno_aligned$batch

if ("is_melas" %in% names(pheno_aligned)) {
  melas_raw <- pheno_aligned$is_melas
  melas_log <- if (is.logical(melas_raw)) {
    melas_raw
  } else {
    v <- tolower(trimws(as.character(melas_raw)))
    v %in% c("1","true","t","yes","y")
  }
  melas_log[is.na(melas_log)] <- FALSE  # boşları "No" say
  ann_col$MELAS <- ifelse(melas_log, "Yes", "No")
}

ord <- order(ann_col$Group, ann_col$Batch, na.last = TRUE)
ann_col <- ann_col[ord, , drop = FALSE]
sub_expr <- sub_expr[, rownames(ann_col), drop = FALSE]

if (TOP_N >= 1) {  
  sub_expr_z <- t(scale(t(sub_expr)))  # satır bazlı z-score
  
  if (nrow(sub_expr_z) >= 2) {
    ph <- pheatmap::pheatmap(sub_expr_z,
                             annotation_col = ann_col,
                             show_rownames = TRUE, show_colnames = TRUE,
                             clustering_distance_rows = "correlation",
                             clustering_distance_cols = "correlation",
                             clustering_method = "average",
                             main = sprintf("Top %d DEG (z-score by gene) – Pompe vs Control", nrow(sub_expr_z)),
                       filename = file.path(plot_dir, sprintf("heatmap_top_FDR%s_main.png", FDR_THRESH))
    )
    
  } else if (nrow(sub_expr_z) == 1) {
    ph <- pheatmap::pheatmap(sub_expr_z,
                             annotation_col = ann_col,
                             show_rownames = TRUE, show_colnames = TRUE,
                             cluster_rows = FALSE, cluster_cols = TRUE,
                             main = "Only 1 gene passed filter (no row clustering)",
                       filename = file.path(plot_dir, sprintf("heatmap_top_FDR%s_main.png", FDR_THRESH))
    )
    
  } else {
    message("Heatmap atlandı: seçilen genlerden hiçbiri ifade matrisinde bulunamadı (probe vs gene adı?).")
  }
  
} else {
  message(sprintf("Heatmap atlandı: FDR<%s ile yeterli gen yok.", FDR_THRESH))
}


if (exists("tt_gene_s") && exists("expr_sens") && exists("ph_sens")) {
  top_gene_tbl_s <- head(tt_gene_s[order(tt_gene_s$adj.P.Val, -abs(tt_gene_s$logFC)), ], TOP_N)
  top_probes_s <- rownames(top_gene_tbl_s)
  sub_expr_s <- expr_sens[intersect(rownames(expr_sens), top_probes_s), , drop = FALSE]
  sub_expr_s_z <- t(scale(t(sub_expr_s)))
  
  ann_col_s <- data.frame(
    Group = ph_sens$group
  )
  if ("batch" %in% names(ph_sens)) ann_col_s$Batch <- ph_sens$batch
  rownames(ann_col_s) <- rownames(ph_sens)
  
  ord_s <- order(ann_col_s$Group, ann_col_s$Batch, decreasing = FALSE, na.last = TRUE)
  sub_expr_s_z <- sub_expr_s_z[, ord_s, drop = FALSE]
  ann_col_s    <- ann_col_s[ord_s, , drop = FALSE]
  
  png(file.path(plot_dir, sprintf("heatmap_top%d_genes_noMELAS.png", TOP_N)),
      width = 1100, height = 1400, res = 150)
  pheatmap(sub_expr_s_z,
           annotation_col = ann_col_s,
           show_rownames = TRUE, show_colnames = TRUE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "average",
           main = sprintf("Top %d genes (z-score) — Pompe vs Control (no MELAS)", TOP_N))
  dev.off()
}

message("Plots written to: ", plot_dir)

gsea_dir <- file.path(root_dir, "results", "deg")
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1234)

stopifnot(exists("tt_gene"))
ranks <- tt_gene$t
names(ranks) <- tt_gene$SYMBOL
ranks <- ranks[!is.na(ranks) & nzchar(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)

if (length(ranks) < 50) {
  warning("GSEA için yeterli sayıda gen yok gibi (rank < 50). Sonuçlar güvenilmez olabilir.")
}

msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
#msig <- try(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"), silent=TRUE)
if (inherits(msig, "try-error") || nrow(msig)==0) {
  msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
}
term2gene <- unique(msig[, c("gs_name","gene_symbol")])

term2gene <- unique(msig[, c("gs_name","gene_symbol")])
gsea_go <- GSEA(geneList = ranks, TERM2GENE = term2gene,
                pvalueCutoff = FDR_THRESH, verbose = FALSE)

write.csv(as.data.frame(gsea_go),
          file.path(gsea_dir, "GSEA_GO_BP_msigdb.csv"), row.names = FALSE)

png(file.path(gsea_dir, "GSEA_GO_BP_dotplot.png"), width=1000, height=800, res=130)
print(dotplot(gsea_go, showCategory = 20))
dev.off()

if (nrow(as.data.frame(gsea_go)) > 0) {
  top_terms <- head(gsea_go@result$ID, 2)
  for (id in top_terms) {
    png(file.path(gsea_dir, paste0("GSEA_GO_BP_runningplot_", gsub("[^A-Za-z0-9]+","_", id), ".png")),
        width=1100, height=700, res=130)
    print(gseaplot2(gsea_go, geneSetID = id))
    dev.off()
  }
}

sym2ent <- bitr(tt_gene$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
sym2ent <- sym2ent[!is.na(sym2ent$ENTREZID) & nzchar(sym2ent$SYMBOL), ]
sym2ent$abs_rank <- abs(ranks[sym2ent$SYMBOL])
sym2ent <- sym2ent[order(sym2ent$abs_rank, decreasing = TRUE), ]
sym2ent_1to1 <- sym2ent[!duplicated(sym2ent$SYMBOL), c("SYMBOL","ENTREZID")]

ranks_kegg <- ranks[sym2ent_1to1$SYMBOL]
names(ranks_kegg) <- sym2ent_1to1$ENTREZID
ranks_kegg <- ranks_kegg[!duplicated(names(ranks_kegg))]
ranks_kegg <- sort(ranks_kegg, decreasing = TRUE)

if (length(ranks_kegg) >= 50) {
  gsea_kegg <- gseKEGG(geneList = ranks_kegg,
                       organism = "hsa",
                       pvalueCutoff = 0.1,
                       verbose = FALSE)
  write.csv(as.data.frame(gsea_kegg),
            file.path(gsea_dir, "GSEA_KEGG.csv"), row.names = FALSE)
  
  png(file.path(gsea_dir, "GSEA_KEGG_dotplot.png"), width=1000, height=800, res=130)
  print(dotplot(gsea_kegg, showCategory = 20))
  dev.off()
} else {
  message("KEGG GSEA atlandı: ENTREZ eşleşmesi az ( < 50 ).")
}

if (exists("tt_gene_s")) {
  ranks_s <- tt_gene_s$t
  names(ranks_s) <- tt_gene_s$SYMBOL
  ranks_s <- ranks_s[!is.na(ranks_s) & nzchar(names(ranks_s))]
  ranks_s <- sort(ranks_s, decreasing = TRUE)
  
  if (length(ranks_s) >= 50) {
    gsea_go_s <- GSEA(geneList = ranks_s, TERM2GENE = term2gene,
                      pvalueCutoff = FDR_THRESH, verbose = FALSE)
    write.csv(as.data.frame(gsea_go_s),
              file.path(gsea_dir, "GSEA_GO_BP_msigdb_noMELAS.csv"), row.names = FALSE)
    
    png(file.path(gsea_dir, "GSEA_GO_BP_dotplot_noMELAS.png"), width=1000, height=800, res=130)
    print(dotplot(gsea_go_s, showCategory = 20))
    dev.off()
  } else {
    message("no-MELAS GSEA atlandı: yeterli rank yok.")
  }
}

if (exists("ph")) {
  pdf(file.path(root_dir,"results","plots","volcano_heatmap_Fig1.pdf"), width=16, height=8, family="Times")
  grid.arrange(
    arrangeGrob(volc, top = textGrob("A", gp=gpar(fontsize=22, fontface="bold"), x=unit(0,"npc"), hjust=0)),
    arrangeGrob(ph$gtable, top = textGrob("B", gp=gpar(fontsize=22, fontface="bold"), x=unit(0,"npc"), hjust=0)),
    ncol = 2
  )
  dev.off()
}

# =======================    ENRICHMENT: GO / KEGG / GSEA + MELAS       =======================

enrich_dir <- file.path(root_dir, "results", "enrichment"); dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir   <- file.path(root_dir, "results", "plots");       dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)

stopifnot(exists("tt_gene"), exists("deg"))
if (nrow(deg) == 0) stop("No DEGs detected; adjust thresholds.")
sig_symbols      <- unique(na.omit(as.character(deg$SYMBOL)))
universe_symbols <- unique(na.omit(as.character(tt_gene$SYMBOL)))

writeLines(sig_symbols,      file.path(enrich_dir, "_sig_gene_symbols.txt"))
writeLines(universe_symbols, file.path(enrich_dir, "_universe_gene_symbols.txt"))

sym2entrez      <- bitr(universe_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_universe <- unique(na.omit(sym2entrez$ENTREZID))
sig2entrez      <- bitr(sig_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_sig      <- unique(na.omit(sig2entrez$ENTREZID))

meta <- list(time=as.character(Sys.time()),
             n_sig_genes=length(sig_symbols), n_universe_genes=length(universe_symbols),
             n_sig_entrez=length(entrez_sig), n_universe_entrez=length(entrez_universe),
             adjust_method="BH", pvalueCutoff=FDR_THRESH, qvalueCutoff=FDR_THRESH)
writeLines(jsonlite::toJSON(meta, pretty=TRUE), file.path(enrich_dir, "_params_enrichment_main.json"))

go_bp <- enrichGO(gene = sig_symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "BP", universe = universe_symbols,
                  pAdjustMethod = "BH", pvalueCutoff = FDR_THRESH, qvalueCutoff = FDR_THRESH,
                  readable = TRUE)
write.csv(as.data.frame(go_bp), file.path(enrich_dir,"GO_BP_main.csv"), row.names = FALSE)

ek_args <- list(gene = entrez_sig, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = FDR_THRESH)
if ("universe" %in% names(formals(clusterProfiler::enrichKEGG))) ek_args$universe <- entrez_universe
kegg <- do.call(clusterProfiler::enrichKEGG, ek_args)
write.csv(as.data.frame(kegg), file.path(enrich_dir,"KEGG_main.csv"), row.names = FALSE)

gl_df <- tt_gene[, c("SYMBOL","t","logFC")]
gl_df$t[is.na(gl_df$t)] <- gl_df$logFC[is.na(gl_df$t)]
gl_df <- merge(gl_df, sym2entrez, by.x="SYMBOL", by.y="SYMBOL", all.x=TRUE)
gl_df <- gl_df[!is.na(gl_df$ENTREZID), c("ENTREZID","t")]
geneList <- gl_df$t; names(geneList) <- gl_df$ENTREZID; geneList <- sort(geneList, decreasing = TRUE)

gsea_go   <- NULL
gsea_kegg <- NULL
if (length(geneList) >= 20) {
  gsea_go   <- try(gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "ENTREZID",
                         ont = "BP", pAdjustMethod = "BH", verbose = FALSE), silent = TRUE)
  gsea_kegg <- try(gseKEGG(geneList = geneList, organism = "hsa", pAdjustMethod = "BH", verbose = FALSE), silent = TRUE)
  if (!inherits(gsea_go, "try-error"))   write.csv(as.data.frame(gsea_go),   file.path(enrich_dir,"GSEA_GO_BP_main.csv"),   row.names = FALSE)
  if (!inherits(gsea_kegg, "try-error")) write.csv(as.data.frame(gsea_kegg), file.path(enrich_dir,"GSEA_KEGG_main.csv"),    row.names = FALSE)
}
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  p_go_dot <- dotplot(go_bp, showCategory = 15, title = "GO BP (BH, explicit universe)")
  ggsave(file.path(plot_dir,"GO_BP_dotplot_main.png"), p_go_dot, width=10, height=9, dpi=300)
  p_go_bar <- barplot(go_bp, showCategory = 15, title = "GO BP (BH, explicit universe)")
  ggsave(file.path(plot_dir,"GO_BP_barplot_main.png"), p_go_bar, width=10, height=9, dpi=300)
  
  hsGO <- try(godata('org.Hs.eg.db', ont="BP"), silent=TRUE)
  if (!inherits(hsGO,"try-error")) {
    s <- try(pairwise_termsim(go_bp, semData=hsGO), silent=TRUE)
    if (!inherits(s,"try-error")) {
      p_emap <- try(emapplot(s, showCategory=15), silent=TRUE)
      if (!inherits(p_emap,"try-error")) ggsave(file.path(plot_dir,"GO_BP_emapplot_main.png"), p_emap, width=10, height=10, dpi=300)
    }
  }
}
if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
  p_kegg <- barplot(kegg, showCategory = 15, title = "KEGG (BH)")
  ggsave(file.path(plot_dir,"KEGG_barplot_main.png"), p_kegg, width=10, height=9, dpi=300)
}
if (!is.null(gsea_go) && !inherits(gsea_go,"try-error") && nrow(as.data.frame(gsea_go))>0) {
  p_gsea1 <- gseaplot2(gsea_go, geneSetID = head(gsea_go@result$ID,1), title = "GSEA GO BP (top set)")
  ggsave(file.path(plot_dir,"GSEA_GO_BP_topset_main.png"), p_gsea1, width=10, height=6, dpi=300)
}
if (!is.null(gsea_kegg) && !inherits(gsea_kegg,"try-error") && nrow(as.data.frame(gsea_kegg))>0) {
  p_gsea2 <- gseaplot2(gsea_kegg, geneSetID = head(gsea_kegg@result$ID,1), title = "GSEA KEGG (top pathway)")
  ggsave(file.path(plot_dir,"GSEA_KEGG_topset_main.png"), p_gsea2, width=10, height=6, dpi=300)
}

enrich_melas <- function(deg_tbl, tag) {
  sig_syms <- unique(na.omit(as.character(deg_tbl$SYMBOL)))
  if (length(sig_syms) < 3) return(invisible(NULL))
  sig_ent  <- unique(na.omit(bitr(sig_syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID))
  
  go <- enrichGO(gene=sig_syms, OrgDb=org.Hs.eg.db, keyType="SYMBOL",
                 ont="BP", universe=universe_symbols,
                 pAdjustMethod="BH", pvalueCutoff=FDR_THRESH, qvalueCutoff=FDR_THRESH, readable=TRUE)
  write.csv(as.data.frame(go), file.path(enrich_dir, paste0("GO_BP_",tag,".csv")), row.names=FALSE)
  
  ek_args <- list(gene=sig_ent, organism="hsa", pAdjustMethod="BH", pvalueCutoff=FDR_THRESH)
  if ("universe" %in% names(formals(clusterProfiler::enrichKEGG))) ek_args$universe <- entrez_universe
  kk <- do.call(clusterProfiler::enrichKEGG, ek_args)
  write.csv(as.data.frame(kk), file.path(enrich_dir, paste0("KEGG_",tag,".csv")), row.names=FALSE)
  
  if (!is.null(go) && nrow(as.data.frame(go))>0) {
    ggsave(file.path(plot_dir, paste0("GO_BP_dotplot_",tag,".png")),
           dotplot(go, showCategory=15, title=paste0("GO BP (",tag,")")), width=10, height=9, dpi=300)
  }
  if (!is.null(kk) && nrow(as.data.frame(kk))>0) {
    ggsave(file.path(plot_dir, paste0("KEGG_barplot_",tag,".png")),
           barplot(kk, showCategory=15, title=paste0("KEGG (",tag,")")), width=10, height=9, dpi=300)
  }
  invisible(list(go=go, kegg=kk))
}

if (exists("deg_sens")) {
  invisible(enrich_melas(deg_sens, "noMELAS"))
  
  has_melas_col <- "is_melas" %in% names(pheno_aligned)
  melas_log <- if (has_melas_col) {
    x <- pheno_aligned$is_melas
    if (is.logical(x)) x else {
      v <- tolower(trimws(as.character(x)))
      v %in% c("1","true","t","yes","y")
    }
  } else rep(FALSE, nrow(pheno_aligned))
  
  
} else if (exists("pheno_aligned") && any(melas_log, na.rm = TRUE)) {
  
  expr_base <- if (exists("expr_use")) expr_use else expr_cb
  
  keep <- !melas_log
  expr_s <- expr_base[, keep, drop = FALSE]
  ph_s   <- droplevels(pheno_aligned[keep, , drop = FALSE])
  
  if (nlevels(ph_s$group) >= 2 && all(table(ph_s$group) > 0)) {
    design_s <- model.matrix(~ 0 + group, data = ph_s)
    colnames(design_s) <- gsub("^group", "", colnames(design_s))  # "Control","Pompe"
    
    fit_s  <- limma::lmFit(expr_s, design_s)
    cont_s <- limma::makeContrasts(Pompe - Control, levels = design_s)
    fit2_s <- limma::eBayes(limma::contrasts.fit(fit_s, cont_s))
    
    tt_all_s <- limma::topTable(fit2_s, coef = 1, number = Inf, adjust.method = "BH")
  
  
  expr_s   <- expr_cb[, keep_idx, drop=FALSE]
  ph_s     <- droplevels(pheno_aligned[keep_idx, ])
  design_s <- model.matrix(~ 0 + group, data = ph_s); colnames(design_s) <- levels(ph_s$group)
  fit_s2   <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr_s, design_s),
                                                 limma::makeContrasts(Pompe - Control, levels = design_s)))

  tt_all_s$abs_logFC <- abs(tt_all_s$logFC)
  tt_gene_s <- tt_all_s[!is.na(tt_all_s$SYMBOL) & nzchar(tt_all_s$SYMBOL), ]
  tt_gene_s <- tt_gene_s[order(tt_gene_s$SYMBOL, -tt_gene_s$abs_logFC), ]
  tt_gene_s <- tt_gene_s[!duplicated(tt_gene_s$SYMBOL), ]
  deg_sens  <- subset(tt_gene_s, adj.P.Val < FDR_THRESH)
  write.csv(deg_sens,
            file.path(enrich_dir,
                      sprintf("DEG_noMELAS_FDR_lt_%s.csv", FDR_THRESH)),
            row.names=FALSE)
  invisible(enrich_melas(deg_sens, "noMELAS"))
  
} else {
  message("MELAS hariçte iki grup kalmadı; Pompe-Control karşılaştırması atlandı.")
  tt_all_s <- NULL
  deg_sens <- NULL
}

sum_lines <- c(
  sprintf("Main GO terms (q<%s): %d", FDR_THRESH, ifelse(is.null(go_bp), 0, nrow(as.data.frame(go_bp)))),
  sprintf("Main KEGG pathways (q<%s): %d", FDR_THRESH, ifelse(is.null(kegg), 0, nrow(as.data.frame(kegg)))),
  sprintf("Universe size: SYMBOL=%d | ENTREZ=%d", length(universe_symbols), length(entrez_universe)),
  "Notes: BH correction everywhere; GO uses explicit universe (all tested gene-level).",
  "Sensitivity: If MELAS present, enrichment repeated without MELAS; GSEA added to reduce threshold dependence."
)
writeLines(sum_lines, file.path(enrich_dir,"_ENRICHMENT_SUMMARY.txt"))

# ====================       miRNA Target Enrichment (FAIR, offline)         ====================

out_dir  <- file.path(root_dir, "results", "miRNA_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20240724)

stopifnot(exists("tt_gene"), "SYMBOL" %in% names(tt_gene))
stopifnot(exists("deg"),     "SYMBOL" %in% names(deg))
if (nrow(deg) == 0) stop("No DEGs detected; adjust thresholds.")
universe_symbols <- unique(tt_gene$SYMBOL[!is.na(tt_gene$SYMBOL) & nzchar(tt_gene$SYMBOL)])
sig_genes        <- unique(deg$SYMBOL[!is.na(deg$SYMBOL) & nzchar(deg$SYMBOL)])

meta <- list(
  dataset          = "GSE38680",
  n_control        = if (exists("groups")) sum(groups=="Control") else NA_integer_,
  n_pompe          = if (exists("groups")) sum(groups=="Pompe") else NA_integer_,
  n_universe_genes = length(universe_symbols),
  n_sig_genes      = length(sig_genes),
  fdr_threshold    = FDR_THRESH,
  adjust_method    = "BH",
  note_small_n     = "Small sample size; interpret enrichment cautiously."
)
jsonlite::write_json(meta, file.path(out_dir, "_mirna_meta.json"), pretty = TRUE)

enrich_mirna <- function(gene_set, universe_set, which_table=c("validated","predicted")) {
  which_table <- match.arg(which_table)
  # multiMiR'den hedef map (yerel kaynaklar)
  mm <- tryCatch(
    multiMiR::get_multimir(org="hsa", target = universe_set, table = which_table, summary = FALSE),
    error = function(e) NULL
  )
  if (is.null(mm) || nrow(mm@data)==0) return(NULL)
  
  df <- mm@data
  target_col <- c("target_symbol","target.gene","target")[
    c("target_symbol","target.gene","target") %in% names(df)][1]
  if (is.na(target_col)) return(NULL)
  
  df <- df[!is.na(df[[target_col]]) & nzchar(df[[target_col]]), c("mature_mirna_id", target_col, "database")]
  colnames(df) <- c("miRNA","SYMBOL","database")
  
  df <- df[df$SYMBOL %in% universe_set, , drop = FALSE]
  if (!nrow(df)) return(NULL)
  
  map_list <- split(df$SYMBOL, df$miRNA)
  res <- lapply(names(map_list), function(m) {
    tgt <- unique(map_list[[m]])
    a <- length(intersect(tgt, gene_set))
    b <- length(setdiff(tgt, gene_set))
    c <- length(setdiff(gene_set, tgt))
    d <- length(setdiff(universe_set, union(tgt, gene_set)))
    if ((a+b)==0 || (c+d)==0) return(NULL)
    ft <- suppressWarnings(fisher.test(matrix(c(a,b,c,d), nrow=2), alternative="greater"))
    data.frame(miRNA=m, a=a, b=b, c=c, d=d, p.value=ft$p.value, targets_in_deg=a,
               targets_total=a+b, stringsAsFactors = FALSE)
  })
  res <- dplyr::bind_rows(res)
  if (is.null(res) || !nrow(res)) return(NULL)
  res$FDR <- p.adjust(res$p.value, method="BH")
  res$method <- which_table
  res <- res[order(res$FDR, -res$targets_in_deg), ]
  res
}

tab_val <- enrich_mirna(sig_genes, universe_symbols, "validated")
tab_pre <- enrich_mirna(sig_genes, universe_symbols, "predicted")
enr_main <- dplyr::bind_rows(tab_val, tab_pre)
if (!is.null(enr_main) && nrow(enr_main)) {
  write.csv(enr_main, file.path(out_dir,"miRNA_enrichment_main.csv"), row.names = FALSE)
  top15 <- head(enr_main[order(enr_main$FDR, -enr_main$targets_in_deg), ], 15)
  if (nrow(top15)) {
    p <- ggplot(top15, aes(x=reorder(miRNA, -log10(FDR)), y=-log10(FDR), fill=method)) +
      geom_col() +
      coord_flip() +
      labs(title = "miRNA Target Enrichment (BH, universe = all tested genes)",
           x = "miRNA", y = expression(-log[10]*"FDR")) +
      theme_minimal(base_size = 12)
    ggsave(file.path(out_dir,"miRNA_enrichment_top15.png"), p, width=8, height=6, dpi=300)
  }
} else {
  message("miRNA enrichment: no signal (check DEG size / universe).")
}

if (!is.null(tab_val) && nrow(tab_val)) {
  best <- head(tab_val$miRNA[order(tab_val$FDR)], 5)
  mm_net <- tryCatch(
    multiMiR::get_multimir(org="hsa", target = sig_genes, table = "validated", summary = FALSE),
    error = function(e) NULL
  )
  if (!is.null(mm_net) && nrow(mm_net@data)) {
    df <- mm_net@data
    target_col <- c("target_symbol","target.gene","target")[
      c("target_symbol","target.gene","target") %in% names(df)][1]
    if (!is.na(target_col)) {
      net <- df[df$mature_mirna_id %in% best & df[[target_col]] %in% sig_genes,
                c("mature_mirna_id", target_col, "database")]
      colnames(net) <- c("miRNA","SYMBOL","database")
      if (nrow(net)) write.csv(net, file.path(out_dir,"network_top5_miRNA_to_DEG.csv"), row.names = FALSE)
    }
  }
}

sig_genes_s <- NULL
if (exists("pheno_aligned") && (("is_melas" %in% names(pheno_aligned)) || ("specialcase" %in% names(pheno_aligned)))) {
  mel_flag <- if ("is_melas" %in% names(pheno_aligned)) pheno_aligned$is_melas else pheno_aligned$specialcase
  mel_flag <- tolower(as.character(mel_flag)) %in% c("1","true","t","yes","y")
  if (any(mel_flag)) {
    if (exists("deg_sens") && "SYMBOL" %in% names(deg_sens)) {
      sig_genes_s <- unique(deg_sens$SYMBOL)
    } else if (exists("expr_cb")) {
      keep_idx <- which(!mel_flag)
      if (length(keep_idx) >= 4) {
        expr_s <- expr_cb[, keep_idx, drop=FALSE]
        ph_s   <- droplevels(pheno_aligned[keep_idx, ])
        des_s  <- model.matrix(~ 0 + group, data = transform(ph_s, group=factor(group, levels=c("Control","Pompe"))))
        colnames(des_s) <- levels(factor(ph_s$group, levels=c("Control","Pompe")))
        fit_s  <- limma::lmFit(expr_s, des_s)
        cont_s <- limma::makeContrasts(Pompe - Control, levels = des_s)
        fit_s2 <- limma::eBayes(limma::contrasts.fit(fit_s, cont_s))
        tt_s   <- limma::topTable(fit_s2, number=Inf, adjust.method="BH")
        tt_s$SYMBOL <- tt_s$SYMBOL %||% NA_character_
        if (!"SYMBOL" %in% names(tt_s)) {
        }
        if ("SYMBOL" %in% names(tt_s)) {
          tt_s$abs_logFC <- abs(tt_s$logFC)
          tt_sg <- tt_s[!is.na(tt_s$SYMBOL) & nzchar(tt_s$SYMBOL), ]
          tt_sg <- tt_sg[order(tt_sg$SYMBOL, -tt_sg$abs_logFC), ]
          tt_sg <- tt_sg[!duplicated(tt_sg$SYMBOL), ]
          sig_genes_s <- unique(tt_sg$SYMBOL[tt_sg$adj.P.Val < FDR_THRESH])
        }
      }
    }
  }
}
`%||%` <- function(x, y) if (is.null(x)) y else x

if (!is.null(sig_genes_s) && length(sig_genes_s)) {
  tab_val_s <- enrich_mirna(sig_genes_s, universe_symbols, "validated")
  tab_pre_s <- enrich_mirna(sig_genes_s, universe_symbols, "predicted")
  enr_sens  <- dplyr::bind_rows(tab_val_s, tab_pre_s)
  if (!is.null(enr_sens) && nrow(enr_sens)) {
    write.csv(enr_sens, file.path(out_dir,"miRNA_enrichment_noMELAS.csv"), row.names = FALSE)

    top_main <- head(enr_main$miRNA[order(enr_main$FDR)], 10)
    top_sens <- head(enr_sens$miRNA[order(enr_sens$FDR)], 10)
    ovlp     <- intersect(top_main, top_sens)
    writeLines(c(
      sprintf("Top10 overlap (validated+predicted): %d", length(ovlp)),
      paste("Common:", paste(ovlp, collapse=", "))
    ), file.path(out_dir,"_sensitivity_overlap.txt"))
  }
}

writeLines(c(capture.output(sessionInfo())), file.path(out_dir,"session_info.txt"))

message("miRNA enrichment (FAIR/offline) tamamlandı. Çıktılar: ", out_dir)

# ===================    MELAS SENSITIVITY VIS (opsiyonel; FAIR) GÖRSELLEŞTİRMYÖNE AL  ===================

has_melas <- "is_melas" %in% names(pheno_aligned) && any(pheno_aligned$is_melas)
if (has_melas) {
  keep_idx <- which(!pheno_aligned$is_melas)
  if (length(keep_idx) >= 4) {  # küçük n koruması
    expr_sens <- expr_use[, keep_idx, drop = FALSE]
    ph_sens   <- droplevels(pheno_aligned[keep_idx, ])
    des_s     <- model.matrix(~0 + group, data = ph_sens)
    colnames(des_s) <- levels(ph_sens$group)
    
    fit_s  <- limma::lmFit(expr_sens, des_s)
    cont_s <- limma::makeContrasts(Pompe - Control, levels = des_s)
    fit_s2 <- limma::eBayes(limma::contrasts.fit(fit_s, cont_s))
    
    tt_all_s <- limma::topTable(fit_s2, coef = 1, number = Inf, adjust.method = "BH")
    stopifnot("SYMBOL" %in% names(tt_all)) # ana tabloda anotasyon vardı
    ann_s <- AnnotationDbi::select(hgu133plus2.db,
                                   keys = rownames(tt_all_s),
                                   columns = c("SYMBOL"),
                                   keytype = "PROBEID")
    tt_all_s$SYMBOL <- ann_s$SYMBOL[match(rownames(tt_all_s), ann_s$PROBEID)]
    tt_all_s$abs_logFC <- abs(tt_all_s$logFC)
    tt_gene_s <- tt_all_s[!is.na(tt_all_s$SYMBOL) & nzchar(tt_all_s$SYMBOL), ]
    tt_gene_s <- tt_gene_s[order(tt_gene_s$SYMBOL, -tt_gene_s$abs_logFC), ]
    tt_gene_s <- tt_gene_s[!duplicated(tt_gene_s$SYMBOL), ]
    deg_sens  <- subset(tt_gene_s, adj.P.Val < FDR_THRESH)
    
    volc_sens <- ggplot(tt_gene_s, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(color = adj.P.Val < FDR_THRESH), alpha = 0.75, size = 1.6) +
      geom_vline(xintercept = c(-LOGFC_THRESH, LOGFC_THRESH), lty = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(FDR_THRESH), lty = "dashed", color = "grey50") +
      scale_color_manual(values = c("grey70","steelblue")) +
      labs(title = "Volcano (No‑MELAS; BH‑FDR; gene‑level)",
           x = expression(log[2]*"FC"), y = expression(-log[10]*"FDR"), color = NULL) +
      theme_bw(base_size = 12)
    
    ggsave(file.path(plots_dir, sprintf("volcano_FDR%s_geneLevel_noMELAS.png", FDR_THRESH)),
           volc_sens, width = 8, height = 6, dpi = 300)
    
    ov   <- length(intersect(deg$SYMBOL, deg_sens$SYMBOL))
    msg  <- sprintf("MELAS dahil DEG: %d | MELAS hariç DEG: %d | ortak: %d",
                    nrow(deg), nrow(deg_sens), ov)
    writeLines(msg, file.path(root_dir,"results","deg","MELAS_sensitivity_summary.txt"))
    message("[MELAS] ", msg)
  } else {
    message("[MELAS] Hariç sonrası örnek sayısı çok az; duyarlılık görselleri atlandı.")
  }
}

# ======================        miRNA VIS (multiMiR; FAIR; ayrı DB'ler)       ======================

mir_dir <- file.path(root_dir, "results", "miRNA_analysis")
dir.create(mir_dir, recursive = TRUE, showWarnings = FALSE)

writeLines("WARNING: small sample size (9 Pompe vs 10 Control). Interpret with caution.",
           file.path(mir_dir, "_NOTE_small_n.txt"))

sig_syms <- if (exists("deg")) head(unique(deg$SYMBOL), 50) else character(0)
dfv <- data.frame()
dfp <- data.frame()
if (length(sig_syms) < 2) {
  message("[miRNA] Yeterli anlamlı gen yok; multiMiR adımları atlandı.")
} else {
          
    mm_val <- multiMiR::get_multimir(org = "hsa",
                                     target = sig_syms,
                                     table = "validated",
                                     summary = FALSE)
    
        class(mm_val); isS4(mm_val); slotNames(mm_val)
    
    mm_to_df <- function(x) {
      if (is.null(x)) return(NULL)
      if (isS4(x) && "data" %in% slotNames(x)) as.data.frame(x@data) else
        if (is.data.frame(x)) x else NULL
    }
      
    dfv <- mm_val@data
    pcol_v <- intersect(c("p.value","p_value"), names(dfv))[1]
    if (!is.na(pcol_v)) {
      dfv$FDR <- p.adjust(dfv[[pcol_v]], method = "BH")
    }

    pltv <- dfv |>
      dplyr::count(database, mature_mirna_id, name="Target_Count") |>
      dplyr::group_by(database) |>
      dplyr::slice_max(Target_Count, n = 10, with_ties = FALSE) |>
      dplyr::ungroup()
    
    p_val <- ggplot(pltv, aes(x = mature_mirna_id, y = Target_Count, fill = database)) +
      geom_bar(stat="identity", position = position_dodge(width=.7), width=.6) +
      coord_flip() + theme_minimal(base_size = 12) +
      labs(title="Validated miRNA Targets (multiMiR; BH‑FDR applied)",
           x="miRNA", y="Target count", fill="DB")
    ggsave(file.path(mir_dir, "Validated_miRNA_targets.png"), p_val, width=8, height=6, dpi=300)
    write.csv(dfv, file.path(mir_dir, "validated_full_with_FDR.csv"), row.names = FALSE)
    }  
    else {
      message("[miRNA] Validated sonuç yok ya da erişilemedi.")
    }
    
  mm_pred <- tryCatch(
    multiMiR::get_multimir(org="hsa", target=sig_syms, table="predicted", summary=TRUE),
    error = function(e) NULL
  )
  if (!is.null(mm_pred) && nrow(mm_pred@data) > 0) {
    dfp <- mm_pred@data
   dfp$FDR <- p.adjust(dfp$p_value, method = "BH")
    
    pltp <- dfp |>
      dplyr::count(database, mature_mirna_id, name="Target_Count") |>
      dplyr::group_by(database) |>
      dplyr::slice_max(Target_Count, n = 10, with_ties = FALSE) |>
      dplyr::ungroup()
    
    p_pred <- ggplot(pltp, aes(x = mature_mirna_id, y = Target_Count, fill = database)) +
      geom_bar(stat="identity", position=position_dodge(width=.7), width=.6) +
      coord_flip() + theme_minimal(base_size = 12) +
      labs(title="Predicted miRNA Targets (multiMiR; BH‑FDR applied)",
           x="miRNA", y="Target count", fill="DB")
    ggsave(file.path(mir_dir, "Predicted_miRNA_targets.png"), p_pred, width=8, height=6, dpi=300)
    write.csv(dfp, file.path(mir_dir, "predicted_full_with_FDR.csv"), row.names = FALSE)
  } else {
    message("[miRNA] Predicted sonuç yok ya da erişilemedi.")
  }

  # ==================== Drug discovery from miRNA targets ====================
  gene_targets <- unique(c(
    if (exists("dfv") && "target_symbol" %in% names(dfv)) dfv$target_symbol else character(),
    if (exists("dfp") && "target_symbol" %in% names(dfp)) dfp$target_symbol else character()
  ))
  if (length(gene_targets) > 0) {
    drug_dir <- file.path(mir_dir, "drug_discovery")
    dir.create(drug_dir, recursive = TRUE, showWarnings = FALSE)
    try({
      dg <- rDGIdb::queryDGIdb(gene_targets)
      drug_tbl <- as.data.frame(dg@matches)
      write.csv(drug_tbl, file.path(drug_dir, "candidate_drugs.csv"), row.names = FALSE)
    }, silent = TRUE)
  }
}

#================= Drug discovery for miRNA targets (rDGIdb) ==================
genes_for_drug <- unique(c(dfv$target_symbol,
                           if (exists("dfp")) dfp$target_symbol else NULL))
genes_for_drug <- genes_for_drug[!is.na(genes_for_drug)]
if (length(genes_for_drug) > 0) {
  dg <- tryCatch(rDGIdb::queryDGIdb(genes_for_drug), error = function(e) NULL)
  if (!is.null(dg) && nrow(dg$matchedTerms) > 0) {
    write.csv(dg$matchedTerms,
              file.path(mir_dir, "miRNA_target_drugs.csv"), row.names = FALSE)
  }
}

expr2 <- na.omit(expr_use)
expr2 <- expr2[!duplicated(rownames(expr2)), ]
set.seed(123)
um <- umap::umap(t(expr2), n_neighbors = 8, random_state = 123)
ump_df <- data.frame(Sample = rownames(um$layout),
                     UMAP1 = um$layout[,1], UMAP2 = um$layout[,2],
                     Group = groups)
p_umap <- ggplot(ump_df, aes(UMAP1, UMAP2, color = Group, label = Sample)) +
  geom_point(size=3) + ggrepel::geom_text_repel(size=3, max.overlaps = 15) +
  theme_minimal(base_size = 13) + labs(title="UMAP QC (batch/SVA sonrası ifade)")
ggsave(file.path(mir_dir, "QC_UMAP.png"), p_umap, width=7, height=5, dpi=300)

meta_mir <- list(dataset="GSE38680", n_pompe=sum(groups=="Pompe"),
                 n_control=sum(groups=="Control"), adj="BH",
                 note_small_n=TRUE, date=as.character(Sys.Date()))
jsonlite::write_json(meta_mir, file.path(mir_dir, "_analysis_meta.json"), pretty=TRUE)







