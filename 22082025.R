# =============================================================================
# Pompe Disease Microarray Transcriptomics Analysis Pipeline
# Author: Fatma Tosun- Harun Bayrak
# Date: 202-07-24
# Affiliation: [Your Affiliation]
# Description: Full pipeline for normalization, DEG analysis, annotation, functional enrichment (GO/KEGG), miRNA analysis, visualization
# License: MIT
# =============================================================================

# ==========================AAAAAAAAAAAAAAA===================================

if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild")
pkgbuild::has_build_tools(debug = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery","affy","limma","Biobase","hgu133plus2.db","R.utils"),
                     ask = FALSE, update = TRUE)

library(GEOquery)
library(affy)
library(limma)
library(Biobase)
library(hgu133plus2.db)
library(R.utils)

base_dir  <- "C:\\Users\\User\\Desktop\\Pompe_Fair_Project_FULL\\Pompe_Fair_Project_FULL"
meta_dir  <- file.path(base_dir, "metadata")
cel_files <- list.files(meta_dir, pattern = "\\.CEL$", 
                        full.names = TRUE, ignore.case = TRUE)

stopifnot(length(cel_files) > 0)

length(cel_files)
head(basename(cel_files))

pheno_path <- file.path(meta_dir, "pheno.csv")
stopifnot(file.exists(pheno_path))
pheno <- read.csv(pheno_path, stringsAsFactors=FALSE)

cel_ids <- tools::file_path_sans_ext(basename(cel_files))    
gsm_from_cel <- sub("^(GSM[0-9]+).*", "\\1", cel_ids)         

id_candidates <- c("sample","Sample","geo_accession","GSM","gsm","SampleID","filename")
id_col <- id_candidates[id_candidates %in% names(pheno)][1]
stopifnot(!is.na(id_col))  

ph_ids <- pheno[[id_col]]
ph_ids <- tools::file_path_sans_ext(basename(trimws(as.character(ph_ids))))
ph_ids <- sub("^(GSM[0-9]+).*", "\\1", ph_ids)               

mi <- match(gsm_from_cel, ph_ids)
if (any(is.na(mi))) {
  cat("Eşleşmeyen GSM'ler (CEL tarafında olup pheno'da bulunamayanlar):\n")
  print(setdiff(gsm_from_cel, ph_ids))
  stop("GSM eşleşmesi bulunamadı; pheno.csv'deki ID kolonunu kontrol et.")
}
pheno <- pheno[mi, , drop = FALSE]

## Analiz boyunca tekil ve tutarlı olsun diye pheno'ya kesin sample etiketi yaz
pheno$sample <- cel_ids
stopifnot(nrow(pheno) == length(cel_ids))

pheno <- pheno[match(cel_ids, pheno$sample), , drop = FALSE]
stopifnot(identical(pheno$sample, cel_ids))

pheno$group <- trimws(pheno$group)
pheno$group <- ifelse(tolower(pheno$group) %in% c("control","healthy","normal"), "Control",
                      ifelse(grepl("pompe", tolower(pheno$group)), "Pompe", pheno$group))
stopifnot(all(pheno$group %in% c("Control","Pompe")))
print(table(pheno$group))

group <- factor(pheno$group, levels = c("Control","Pompe"))
design_group <- model.matrix(~ 0 + group)
colnames(design_group) <- levels(group)

# ========================== BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB ===============

sn <- sampleNames(eset)
gsm_from_sn <- sub(".*(GSM[0-9]+).*", "\\1", sn, ignore.case = TRUE)

rown_gsm <- rownames(pheno)
if (!all(grepl("^GSM[0-9]+$", rown_gsm))) {
  if ("geo_accession" %in% colnames(pheno)) rown_gsm <- pheno$geo_accession
}

expr_matrix <- exprs(eset)
sn_files <- basename(sampleNames(eset))        
pheno$filename <- trimws(pheno$filename)

norm <- function(x) toupper(gsub("\\\\", "/", basename(x)))
sn_norm    <- norm(sn_files)
pheno_norm <- norm(pheno$filename)

idx <- match(sn_norm, pheno_norm)
if (any(is.na(idx))) {
  missing <- sn_files[is.na(idx)]
  stop("Sample(s) in expression but not in pheno$filename: ",
       paste(missing, collapse = ", "))
}
pheno_aligned <- pheno[idx, , drop = FALSE]

stopifnot(ncol(expr_matrix) == nrow(pheno_aligned))

pheno_aligned$group <- trimws(pheno_aligned$group)
pheno_aligned$group <- ifelse(tolower(pheno_aligned$group) %in% c("control","healthy","normal"), "Control",
                              ifelse(grepl("pompe", tolower(pheno_aligned$group)), "Pompe", pheno_aligned$group))
stopifnot(all(pheno_aligned$group %in% c("Control","Pompe")))
groups <- factor(pheno_aligned$group, levels = c("Control","Pompe"))
design_group <- model.matrix(~ 0 + groups); colnames(design_group) <- levels(groups)

# ========================== CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC ===============

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("affy","limma","Biobase","hgu133plus2.db","annotate","genefilter","sva",
          "ggplot2","ggrepel","pheatmap","data.table")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
invisible(lapply(pkgs, library, character.only = TRUE))

set.seed(20240724)

root     <- normalizePath("C:\\Users\\User\\Desktop\\Pompe_Fair_Project_FULL\\Pompe_Fair_Project_FULL",
                          winslash = "/")
meta_dir <- file.path(root, "metadata")
stopifnot(dir.exists(meta_dir))

dir.create(file.path(root, "results", "qc"),    recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results", "deg"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "results", "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root, "logs"),             recursive = TRUE, showWarnings = FALSE)

writeLines(c(capture.output(sessionInfo())), file.path(root, "logs", "session_info.txt"))

pheno <- data.table::fread(file.path(meta_dir, "pheno.csv"), na.strings = c("", "NA", "NaN"))

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

# Batch (optional)
if (!"batch" %in% names(pheno)) pheno$batch <- NA


raw <- ReadAffy(filenames = cel_files)
pheno_df <- as.data.frame(pheno, stringsAsFactors = FALSE)
sn0      <- sampleNames(raw)                      
sn_key   <- toupper(basename(sn0))                
file_key <- toupper(basename(trimws(pheno_df$filename)))

idx <- match(sn_key, file_key)
if (any(is.na(idx))) {
  stop("pheno$filename içinde bulunamayan örnek(ler): ",
       paste(sn0[is.na(idx)], collapse = ", "))
}
pheno_aligned <- pheno_df[idx, , drop = FALSE]

rownames(pheno_aligned) <- sn0
eset <- rma(raw)
expr <- exprs(eset)   # QC bloğu 'expr' kullanıyor; bunu tanımla


## (İsteğe bağlı kontrol: isimler birebir mi?)
stopifnot(identical(sampleNames(eset), rownames(pheno_aligned)))

Biobase::pData(eset) <- pheno_aligned
validObject(eset)

iqr <- apply(expr, 1, IQR)
thr <- quantile(iqr, 0.20, na.rm = TRUE)   # alt %20 IQR dışarı
keep <- iqr > thr
expr_f <- expr[keep, ]
message("Low-expression removed: ", sum(!keep), " probes (", round(mean(!keep)*100,2), "%)")
write.csv(data.frame(probe=rownames(expr)[!keep], IQR=iqr[!keep]),
          file.path(root,"results","qc","low_expression_removed.csv"), row.names = FALSE)

# PCA ile QC (FAIR: şeffaflık)
pca <- prcomp(t(expr_f), scale. = TRUE)
pca_df <- data.frame(pca$x[,1:2], group = pheno_aligned$group,
                     batch = pheno_aligned$batch, sample = pheno_aligned$sample)
gg <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = factor(batch), label = sample)) +
  geom_point(size=3) + ggrepel::geom_text_repel(size=2.8) +
  theme_bw(14) + labs(title="PCA – RMA + IQR filtre", shape = "batch")
ggplot2::ggsave(file.path(root,"results","qc","pca_groups_batches.png"), gg, width=9, height=7)

# ========================== DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD ===============

stopifnot(exists("pheno_aligned"))
expr_cb <- if (exists("expr_f")) expr_f else expr_matrix
stopifnot(ncol(expr_cb) == nrow(pheno_aligned))

dir.create(file.path(root,"results","deg"), showWarnings = FALSE, recursive = TRUE)

do_combat <- "batch" %in% names(pheno_aligned) &&
  length(unique(na.omit(pheno_aligned$batch))) > 1
if (do_combat) {
  batch_vec <- factor(pheno_aligned$batch)
  mod_cb    <- model.matrix(~ pheno_aligned$group)   # grup kovaryatıyla düzelt
  expr_cb   <- sva::ComBat(dat = expr_cb, batch = batch_vec, mod = mod_cb,
                           par.prior = TRUE, prior.plots = FALSE)
} else {
  message("ComBat skipped: batch has <2 levels or all NA.")
}

groups <- droplevels(factor(pheno_aligned$group, levels = c("Control","Pompe")))
mod  <- model.matrix(~ groups)                   # tam model
mod0 <- model.matrix(~ 1, data = pheno_aligned)  # null model
k_sv <- tryCatch(sva::num.sv(expr_cb, mod, method = "leek"), error = function(e) 0)
k_sv <- if (is.na(k_sv) || k_sv < 1) 0 else k_sv
svobj <- if (k_sv > 0) sva::sva(expr_cb, mod, mod0, n.sv = k_sv) else NULL

design_df <- data.frame(group = groups, stringsAsFactors = FALSE)

add_cov <- function(nm) {
  if (!nm %in% names(pheno_aligned)) return(invisible(NULL))
  v <- pheno_aligned[[nm]]
  if (is.numeric(v)) {
    if (sum(!is.na(v)) >= 2 && stats::sd(v, na.rm = TRUE) > 0)
      design_df[[nm]] <<- v else message("Dropped numeric covariate (no variance): ", nm)
  } else {
    f <- droplevels(factor(v))
    if (nlevels(f) >= 2) design_df[[nm]] <<- f
    else message("Dropped factor covariate (only 1 level): ", nm)
  }
}

invisible(lapply(c("age","sex","rin"), add_cov))

if (!is.null(svobj) && !is.null(svobj$sv)) {
  sva_mat <- scale(svobj$sv, center = TRUE, scale = FALSE)       # merkezle
  colnames(sva_mat) <- paste0("SV", seq_len(ncol(sva_mat)))
  design_df <- cbind(design_df, sva_mat)
}

design <- model.matrix(~ 0 + ., data = design_df)

cn <- colnames(design)
cn[cn == "groupControl"] <- "Control"
cn[cn == "groupPompe"]   <- "Pompe"
colnames(design) <- cn

stopifnot(all(c("Control","Pompe") %in% colnames(design)))
stopifnot(nrow(design) == ncol(expr_cb))

cat("\n>> DESIGN sütunları:\n"); print(colnames(design))
cat(">> Num. SV (SVA):", k_sv, "\n")

fit  <- limma::lmFit(expr_cb, design)
cont <- limma::makeContrasts(PompeVsControl = Pompe - Control, levels = design)
fit2 <- limma::eBayes(limma::contrasts.fit(fit, cont))

tt_all <- limma::topTable(fit2, coef = "PompeVsControl", number = Inf, adjust.method = "BH")

probes <- rownames(expr_cb)
ann <- AnnotationDbi::select(hgu133plus2.db, keys = probes,
                             columns = c("SYMBOL","GENENAME"),
                             keytype = "PROBEID")
sym    <- ann$SYMBOL[match(probes, ann$PROBEID)]
genenm <- ann$GENENAME[match(probes, ann$PROBEID)]

tt_all$SYMBOL   <- sym[rownames(tt_all)]
tt_all$GENENAME <- genenm[rownames(tt_all)]

tt_all$abs_logFC <- abs(tt_all$logFC)
tt_gene <- tt_all[!is.na(tt_all$SYMBOL) & nzchar(tt_all$SYMBOL), ]
tt_gene <- tt_gene[order(tt_gene$SYMBOL, -tt_gene$abs_logFC), ]
tt_gene <- tt_gene[!duplicated(tt_gene$SYMBOL), ]

deg <- subset(tt_gene, adj.P.Val < 0.05)
up_n   <- sum(deg$logFC > 0, na.rm = TRUE)
down_n <- sum(deg$logFC < 0, na.rm = TRUE)

write.csv(tt_all,  file.path(root,"results","deg","all_probes_BH.csv"))
write.csv(tt_gene, file.path(root,"results","deg","collapsed_by_gene.csv"), row.names = TRUE)
write.csv(deg,     file.path(root,"results","deg","DEG_FDR_lt_0.05.csv"), row.names = TRUE)

cat("\n=== ÖZET ===\n")
cat("ComBat: ", if (do_combat) "APPLIED" else "SKIPPED", "\n")
cat("SVA (n.sv): ", k_sv, if (k_sv > 0) "  [SV'ler design'a eklendi]\n" else "  [kullanılmadı]\n", sep = "")
cat("Design boyutu: ", nrow(design), " örnek x ", ncol(design), " değişken\n", sep = "")
cat("Toplam gen (probe-level): ", nrow(tt_all), "\n", sep = "")
cat("Tekilleştirilmiş gen sayısı: ", nrow(tt_gene), "\n", sep = "")
cat("DEG (FDR<0.05): ", nrow(deg), "  |  Up: ", up_n, "  Down: ", down_n, "\n", sep = "")
cat("İlk 5 DEG:\n"); print(utils::head(deg[, c("SYMBOL","logFC","adj.P.Val","GENENAME")], 5))

de_all <- tt_all          # tüm sonuçlar (probe-level, anotasyonlu)
de_sig <- deg             # FDR<0.05 gene-level
expr_use <- expr_cb       # batch/SVA sonrası matris (grafikler için)

# ========================== EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ===============

if (exists("expr_use") && "is_melas" %in% names(pheno_aligned)) {
  keep_idx <- which(!pheno_aligned$is_melas)
  if (length(keep_idx) >= 4 && length(unique(pheno_aligned$group[keep_idx])) == 2) {
    expr_sens <- expr_use[, keep_idx, drop = FALSE]
    ph_sens   <- droplevels(pheno_aligned[keep_idx, ])
    
    n.sv.s <- tryCatch(sva::num.sv(expr_sens, model.matrix(~ ph_sens$group), method="leek"),
                       error=function(e) 0)
    if (is.na(n.sv.s) || n.sv.s < 1) n.sv.s <- 0
    sva_mat_s <- NULL
    if (n.sv.s > 0) {
      mod_s  <- model.matrix(~ ph_sens$group)
      mod0_s <- model.matrix(~ 1)
      sv_s   <- sva::sva(expr_sens, mod_s, mod0_s, n.sv = n.sv.s)
      sva_mat_s <- scale(sv_s$sv, center = TRUE, scale = FALSE)
      colnames(sva_mat_s) <- paste0("SVs", seq_len(ncol(sva_mat_s)))
    }
    
    design_s_df <- data.frame(group = ph_sens$group)
    lapply(c("age","sex","rin"), function(z) {
      if (z %in% names(ph_sens)) {
        v <- ph_sens[[z]]
        if (is.numeric(v) && sd(v, na.rm=TRUE) > 0) design_s_df[[z]] <<- v
        if (!is.numeric(v)) {
          f <- droplevels(factor(v)); if (nlevels(f) >= 2) design_s_df[[z]] <<- f
        }
      }
    })
    if (!is.null(sva_mat_s)) design_s_df <- cbind(design_s_df, sva_mat_s)
    
    design_s <- model.matrix(~ 0 + ., data = design_s_df)

    idxg <- which(colnames(design_s) %in% paste0("group", levels(ph_sens$group)))
    if (length(idxg) == 2) {
      colnames(design_s)[idxg] <- levels(ph_sens$group)
    } else {
      colnames(design_s)[1:2] <- levels(ph_sens$group)
    }
    
    fitS  <- limma::lmFit(expr_sens, design_s)
    contS <- limma::makeContrasts(Pompe - Control, levels = design_s)
    fitS2 <- limma::eBayes(limma::contrasts.fit(fitS, contS))
    
    tt_all_s <- limma::topTable(fitS2, coef = 1, number = Inf, adjust.method = "BH")
    
    tt_all_s$SYMBOL <- sym[match(rownames(tt_all_s), probes)]
    tt_all_s$GENENAME <- genenm[match(rownames(tt_all_s), probes)]
    tt_all_s$abs_logFC <- abs(tt_all_s$logFC)
    tt_gene_s <- tt_all_s[!is.na(tt_all_s$SYMBOL) & nzchar(tt_all_s$SYMBOL), ]
    tt_gene_s <- tt_gene_s[order(tt_gene_s$SYMBOL, -tt_gene_s$abs_logFC), ]
    tt_gene_s <- tt_gene_s[!duplicated(tt_gene_s$SYMBOL), ]
    deg_sens  <- subset(tt_gene_s, adj.P.Val < 0.05)
    
    overlap <- length(intersect(deg$SYMBOL, deg_sens$SYMBOL))
    msg <- sprintf("MELAS dahil DEG: %d | MELAS hariç DEG: %d | Ortak gen: %d",
                   nrow(deg), nrow(deg_sens), overlap)
    message(msg)
    writeLines(msg, file.path(root,"results","deg","MELAS_sensitivity_summary.txt"))
    write.csv(deg_sens, file.path(root,"results","deg","DEG_noMELAS_FDR_lt_0.05.csv"),
              row.names = TRUE)
  } else {
    message("MELAS duyarlılık atlandı: yeterli örnek/iki grup yok.")
  }
}

# =================== FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF ===================

plot_dir <- file.path(root, "results", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

expr_plot <- if (exists("expr_use")) expr_use else if (exists("expr_cb")) expr_cb else expr_f
stopifnot(ncol(expr_plot) == nrow(pheno_aligned))

stopifnot(exists("tt_all"), exists("tt_gene"), exists("deg"))

library(ggplot2); library(ggrepel)

volc_main <- within(tt_all, {
  negLog10FDR <- -log10(adj.P.Val)
  sig <- adj.P.Val < 0.05
  dir <- ifelse(sig & logFC > 0, "Up",
                ifelse(sig & logFC < 0, "Down", "NS"))
})

lab_main <- head(volc_main[order(volc_main$adj.P.Val, -abs(volc_main$logFC)), ], 15)
lab_main$SYMBOL[is.na(lab_main$SYMBOL) | lab_main$SYMBOL == ""] <- rownames(lab_main)[is.na(lab_main$SYMBOL) | lab_main$SYMBOL == ""]

p_main <- ggplot(volc_main, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(shape = dir), alpha = 0.6) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = 3) +
  ggrepel::geom_text_repel(
    data = lab_main,
    aes(label = SYMBOL),
    size = 3, max.overlaps = 100
  ) +
  theme_bw(14) +
  labs(title = "Volcano (Pompe vs Control) — Main",
       x = "log2 Fold Change", y = "-log10(FDR)")
ggsave(file.path(plot_dir, "volcano_main.png"), p_main, width = 8, height = 6, dpi = 300)

if (exists("tt_all_s") && exists("deg_sens")) {
  volc_s <- within(tt_all_s, {
    negLog10FDR <- -log10(adj.P.Val)
    sig <- adj.P.Val < 0.05
    dir <- ifelse(sig & logFC > 0, "Up",
                  ifelse(sig & logFC < 0, "Down", "NS"))
  })
  lab_s <- head(volc_s[order(volc_s$adj.P.Val, -abs(volc_s$logFC)), ], 15)
  lab_s$SYMBOL[is.na(lab_s$SYMBOL) | lab_s$SYMBOL == ""] <- rownames(lab_s)[is.na(lab_s$SYMBOL) | lab_s$SYMBOL == ""]
  
  p_s <- ggplot(volc_s, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(shape = dir), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_vline(xintercept = c(-1, 1), linetype = 3) +
    ggrepel::geom_text_repel(
      data = lab_s,
      aes(label = SYMBOL),
      size = 3, max.overlaps = 100
    ) +
    theme_bw(14) +
    labs(title = "Volcano (Pompe vs Control) — MELAS excluded",
         x = "log2 Fold Change", y = "-log10(FDR)")
  ggsave(file.path(plot_dir, "volcano_noMELAS.png"), p_s, width = 8, height = 6, dpi = 300)
}

library(pheatmap)

stopifnot(exists("tt_gene"), exists("deg"))
expr_plot <- if (exists("expr_use")) expr_use else if (exists("expr_cb")) expr_cb else expr_f  # batch/SVA sonrası

TOP_N <- min(50, nrow(deg))
if (TOP_N >= 1) {

    top_syms <- head(deg$SYMBOL[order(deg$adj.P.Val)], TOP_N)
  
  reps <- rownames(tt_gene)[match(top_syms, tt_gene$SYMBOL)]
  reps <- reps[!is.na(reps)]
  
  sub_expr <- expr_plot[intersect(rownames(expr_plot), reps), , drop = FALSE]
  
  rownames(sub_expr) <- tt_gene$SYMBOL[match(rownames(sub_expr), rownames(tt_gene))]
  
  ann_col <- data.frame(Group = pheno_aligned$group, row.names = rownames(pheno_aligned))
  if ("batch" %in% names(pheno_aligned))   ann_col$Batch  <- pheno_aligned$batch
  if ("is_melas" %in% names(pheno_aligned)) ann_col$MELAS <- ifelse(pheno_aligned$is_melas, "Yes","No")
  
  ord <- order(ann_col$Group, ann_col$Batch, na.last = TRUE)
  ann_col <- ann_col[ord, , drop = FALSE]
  sub_expr <- sub_expr[, rownames(ann_col), drop = FALSE]
  
  sub_expr_z <- t(scale(t(sub_expr)))  # satır bazlı z-score
  
  if (nrow(sub_expr_z) >= 2) {
    pheatmap::pheatmap(sub_expr_z,
                       annotation_col = ann_col,
                       show_rownames = TRUE, show_colnames = TRUE,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation",
                       clustering_method = "average",
                       main = sprintf("Top %d DEG (z-score by gene) – Pompe vs Control", nrow(sub_expr_z)),
                       filename = file.path(plot_dir, "heatmap_top_FDR.png")
    )
  } else if (nrow(sub_expr_z) == 1) {
    pheatmap::pheatmap(sub_expr_z,
                       annotation_col = ann_col,
                       show_rownames = TRUE, show_colnames = TRUE,
                       cluster_rows = FALSE, cluster_cols = TRUE,
                       main = "Only 1 gene passed filter (no row clustering)",
                       filename = file.path(plot_dir, "heatmap_top_FDR.png")
    )
  } else {
    message("Heatmap atlandı: seçilen genlerden hiçbiri ifade matrisinde bulunamadı (probe vs gene adı?).")
  }
} else {
  message("Heatmap atlandı: FDR<0.05 ile yeterli gen yok.")
}


if (exists("tt_gene_s") && exists("expr_sens") && exists("ph_sens")) {
  top_gene_tbl_s <- head(tt_gene_s[order(tt_gene_s$adj.P.Val, -abs(tt_gene_s$logFC)), ], topN)
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
  
  png(file.path(plot_dir, sprintf("heatmap_top%d_genes_noMELAS.png", topN)),
      width = 1100, height = 1400, res = 150)
  pheatmap(sub_expr_s_z,
           annotation_col = ann_col_s,
           show_rownames = TRUE, show_colnames = TRUE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "average",
           main = sprintf("Top %d genes (z-score) — Pompe vs Control (no MELAS)", topN))
  dev.off()
}

message("Plots written to: ", plot_dir)


# =================== GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG =====================

plot_dir <- file.path(root, "results", "plots")
deg_dir  <- file.path(root, "results", "deg")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir,  recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(hgu133plus2.db)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

if (!exists("tt_all")) {
  stopifnot(exists("fit2"))
  tt_all <- limma::topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
}

if (!"SYMBOL" %in% names(tt_all) || any(is.na(tt_all$SYMBOL))) {
  probes <- rownames(tt_all)
  sym  <- AnnotationDbi::mapIds(hgu133plus2.db, keys = probes,
                                keytype = "PROBEID", column = "SYMBOL",
                                multiVals = "first")
  gnm  <- AnnotationDbi::mapIds(hgu133plus2.db, keys = probes,
                                keytype = "PROBEID", column = "GENENAME",
                                multiVals = "first")
  
  #if (!"SYMBOL" %in% names(tt_all))   tt_all$SYMBOL   <- NA_character_
  #if (!"GENENAME" %in% names(tt_all)) tt_all$GENENAME <- NA_character_
  
  tt_all$SYMBOL[rownames(tt_all)  ]   <- sym[rownames(tt_all)]
  tt_all$GENENAME[rownames(tt_all)]   <- gnm[rownames(tt_all)]
}

if (!exists("tt_gene")) {
  tt_all$abs_logFC <- abs(tt_all$logFC)
  tt_gene <- tt_all[!is.na(tt_all$SYMBOL) & nzchar(tt_all$SYMBOL), ]
  tt_gene <- tt_gene[order(tt_gene$SYMBOL, -tt_gene$abs_logFC), ]
  tt_gene <- tt_gene[!duplicated(tt_gene$SYMBOL), ]
}

DEG_FDR_CUTOFF <- 0.05
if (!exists("deg")) {
  deg <- subset(tt_gene, adj.P.Val < DEG_FDR_CUTOFF)
}

groups <- droplevels(factor(pheno_aligned$group, levels = c("Control","Pompe")))
note <- sprintf("DEG (FDR<%.2g): %d gene(s). Sample sizes — Pompe: %d, Control: %d.",
                DEG_FDR_CUTOFF, nrow(deg),
                sum(groups == "Pompe"), sum(groups == "Control"))
writeLines(note, file.path(deg_dir, "_DEG_summary.txt"))

LFC_VISUAL_THR <- 1
stopifnot(all(c("logFC","adj.P.Val","SYMBOL") %in% names(tt_gene)))
labelN <- head(tt_gene[order(tt_gene$adj.P.Val), c("SYMBOL","logFC","adj.P.Val")], 20)
labelN <- labelN[!is.na(labelN$SYMBOL) & nzchar(labelN$SYMBOL), ]

volc <- ggplot(tt_gene, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = adj.P.Val < DEG_FDR_CUTOFF), alpha = 0.75, size = 1.8) +
  geom_vline(xintercept = c(-LFC_VISUAL_THR, LFC_VISUAL_THR), linetype = "dashed") +
  geom_hline(yintercept = -log10(DEG_FDR_CUTOFF), linetype = "dashed") +
  scale_color_manual(values = c("grey70","tomato"), labels = c("FDR ≥ cutoff","FDR < cutoff")) +
  labs(title = "Volcano (FDR-based)", x = expression(log[2]*"FC"), y = expression(-log[10]*"FDR"), color = NULL) +
  theme_bw(base_size = 12) +
  ggrepel::geom_text_repel(
    data = labelN,
    aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL),
    size = 3, max.overlaps = 200
  )
ggsave(file.path(plot_dir, "volcano_FDR_consistent.png"), volc, width = 8, height = 6, dpi = 300)

TOP_N <- min(50, nrow(deg))
if (TOP_N >= 2) {
  # temsilci probe = tt_gene’deki satır isimleri
  reps <- rownames(tt_gene)[match(head(deg$SYMBOL[order(deg$adj.P.Val)], TOP_N), tt_gene$SYMBOL)]
  # ifade matrisi: batch/SVA sonrası olanı tercih et
  expr_plot <- if (exists("expr_use")) expr_use else if (exists("expr_cb")) expr_cb else expr_f
  sub_expr  <- expr_plot[intersect(rownames(expr_plot), reps), , drop = FALSE]
  rownames(sub_expr) <- tt_gene$SYMBOL[match(rownames(sub_expr), rownames(tt_gene))]
  
  ann_col <- data.frame(Group = pheno_aligned$group, row.names = rownames(pheno_aligned))
  if ("batch" %in% names(pheno_aligned)) ann_col$Batch <- pheno_aligned$batch
  if ("is_melas" %in% names(pheno_aligned)) ann_col$MELAS <- ifelse(pheno_aligned$is_melas, "Yes","No")
  
  pheatmap(sub_expr, scale = "row",
           annotation_col = ann_col,
           show_rownames = TRUE, show_colnames = TRUE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           clustering_method = "average",
           main = sprintf("Top %d DEG (FDR<0.05)", TOP_N),
           filename = file.path(plot_dir, "heatmap_top_FDR.png"))
}

if (exists("tt_all_s") && exists("tt_gene_s") && exists("deg_sens")) {
  # Volcano no-MELAS
  labelNs <- head(tt_gene_s[order(tt_gene_s$adj.P.Val), c("SYMBOL","logFC","adj.P.Val")], 20)
  labelNs <- labelNs[!is.na(labelNs$SYMBOL) & nzchar(labelNs$SYMBOL), ]
  volc_s <- ggplot(tt_gene_s, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = adj.P.Val < DEG_FDR_CUTOFF), alpha = 0.75, size = 1.8) +
    geom_vline(xintercept = c(-LFC_VISUAL_THR, LFC_VISUAL_THR), linetype = "dashed") +
    geom_hline(yintercept = -log10(DEG_FDR_CUTOFF), linetype = "dashed") +
    scale_color_manual(values = c("grey70","tomato"), labels = c("FDR ≥ cutoff","FDR < cutoff")) +
    labs(title = "Volcano (no-MELAS)", x = expression(log[2]*"FC"), y = expression(-log[10]*"FDR"), color = NULL) +
    theme_bw(base_size = 12) +
    ggrepel::geom_text_repel(
      data = labelNs,
      aes(x = logFC, y = -log10(adj.P.Val), label = SYMBOL), size = 3, max.overlaps = 200)
  ggsave(file.path(plot_dir, "volcano_noMELAS.png"), volc_s, width = 8, height = 6, dpi = 300)
  
  TOP_Ns <- min(50, nrow(deg_sens))
  if (TOP_Ns >= 2 && exists("expr_sens") && exists("ph_sens")) {
    reps_s <- rownames(tt_gene_s)[match(head(deg_sens$SYMBOL[order(deg_sens$adj.P.Val)], TOP_Ns), tt_gene_s$SYMBOL)]
    sub_expr_s <- expr_sens[intersect(rownames(expr_sens), reps_s), , drop = FALSE]
    rownames(sub_expr_s) <- tt_gene_s$SYMBOL[match(rownames(sub_expr_s), rownames(tt_gene_s))]
    
    ann_col_s <- data.frame(Group = ph_sens$group, row.names = rownames(ph_sens))
    if ("batch" %in% names(ph_sens)) ann_col_s$Batch <- ph_sens$batch
    
    pheatmap(sub_expr_s, scale = "row",
             annotation_col = ann_col_s,
             show_rownames = TRUE, show_colnames = TRUE,
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             clustering_method = "average",
             main = sprintf("Top %d DEG (no-MELAS)", TOP_Ns),
             filename = file.path(plot_dir, "heatmap_top_FDR_noMELAS.png"))
  }
}

message("Volcano & Heatmap dosyaları yazıldı: ", plot_dir)

# ========================= HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH=========================

if (!requireNamespace("msigdbr", quietly=TRUE)) install.packages("msigdbr")
if (!requireNamespace("clusterProfiler", quietly=TRUE)) install.packages("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly=TRUE)) BiocManager::install("org.Hs.eg.db")

library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

gsea_dir <- file.path(root, "results", "deg")
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

msig <- try(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP"), silent=TRUE)
if (inherits(msig, "try-error") || nrow(msig)==0) {
  msig <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
}
term2gene <- unique(msig[, c("gs_name","gene_symbol")])

gsea_go <- GSEA(geneList = ranks,
                TERM2GENE = term2gene,
                pvalueCutoff = 0.05,
                verbose = FALSE)

write.csv(as.data.frame(gsea_go),
          file.path(gsea_dir, "GSEA_GO_BP_msigdb.csv"), row.names = FALSE)

png(file.path(gsea_dir, "GSEA_GO_BP_dotplot.png"), width=1000, height=800, res=130)
print(dotplot(gsea_go, showCategory = 20))
dev.off()

# En anlamlı 1–2 yolak için koşu grafiği
if (nrow(as.data.frame(gsea_go)) > 0) {
  top_terms <- head(gsea_go@result$ID, 2)
  for (id in top_terms) {
    png(file.path(gsea_dir, paste0("GSEA_GO_BP_runningplot_", gsub("[^A-Za-z0-9]+","_", id), ".png")),
        width=1100, height=700, res=130)
    print(gseaplot2(gsea_go, geneSetID = id))
    dev.off()
  }
}

sym2ent <- AnnotationDbi::select(org.Hs.eg.db,
                                 keys = unique(names(ranks)),
                                 keytype = "SYMBOL",
                                 columns = "ENTREZID")
sym2ent <- sym2ent[!is.na(sym2ent$ENTREZID) & nzchar(sym2ent$SYMBOL), ]

ranks_kegg <- ranks[names(ranks) %in% sym2ent$SYMBOL]
names(ranks_kegg) <- sym2ent$ENTREZID[match(names(ranks_kegg), sym2ent$SYMBOL)]
ranks_kegg <- sort(ranks_kegg, decreasing = TRUE)

if (length(ranks_kegg) >= 50) {
  gsea_kegg <- gseKEGG(geneList = ranks_kegg,
                       organism = "hsa",
                       pvalueCutoff = 0.05,
                       verbose = FALSE)
  write.csv(as.data.frame(gsea_kegg),
            file.path(gsea_dir, "GSEA_KEGG.csv"), row.names = FALSE)
  
  png(file.path(gsea_dir, "GSEA_KEGG_dotplot.png"), width=1000, height=800, res=130)
  print(dotplot(gsea_kegg, showCategory = 20))
  dev.off()
} else {
  message("KEGG GSEA atlandı: ENTREZ eşleşmesi az ( < 50 ).")
}

# ========================= IIIIIIIIIIIIIIIIIIIIIIIIIII =========================

if (exists("tt_gene_s")) {
  ranks_s <- tt_gene_s$t
  names(ranks_s) <- tt_gene_s$SYMBOL
  ranks_s <- ranks_s[!is.na(ranks_s) & nzchar(names(ranks_s))]
  ranks_s <- sort(ranks_s, decreasing = TRUE)
  
  if (length(ranks_s) >= 50) {
    gsea_go_s <- GSEA(geneList = ranks_s, TERM2GENE = term2gene,
                      pvalueCutoff = 0.05, verbose = FALSE)
    write.csv(as.data.frame(gsea_go_s),
              file.path(gsea_dir, "GSEA_GO_BP_msigdb_noMELAS.csv"), row.names = FALSE)
    
    png(file.path(gsea_dir, "GSEA_GO_BP_dotplot_noMELAS.png"), width=1000, height=800, res=130)
    print(dotplot(gsea_go_s, showCategory = 20))
    dev.off()
  } else {
    message("no-MELAS GSEA atlandı: yeterli rank yok.")
  }
}

if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra); library(grid)

if (exists("ph")) {
  pdf(file.path(root,"results","plots","volcano_heatmap_Fig1.pdf"), width=16, height=8, family="Times")
  grid.arrange(
    arrangeGrob(volc, top = textGrob("A", gp=gpar(fontsize=22, fontface="bold"), x=unit(0,"npc"), hjust=0)),
    arrangeGrob(ph$gtable, top = textGrob("B", gp=gpar(fontsize=22, fontface="bold"), x=unit(0,"npc"), hjust=0)),
    ncol = 2
  )
  dev.off()
}

# ======================= ENRICHMENT: GO / KEGG / GSEA + MELAS =======================
# Ön-koşullar: root, tt_all, tt_gene, deg, groups, expr_cb/expr_matrix
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("clusterProfiler","org.Hs.eg.db","enrichplot","DOSE","GOSemSim")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}
for (p in c("ggplot2","ggrepel","jsonlite")) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
library(DOSE); library(GOSemSim); library(ggplot2); library(ggrepel); library(jsonlite)

## 0) Klasörler
enrich_dir <- file.path(root, "results", "enrichment"); dir.create(enrich_dir, recursive = TRUE, showWarnings = FALSE)
plot_dir   <- file.path(root, "results", "plots");       dir.create(plot_dir,   recursive = TRUE, showWarnings = FALSE)

## 1) Girdi listeleri (main analiz)
stopifnot(exists("tt_gene"), exists("deg"))
sig_symbols      <- unique(na.omit(as.character(deg$SYMBOL)))
universe_symbols <- unique(na.omit(as.character(tt_gene$SYMBOL)))

# FAIR: evren/sig listelerini kaydet
writeLines(sig_symbols,      file.path(enrich_dir, "_sig_gene_symbols.txt"))
writeLines(universe_symbols, file.path(enrich_dir, "_universe_gene_symbols.txt"))

## 2) ID dönüşümü (Entrez)
sym2entrez      <- bitr(universe_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_universe <- unique(na.omit(sym2entrez$ENTREZID))
sig2entrez      <- bitr(sig_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
entrez_sig      <- unique(na.omit(sig2entrez$ENTREZID))

## 3) Provenans/param yaz
meta <- list(time=as.character(Sys.time()),
             n_sig_genes=length(sig_symbols), n_universe_genes=length(universe_symbols),
             n_sig_entrez=length(entrez_sig), n_universe_entrez=length(entrez_universe),
             adjust_method="BH", pvalueCutoff=0.05, qvalueCutoff=0.05)
writeLines(jsonlite::toJSON(meta, pretty=TRUE), file.path(enrich_dir, "_params_enrichment_main.json"))

## 4) GO: BP (SYMBOL + explicit universe)
go_bp <- enrichGO(gene = sig_symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "BP", universe = universe_symbols,
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                  readable = TRUE)
write.csv(as.data.frame(go_bp), file.path(enrich_dir,"GO_BP_main.csv"), row.names = FALSE)

## 5) KEGG (ENTREZ) – bazı sürümlerde universe parametresi yok → koşullu ekle
ek_args <- list(gene = entrez_sig, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)
if ("universe" %in% names(formals(clusterProfiler::enrichKEGG))) ek_args$universe <- entrez_universe
kegg <- do.call(clusterProfiler::enrichKEGG, ek_args)
write.csv(as.data.frame(kegg), file.path(enrich_dir,"KEGG_main.csv"), row.names = FALSE)

## 6) GSEA (eşik bağımlılığını azaltır) — geneList = t-istatistiği (daha güçlü)
# tt_all içinde t sütunu olmalı; yoksa logFC kullan
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

## 7) Görseller (statik, makale‑hazır)
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

## 8) MELAS duyarlılık zenginleştirme (varsa)
enrich_melas <- function(deg_tbl, tag) {
  sig_syms <- unique(na.omit(as.character(deg_tbl$SYMBOL)))
  if (length(sig_syms) < 3) return(invisible(NULL))
  sig_ent  <- unique(na.omit(bitr(sig_syms, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID))

  go <- enrichGO(gene=sig_syms, OrgDb=org.Hs.eg.db, keyType="SYMBOL",
                 ont="BP", universe=universe_symbols,
                 pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.05, readable=TRUE)
  write.csv(as.data.frame(go), file.path(enrich_dir, paste0("GO_BP_",tag,".csv")), row.names=FALSE)

  ek_args <- list(gene=sig_ent, organism="hsa", pAdjustMethod="BH", pvalueCutoff=0.05)
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
} else if (exists("pheno_aligned") && "is_melas" %in% names(pheno_aligned) && any(pheno_aligned$is_melas)) {
  # MELAS'lar hariç yeni DEG üret (ana ayarlarla tutarlı)
  keep_idx <- which(!pheno_aligned$is_melas)
  expr_s   <- expr_cb[, keep_idx, drop=FALSE]
  ph_s     <- droplevels(pheno_aligned[keep_idx, ])
  design_s <- model.matrix(~ 0 + group, data = ph_s); colnames(design_s) <- levels(ph_s$group)
  fit_s2   <- limma::eBayes(limma::contrasts.fit(limma::lmFit(expr_s, design_s),
                                                 limma::makeContrasts(Pompe - Control, levels = design_s)))
  tt_all_s <- limma::topTable(fit_s2, coef=1, number=Inf, adjust.method="BH")
  # gene-level tekilleştirme
  tt_all_s$abs_logFC <- abs(tt_all_s$logFC)
  tt_gene_s <- tt_all_s[!is.na(tt_all_s$SYMBOL) & nzchar(tt_all_s$SYMBOL), ]
  tt_gene_s <- tt_gene_s[order(tt_gene_s$SYMBOL, -tt_gene_s$abs_logFC), ]
  tt_gene_s <- tt_gene_s[!duplicated(tt_gene_s$SYMBOL), ]
  deg_sens  <- subset(tt_gene_s, adj.P.Val < 0.05)
  write.csv(deg_sens, file.path(enrich_dir,"DEG_noMELAS_FDR_lt_0.05.csv"), row.names=FALSE)
  invisible(enrich_melas(deg_sens, "noMELAS"))
}

## 9) Kısa özet
sum_lines <- c(
  sprintf("Main GO terms (q<0.05): %d", ifelse(is.null(go_bp), 0, nrow(as.data.frame(go_bp)))),
  sprintf("Main KEGG pathways (q<0.05): %d", ifelse(is.null(kegg), 0, nrow(as.data.frame(kegg)))),
  sprintf("Universe size: SYMBOL=%d | ENTREZ=%d", length(universe_symbols), length(entrez_universe)),
  "Notes: BH correction everywhere; GO uses explicit universe (all tested gene-level).",
  "Sensitivity: If MELAS present, enrichment repeated without MELAS; GSEA added to reduce threshold dependence."
)
writeLines(sum_lines, file.path(enrich_dir,"_ENRICHMENT_SUMMARY.txt"))
# =============================================================================== 



# ==================== miRNA Target Enrichment (FAIR, offline) ====================
# Koşullar: tt_gene (gene-level tablo) ve deg (FDR<0.05) mevcut olmalı
# Kullanılan nesneler: tt_gene, deg, pheno_aligned (MELAS için), groups
# Çıktılar: results/miRNA_analysis/ altına CSV/PNG/JSON/TXT

# ---- Paketler
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
need_bioc <- c("multiMiR")
for (p in need_bioc) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE)
need_cran <- c("dplyr","ggplot2","ggrepel","jsonlite")
for (p in need_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
lapply(c("multiMiR","dplyr","ggplot2","ggrepel","jsonlite"), library, character.only = TRUE)

# ---- Klasörler + meta
root_dir <- if (exists("root")) root else getwd()
out_dir  <- file.path(root_dir, "results", "miRNA_analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20240724)

# ---- Girdiler (FAIR kontroller)
stopifnot(exists("tt_gene"), "SYMBOL" %in% names(tt_gene))
stopifnot(exists("deg"),     "SYMBOL" %in% names(deg))
universe_symbols <- unique(tt_gene$SYMBOL[!is.na(tt_gene$SYMBOL) & nzchar(tt_gene$SYMBOL)])
sig_genes        <- unique(deg$SYMBOL[!is.na(deg$SYMBOL) & nzchar(deg$SYMBOL)])

meta <- list(
  dataset              = "GSE38680",
  n_control            = if (exists("groups")) sum(groups=="Control") else NA_integer_,
  n_pompe              = if (exists("groups")) sum(groups=="Pompe") else NA_integer_,
  n_universe_genes     = length(universe_symbols),
  n_sig_genes_FDR_lt05 = length(sig_genes),
  adjust_method        = "BH",
  note_small_n         = "Small sample size; interpret enrichment cautiously."
)
jsonlite::write_json(meta, file.path(out_dir, "_mirna_meta.json"), pretty = TRUE)

# ---- Yardımcı: miRNA hedef zenginleştirme (Fisher + BH)
enrich_mirna <- function(gene_set, universe_set, which_table=c("validated","predicted")) {
  which_table <- match.arg(which_table)
  # multiMiR'den hedef map (yerel kaynaklar)
  mm <- tryCatch(
    multiMiR::get_multimir(org="hsa", target = universe_set, table = which_table, summary = FALSE),
    error = function(e) NULL
  )
  if (is.null(mm) || nrow(mm@data)==0) return(NULL)

  df <- mm@data
  # multiMiR sütunları değişebilir: hedef sembol için esnek alan seçelim
  target_col <- c("target_symbol","target.gene","target")[
    c("target_symbol","target.gene","target") %in% names(df)][1]
  if (is.na(target_col)) return(NULL)

  df <- df[!is.na(df[[target_col]]) & nzchar(df[[target_col]]), c("mature_mirna_id", target_col, "database")]
  colnames(df) <- c("miRNA","SYMBOL","database")

  # Evren dışındaki satırları at
  df <- df[df$SYMBOL %in% universe_set, , drop = FALSE]
  if (!nrow(df)) return(NULL)

  # Her miRNA için kontenjans tablosu
  #            |  DEG (sig)  |  non-DEG  |
  # --------------------------------------
  # target     |    a        |    b      |
  # notarget   |    c        |    d      |
  # a = miRNA'nın hedeflediği DEGsayısı; b = miRNA hedefleri - a; c = (DEG - a); d = (Universe - a - b - c)
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

# ---- Validated ve Predicted için çalıştır
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

# ---- Ağ verisi (top miRNA → hedef DEGs) – referans için CSV
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

# ---- MELAS duyarlılık: varsa, aynı zenginleştirmeyi tekrarla ve karşılaştır
sig_genes_s <- NULL
if (exists("pheno_aligned") && (("is_melas" %in% names(pheno_aligned)) || ("specialcase" %in% names(pheno_aligned)))) {
  mel_flag <- if ("is_melas" %in% names(pheno_aligned)) pheno_aligned$is_melas else pheno_aligned$specialcase
  mel_flag <- tolower(as.character(mel_flag)) %in% c("1","true","t","yes","y")
  if (any(mel_flag)) {
    # DEG_sens zaten hesaplandıysa kullan; yoksa hızlıca hesapla
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
        # gene-level tekilleştirme
        tt_s$SYMBOL <- tt_s$SYMBOL %||% NA_character_
        if (!"SYMBOL" %in% names(tt_s)) {
          # eğer anotasyon nesneleri hazırsa tt_all/ann'den eşleştirilebilirdi; basitçe geçelim
        }
        if ("SYMBOL" %in% names(tt_s)) {
          tt_s$abs_logFC <- abs(tt_s$logFC)
          tt_sg <- tt_s[!is.na(tt_s$SYMBOL) & nzchar(tt_s$SYMBOL), ]
          tt_sg <- tt_sg[order(tt_sg$SYMBOL, -tt_sg$abs_logFC), ]
          tt_sg <- tt_sg[!duplicated(tt_sg$SYMBOL), ]
          sig_genes_s <- unique(tt_sg$SYMBOL[tt_sg$adj.P.Val < 0.05])
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
    # Kısa karşılaştırma özeti
    top_main <- head(enr_main$miRNA[order(enr_main$FDR)], 10)
    top_sens <- head(enr_sens$miRNA[order(enr_sens$FDR)], 10)
    ovlp     <- intersect(top_main, top_sens)
    writeLines(c(
      sprintf("Top10 overlap (validated+predicted): %d", length(ovlp)),
      paste("Common:", paste(ovlp, collapse=", "))
    ), file.path(out_dir,"_sensitivity_overlap.txt"))
  }
}

# ---- Oturum bilgisi
writeLines(c(capture.output(sessionInfo())), file.path(out_dir,"session_info.txt"))

message("miRNA enrichment (FAIR/offline) tamamlandı. Çıktılar: ", out_dir)


