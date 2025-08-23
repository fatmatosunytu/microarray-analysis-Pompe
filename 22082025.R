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

