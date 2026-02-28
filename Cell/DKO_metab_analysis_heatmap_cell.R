# =========================
# Derlin KO: end-to-end run
# =========================
setwd("~/Dropbox/Derlin KO Metabolomics")

# ---- packages ----
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(tibble)
library(UpSetR)
library(ggrepel)
library(pheatmap)
library(svglite)  # for SVG output
library(grid)     # for grid.draw

# ---- output dir ----
out_dir <- "figures_svg"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
svg_path <- function(filename) file.path(out_dir, filename)

# Helper: save a pheatmap object to SVG
save_pheatmap_svg <- function(ph, file, width = 8, height = 6) {
  svglite::svglite(file, width = width, height = height)
  on.exit(dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  invisible(file)
}

# =================
# read in the data
# =================
data <- read_tsv(
  "Derlin_KO_Data_Clean.txt",
  quote = "",
  comment = "",
  name_repair = "minimal"
)

# subset to only have the cellular extract data
data <- data %>% select(1, contains("Cell", ignore.case = FALSE))

# ================
# data processing
# ================
cell_cols <- names(data)[-1]                 # get the cols we need
data[cell_cols] <- lapply(data[cell_cols], as.double) # ensure numeric

# create matrix with data
X <- as.matrix(data[, cell_cols, drop = FALSE])

# replace negative values with NA
X[X < 0] <- NA_real_

# get half-min for each row
hmin <- apply(X, 1, function(v) {
  pv <- v[!is.na(v) & v > 0]
  if (length(pv)) 0.5 * min(pv) else NA_real_
})

# keep only rows where half-min is defined (>=1 positive)
keep <- !is.na(hmin)
X    <- X[keep, , drop = FALSE]
hmin <- hmin[keep]

# replace zeros with half-min
zero_idx <- which(X == 0, arr.ind = TRUE)
if (nrow(zero_idx)) X[zero_idx] <- hmin[zero_idx[, 1]]

# replace NA with half-min
na_idx <- which(is.na(X), arr.ind = TRUE)
if (nrow(na_idx)) X[na_idx] <- hmin[na_idx[, 1]]

# log2 transform
X_log2 <- log2(X)

# put it all together
data_log2 <- cbind(data[keep, 1, drop = FALSE], as.data.frame(X_log2))

# deal with duplicates: keep the one with highest row mean
data_log2 <- data_log2 %>%
  mutate(.mean = rowMeans(across(-1), na.rm = TRUE)) %>%
  arrange(compound, desc(.mean)) %>%
  distinct(compound, .keep_all = TRUE) %>%
  select(-.mean)

# ===========
# sample info
# ===========
samples <- colnames(data_log2)[-1]

# Extract the line prefix at the start of each name: D1, D2, D3, DT, D_WT, DWT
grp_raw <- sub("^((D_WT|DWT|DT|D1|D2|D3)).*", "\\1", samples, perl = TRUE)
grp <- gsub("^(D_WT|DWT)$", "WT", grp_raw) # make D_WT/DWT -> WT

sample_info <- data.frame(
  Sample = samples,
  Group  = factor(grp, levels = c("WT", "D1", "D2", "D3", "DT")),
  stringsAsFactors = FALSE
)

# ==========================
# differential (limma/ebayes)
# ==========================
expr <- as.matrix(data_log2[, -1, drop = FALSE])
rownames(expr) <- data_log2[[1]]

Group <- sample_info$Group
design <- model.matrix(~ 0 + Group)
colnames(design) <- levels(Group)

fit  <- lmFit(expr, design)
contr <- makeContrasts(
  D1vsWT = D1 - WT,
  D2vsWT = D2 - WT,
  D3vsWT = D3 - WT,
  DTvsWT = DT - WT,
  KOvsWT = (D1 + D2 + D3 + DT)/4 - WT,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, contr))

res_D1 <- topTable(fit2, coef = "D1vsWT", number = Inf, sort.by = "P")
res_D2 <- topTable(fit2, coef = "D2vsWT", number = Inf, sort.by = "P")
res_D3 <- topTable(fit2, coef = "D3vsWT", number = Inf, sort.by = "P")
res_DT <- topTable(fit2, coef = "DTvsWT", number = Inf, sort.by = "P")
res_KO <- topTable(fit2, coef = "KOvsWT", number = Inf, sort.by = "P")

add_id <- function(df, name) rownames_to_column(df, "ID") %>%
  dplyr::mutate(contrast = name, .before = 1)

results_long <- dplyr::bind_rows(
  add_id(res_D1, "D1vsWT"),
  add_id(res_D2, "D2vsWT"),
  add_id(res_D3, "D3vsWT"),
  add_id(res_DT, "DTvsWT"),
  add_id(res_KO, "KOvsWT")
)

# thresholds / order
ord    <- rev(c("D1vsWT","D2vsWT","D3vsWT","DTvsWT","KOvsWT"))
alpha  <- 0.05
lfc_min <- 0

# ================================
# Heatmap: all significant overall
# ================================
sig_ids <- results_long %>%
  filter(adj.P.Val < alpha, abs(logFC) >= lfc_min) %>%
  pull(ID) %>% unique()

if (length(sig_ids) == 0) stop("No significant metabolites at the chosen thresholds.")

sub <- data_log2[data_log2[[1]] %in% sig_ids, ]
mat <- as.matrix(sub[, -1, drop = FALSE])
rownames(mat) <- sub[[1]]

# order columns by Group (WT, D1, D2, D3, DT)
sample_order <- sample_info$Sample[order(sample_info$Group)]
mat <- mat[, sample_order, drop = FALSE]

# row Z-scores; drop zero-variance rows
keep_rows <- apply(mat, 1, sd, na.rm = TRUE) > 0
mat_z <- t(scale(t(mat[keep_rows, , drop = FALSE])))

# column annotation
ann_col <- data.frame(Group = droplevels(sample_info$Group[match(colnames(mat_z), sample_info$Sample)]))
rownames(ann_col) <- colnames(mat_z)

# draw & save SVG
ph_all <- pheatmap(
  mat_z,
  annotation_col = ann_col,
  show_rownames = TRUE,      # set FALSE if too many rows
  show_colnames = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = sprintf("All significant metabolites (FDR < %.2f%s)",
                 alpha, if (lfc_min > 0) paste0(", |log2FC| ≥ ", lfc_min) else ""),
  silent = TRUE
)
save_pheatmap_svg(ph_all, svg_path("heatmap_all_significant_cell.svg"), width = 9, height = 10)

# ==========================================
# Per-contrast heatmaps (all groups displayed)
# ==========================================
results_long$contrast <- factor(as.character(results_long$contrast), levels = ord)

plot_contrast_heatmap_allgroups <- function(cn, show_rownames = TRUE) {
  hits <- results_long %>%
    filter(contrast == cn, adj.P.Val < alpha, abs(logFC) >= lfc_min)
  
  ids <- unique(hits$ID)
  if (length(ids) == 0) {
    message("No significant metabolites for ", cn)
    return(invisible(NULL))
  }
  
  sub <- data_log2[data_log2[[1]] %in% ids, ]
  mat <- as.matrix(sub[, -1, drop = FALSE])
  rownames(mat) <- sub[[1]]
  
  # order columns by full Group
  sample_order <- sample_info$Sample[order(sample_info$Group)]
  mat <- mat[, sample_order, drop = FALSE]
  
  # row Z-scores; drop zero-variance rows
  keep_rows <- apply(mat, 1, sd, na.rm = TRUE) > 0
  mat_z <- t(scale(t(mat[keep_rows, , drop = FALSE])))
  
  if (nrow(mat_z) == 0) {
    message("All significant metabolites for ", cn, " had zero variance; skipping.")
    return(invisible(NULL))
  }
  
  # annotations
  ann_col <- data.frame(Group = droplevels(sample_info$Group[match(colnames(mat_z), sample_info$Sample)]))
  rownames(ann_col) <- colnames(mat_z)
  
  dir_tbl <- hits %>%
    select(ID, logFC) %>%
    distinct() %>%
    mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%
    tibble::column_to_rownames("ID")
  ann_row <- data.frame(Direction = dir_tbl[rownames(mat_z), "Direction", drop = TRUE])
  rownames(ann_row) <- rownames(mat_z)
  
  cl_rows <- nrow(mat_z) > 1
  cl_cols <- ncol(mat_z) > 1
  
  ph <- pheatmap(
    mat_z,
    annotation_col = ann_col,
    annotation_row = ann_row,
    annotation_names_row = FALSE,
    show_rownames = show_rownames,
    show_colnames = FALSE,
    cluster_rows = cl_rows,
    cluster_cols = cl_cols,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = sprintf("%s: significant metabolites (FDR < %.2f%s)",
                   cn, alpha, if (lfc_min > 0) paste0(", |log2FC| ≥ ", lfc_min) else ""),
    silent = TRUE
  )
  
  save_pheatmap_svg(ph, svg_path(paste0("heatmap_", cn, "_all-groups_cell.svg")), width = 9, height = 10)
  invisible(NULL)
}

for (cn in ord) plot_contrast_heatmap_allgroups(cn)

# ===================================================
# Per-contrast heatmaps (only the relevant groups)
# ===================================================
groups_for <- function(cn) switch(
  cn,
  "D1vsWT" = c("WT","D1"),
  "D2vsWT" = c("WT","D2"),
  "D3vsWT" = c("WT","D3"),
  "DTvsWT" = c("WT","DT"),
  "KOvsWT" = c("WT","D1","D2","D3","DT"),
  c("WT","D1","D2","D3","DT") # fallback
)

plot_contrast_heatmap_relevant <- function(cn, show_rownames = TRUE) {
  hits <- results_long %>%
    filter(contrast == cn, adj.P.Val < alpha, abs(logFC) >= lfc_min)
  
  ids <- unique(hits$ID)
  if (length(ids) == 0) {
    message("No significant metabolites for ", cn)
    return(invisible(NULL))
  }
  
  sub <- data_log2[data_log2[[1]] %in% ids, ]
  mat <- as.matrix(sub[, -1, drop = FALSE])
  rownames(mat) <- sub[[1]]
  
  keep_groups <- groups_for(cn)
  keep_samples <- sample_info %>%
    filter(Group %in% keep_groups) %>%
    mutate(Group = factor(Group, levels = keep_groups)) %>%
    arrange(Group) %>%
    pull(Sample)
  
  mat <- mat[, keep_samples, drop = FALSE]
  
  keep_rows <- apply(mat, 1, sd, na.rm = TRUE) > 0
  mat_z <- t(scale(t(mat[keep_rows, , drop = FALSE])))
  
  if (nrow(mat_z) == 0) {
    message("All significant metabolites for ", cn, " had zero variance; skipping.")
    return(invisible(NULL))
  }
  
  ann_col <- data.frame(Group = droplevels(sample_info$Group[match(colnames(mat_z), sample_info$Sample)]))
  rownames(ann_col) <- colnames(mat_z)
  
  dir_tbl <- hits %>%
    select(ID, logFC) %>%
    distinct() %>%
    mutate(Direction = ifelse(logFC > 0, "Up", "Down")) %>%
    tibble::column_to_rownames("ID")
  ann_row <- data.frame(Direction = dir_tbl[rownames(mat_z), "Direction", drop = TRUE])
  rownames(ann_row) <- rownames(mat_z)
  
  cl_rows <- nrow(mat_z) > 1
  cl_cols <- ncol(mat_z) > 1
  
  ph <- pheatmap(
    mat_z,
    annotation_col = ann_col,
    annotation_row = ann_row,
    annotation_names_row = FALSE,
    show_rownames = show_rownames,
    show_colnames = FALSE,
    cluster_rows = cl_rows,
    cluster_cols = cl_cols,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    main = sprintf("%s: significant metabolites (FDR < %.2f%s)",
                   cn, alpha, if (lfc_min > 0) paste0(", |log2FC| ≥ ", lfc_min) else ""),
    silent = TRUE
  )
  
  save_pheatmap_svg(ph, svg_path(paste0("heatmap_", cn, "_relevant-groups_cell.svg")), width = 7.5, height = 8.5)
  invisible(NULL)
}

for (cn in ord) plot_contrast_heatmap_relevant(cn)

# =========
# Finished!
# =========
message("SVGs written to: ", normalizePath(out_dir))

