

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
library(svglite)

# ---- read + prep ----
data <- read_tsv(
  "Derlin_KO_Data_Clean.txt",
  quote = "",
  comment = "",
  name_repair = "minimal"
)

# media only
data <- data %>% select(1, contains("Media", ignore.case = FALSE))

media_cols <- names(data)[-1]
data[media_cols] <- lapply(data[media_cols], as.double)

X <- as.matrix(data[, media_cols, drop = FALSE])
X[X < 0] <- NA_real_

hmin <- apply(X, 1, function(v) {
  pv <- v[!is.na(v) & v > 0]
  if (length(pv)) 0.5 * min(pv) else NA_real_
})

keep <- !is.na(hmin)
X    <- X[keep, , drop = FALSE]
hmin <- hmin[keep]

zero_idx <- which(X == 0, arr.ind = TRUE)
if (nrow(zero_idx)) X[zero_idx] <- hmin[zero_idx[, 1]]

na_idx <- which(is.na(X), arr.ind = TRUE)
if (nrow(na_idx)) X[na_idx] <- hmin[na_idx[, 1]]

X_log2 <- log2(X)

data_log2 <- cbind(data[keep, 1, drop = FALSE], as.data.frame(X_log2))

data_log2 <- data_log2 %>%
  mutate(.mean = rowMeans(across(-1), na.rm = TRUE)) %>%
  arrange(compound, desc(.mean)) %>%
  distinct(compound, .keep_all = TRUE) %>%
  select(-.mean)

# ---- PCA ----
samples <- colnames(data_log2)[-1]
grp_raw <- sub("^((D_WT|DWT|DT|D1|D2|D3)).*", "\\1", samples, perl = TRUE)
grp <- gsub("^(D_WT|DWT)$", "WT", grp_raw)

sample_info <- data.frame(
  Sample = samples,
  Group  = factor(grp, levels = c("WT", "D1", "D2", "D3", "DT")),
  stringsAsFactors = FALSE
)

mat <- t(as.matrix(data_log2[, -1, drop = FALSE]))

if (ncol(mat) > 1) {
  keep_feat <- apply(mat, 2, var) > 0
  mat <- mat[, keep_feat, drop = FALSE]
}

pca <- prcomp(mat, center = TRUE, scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)

scores <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  row.names = NULL
)
scores <- merge(scores, sample_info, by = "Sample", sort = FALSE)

p_pca <- ggplot(scores, aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) +
  labs(
    x = paste0("PC1 (", round(100 * var_expl[1], 1), "%)"),
    y = paste0("PC2 (", round(100 * var_expl[2], 1), "%)"),
    color = "Line"
  ) +
  theme_classic(base_size = 14)
print(p_pca)
ggsave("PCA_Media.svg", plot = p_pca, device = "svg", width = 8, height = 6, dpi = 300)

# ---- Differential analysis (limma) ----
expr <- as.matrix(data_log2[, -1, drop = FALSE])
rownames(expr) <- data_log2[[1]]
Group <- sample_info$Group

design <- model.matrix(~ 0 + Group)
colnames(design) <- levels(Group)

fit <- lmFit(expr, design)

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

# ---- UpSet plots (save as SVG) ----
ord <- rev(c("D1vsWT","D2vsWT","D3vsWT","DTvsWT","KOvsWT"))
alpha <- 0.05; lfc_min <- 0

ids_for <- function(cn, dir = c("up","down")) {
  dir <- match.arg(dir)
  ix <- results_long$contrast == cn &
    results_long$adj.P.Val < alpha &
    abs(results_long$logFC) >= lfc_min
  ix <- if (dir == "up") ix & results_long$logFC > 0 else ix & results_long$logFC < 0
  unique(results_long$ID[ix])
}

sets_up   <- setNames(lapply(ord, \(cn) ids_for(cn, "up")),   ord)
sets_down <- setNames(lapply(ord, \(cn) ids_for(cn, "down")), ord)

# Upregulated
svglite::svglite("UpSet_Upregulated_Media.svg", width = 10, height = 6, pointsize = 12)
upset(UpSetR::fromList(sets_up), keep.order = TRUE, order.by = "freq",
      mainbar.y.label = "Intersection size (Up)",
      sets.x.label    = "Set size (Up)")
dev.off()

# Downregulated
svglite::svglite("UpSet_Downregulated_Media.svg", width = 10, height = 6, pointsize = 12)
upset(UpSetR::fromList(sets_down), keep.order = TRUE, order.by = "freq",
      mainbar.y.label = "Intersection size (Down)",
      sets.x.label    = "Set size (Down)")
dev.off()

# ---- Volcano (faceted) ----
volc_df <- results_long %>%
  transmute(
    contrast = factor(contrast, levels = ord),
    ID, logFC, P.Value, adj.P.Val,
    neglog10P = -log10(P.Value),
    sig = case_when(
      adj.P.Val < alpha & logFC >=  lfc_min ~ "Up",
      adj.P.Val < alpha & logFC <= -lfc_min ~ "Down",
      TRUE                                  ~ "NS"
    )
  )

n_label <- 10
label_df <- volc_df %>%
  filter(sig != "NS") %>%
  group_by(contrast) %>%
  slice_max(order_by = neglog10P, n = n_label, with_ties = FALSE) %>%
  ungroup()

p_volc <- ggplot(volc_df, aes(x = logFC, y = neglog10P, color = sig)) +
  geom_point(size = 1.6, alpha = 0.85) +
  geom_vline(xintercept = c(-lfc_min, lfc_min), linetype = "dashed") +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
  facet_wrap(~ contrast, ncol = 3) +
  labs(
    x = "log2 fold-change (KO − WT)",
    y = expression(-log[10](p)),
    color = ""
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")
print(p_volc)
ggsave("Volcano_AllContrasts_Media.svg", plot = p_volc, device = "svg",
       width = 10, height = 8, dpi = 300)
