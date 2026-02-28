

#load packages
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(limma)
library(tibble)
library(UpSetR)
library(ggrepel)
library(pheatmap)
library(purrr)
library(fgsea)

#read in data
data <- read_tsv(
  "Derlin_KO_Data_Clean.txt",
  quote = "",
  comment = "",
  name_repair = "minimal"
)

#subset to only have the cellular extract data
data <- data %>% select(1, contains("Cell", ignore.case = FALSE))

### Data Processing ###
cell_cols <- names(data)[-1] #get the cols we need

data[cell_cols] <- lapply(data[cell_cols], as.double) #make sure things are numbers

# create matrix with data
X <- as.matrix(data[, cell_cols, drop = FALSE])

# replace neg values with NA
X[X < 0] <- NA_real_

# get half min for each row
hmin <- apply(X, 1, function(v) {
  pv <- v[!is.na(v) & v > 0]
  if (length(pv)) 0.5 * min(pv) else NA_real_
})

# keep only rows where half-min is defined (i.e., at least one positive value exists)
keep <- !is.na(hmin)
X    <- X[keep, , drop = FALSE]
hmin <- hmin[keep]

# replace 0s with half min for that row
zero_idx <- which(X == 0, arr.ind = TRUE)
if (nrow(zero_idx)) X[zero_idx] <- hmin[zero_idx[, 1]]

# replace NA with half min
na_idx <- which(is.na(X), arr.ind = TRUE)
if (nrow(na_idx)) X[na_idx] <- hmin[na_idx[, 1]]

# log transform
X_log2 <- log2(X)

#put it all together
data_log2 <- cbind(data[keep, 1, drop = FALSE], as.data.frame(X_log2))

#deal with duplicates
data_log2 <- data_log2 %>%
  mutate(.mean = rowMeans(across(-1), na.rm = TRUE)) %>%
  arrange(compound, desc(.mean)) %>%
  distinct(compound, .keep_all = TRUE) %>%
  select(-.mean)

# Extract the line prefix at the start of each name: D1, D2, D3, DT, D_WT 
samples <- colnames(data_log2)[-1]
grp_raw <- sub("^((D_WT|DWT|DT|D1|D2|D3)).*", "\\1", samples, perl = TRUE)
grp <- gsub("^(D_WT|DWT)$", "WT", grp_raw) #make D_WT WT

sample_info <- data.frame(
  Sample = samples,
  Group  = factor(grp, levels = c("WT", "D1", "D2", "D3", "DT")),
  stringsAsFactors = FALSE
)


#### Differentiation Analysis ####
# Expression matrix: rows = metabolites, cols = samples
expr <- as.matrix(data_log2[, -1, drop = FALSE])
rownames(expr) <- data_log2[[1]]

# Use the Group you already parsed for PCA
Group <- sample_info$Group  # factor with levels WT, D1, D2, D3, DT

# Design (no intercept so each group is a column)
design <- model.matrix(~ 0 + Group)
colnames(design) <- levels(Group)

# Fit per-metabolite linear model and apply empirical-Bayes moderation
fit <- lmFit(expr, design)

# Contrasts: KOs vs WT, and overall KO mean vs WT
contr <- makeContrasts(
  D1vsWT = D1 - WT,
  D2vsWT = D2 - WT,
  D3vsWT = D3 - WT,
  DTvsWT = DT - WT,
  KOvsWT = (D1 + D2 + D3 + DT)/4 - WT,   # equal-weight average of KOs
  levels = design
)

fit2 <- eBayes(contrasts.fit(fit, contr))

# Results tables (logFC = KO - WT on log2 scale)
res_D1 <- topTable(fit2, coef = "D1vsWT", number = Inf, sort.by = "P")
res_D2 <- topTable(fit2, coef = "D2vsWT", number = Inf, sort.by = "P")
res_D3 <- topTable(fit2, coef = "D3vsWT", number = Inf, sort.by = "P")
res_DT <- topTable(fit2, coef = "DTvsWT", number = Inf, sort.by = "P")
res_KO <- topTable(fit2, coef = "KOvsWT", number = Inf, sort.by = "P")

# Combine into one long table with a contrast label
tag <- function(df, name) df %>% mutate(contrast = name, .before = 1)

add_id <- function(df, name) rownames_to_column(df, "ID") %>%
  dplyr::mutate(contrast = name, .before = 1)

results_long <- dplyr::bind_rows(
  add_id(res_D1, "D1vsWT"),
  add_id(res_D2, "D2vsWT"),
  add_id(res_D3, "D3vsWT"),
  add_id(res_DT, "DTvsWT"),
  add_id(res_KO, "KOvsWT")
)


#name conversion
metab_analyst_names <- read.csv("Covnerted_metabolite_all.csv", header = T)[,c(1,6)]

names_map <- metab_analyst_names %>% distinct(ID, KEGG)

results_long <- results_long %>%
  left_join(names_map, by = "ID") %>%
  mutate(ID = coalesce(KEGG, ID)) %>%  # replace with KEGG if present
  select(-KEGG)

#MSEA
# ---- Make MetaboAnalyst GUI QEA files (one per contrast) ----
# Assumes you already have: data_log2 (compound + sample columns), sample_info (Sample, Group)

# 1) Matrix: rows = samples, cols = metabolites
met_names <- data_log2$compound
mat <- t(as.matrix(data_log2[, -1, drop = FALSE]))
colnames(mat) <- met_names

# 2) Align sample metadata
samples_df <- data.frame(Sample = rownames(mat), stringsAsFactors = FALSE) |>
  dplyr::left_join(sample_info, by = "Sample")  # adds Group

# 3) Helper to write one contrast in GUI format
write_qea_contrast <- function(contrast, out_file) {
  if (contrast == "KOvsWT") {
    keep <- samples_df$Group %in% c("WT","D1","D2","D3","DT")
    cls  <- ifelse(samples_df$Group[keep] == "WT", "WT", "KO")
  } else {
    g <- sub("vsWT$", "", contrast)                 # e.g., "D1"
    keep <- samples_df$Group %in% c("WT", g)
    cls  <- ifelse(samples_df$Group[keep] == g, g, "WT")
  }
  X <- mat[keep, , drop = FALSE]
  out <- data.frame(Sample = rownames(X), Class = cls, X, check.names = FALSE)
  write.csv(out, out_file, row.names = FALSE)
}

# 4) Write all 5 files
dir.create("metaboanalyst_QEA_Cell", showWarnings = FALSE)
contrasts <- c("D1vsWT","D2vsWT","D3vsWT","DTvsWT","KOvsWT")
for (cn in contrasts) {
  write_qea_contrast(cn, file.path("metaboanalyst_QEA_Cell", paste0("QEA_", cn, ".csv")))
}


