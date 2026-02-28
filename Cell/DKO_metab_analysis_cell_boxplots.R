


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
library(svglite)
library(ggh4x)


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


#Build sample metadata from column names (e.g., "D1_Cell_1", "D_WT_Cell_2", etc.)
samples <- colnames(data_log2)[-1]

# Extract the line prefix at the start of each name: D1, D2, D3, DT, D_WT 
grp_raw <- sub("^((D_WT|DWT|DT|D1|D2|D3)).*", "\\1", samples, perl = TRUE)
grp <- gsub("^(D_WT|DWT)$", "WT", grp_raw) #make D_WT WT

sample_info <- data.frame(
  Sample = samples,
  Group  = factor(grp, levels = c("WT", "D1", "D2", "D3", "DT")),
  stringsAsFactors = FALSE
)

# NAD+ salvage pathway plot
targets <- c("NAD+", "NADH", "nicotinamide mononucleotide",
             "nicotinamide", "methylnicotinamide")

# Subset rows for targets (case-insensitive)
comp_lower <- tolower(data_log2$compound)
pos <- match(tolower(targets), comp_lower)
if (anyNA(pos)) message("Not found in data: ", paste(targets[is.na(pos)], collapse = ", "))

targets_found <- targets[!is.na(pos)]
plot_rows <- data_log2[pos[!is.na(pos)], , drop = FALSE] %>%
  dplyr::rename(Metabolite = compound)

# Long format + join group labels
long_df <- plot_rows %>%
  tidyr::pivot_longer(cols = -Metabolite, names_to = "Sample", values_to = "Abundance_log2") %>%
  dplyr::left_join(sample_info, by = "Sample") %>%
  dplyr::mutate(
    Group = factor(Group, levels = c("WT","D1","D2","D3","DT")),
    Metabolite = factor(Metabolite, levels = targets_found)  # facet order
  )

# Okabe–Ito palette
pal <- c(WT="#000000", D1="#E69F00", D2="#56B4E9", D3="#009E73", DT="#D55E00")

# Build plot
p <- ggplot(long_df, aes(x = Group, y = Abundance_log2)) +
  geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.7, alpha = 0.35) +
  geom_point(aes(color = Group),
             position = position_jitter(width = 0.12, height = 0, seed = 1),
             size = 2, alpha = 0.9) +
  scale_fill_manual(values = pal, name = "Group") +
  scale_color_manual(values = pal, guide = "none") +
  ggh4x::facet_wrap2(~ Metabolite, ncol = 3, scales = "free_y", axes = "all") +
  labs(title = "NAD+ Salvage Pathway", x = NULL, y = "log2 abundance") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(vjust = 0.9)
  )

# Show plot in RStudio
print(p)

# Save as SVG with "MEDIA" in the file name
ggsave("NAD_Salvage_Pathway_CELL.svg", plot = p, device = "svg", width = 10, height = 7, dpi = 300)
