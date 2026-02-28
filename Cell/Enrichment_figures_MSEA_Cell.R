setwd("~/Dropbox/Derlin KO Metabolomics/metaboanalyst_QEA_Cell/")

# ---- Plot MetaboAnalyst QEA results for all contrasts (top 30) ----
library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(tools)
library(svglite)

# Settings
pattern <- "^msea_qea_result_.*\\.csv$"
topN    <- 30

# Helper: read + tidy one QEA CSV
read_qea <- function(infile) {
  df <- read_csv(infile, show_col_types = FALSE)
  
  # First column is the pathway name
  if ("...1" %in% names(df)) {
    df <- df %>% rename(Pathway = `...1`)
  } else {
    names(df)[1] <- "Pathway"
  }
  
  # Numeric columns we expect
  num_cols <- c("Total Cmpd","Hits","Statistic Q","Expected Q","Raw p","Holm p","FDR")
  present_num <- intersect(num_cols, names(df))
  df <- df %>% mutate(across(all_of(present_num), ~ suppressWarnings(as.numeric(.x))))
  
  df %>%
    mutate(
      hit_ratio   = Hits / `Total Cmpd`,
      neglog10FDR = -log10(pmax(FDR, 1e-300))
    ) %>%
    arrange(FDR, desc(hit_ratio))
}

# Helper: make, print, and save the plot
plot_qea <- function(df, infile, topN = 30) {
  df_top <- df %>% slice_head(n = topN)
  if (nrow(df_top) == 0) {
    message(sprintf("No rows to plot for %s", infile))
    return(invisible(NULL))
  }
  
  # Clean title: e.g. msea_qea_result_D1vsWT.csv -> "MSEA D1 vs WT"
  contrast <- str_replace(basename(infile), "msea_qea_result_", "")
  contrast <- str_replace(contrast, "\\.csv$", "")
  contrast <- str_replace_all(contrast, "(?<=\\D)(vs)(?=\\D)", " vs ")
  title_txt <- paste("MSEA", contrast)
  
  # Create plot
  p <- ggplot(df_top, aes(x = neglog10FDR, y = reorder(Pathway, neglog10FDR))) +
    geom_segment(aes(x = 0, xend = neglog10FDR, yend = Pathway),
                 linewidth = 0.4, alpha = 0.5) +
    geom_point(aes(size = hit_ratio, color = `Statistic Q`)) +
    geom_vline(xintercept = -log10(0.05), linetype = 2, linewidth = 0.4) +
    scale_size_continuous(name = "Hits / Total", range = c(2, 8)) +
    scale_color_viridis_c(name = "Statistic Q") +
    labs(x = expression(-log[10]~FDR), y = NULL, title = title_txt) +
    theme_minimal(base_size = 12) +
    theme(axis.text.y = element_text(size = 10),
          plot.title  = element_text(face = "bold"))
  
  print(p)
  
  # Save as SVG (file name with CELL added)
  out_file <- paste0(str_replace(tools::file_path_sans_ext(basename(infile)),
                                 "msea_qea_result_", ""),
                     "_CELL.svg")
  ggsave(out_file, plot = p, device = "svg", width = 8, height = 6, dpi = 300)
  
  invisible(p)
}

# Run for all files in your working directory
files <- list.files(pattern = pattern, full.names = TRUE)
if (length(files) == 0) stop("No files found matching pattern: ", pattern)

invisible(lapply(files, function(f) {
  df <- read_qea(f)
  plot_qea(df, f, topN = topN)
}))
