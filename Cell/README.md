# Cell Metabolomics Analysis

Scripts and data used for analysis of cellular metabolomics in the manuscript  
*“Derlin-mediated ERAD of lipid regulator ORMDL3 safeguards mitochondrial function”* (DOI: placeholder).

This directory contains all R scripts used to generate the major metabolomic analyses and figures derived from **cell extracts**, including:

- Data preprocessing and normalization
- PCA analysis
- Differential metabolite testing (limma)
- Volcano plots and UpSet analyses
- Heatmaps of significant metabolites
- Targeted metabolite boxplots
- MetaboAnalyst QEA file generation
- Pathway enrichment figure generation

## Contents

- `DKO_metab_analysis_pca_volcano_upset_cell.R`  
  PCA, differential analysis, volcano plots, and overlap analyses.

- `DKO_metab_analysis_heatmap_cell.R`  
  Heatmaps of significant metabolites.

- `DKO_metab_analysis_cell_boxplots.R`  
  Targeted metabolite visualization.

- `DKO_metab_enrichment_cell_file_generation.R`  
  Generates MetaboAnalyst QEA input files.

- `Enrichment_figures_MSEA_Cell.R`  
  Generates enrichment figures from MetaboAnalyst output.

- Data files included here are those used directly by the scripts.

## Usage

Scripts can be run independently but are typically executed in the following order:

1. PCA and differential analysis
2. Heatmap generation
3. Targeted metabolite plots
4. MetaboAnalyst QEA file generation
5. Enrichment figure generation

## Contact

Mark Phillips — philmark@oregonstate.edu
