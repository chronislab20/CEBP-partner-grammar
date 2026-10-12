library(Seurat)
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(readxl)
library(readr)
library(stringr)
library(RColorBrewer)
library(viridisLite)

root_dir <- here::here()
data_dir <- file.path(root_dir, "data")

fig7_dir <- file.path(root_dir, "figures", "Fig7")
figS7_dir <- file.path(root_dir, "figures", "FigS7")
fig8_dir <- file.path(root_dir, "figures", "Fig8")
figS8_dir <- file.path(root_dir, "figures", "FigS8")

for (d in c(fig7_dir, figS7_dir, fig8_dir, figS8_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

metascape_dir <- file.path(data_dir, "metascape")
chromhmm_dir <- file.path(data_dir, "chromhmm")
chrombpnet_dir <- file.path(data_dir, "chrombpnet")

paths <- list(
  seurat = file.path(data_dir, "CITEseq_final_RNA_object.rds"),
  atac = file.path(data_dir, "multiome_final_ATACusedonly.rds"),

  metascape_adipo = file.path(metascape_dir, "common_adipocyte", "metascape_result.xlsx"),
  metascape_myeloid = file.path(metascape_dir, "common_myeloid", "metascape_result.xlsx"),
  metascape_fig7f = file.path(
    metascape_dir, "Summary_enrichment_plots",
    "Metascape_seven_population_final_selected_terms_red3_blue3_data.xlsx"
  ),
  metascape_figS7i = file.path(
    metascape_dir, "Summary_enrichment_plots",
    "Metascape_supplemental_selected_terms_red3_blue3_data.xlsx"
  ),

  heatmap = file.path(data_dir, "heatmap", "FINAL_HEATMAP_REPRODUCTION_ALL_VALUES.rds"),

  human_mouse = file.path(data_dir, "transgenes", "0.1.1.Human_Mouse_Cells_matched_cebpspisamples_lite.rds"),
  zsgreen = file.path(data_dir, "transgenes", "0.1.2.mm10_ZsGreen_aligned_obj.rds"),
  flag = file.path(data_dir, "transgenes", "0.1.2.mm10_flag_aligned_obj.rds"),

  chromhmm_bed = file.path(chromhmm_dir, "Fig8_scATAC_ChromHMM_20state_segments.bed"),
  chromhmm_rpkm = file.path(chromhmm_dir, "Fig8_scATAC_ChromHMM_St1_St20_log2RPKM_by_cell.rds"),

  fig8h_values = file.path(chrombpnet_dir, "Fig8j_FINAL_highlighted_50_plot_values_pct_peaks_containing_motif.csv"),
  modisco_html = file.path(chrombpnet_dir, "merged_modisco.html", "motifs.html"),
  modisco_cwm = file.path(chrombpnet_dir, "finemo", "merged", "CWMs"),
  motifs_meme = file.path(chrombpnet_dir, "motifs.meme.txt")
)

check_inputs <- function(keys) {
  missing <- keys[!file.exists(unlist(paths[keys]))]
  if (length(missing) > 0) {
    stop(
      "Missing input files (see README.md for the expected data/ layout):\n",
      paste0("  - ", unlist(paths[missing]), collapse = "\n"),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

check_features <- function(features, mat, what = "features") {
  missing <- setdiff(features, rownames(mat))
  if (length(missing) > 0) {
    stop("Missing ", what, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

warn_unmatched <- function(raw, levels, what) {
  bad <- setdiff(unique(raw), levels)
  if (length(bad) > 0) {
    warning(
      length(bad), " ", what, " value(s) not found in the plotting order and ",
      "will be shown as NA/dropped: ", paste(head(bad, 5), collapse = "; "),
      call. = FALSE
    )
  }
  invisible(bad)
}

sample_order <- c("MEF_CEBPA", "MEF_CEBPA_SPI1", "MEF_SPI1")

sample_cols <- c(
  "MEF_CEBPA" = "#1b4332",
  "MEF_CEBPA_SPI1" = "#4e3594",
  "MEF_SPI1" = "#c75000"
)

sample_labels <- c(
  "MEF_CEBPA" = "CEBPA",
  "MEF_CEBPA_SPI1" = "CEBPA+SPI1",
  "MEF_SPI1" = "SPI1"
)

identity_cols <- c(
  "Residual Fibroblast" = "#7FA37D",
  "Adhesion-remodel mesenchyme" = "#88EECC",
  "Cycling" = "#D77AA6",
  "Adipogenic transition (early)" = "#F1AE5B",
  "Sca1+Cd34+CD140+ Fibro-adipogenic Progenitor" = "#BBA14F",
  "Dual-prog transition" = "#B8DE6F",
  "Chromatin remodeling-high adipogenic" = "#EDDEA4",
  "Fap+ activated adipogenic progenitor" = "#C96B1B",
  "Early adipocyte" = "#9E4F12",
  "Pparg+ Adipo-mac hybrid" = "#B8A1D9",
  "Cycling myeloid-like" = "#A8D2F0",
  "Inflammatory Mac-like" = "#2B5C9A",
  "CMP-like progenitor" = "#6FAFE0",
  "APC Mac-like" = "#163E73",
  "Myeloid-like" = "#6C8FC5",
  "Skeletal/mesenchyme" = "#ADB5BD",
  "Low-signature-SPI1" = "#D7D9CE"
)

umap_theme <- theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank()
  )

centered_title <- theme(plot.title = element_text(hjust = 0.5, face = "plain"))
