# Single-cell CEBPA/SPI1 reprogramming analysis

This directory contains code used for the **CITE-seq and single-cell multiome ATAC analyses** of CEBPA, SPI1, and CEBPA+SPI1 reprogramming of mouse embryonic fibroblasts, supporting **Figures 7 and 8** and **Supplementary Figures S7 and S8**.

For manuscript-level information, data availability, citation, and contacts, see the [main repository README](../../README.md).

## Input data

The scripts expect the required input files in a local `data/` directory:

```text
data/
├── CITEseq_final_RNA_object.rds
├── multiome_final_ATACusedonly.rds
├── chromhmm/
│   ├── Fig8_scATAC_ChromHMM_20state_segments.bed
│   └── Fig8_scATAC_ChromHMM_St1_St20_log2RPKM_by_cell.rds
├── chrombpnet/
│   ├── Fig8j_FINAL_highlighted_50_plot_values_pct_peaks_containing_motif.csv
│   ├── merged_modisco.html/motifs.html
│   ├── finemo/merged/CWMs/
│   └── motifs.meme.txt
├── heatmap/
│   └── FINAL_HEATMAP_REPRODUCTION_ALL_VALUES.rds
├── transgenes/
│   ├── 0.1.1.Human_Mouse_Cells_matched_cebpspisamples_lite.rds
│   ├── 0.1.2.mm10_ZsGreen_aligned_obj.rds
│   └── 0.1.2.mm10_flag_aligned_obj.rds
└── metascape/
    ├── common_adipocyte/metascape_result.xlsx
    ├── common_myeloid/metascape_result.xlsx
    └── Summary_enrichment_plots/
        ├── Metascape_seven_population_final_selected_terms_red3_blue3_data.xlsx
        └── Metascape_supplemental_selected_terms_red3_blue3_data.xlsx
```

If the per-cell ChromHMM state scores (`..._log2RPKM_by_cell.rds`) are not provided, they are recomputed from the scATAC fragment files referenced in the ATAC object.

## Scripts

```text
R/
├── setup.R              # packages, file paths, and shared colours
├── citeseq_helpers.R    # loads the CITE-seq object and plotting functions
├── atac_helpers.R       # loads the scATAC object and ChromHMM functions
├── fig7.R               # Figure 7a–g
├── figS7.R              # Supplementary Figure S7a–j
├── fig8.R               # Figure 8a and 8c
├── fig8h.R              # Figure 8h
└── figS8.R              # Supplementary Figure S8a–b
```

## Software requirements

Requires R with `Seurat` (v5), `Signac`, `GenomicRanges`, `ggplot2`, `dplyr`, `tidyr`, `patchwork`, `ComplexHeatmap`, `circlize`, `readxl`, `readr`, `stringr`, `RColorBrewer`, `viridisLite`, `rvest`, `ggseqlogo`, `universalmotif`, `scales`, and `here`.

## Usage note

These scripts are provided as analysis and figure-generation code. File paths and input objects may need to be configured for the local environment before execution.
