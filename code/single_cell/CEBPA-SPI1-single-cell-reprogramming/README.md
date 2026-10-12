# Sequence context and partner availability specify divergent cell fates from a single transcription factor

Code to reproduce Figures 7 and 8 and Supplementary Figures 7 and 8 (CITE-seq and single-cell multiome ATAC analysis of CEBPA, SPI1 and CEBPA+SPI1 reprogramming of mouse embryonic fibroblasts).

Hammad Farooq<sup>#</sup>, Kenneth Joves<sup>#</sup>, Aydolun Petenkaya<sup>#</sup>, Itay Solomon<sup>#</sup>, Eva H.W. Leung, Sen Zhang, Mary Menhart, Ali Farhat, Lin Du, Nissim Hay, Sandra Pinho\*, and Constantinos Chronis\*

<sup>#</sup> These authors contributed equally. \* Correspondence

## Data

Sequencing data have been deposited in the Gene Expression Omnibus under SuperSeries accession **GSE344075**, which will be made public upon publication.

The scripts expect the input files in a `data/` folder:

```
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

## Code

```
R/
├── setup.R              packages, file paths, shared colours
├── citeseq_helpers.R    loads the CITE-seq object and plotting functions
├── atac_helpers.R       loads the scATAC object and ChromHMM functions
├── fig7.R               Figure 7a–g
├── figS7.R              Supplementary Figure 7a–j
├── fig8.R               Figure 8a, 8c
├── fig8h.R              Figure 8h
└── figS8.R              Supplementary Figure 8a–b
```


Requires R with Seurat (v5), Signac, GenomicRanges, ggplot2, dplyr, tidyr, patchwork, ComplexHeatmap, circlize, readxl, readr, stringr, RColorBrewer, viridisLite, rvest, ggseqlogo, universalmotif, scales and here.

## Contact

Please direct any requests for further information to the lead contact, Dr Constantinos Chronis (chronis@uic.edu), Department of Biochemistry and Molecular Genetics, University of Illinois Chicago, Chicago, IL 60607, USA.

Questions about the code can be directed to Aydolun Petenkaya (apeten2@uic.edu), Department of Biomedical Engineering, University of Illinois Chicago, Chicago, IL 60607, USA.
