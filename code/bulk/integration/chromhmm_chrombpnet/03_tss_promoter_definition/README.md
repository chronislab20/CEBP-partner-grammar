# Starting-cell promoter/TSS definition

## Purpose

Define cell-type-associated promoter/TSS groups for the promoter side of Fig. S5e.

## Workflow

### 1. RNA-based starting-cell groups

`01_cluster_starting_cell_rna_k4.ipynb` performs pairwise differential-expression analysis across the four starting cell types. Differential genes are combined, VST expression is averaged across replicates and standardized by gene, and the final expression solution uses K=4 clusters corresponding to MEF-, NSC-, ESC-, and B-cell-associated expression groups. Genes are mapped to mm10 TSS coordinates and expanded to +/-1 kb intervals.

### 2. Promoter-associated histone filtering

`02_filter_tss_by_histone_marks.ipynb` intersects the RNA-defined TSS intervals with matched H3K27ac and H3K4me3 information to retain promoter regions supported by the expected starting-cell chromatin state.

### 3. Cell-type-selective H3K27ac filtering

`03_filter_tss_by_celltype_h3k27ac_signal.ipynb` further refines each RNA/TSS group using cell-type H3K27ac signal. The final BED sets are `RNA_cluster_1.bed` through `RNA_cluster_4.bed`, corresponding to MEF, NSC, ESC, and B cell.

## Naming note

The original analysis directory contained `k5` in its path name, but this was only a legacy naming error. The final analysis uses K=4. Public paths use K=4 terminology.
