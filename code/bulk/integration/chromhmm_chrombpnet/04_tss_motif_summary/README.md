# TSS-group motif summary

## Purpose

Summarize ChromBPNet/FiNeMo motif activity across the four RNA-defined starting-cell promoter/TSS groups used in Fig. S5e.

## Method

The final `RNA_cluster_1`-`RNA_cluster_4` BED sets are intersected with ChromBPNet interpreted ATAC regions from the matching starting cell type. FiNeMo hits are restricted to these overlapping regions. For each motif and promoter/TSS group, `hit_importance` values are summed and normalized by the number of TSS regions in that group:

`normalized promoter score = sum(hit_importance) / group size`

The motif-by-group matrix is row-standardized within motif for visualization. The script also computes the fraction of regions in each promoter group containing each motif and records the maximum fraction across groups.

## Main outputs

- `Finemo_merged_sum_imp_per_size_cleaned_Starting_cell_types_clusters_counts.csv`
- `TSS_cluster_motif_max_frac_counts.csv`

## Script

`summarize_finemo_by_tss_group.R`
