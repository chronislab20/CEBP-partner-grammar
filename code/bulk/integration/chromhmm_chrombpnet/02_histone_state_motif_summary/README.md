# Histone-state motif summary

## Purpose

Generate motif-by-histone-ChromHMM-state summaries used for the enhancer side of Fig. S5e.

## Method

Merged FiNeMo motif instances are overlapped with the refined `Histones_22_filtered_segments.bed` segmentation at 200-bp resolution. Histone-state/sample combinations are restricted according to `../../configs/histone_state_sample_map.tsv`.

For each motif/state combination, FiNeMo `hit_importance` is summed and divided by the genomic fraction occupied by the corresponding histone state. Motif occupancy proportions are calculated as the fraction of unique 200-bp bins in a state that overlap the motif. The matrices used for plotting are row-standardized within motif.

The enhancer panel used in Fig. S5e focuses on the cell-type-associated enhancer states H4-H8.

## Main outputs

- `Finemo_merged_Pr_Histone_sum_importance_over_state_size_sep_counts.csv`
- `Finemo_merged_Pr_histone_motif_prop_counts.csv`

## Script

`summarize_finemo_by_histone_state.R`
