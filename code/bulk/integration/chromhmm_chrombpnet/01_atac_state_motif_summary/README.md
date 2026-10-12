# ATAC-state motif summary

## Purpose

Generate motif-by-sample and motif-by-ATAC-ChromHMM-state matrices used in Fig. 5c.

## Method

The script reads merged FiNeMo `hits.tsv` motif instances and the refined `ATAC_23_filtered_segments.bed` segmentation. ChromHMM intervals are represented at 200-bp resolution and overlapped with FiNeMo motif instances. State/sample combinations are restricted to the biologically assigned conditions recorded in `../../configs/atac_state_sample_map.tsv`.

For each motif and ATAC state, FiNeMo `hit_importance` values are summed. The state-level score is normalized by the genomic fraction occupied by that state:

`normalized state score = sum(hit_importance) / state genomic fraction`

The resulting motif-by-state and motif-by-sample matrices are row-standardized within motif for visualization. Motif occupancy proportions are also calculated from the fraction of unique 200-bp state bins overlapping a motif instance.

## Main output

`Finemo_merged_Pr_ATAC_sum_importance_over_state_size_sep_counts.csv` is the principal matrix consumed by the Fig. 5c plotting notebook.

## Script

`summarize_finemo_by_atac_state.R`
