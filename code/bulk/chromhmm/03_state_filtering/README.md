# ChromHMM state refinement by peak summits

## Purpose

Refine the cell-type/condition-associated ChromHMM states before integrating them with ChromBPNet/FiNeMo motif instances. This step generates the filtered segmentations consumed by the downstream motif-summary scripts.

## ATAC state refinement

`filter_atac_states_by_summits.ipynb` intersects selected ATAC ChromHMM states with ATAC-seq MACS2 summits from the conditions assigned to each state. For these selected states, intervals supported by the expected condition(s) are retained while intervals also supported by conditions outside the state assignment are excluded. The refined selected states are then combined with the remaining unfiltered states to produce `ATAC_23_filtered_segments.bed`.

## Histone state refinement

`filter_histone_states_by_h3k27ac_summits.ipynb` applies the analogous procedure to enhancer-associated histone states using starting-cell H3K27ac summits. The refined enhancer states are combined with the remaining histone states to produce `Histones_22_filtered_segments.bed`.

## Downstream use

- `ATAC_23_filtered_segments.bed` -> `integration/chromhmm_chrombpnet/01_atac_state_motif_summary/`
- `Histones_22_filtered_segments.bed` -> `integration/chromhmm_chrombpnet/02_histone_state_motif_summary/`

State-to-sample assignments used for downstream interpretation are recorded in the integration `configs/` directory.
