# ATAC stacked ChromHMM model

## Purpose

Build the ATAC-based ChromHMM segmentation across the starting and CEBP/OSKM perturbation conditions. The selected 23-state model supports Fig. S5a and the motif-by-state analysis in Fig. 5c.

## Method

Each of the 20 ATAC-seq conditions is represented as an individual mark in one stacked/virtual ChromHMM cell type. Reads are binarized in non-overlapping 200-bp bins using `BinarizeBam -mixed`, and `LearnModel` is run at 200-bp resolution while preserving input-column order. The selected downstream model contains 23 states.

`reorder_23_state_model.ipynb` applies the final raw-state to A1-A23 mapping. `plot_23_state_model.ipynb` visualizes the reordered emission matrix used for Fig. S5a.

## Files

- `run_stacked_atac_chromhmm.sh` - ChromHMM binarization/model fitting and state BED generation.
- `reorder_23_state_model.ipynb` - applies the final raw-state to A1-A23 mapping.
- `plot_23_state_model.ipynb` - plots the reordered emission matrix.

## Inputs and outputs

The stacked condition table is in `../configs/atac_cellmarkfiletable.txt`. The reordered 23-state segmentation is refined for selected open/changing states in `../03_state_filtering/` before integration with FiNeMo motif instances.
