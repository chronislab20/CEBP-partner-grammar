# Histone stacked ChromHMM model

## Purpose

Build the histone-mark ChromHMM segmentation used to define chromatin states across the four starting cell types and to support Fig. S5d and the enhancer-state analysis in Fig. S5e.

## Method

Nine histone marks from MEF, NSC, ESC, and B-cell starting states are represented in a single stacked/virtual ChromHMM cell type. Reads are binarized in non-overlapping 200-bp bins with `BinarizeBam -mixed`. `LearnModel` is run with 200-bp bins and the original mark order preserved. The selected downstream model contains 22 states.

The raw 22-state solution is reordered into the H1-H22 display/analysis order in `reorder_22_state_model.ipynb`. `plot_22_state_model.ipynb` visualizes the reordered emission matrix used in Fig. S5d.

## Files

- `run_stacked_histone_chromhmm.sh` - ChromHMM binarization/model fitting and state BED generation.
- `reorder_22_state_model.ipynb` - applies the final raw-state to H1-H22 mapping.
- `plot_22_state_model.ipynb` - plots the reordered emission matrix.

## Inputs

The stacked cell/mark table is in `../configs/histone_cellmarkfiletable.txt`. Paths to BAMs, ChromHMM, genome resources, and output locations should be supplied through the public path configuration rather than hard-coded server paths.

## Outputs

The key downstream output is the reordered 22-state genomic segmentation. Selected enhancer-associated states are further refined in `../03_state_filtering/` before motif integration.
