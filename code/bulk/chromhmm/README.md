# Stacked ChromHMM workflow

Two stacked ChromHMM models were used at 200-bp resolution: a 22-state histone model across four starting cell types and nine histone marks, and a 23-state ATAC model across 20 starting/perturbed conditions.

## Stages

- `01_histone_model/`: binarization/model fitting, state reordering, and Fig. S5d plotting.
- `02_atac_model/`: binarization/model fitting, state reordering, and Fig. S5a plotting.
- `03_state_filtering/`: refinement of selected open/enhancer states using matched peak summits before motif integration.

State reordering and biological annotations are preserved in `configs/`. The exact lower bound of the exploratory K grid is not scientifically material to downstream figures; the selected final models are K=22 and K=23.
