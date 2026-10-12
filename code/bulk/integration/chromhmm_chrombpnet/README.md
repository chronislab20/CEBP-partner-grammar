# ChromHMM x ChromBPNet / FiNeMo integration

This workflow links motif instances learned from ChromBPNet/TF-MoDISco/FiNeMo to refined ChromHMM states and to RNA-defined starting-cell promoter/TSS groups. It contains the quantitative preprocessing and final plotting code for Fig. 5c and Fig. S5e.

## Workflow

1. `01_atac_state_motif_summary/` - overlap FiNeMo hits with refined ATAC ChromHMM states and summarize motif importance.
2. `02_histone_state_motif_summary/` - overlap FiNeMo hits with refined histone ChromHMM states and summarize motif importance.
3. `03_tss_promoter_definition/` - define K=4 RNA-associated starting-cell TSS groups and refine them with promoter-associated histone signals.
4. `04_tss_motif_summary/` - summarize FiNeMo motif importance across the four final TSS groups.
5. `05_fig5c/` - assemble Fig. 5c motif/state, logo, and TF-expression panels.
6. `06_figS5e/` - assemble Fig. S5e enhancer-state and promoter/TSS motif panels.

## ChromHMM state scores

For each motif and state, FiNeMo `hit_importance` values are summed and normalized by the genomic fraction occupied by the state. Matrices displayed in the figures are subsequently row-standardized within motif.

## Promoter/TSS scores

For each motif and RNA-defined promoter/TSS group, FiNeMo `hit_importance` is summed and normalized by the number of regions in the group. The resulting motif-by-group matrix is row-standardized within motif for display.

## Notebook-output policy

Most notebooks are distributed without outputs because they depend on large external intermediate files. Compact publication-relevant outputs may be retained after execution, provided machine-specific information is removed.
