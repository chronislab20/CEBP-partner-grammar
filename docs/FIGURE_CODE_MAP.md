# Figure-to-code map

This document maps the figure panels represented by the current public code release to their primary analysis directories.

## Figure 5a — ChromBPNet analysis workflow

The schematic summarizes the pipeline implemented in:

1. `code/bulk/chrombpnet/01_prepare_inputs/`
2. `code/bulk/chrombpnet/02_train_bias_models/`
3. `code/bulk/chrombpnet/03_select_bias_models/`
4. `code/bulk/chrombpnet/04_train_models/`
5. `code/bulk/chrombpnet/06_contribution_scores/`
6. `code/bulk/chrombpnet/08_tfmodisco/`
7. `code/bulk/chrombpnet/09_merge_motif_lexicon/`
8. `code/bulk/chrombpnet/10_finemo/`

The panel itself is a conceptual schematic rather than a direct program output.

## Figure 5b — ChromBPNet model performance

- `code/bulk/chrombpnet/05_model_evaluation/ChromBPNet_pearsonr.ipynb`

The notebook summarizes fold-specific ChromBPNet evaluation metrics and generates the model-performance panel.

## Figure 5c — cell-type-specific regulatory grammar

1. ChromBPNet / TF-MoDISco / FiNeMo: `code/bulk/chrombpnet/`
2. ATAC ChromHMM 23-state model: `code/bulk/chromhmm/02_atac_model/`
3. ATAC-state refinement: `code/bulk/chromhmm/03_state_filtering/filter_atac_states_by_summits.ipynb`
4. Motif importance summarized by ATAC ChromHMM state: `code/bulk/integration/chromhmm_chrombpnet/01_atac_state_motif_summary/`
5. Final motif/state, motif-logo, and TF-expression assembly: `code/bulk/integration/chromhmm_chrombpnet/05_fig5c/`

## Supplementary Figure S5a — stacked ATAC ChromHMM states

- `code/bulk/chromhmm/02_atac_model/`

This directory contains the 23-state ATAC ChromHMM workflow, state reordering, and plotting notebook.

## Supplementary Figure S5d — stacked histone ChromHMM states

- `code/bulk/chromhmm/01_histone_model/`

This directory contains the 22-state histone ChromHMM workflow, state reordering, and plotting notebook.

## Supplementary Figure S5e — enhancer/promoter motif grammar

### Enhancer side

1. Histone ChromHMM 22-state model.
2. H3K27ac summit refinement of enhancer-associated states.
3. FiNeMo motif importance summarized across histone states.
4. Final panel restricted to the selected enhancer-associated states H4–H8.

Primary code:

- `code/bulk/chromhmm/01_histone_model/`
- `code/bulk/chromhmm/03_state_filtering/filter_histone_states_by_h3k27ac_summits.ipynb`
- `code/bulk/integration/chromhmm_chrombpnet/02_histone_state_motif_summary/`
- `code/bulk/integration/chromhmm_chrombpnet/06_figS5e/`

### Promoter/TSS side

1. Starting-cell RNA-seq differential-expression union.
2. K=4 expression clustering.
3. Gene-to-TSS mapping and ±1 kb TSS intervals.
4. H3K4me3/H3K27ac-supported promoter filtering.
5. Cell-type-selective H3K27ac refinement.
6. FiNeMo motif importance normalized by TSS-group size.
7. Final S5e plotting notebook.

Primary code:

- `code/bulk/integration/chromhmm_chrombpnet/03_tss_promoter_definition/`
- `code/bulk/integration/chromhmm_chrombpnet/04_tss_motif_summary/`
- `code/bulk/integration/chromhmm_chrombpnet/06_figS5e/`

## Figures 7–8 and Supplementary Figures S7–S8 — single-cell CEBPA/SPI1 reprogramming

The authors' existing single-cell analysis package is included unchanged at:

- `code/single_cell/CEBPA-SPI1-single-cell-reprogramming-main.zip`

Its internal README maps the included R scripts to Figures 7, 8, S7, and S8 and documents the required CITE-seq, single-cell multiome ATAC, ChromHMM, ChromBPNet/FiNeMo, and supporting input files. The archive was not modified during repository cleanup.

## Panels not covered by this release

This release includes the organized Fig. 5/S5 bulk regulatory-genomics workflows and the unmodified single-cell archive for Figs. 7–8/S7–S8. Other manuscript panels should only be added to this map when their corresponding analysis code has been cleaned and included in the repository.
