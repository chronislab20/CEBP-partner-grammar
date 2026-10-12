# Figure 5c assembly

## Purpose

Assemble the motif-level panels shown in Fig. 5c from the processed ChromBPNet/FiNeMo and ATAC-ChromHMM outputs.

## Inputs

The main quantitative input is the ATAC-state motif summary produced by `../01_atac_state_motif_summary/`. Additional inputs include the merged TF-MoDISco contribution-weighted motifs, matched known motif PWMs, FiNeMo motif instances, and RNA-seq expression used for the TF-expression panel.

## Motif annotation table

`fig5c_motif_annotations.tsv` is a compact publication-facing annotation table used by all three notebooks in this directory. It consolidates the final motif order, motif IDs, cluster assignments, representative TF labels, and matched known motifs required for figure generation. The two larger historical working annotation files used during manual triage are intentionally not distributed.

## Notebooks

- `plot_motif_state_summary.ipynb` - motif instance counts/percentages and row-standardized FiNeMo importance across samples and selected ATAC ChromHMM states.
- `plot_motif_logos.ipynb` - de novo contribution-weighted motif logos and matched known PWM logos.
- `plot_tf_expression.ipynb` - expression panel for representative TFs associated with the displayed motifs.

## Notebook outputs

These plotting notebooks are distributed without execution outputs because they depend on large intermediate files that are not stored in Git. After local execution, retain only compact scientific outputs and remove machine-specific paths before committing.
