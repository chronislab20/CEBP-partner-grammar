# Figure S5e assembly

## Purpose

Assemble the enhancer-state and promoter/TSS motif panels in Fig. S5e.

## Enhancer-side input

The histone-state motif summary from `../02_histone_state_motif_summary/` provides motif activity across histone ChromHMM states. Fig. S5e uses the cell-type-associated enhancer states H4-H8 (MEF, NSC, two ESC enhancer states, and B-cell enhancer state).

## Promoter-side input

The promoter/TSS summary from `../04_tss_motif_summary/` provides motif activity for the four RNA-defined K=4 starting-cell TSS groups. These promoter groups are not ChromHMM promoter states; they are RNA-defined TSS groups refined using H3K4me3/H3K27ac information.

## Notebooks and annotations

- `plot_enhancer_promoter_heatmap.ipynb` - combines the enhancer-state and promoter/TSS matrices and row-standardizes motif scores for visualization.
- `plot_motif_logos.ipynb` - plots selected motif logos in the final S5e order.
- `5e_EP_motifs_min_prop_0.01_with_KC_names.csv` - final S5e motif selection/order and display-name annotations used by the plotting notebooks.

## Notebook outputs

These plotting notebooks are distributed without execution outputs because they depend on large intermediate files that are not stored in Git. After local execution, retain only compact scientific outputs and remove machine-specific paths before committing.
