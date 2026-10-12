# 09 — Merge TF-MoDISco motif lexicon and annotate merged motifs

This stage merges de novo TF-MoDISco motifs discovered independently across the **20 bulk C/EBP manuscript ATAC-seq conditions** into a non-redundant motif lexicon.

## Files

- `CEBPB_All_Merge_modisco.py` — motif-merging workflow.
- `Merge_modisco_functions.py` — helper functions, including TOMTOM annotation functions.

The public version retains the original sample order, motif-merging parameters, and analysis logic while replacing server-specific paths with configurable environment variables.

## Required path variables

Set these before running:

```bash
export MODEL_ROOT=/path/to/chrombpnet_model_selected_bias
export MERGED_MOTIF_DIR=/path/to/merged_motif_lexicon
export JASPAR_MEME=/path/to/curated_JASPAR_motifs.meme
export KNOWN_MOTIF_LOGO_DIR=/path/to/known_motif_logos   # optional
```

For each condition the script associates:

- the 1,000-bp / 5,000,000-seqlet TF-MoDISco HDF5 result;
- the fold-averaged counts contribution HDF5;
- the fold-0 interpreted-regions BED used to recover peak/summit coordinates;
- the manuscript sample label.

## TOMTOM output

The workflow creates `tomtom_cwm/` and generates, for each merged pattern:

```text
pos_patterns.pattern_X_cwm.txt
pos_patterns.pattern_X_tomtom.csv
```

TOMTOM is therefore part of this stage rather than a separate pipeline.
