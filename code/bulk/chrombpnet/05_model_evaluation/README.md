# 05 — Model evaluation / Figure 5b

This directory contains the **executed notebook with outputs** used to summarize ChromBPNet model performance for Figure 5b.

## File

- `ChromBPNet_pearsonr.ipynb` — reads fold-specific `chrombpnet_metrics.json` files, extracts Pearson correlation between predicted and observed counts in peak regions, summarizes the five chromosome hold-out folds for each model, and generates the Figure 5b performance panel.

The notebook is retained with outputs so that the model-level correlations and rendered panels can be inspected without access to the trained-model directory.

## Inputs

```text
<MODEL_ROOT>/<sample>/fold_<N>/evaluation/chrombpnet_metrics.json
```

For the representative scatter panel the notebook reads the ChromBPNet evaluation output:

```text
<MODEL_ROOT>/<sample>/<fold>/evaluation/chrombpnet_only_peaks.counts_pearsonr.png
```

## Outputs

The notebook writes:

```text
fig5b_outputs/
├── chrombpnet_fold_level_pearson.tsv
├── chrombpnet_sample_mean_pearson.tsv
├── chrombpnet_fig5b_left_panel.pdf
├── chrombpnet_fig5b_left_panel.png
├── chrombpnet_fig5b_right_panel.pdf
└── chrombpnet_fig5b_right_panel.png
```

In the executed analysis, 100 fold-level evaluation files corresponding to 20 sample/condition models × 5 folds were summarized.
