# 03 — Inspect and select bias models

`Plot_Bias_model_report.R` is the script used to **inspect and select bias models** for the bulk C/EBP ChromBPNet analysis.

Candidate bias models were trained at `bias_threshold_factor` values of 0.5, 0.6, and 0.7 across five chromosome folds. The purpose of this hyperparameter is to control which non-peak regions are retained for bias-model training. ChromBPNet filters non-peak regions based on their total ATAC counts relative to a low-count quantile of peak regions, scaled by `bias_threshold_factor`.

## Rationale for choosing `bias_threshold_factor`

The selected value should balance two competing considerations:

- **If `bias_threshold_factor` is too low**, too many non-peak regions can be filtered out. This can make the retained non-peak training set differ from peak regions in sequence composition, particularly GC content, and may cause the bias model to learn an AT-rich sequence preference that does not transfer appropriately to peak regions.
- **If `bias_threshold_factor` is too high**, non-peak regions with relatively high accessibility counts can remain in the training set. These regions may contain cell-type-specific transcription-factor signal, causing the bias model to learn biological TF motifs that should instead remain in the bias-corrected ChromBPNet model.

Following the ChromBPNet guidance for ATAC-seq, candidate models were therefore inspected for both **bias-only motif content** and **performance in peak regions**. A suitable model should capture Tn5-associated sequence preference without strong cell-type-specific TF motifs, while avoiding strongly negative peak-region Pearson correlation. If TF motifs are captured, ChromBPNet recommends reducing `bias_threshold_factor`; if no TF motifs are captured but peak-region Pearson correlation remains too negative (approximately below -0.3), the threshold can be increased while continuing to check for TF leakage.

For this study, `Plot_Bias_model_report.R` was used to inspect the candidate bias-model reports and select the final threshold/model for each sample. The selected bias threshold/model for each sample is recorded in `../configs/samples.tsv` and is used during full ChromBPNet model training.

The script preserves the QC logic used for the manuscript. Set `project_path` and `report_path` at the top of the script (or adapt them to your local configuration) before execution.

For the rationale and recommended procedure for selecting `bias_threshold_factor`, see the official ChromBPNet FAQ:
https://github.com/kundajelab/chrombpnet/wiki/FAQ

For general bias-model training and QC information, see:
https://github.com/kundajelab/chrombpnet/wiki/ChromBPNet-training
