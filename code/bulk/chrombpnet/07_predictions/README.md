# 07 — Predicted accessibility tracks

`run_predictions.py` runs `chrombpnet pred_bw` using the bias-corrected `chrombpnet_nobias.h5` model for each fold.

`avg_pred_bw.sh` averages the five fold-specific prediction bigWigs for each of the 20 manuscript models. The public helper uses configurable `MODEL_ROOT`, `CHROM_SIZES`, and optional `PREDICTION_OUTPUT_DIR` variables rather than machine-specific paths.
