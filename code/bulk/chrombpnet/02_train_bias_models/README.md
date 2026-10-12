# Train bias models

`run_bias_models.py` contains the Parsl bias-model task and the task-submission logic. The ChromBPNet command implements the manuscript workflow: ATAC mode, sample-specific BAM/peaks/nonpeaks, fold-specific chromosome split, and candidate bias threshold.

Candidate thresholds evaluated for the manuscript were 0.5, 0.6, and 0.7 across five folds.
