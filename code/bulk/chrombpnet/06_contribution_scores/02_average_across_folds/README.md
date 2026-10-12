# Average contribution scores across folds

## Parsl task
`run_average_deepshaps.py` launches the project helper `05-average_deepshaps.py` to combine contribution-score HDF5 files across folds 0–4.

## Role in this study
This is study-specific post-processing rather than a standalone ChromBPNet wiki command. The five chromosome hold-out models provide fold-specific interpretations; this step creates an averaged representation used for subsequent motif discovery and browser-track generation.
