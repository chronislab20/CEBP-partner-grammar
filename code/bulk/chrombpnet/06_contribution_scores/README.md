# Contribution scores

The manuscript workflow uses **counts contribution scores** from the bias-corrected `chrombpnet_nobias.h5` model.

- `01_bias_corrected_contribs/`: run `chrombpnet contribs_bw --profile-or-counts counts` for each sample and fold.
- `02_average_across_folds/`: average fold-specific DeepSHAP HDF5 outputs with the original `05-average_deepshaps.py` helper.
- `03_hdf5_to_bigwig/`: convert averaged HDF5 contribution scores to bigWig using the original helper script.

Profile contribution scores and with-bias contribution-score diagnostics are intentionally not included because they were not used for the manuscript analyses.
