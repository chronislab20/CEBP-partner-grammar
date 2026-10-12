# Bias-corrected contribution scores

## Parsl task
`run_contribution_scores.py` contains the Parsl task and submission logic. It runs `chrombpnet contribs_bw` using `chrombpnet_nobias.h5`.

## What ChromBPNet documents
ChromBPNet has profile and counts prediction heads. The contribution-score documentation specifically recommends using `chrombpnet_nobias.h5` to derive sequence contribution scores. `chrombpnet contribs_bw` can interpret either counts, profile, or both and writes HDF5 contribution scores, bigWig tracks, and the interpreted-regions BED file. The HDF5 stores one-hot sequence, hypothetical contribution scores, and projected contribution scores used for motif discovery and browser tracks.

Official documentation: https://github.com/kundajelab/chrombpnet/wiki/Generate-contribution-score-bigwigs

## This study
The manuscript Methods specify count contribution scores from the bias-corrected model as the primary interpretation used for downstream TF-MoDISco. Fold-specific scores were subsequently averaged across the five chromosome hold-out models.
