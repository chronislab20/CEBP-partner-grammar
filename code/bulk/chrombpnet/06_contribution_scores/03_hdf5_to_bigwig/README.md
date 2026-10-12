# Convert averaged contribution HDF5 to bigWig

## Parsl task
`run_hdf5_to_bigwig.py` invokes the project helper `importance_hdf5_to_bigwig.py`.

## Role in this study
This is project-specific post-processing of the fold-averaged contribution HDF5. ChromBPNet itself can emit contribution-score bigWigs from `contribs_bw`; here the custom helper was used after fold averaging so the browser track corresponds to the averaged contribution scores.

Related ChromBPNet documentation: https://github.com/kundajelab/chrombpnet/wiki/Generate-contribution-score-bigwigs
