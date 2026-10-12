# Chromosome cross-validation folds

These are the five chromosome split JSON files used for the bulk mouse ChromBPNet analysis.

ChromBPNet uses a fold JSON to define three disjoint chromosome sets: `train`, `valid`, and `test`. The same fold definition is supplied to background-region generation and model training so that GC-matched nonpeaks and trained models use the same chromosome partition.

The five files included here are:

- `fold_0.json`
- `fold_1.json`
- `fold_2.json`
- `fold_3.json`
- `fold_4.json`

They are mouse (`mm10`) folds and contain chromosomes chr1–chr19, chrX, and chrY. They should not be replaced by the human (`hg38`) five-fold JSON files distributed with the ChromBPNet Zenodo example data.

ChromBPNet preprocessing documentation describes the fold-file format and how train/validation/test chromosomes are used:
https://github.com/kundajelab/chrombpnet/wiki/Preprocessing

The human example folds are available in the ChromBPNet Zenodo data record (10.5281/zenodo.7445373), but those are not the folds used here.
