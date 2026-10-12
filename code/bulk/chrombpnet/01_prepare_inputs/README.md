# 01 — Prepare ChromBPNet input regions

This stage prepares the bulk ATAC-seq inputs used for ChromBPNet. It was executed on a local laboratory compute server before the prepared files were transferred to Aurora for model training and interpretation.

## Script

- `prepare_bulk_chrombpnet_inputs.sh` — public/configurable version of the input-preparation commands used for the manuscript.

Server-specific absolute paths have been removed. Set `WORK_ROOT` and `RAW_BAM_DIR` before running; optional Conda variables are described at the top of the script.

## Processing steps

1. Download the mm10 no-alt ENCODE analysis-set FASTA.
2. Download mm10 chromosome sizes.
3. Download ENCODE blacklist `ENCFF547MET`.
4. Sort and index ATAC-seq BAM files with SAMtools.
5. Filter alignments with `samtools view -F 780 -q 30`.
6. Call relaxed ATAC peaks with MACS2 using `-p 0.01 --shift -75 --extsize 150 --nomodel --keep-dup all --call-summits`.
7. Rank peaks by MACS2 `signalValue` and retain the top 200,000 peaks per sample.
8. Expand the blacklist by 1,057 bp on each side and remove overlapping peaks.
9. Generate GC-matched background/nonpeak regions with `chrombpnet prep nonpeaks` separately for each of the five chromosome folds.

## Fold definitions

The exact five fold JSON files used in the manuscript are stored under `../configs/folds/`. Each fold assigns non-overlapping mm10 chromosomes to training, validation and test sets.

## Outputs passed to model training

The downstream Aurora workflow uses:

- filtered/indexed BAM files;
- blacklist-filtered top-200,000 peak BED files;
- fold-specific GC-matched nonpeak/background regions;
- mm10 FASTA and chromosome sizes;
- the five fold JSON files.
