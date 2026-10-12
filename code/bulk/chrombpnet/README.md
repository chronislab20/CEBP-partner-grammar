# ChromBPNet / TF-MoDISco / FiNeMo workflow

This directory contains the ChromBPNet analysis branch used to learn base-resolution sequence features from ATAC-seq, evaluate the trained models, derive bias-corrected nucleotide contribution scores, discover de novo motifs with TF-MoDISco, consolidate motifs across conditions, and identify predictive motif instances with FiNeMo.

## Workflow order

1. `01_prepare_inputs/` — prepare filtered BAMs, relaxed peaks, blacklist-filtered top-200k peaks, and fold-specific GC-matched nonpeaks.
2. `02_train_bias_models/` — train candidate ChromBPNet bias models across samples, thresholds, and folds.
3. `03_select_bias_models/` — inspect bias-model QC and select the sample-specific threshold.
4. `04_train_models/` — train one full ChromBPNet model per sample and chromosome fold.
5. `05_model_evaluation/` — summarize held-out model performance (Fig. 5b).
6. `06_contribution_scores/` — compute count contribution scores from `chrombpnet_nobias.h5`, average across folds, and convert averaged scores to bigWig.
7. `07_predictions/` — generate and average predicted accessibility tracks across folds.
8. `08_tfmodisco/` — run TF-MoDISco-lite on fold-averaged contribution scores.
9. `09_merge_motif_lexicon/` — consolidate related motifs across conditions while retaining seqlet provenance.
10. `10_finemo/` — annotate the merged motif lexicon and call contribution-aware motif instances with FiNeMo.

## Configuration

- `configs/samples.tsv` records the 20 included cell-state/condition models and selected bias thresholds.
- `configs/workflow_parameters.tsv` records key manuscript parameters.
- `configs/folds/` contains the five chromosome hold-out split JSON files used for the analysis.
- `configs/paths.example.sh` documents the expected local directory layout.

## Execution environment

The original large-scale model training and interpretation jobs were launched on ALCF Aurora with Parsl. The executor configuration is retained under `aurora/` for provenance. The analysis commands can be adapted to another scheduler without changing the scientific parameters.

## Major external tools

ChromBPNet, TensorFlow, TF-MoDISco-lite, FiNeMo, SAMtools, MACS2, BEDTools, Python, R, and Parsl.

See the stage-specific README files for workflow-specific details, inputs, parameters, and outputs.
