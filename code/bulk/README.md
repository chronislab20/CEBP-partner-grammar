# Bulk regulatory-genomics analyses

This directory contains code for the bulk computational analyses included in the manuscript, with a focus on the ChromBPNet, TF-MoDISco/FiNeMo, ChromHMM, and integrated regulatory-grammar analyses.

The code is organized by analysis component rather than as a single end-to-end workflow. Individual subdirectories contain workflow-specific documentation, scripts, configuration files, and figure-generation code.

## Analysis overview

The ChromBPNet and chromatin-state analyses are organized as follows:

```text
ATAC-seq BAMs
   │
   ├── ChromBPNet input preparation
   │      └── 5-fold chromosome hold-out models
   │             ├── bias-model selection
   │             ├── model evaluation
   │             ├── bias-corrected contribution scores
   │             ├── TF-MoDISco motif discovery
   │             ├── cross-condition motif lexicon merging
   │             └── FiNeMo contribution-aware motif-instance calling
   │
   └── stacked ATAC ChromHMM ───────────────┐
                                             ├── motif/state integration ── Fig. 5c
Histone-mark data ── stacked ChromHMM ──────┤
                                             └── enhancer/promoter analysis ── Fig. S5e
RNA-seq ── starting-cell expression groups ── TSS/promoter refinement ────────┘
```

## Directory structure

```text
bulk/
├── README.md
├── chrombpnet/
│   ├── 01_prepare_inputs/
│   ├── 02_train_bias_models/
│   ├── 03_select_bias_models/
│   ├── 04_train_models/
│   ├── 05_model_evaluation/
│   ├── 06_contribution_scores/
│   ├── 07_predictions/
│   ├── 08_tfmodisco/
│   ├── 09_merge_motif_lexicon/
│   ├── 10_finemo/
│   ├── aurora/
│   └── configs/
├── chromhmm/
│   ├── 01_histone_model/
│   ├── 02_atac_model/
│   ├── 03_state_filtering/
│   └── configs/
└── integration/
    └── chromhmm_chrombpnet/
        ├── 01_atac_state_motif_summary/
        ├── 02_histone_state_motif_summary/
        ├── 03_tss_promoter_definition/
        ├── 04_tss_motif_summary/
        ├── 05_fig5c/
        ├── 06_figS5e/
        └── configs/
```

## ChromBPNet analysis

[`chrombpnet/`](chrombpnet/) contains the code used to prepare ATAC-seq inputs, train and select bias models, train five chromosome hold-out ChromBPNet models per condition, evaluate model performance, generate bias-corrected contribution scores and predictions, identify de novo motifs with TF-MoDISco, merge motifs across conditions into a consolidated motif lexicon, and call contribution-aware motif instances with FiNeMo.

See [`chrombpnet/README.md`](chrombpnet/README.md) for the detailed workflow and expected inputs and outputs.

## ChromHMM analysis

[`chromhmm/`](chromhmm/) contains the stacked ChromHMM analyses used to summarize chromatin states across cell types and experimental conditions. Separate models were generated from ATAC-seq and histone-mark datasets, followed by state filtering/refinement for downstream integration.

See [`chromhmm/README.md`](chromhmm/README.md) for details.

## ChromBPNet–ChromHMM integration

[`integration/chromhmm_chrombpnet/`](integration/chromhmm_chrombpnet/) contains analyses integrating FiNeMo motif instances and ChromBPNet contribution information with ATAC- and histone-defined ChromHMM states and RNA-defined promoter/TSS groups. This branch includes code associated with the motif/state summary in Fig. 5c and the enhancer/promoter analysis in Fig. S5e.

See [`integration/chromhmm_chrombpnet/README.md`](integration/chromhmm_chrombpnet/README.md) for details.

## Key analysis settings represented in the code

The included scripts and configuration files preserve the main settings used for these analyses, including:

- mm10 reference genome for the ATAC/ChromBPNet analysis;
- top **200,000** relaxed ATAC peaks per sample after blacklist filtering;
- **2,114-bp** ChromBPNet sequence input with prediction over the central **1,000 bp**;
- **five** chromosome hold-out folds;
- sample-specific ChromBPNet bias-threshold selection from **0.5, 0.6, and 0.7**;
- count contribution scores from the bias-corrected `chrombpnet_nobias.h5` model;
- TF-MoDISco motif discovery from fold-averaged contribution scores;
- a consolidated cross-condition motif lexicon followed by FiNeMo motif-instance calling;
- a **23-state** stacked ATAC ChromHMM model and a **22-state** stacked histone ChromHMM model;
- K=4 RNA-defined starting-cell groups for the promoter/TSS analysis.

Parameter tables, sample definitions, fold files, and example path configurations are provided under the relevant `configs/` directories.

## Input and intermediate files

Raw sequencing files and large intermediate/model outputs are not distributed with the repository. Depending on the analysis step, the scripts may require processed BAM files, peak/summit BED files, ChromBPNet model files, contribution-score HDF5 or bigWig files, FiNeMo outputs, ChromHMM segmentations, or RNA-seq-derived tables.

Example path configuration files are provided at:

```text
chrombpnet/configs/paths.example.sh
chromhmm/configs/paths.example.sh
integration/chromhmm_chrombpnet/configs/paths.example.sh
```

These should be adapted to the local filesystem and computing environment.

## Computing environment

Parts of the ChromBPNet analysis were run on ALCF Aurora, while other analyses were run in separate computing environments. Machine-specific paths have been removed from the public code. The Aurora-specific Parsl configuration is retained under [`chrombpnet/aurora/`](chrombpnet/aurora/) to document the execution setup used for those jobs.

This directory is intended to document and support reproduction of the custom bulk analyses rather than provide a single fully automated pipeline. For repository-wide scope and data availability, see the [main README](../../README.md).
