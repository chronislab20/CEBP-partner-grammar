# Sequence context and partner availability specify divergent cell fates from a single transcription factor

This repository contains computational analysis code accompanying the manuscript:

**Hammad Farooq<sup>#</sup>, Kenneth Joves<sup>#</sup>, Aydolun Petenkaya<sup>#</sup>, Itay Solomon<sup>#</sup>, Eva H.W. Leung, Sen Zhang, Mary Menhart, Ali Farhat, Lin Du, Nissim Hay, Sandra Pinho\*, and Constantinos Chronis\***  
<sup>#</sup> These authors contributed equally. \* Correspondence

**Sequence context and partner availability specify divergent cell fates from a single transcription factor**

The repository includes the organized **ChromBPNet**, **TF-MoDISco/FiNeMo**, **ChromHMM**, and ChromBPNet–ChromHMM integration workflows supporting **Fig. 5** and **Supplementary Fig. S5**, together with the **CITE-seq and single-cell multiome** analysis supporting **Figs. 7–8** and **Supplementary Figs. S7–S8**.

> **Repository scope.** This is a manuscript-analysis repository, not a general-purpose software package. Large sequencing files, trained neural-network models, intermediate bigWig/HDF5 files, and other generated data products are intentionally not stored in Git.

## Analysis overview

The workflow used here is:

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

## Repository structure

```text
code/
├── bulk/
│   ├── chrombpnet/                  # ChromBPNet → TF-MoDISco → FiNeMo workflow
│   ├── chromhmm/                    # stacked ATAC and histone ChromHMM workflows
│   └── integration/
│       └── chromhmm_chrombpnet/     # motif/state and promoter/TSS integration
└── single_cell/
    ├── README.md                    # scope of the single-cell analysis
    └── CEBPA-SPI1-single-cell-reprogramming/

docs/
├── FIGURE_CODE_MAP.md               # figure-to-code mapping
└── REPRODUCIBILITY.md               # configuration and reproducibility notes
```

Each analysis directory contains a local `README.md` describing its inputs, processing logic, and expected outputs.

This repository is intended to complement the Materials and Methods section by providing code used for custom analyses described in the manuscript, with the goal of improving transparency and reproducibility of the main results. It is not intended to serve as a fully automated or end-to-end executable workflow, and some analyses may require intermediate files, software environments, or computing resources specific to the original analysis setup.


## Quick start

### 1. Clone the repository

```bash
git clone <repository-url>
cd <repository-name>
```

### 2. Configure local paths

The public code does not contain laboratory- or cluster-specific absolute paths. Example path configuration files are provided at:

```text
code/bulk/chrombpnet/configs/paths.example.sh
code/bulk/chromhmm/configs/paths.example.sh
code/bulk/integration/chromhmm_chrombpnet/configs/paths.example.sh
```

Copy the relevant example to a local file, edit it for your system, and keep the local copy untracked by Git.

### 3. Follow the workflow-specific documentation

Start with:

- [`code/bulk/chrombpnet/README.md`](code/bulk/chrombpnet/README.md) for model training, interpretation, motif discovery, and FiNeMo.
- [`code/bulk/chromhmm/README.md`](code/bulk/chromhmm/README.md) for the stacked ATAC and histone ChromHMM models.
- [`code/bulk/integration/chromhmm_chrombpnet/README.md`](code/bulk/integration/chromhmm_chrombpnet/README.md) for the Fig. 5c / Fig. S5e integration analyses.
- [`code/single_cell/README.md`](code/single_cell/README.md) for the single-cell CEBPA/SPI1 analysis.
- [`docs/FIGURE_CODE_MAP.md`](docs/FIGURE_CODE_MAP.md) to locate the code associated with each included figure panel.

## Key manuscript parameters represented in the repository

The repository preserves the manuscript analysis settings, including:

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

Parameter tables and fold definitions are stored under each workflow’s `configs/` directory.

## Inputs not distributed in this repository

Because of file size and data-distribution constraints, this repository does not include raw FASTQ files, processed BAM files, trained ChromBPNet models, contribution-score HDF5/bigWig files, or full intermediate ChromHMM/FiNeMo output directories. The code expects those inputs to be supplied through the documented configuration paths.

Data accessions and experimental details should be taken from the associated manuscript and its data-availability statement.

## Reproducibility and portability

The analysis was run across more than one computing environment, including local laboratory infrastructure and ALCF Aurora. Machine-specific paths have been removed from the public code and replaced with configurable path variables. Aurora-specific Parsl configuration is retained under `code/bulk/chrombpnet/aurora/` for provenance; users on other systems should adapt the launcher/executor layer while preserving the analysis commands and parameters.

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for details.

## Notebooks

Notebook outputs may be retained when they provide compact, publication-relevant QC or figure-generation provenance. Large or environment-specific outputs are not required for the public release, and machine-specific paths or credentials should not be present in committed notebooks.

## Citation

If you use this code, please cite the associated manuscript:

> Hammad Farooq<sup>#</sup>, Kenneth Joves<sup>#</sup>, Aydolun Petenkaya<sup>#</sup>, Itay Solomon<sup>#</sup>, Eva H.W. Leung, Sen Zhang, Mary Menhart, Ali Farhat, Lin Du, Nissim Hay, Sandra Pinho\*, and Constantinos Chronis\*. **Sequence context and partner availability specify divergent cell fates from a single transcription factor.**

<sup>#</sup> These authors contributed equally. \* Correspondence

A final journal citation/DOI can be added here after publication.

## Data

Sequencing data have been deposited in the Gene Expression Omnibus under SuperSeries accession **GSE344075**, which will be made public upon publication.

## Contact

Please direct any requests for further information to the lead contact, Dr. Constantinos Chronis (chronis@uic.edu).

Questions about the code can be directed to Hammad Farooq (hfarooq4@jh.edu) and Aydolun Petenkaya (apeten2@uic.edu).

