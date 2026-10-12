# Sequence context and partner availability specify divergent cell fates from a single transcription factor

This repository contains computational analysis code accompanying the manuscript:

**Hammad Farooq<sup>#</sup>, Kenneth Joves<sup>#</sup>, Aydolun Petenkaya<sup>#</sup>, Itay Solomon<sup>#</sup>, Eva H.W. Leung, Sen Zhang, Mary Menhart, Ali Farhat, Lin Du, Nissim Hay, Sandra Pinho\*, and Constantinos Chronis\***  
<sup>#</sup> These authors contributed equally. \* Correspondence

**Sequence context and partner availability specify divergent cell fates from a single transcription factor**

This repository provides code for custom computational analyses used in the manuscript, including **ChromBPNet**, **TF-MoDISco/FiNeMo**, **ChromHMM**, ChromBPNet–ChromHMM integration, and the **CITE-seq and single-cell multiome** analyses associated with the included figures.

> **Repository scope.** This repository is intended to complement the Materials and Methods section by providing code used for custom analyses described in the manuscript, with the goal of improving transparency and reproducibility of the main results. It is not intended to serve as a fully automated or end-to-end executable workflow, and some analyses may require intermediate files, software environments, or computing resources specific to the original analysis setup.

## Repository structure

```text
code/
├── bulk/
│   ├── README.md                      # overview of bulk analysis workflows
│   ├── chrombpnet/                  # ChromBPNet → TF-MoDISco → FiNeMo workflow
│   ├── chromhmm/                    # stacked ATAC and histone ChromHMM workflows
│   └── integration/
│       └── chromhmm_chrombpnet/     # motif/state and promoter/TSS integration
└── single_cell/
    ├── README.md                    # single-cell inputs, scripts, and dependencies
    └── R/                           # CITE-seq and single-cell multiome figure scripts

docs/
├── FIGURE_CODE_MAP.md               # figure-to-code mapping
└── REPRODUCIBILITY.md               # configuration and reproducibility notes
```

Each analysis directory contains a local `README.md` describing its inputs, processing logic, and expected outputs.

## Quick start

### 1. Clone the repository

```bash
git clone <repository-url>
cd <repository-name>
```

### 2. Configure local paths

The repository does not contain laboratory- or cluster-specific absolute paths. Example path configuration files are provided at:

```text
code/bulk/chrombpnet/configs/paths.example.sh
code/bulk/chromhmm/configs/paths.example.sh
code/bulk/integration/chromhmm_chrombpnet/configs/paths.example.sh
```

Copy the relevant example to a local file, edit it for your system, and keep the local copy untracked by Git.

### 3. Follow the workflow-specific documentation

Start with:

- [`code/bulk/README.md`](code/bulk/README.md) for an overview of the bulk ChromBPNet, ChromHMM, and integration analyses.
- [`code/single_cell/README.md`](code/single_cell/README.md) for the CITE-seq and single-cell multiome analyses.
- [`docs/FIGURE_CODE_MAP.md`](docs/FIGURE_CODE_MAP.md) to locate code associated with the included figure panels.

Detailed workflow instructions are provided in the README files within the corresponding analysis directories.

## Inputs not distributed in this repository

Because of file size and data-distribution constraints, this repository does not include raw FASTQ files, processed BAM files, trained ChromBPNet models, contribution-score HDF5/bigWig files, or full intermediate ChromHMM/FiNeMo output directories. The code expects those inputs to be supplied through the documented configuration paths.

Data accessions and experimental details should be taken from the associated manuscript and its data-availability statement.

## Reproducibility and portability

The analysis was run across more than one computing environment, including local laboratory infrastructure and ALCF Aurora. Machine-specific paths have been removed from the public code and replaced with configurable path variables. Aurora-specific Parsl configuration is retained under `code/bulk/chrombpnet/aurora/` for provenance; users on other systems should adapt the launcher/executor layer while preserving the analysis commands and parameters.

See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) for details.

## Computing environment

ChromBPNet model training was performed on the Aurora supercomputer at the Argonne Leadership Computing Facility (ALCF).

Aurora system documentation:
https://docs.alcf.anl.gov/aurora/

Workflow orchestration used Parsl on Aurora:
https://docs.alcf.anl.gov/aurora/workflows/parsl/

The scripts in this repository reflect the original analysis environment and may require adaptation for other HPC systems, schedulers, or accelerator backends.

## Notebooks

Notebook outputs may be retained when they provide compact, publication-relevant QC or figure-generation provenance. Large or environment-specific outputs are omitted, and committed notebooks should not contain machine-specific paths or credentials.

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
