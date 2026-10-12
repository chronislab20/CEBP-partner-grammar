# Reproducibility notes

## Design principle

The repository separates **analysis logic** from **machine-specific configuration**. Paths to local sequencing data, reference files, model directories, and cluster scratch space should be supplied through configuration files or environment variables rather than committed as absolute paths.

## Computing environments

The manuscript analyses were run across multiple systems. In particular, ChromBPNet model training and interpretation were run on ALCF Aurora using Parsl launchers, while some preprocessing and downstream integration steps were run on laboratory compute infrastructure.

The Aurora-specific executor definition under `code/bulk/chrombpnet/aurora/` is retained to document the original execution environment. Users on another scheduler or workstation do not need to reproduce the Parsl layer exactly; the scientific commands, inputs, parameters, and fold definitions are the reproducibility-critical components.

## Local path configuration

Example path files are provided under each workflow:

```text
code/bulk/chrombpnet/configs/paths.example.sh
code/bulk/chromhmm/configs/paths.example.sh
code/bulk/integration/chromhmm_chrombpnet/configs/paths.example.sh
```

Recommended practice:

```bash
cp code/bulk/chrombpnet/configs/paths.example.sh \
   code/bulk/chrombpnet/configs/paths.local.sh

# edit paths.local.sh, then
source code/bulk/chrombpnet/configs/paths.local.sh
```

Files named `*.local.sh` are excluded by the repository `.gitignore`.

## ChromBPNet cross-validation

The exact five chromosome-split JSON files distributed with this repository should be used when reproducing the reported models. They are stored at:

```text
code/bulk/chrombpnet/configs/folds/
```

Each fold assigns non-overlapping chromosomes to training, validation, and test partitions.

## Large intermediate files

The following are not suitable for Git and are therefore external inputs/outputs:

- FASTQ/BAM/BAI files;
- bigWig files;
- ChromBPNet `.h5` model files;
- DeepLIFT/DeepSHAP HDF5 files;
- large TF-MoDISco and FiNeMo result directories;
- full ChromHMM emission/segmentation output directories;
- generated figure PDFs/PNGs unless intentionally selected for release.

The `.gitignore` excludes the most common forms of these files.

## Software dependencies

Major external tools used by the included workflows include:

- ChromBPNet
- TensorFlow backend compatible with ChromBPNet
- TF-MoDISco-lite
- FiNeMo
- ChromHMM
- BEDTools
- SAMtools
- MACS2
- Python 3 with NumPy/Pandas/HDF5-related packages used by the helper scripts
- R with the CRAN/Bioconductor packages imported by the analysis scripts
- Parsl for the original Aurora launchers

ChromHMM v1.27 and STAR v2.7.10b are explicitly recorded in the manuscript Methods. Exact versions for other dependencies should be matched to the environments used for the corresponding analyses where available.
