#!/usr/bin/env bash
# Example path configuration for the BULK ChromBPNet workflow.
# Copy to a local, untracked file (for example paths.local.sh) and edit for your system.

export PROJECT_ROOT="/path/to/CEBP_project"
export CHROMBPNET_ROOT="${PROJECT_ROOT}/chrombpnet"

# Prepared bulk ATAC inputs
export BAM_DIR="${CHROMBPNET_ROOT}/bam_filtered_final"
export PEAK_DIR="${CHROMBPNET_ROOT}/macs2_peaks_no_blacklist"
export NONPEAK_DIR="${CHROMBPNET_ROOT}/chrombpnet_nonpeak_background"
export SPLIT_DIR="${CHROMBPNET_ROOT}/splits"

# mm10 reference files
export REFERENCE_DIR="${CHROMBPNET_ROOT}/reference_data"
export GENOME_FASTA="${REFERENCE_DIR}/mm10.fa"
export CHROM_SIZES="${REFERENCE_DIR}/mm10.chrom.sizes"
export BLACKLIST_BED="${REFERENCE_DIR}/mm10.blacklist.bed.gz"

# Model/output roots
export BIAS_MODEL_DIR="${CHROMBPNET_ROOT}/bias_models"
export MODEL_DIR="${CHROMBPNET_ROOT}/models"
export INTERPRETATION_DIR="${CHROMBPNET_ROOT}/interpretation"
export MOTIF_DIR="${CHROMBPNET_ROOT}/motifs"

# Curated known-motif database used for annotation
export JASPAR_MEME="/path/to/curated_JASPAR_motifs.meme"

# Fold-averaged model outputs / downstream motif analysis
export MODEL_ROOT="${CHROMBPNET_ROOT}/chrombpnet_model_selected_bias"
export PREDICTION_OUTPUT_DIR="${CHROMBPNET_ROOT}/avg_pred_bw"
export MERGED_MOTIF_DIR="${CHROMBPNET_ROOT}/merged_motif_lexicon"
export KNOWN_MOTIF_LOGO_DIR="${MERGED_MOTIF_DIR}/known_motif_logos"
