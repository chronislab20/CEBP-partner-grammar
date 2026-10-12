#!/usr/bin/env bash
set -euo pipefail

# Public/repository version of the bulk ChromBPNet input-preparation workflow.
# Server-specific absolute paths have been replaced by configurable variables.
# Analysis parameters are those used for the manuscript.
#
# Usage:
#   export WORK_ROOT=/path/to/chrombpnet_work
#   export RAW_BAM_DIR=/path/to/input_atac_bams
#   export CONDA_SH=/path/to/conda.sh          # optional, if MACS2 is in a conda env
#   export PEAK_ENV=scprinter                  # optional
#   bash prepare_bulk_chrombpnet_inputs.sh

: "${WORK_ROOT:=/path/to/chrombpnet_work}"
: "${RAW_BAM_DIR:=/path/to/input_atac_bams}"

REFERENCE_DIR="${WORK_ROOT}/reference_data"
SORTED_BAM_DIR="${WORK_ROOT}/bam_sorted_indexed"
FILTERED_BAM_DIR="${WORK_ROOT}/bam_filtered_final"
RELAXED_PEAK_DIR="${WORK_ROOT}/macs2_relaxed_peaks"
FINAL_PEAK_DIR="${WORK_ROOT}/macs2_peaks_no_blacklist"
NONPEAK_DIR="${WORK_ROOT}/chrombpnet_nonpeak_background"
SPLIT_DIR="${WORK_ROOT}/splits"

mkdir -p "$REFERENCE_DIR" "$SORTED_BAM_DIR" "$FILTERED_BAM_DIR" \
         "$RELAXED_PEAK_DIR" "$FINAL_PEAK_DIR" "$NONPEAK_DIR"

# -----------------------------------------------------------------------------
# 1. Reference data (mm10)
# -----------------------------------------------------------------------------
wget -O "${REFERENCE_DIR}/mm10.fa.gz" \
  "https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz"
gunzip -f "${REFERENCE_DIR}/mm10.fa.gz"

wget -O "${REFERENCE_DIR}/mm10.chrom.sizes" \
  "http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes"

wget -O "${REFERENCE_DIR}/mm10.blacklist.bed.gz" \
  "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"

# -----------------------------------------------------------------------------
# 2. Sort and index BAM files
# -----------------------------------------------------------------------------
for bam_file in "${RAW_BAM_DIR}"/*.bam; do
    base_name=$(basename "$bam_file" .bam)
    sorted_bam="${SORTED_BAM_DIR}/${base_name}.sorted.bam"
    samtools sort -@40 "$bam_file" -o "$sorted_bam"
    samtools index "$sorted_bam"
done

# -----------------------------------------------------------------------------
# 3. Filter BAM files
#    -F 780 -q 30 are the manuscript filtering settings.
# -----------------------------------------------------------------------------
for bam_file in "${SORTED_BAM_DIR}"/*.sorted.bam; do
    base_name=$(basename "$bam_file" .sorted.bam)
    out_bam="${FILTERED_BAM_DIR}/${base_name}.filtered.bam"

    # This single-end ESC file was handled separately in the original workflow.
    if [[ "$bam_file" == *"ESC_ATAC-Seq_mm10.bam"* ]]; then
        echo "Skipping single-end file in paired-end loop: $bam_file"
        continue
    fi

    samtools view -b -@50 -F 780 -q 30 "$bam_file" > "$out_bam"
    samtools index "$out_bam"
done

# -----------------------------------------------------------------------------
# 4. Relaxed MACS2 peaks; retain top 200,000 by signalValue
# -----------------------------------------------------------------------------
if [[ -n "${CONDA_SH:-}" ]]; then
    # shellcheck disable=SC1090
    source "$CONDA_SH"
    if [[ -n "${PEAK_ENV:-}" ]]; then
        conda activate "$PEAK_ENV"
    fi
fi

pval_thresh=0.01
NPEAKS=200000
smooth_window=150
shiftsize=-75
gensz="mm"

for bam_file in "${FILTERED_BAM_DIR}"/*.filtered.bam; do
    base_name=$(basename "$bam_file" .filtered.bam)

    macs2 callpeak \
        -t "$bam_file" \
        -f BAM \
        -g "$gensz" \
        -p "$pval_thresh" \
        --shift "$shiftsize" \
        --extsize "$smooth_window" \
        --nomodel \
        --keep-dup all \
        --call-summits \
        -n "$base_name" \
        --outdir "$RELAXED_PEAK_DIR"

    peak_file="${RELAXED_PEAK_DIR}/${base_name}_peaks.narrowPeak"
    capped_file="${RELAXED_PEAK_DIR}/${base_name}_peaks.top${NPEAKS}.narrowPeak"
    sort -k7,7nr "$peak_file" | head -n "$NPEAKS" > "$capped_file"
done

# -----------------------------------------------------------------------------
# 5. Expand ENCODE blacklist by 1,057 bp and remove overlapping peaks
# -----------------------------------------------------------------------------
BLACKLIST_EXTENDED="${REFERENCE_DIR}/mm10.blacklist.extended1057.bed"
bedtools slop \
  -i "${REFERENCE_DIR}/mm10.blacklist.bed.gz" \
  -g "${REFERENCE_DIR}/mm10.chrom.sizes" \
  -b 1057 \
  > "$BLACKLIST_EXTENDED"

for capped_file in "${RELAXED_PEAK_DIR}"/*_peaks.top*.narrowPeak; do
    base_name=$(basename "$capped_file" .narrowPeak)
    filtered_file="${FINAL_PEAK_DIR}/${base_name}.no_blacklist.bed"
    bedtools intersect -v -a "$capped_file" -b "$BLACKLIST_EXTENDED" > "$filtered_file"
done

# -----------------------------------------------------------------------------
# 6. Generate GC-matched nonpeak regions for each chromosome fold
# -----------------------------------------------------------------------------
# The exact five fold JSON files used in the manuscript are included under
# configs/folds/. Copy or point SPLIT_DIR to those files before this step.
for peak_file in "${FINAL_PEAK_DIR}"/*.no_blacklist.bed; do
    peak_base=$(basename "$peak_file" .no_blacklist.bed)

    for fold_file in "${SPLIT_DIR}"/fold_*.json; do
        fold_name=$(basename "$fold_file" .json)
        out_dir="${NONPEAK_DIR}/${peak_base}_${fold_name}_nonpeaks"
        mkdir -p "$out_dir"

        chrombpnet prep nonpeaks \
            -g "${REFERENCE_DIR}/mm10.fa" \
            -p "$peak_file" \
            -c "${REFERENCE_DIR}/mm10.chrom.sizes" \
            -fl "$fold_file" \
            -br "${REFERENCE_DIR}/mm10.blacklist.bed.gz" \
            -o "$out_dir"
    done
done
