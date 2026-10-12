#!/usr/bin/env bash
set -euo pipefail

# Average the five fold-specific bias-corrected ChromBPNet prediction bigWigs.
: "${MODEL_ROOT:?Set MODEL_ROOT to the directory containing the trained sample directories}"
: "${CHROM_SIZES:?Set CHROM_SIZES to the mm10 chromosome sizes file}"
: "${PREDICTION_OUTPUT_DIR:=${MODEL_ROOT}/avg_pred_bw}"

samples=(
  "Bcell_Graf_Cebpa18h_merged_mm10_Train_0.7"
  "Bcell_Graf_merged_mm10_Train_0.5"
  "Bcell_Graf_OKSMd1_merged_mm10_Train_0.7"
  "Bcell_Graf_OKSMd2_merged_mm10_Train_0.7"
  "ESC_ATAC-Seq_mm10_Train_0.7"
  "ESC_d1_hCEBPA_ATAC_mm10_Train_0.6"
  "ESC_d2_hCEBPA_ATAC_mm10_Train_0.7"
  "MEF_ATAC_mm10_Train_0.5"
  "MEF_d2_mCEBPAflag_ATAC_mm10_Train_0.5"
  "MEF_d2_mCEBPBflag_ATAC_mm10_Train_0.5"
  "MEF_d2_OSKM_ATAC_L002_mm10_Train_0.5"
  "MEF_d2_OSKM_hCEBPDflag_ATAC_mm10_Train_0.5"
  "MEF_d2_OSKM_mCEBPAflag_ATAC_mm10_Train_0.5"
  "MEF_d2_OSKM_mCEBPBflag_ATAC_mm10_Train_0.5"
  "NSC_ATAC_S23_L005_mm10_Train_0.6"
  "NSC_d2_mCEBPAflag_ATAC_S15_mm10_Train_0.6"
  "NSC_d2_mCEBPBflag_ATAC_mm10_Train_0.6"
  "NSC_d2_OSKM_ATAC_mm10_Train_0.6"
  "NSC_d2_OSKM_mCEBPAflag_ATAC_mm10_Train_0.6"
  "NSC_d2_OSKM_mCEBPBflag_ATAC_mm10_Train_0.6"
)

mkdir -p "$PREDICTION_OUTPUT_DIR"
for sample_name in "${samples[@]}"; do
    bw_files=()
    for fold in {0..4}; do
        input_file="${MODEL_ROOT}/${sample_name}/fold_${fold}/pred_bw_cmb_chrombpnet_nobias.bw"
        [[ -f "$input_file" ]] && bw_files+=("$input_file") || echo "Missing: $input_file" >&2
    done

    if [[ ${#bw_files[@]} -ne 5 ]]; then
        echo "Expected 5 fold bigWigs for ${sample_name}; skipping." >&2
        continue
    fi

    output_bedgraph="${PREDICTION_OUTPUT_DIR}/${sample_name}_chrombpnet_nobias.bedGraph"
    output_bw="${PREDICTION_OUTPUT_DIR}/${sample_name}_chrombpnet_nobias.bw"

    wiggletools write_bg - mean "${bw_files[@]}" | sort -k1,1 -k2,2n > "$output_bedgraph"
    bedGraphToBigWig "$output_bedgraph" "$CHROM_SIZES" "$output_bw"
    rm -f "$output_bedgraph"
done
