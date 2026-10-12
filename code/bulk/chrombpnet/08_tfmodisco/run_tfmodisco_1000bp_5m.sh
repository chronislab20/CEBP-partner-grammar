#!/usr/bin/env bash
set -euo pipefail

# Final TF-MoDISco analysis used in the manuscript.
# Environment-specific GPU/CUDA settings should be configured before launching.
: "${MODEL_ROOT:?Set MODEL_ROOT to the directory containing the trained sample directories}"
: "${JASPAR_MEME:?Set JASPAR_MEME to the curated MEME-format motif database}"

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

shaptype="counts"
window=1000
max_seqlets=5000000

for sample in "${samples[@]}"; do
    out_dir="${MODEL_ROOT}/${sample}/Ave/modisco_${window}_${max_seqlets}/${shaptype}"
    mkdir -p "$out_dir"

    modisco motifs \
      -i "${MODEL_ROOT}/${sample}/Ave/average_shaps.${shaptype}.h5" \
      -n "$max_seqlets" \
      -l 2 \
      -w "$window" \
      --verbose \
      -o "${out_dir}/merged_allchrs.${shaptype}_modisco.h5"

    modisco report \
      -i "${out_dir}/merged_allchrs.${shaptype}_modisco.h5" \
      -o "${out_dir}/merged_allchrs.${shaptype}_modisco_report" \
      -s "${out_dir}/merged_allchrs.${shaptype}_modisco_report" \
      -m "$JASPAR_MEME" \
      -n 10
done
