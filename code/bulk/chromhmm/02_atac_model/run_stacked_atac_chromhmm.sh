#!/usr/bin/env bash
set -euo pipefail
: "${CHROMHMM_JAR:?Set CHROMHMM_JAR}"
: "${MM10_CHROMSIZES:?Set MM10_CHROMSIZES}"
: "${ATAC_BAM_DIR:?Set ATAC_BAM_DIR}"
: "${ATAC_CHROMHMM_OUT:?Set ATAC_CHROMHMM_OUT}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TABLE="${SCRIPT_DIR}/../configs/atac_cellmarkfiletable.txt"
OUTBIN="${ATAC_CHROMHMM_OUT}/bin_mm10_200"
mkdir -p "${ATAC_CHROMHMM_OUT}" "${OUTBIN}"

java -mx100G -jar "${CHROMHMM_JAR}" BinarizeBam -b 200 -mixed \
  "${MM10_CHROMSIZES}" "${ATAC_BAM_DIR}" "${TABLE}" "${OUTBIN}"

# Original final analysis script evaluated K=15..40; K=23 was selected.
for K in $(seq 15 40); do
  KDIR="${ATAC_CHROMHMM_OUT}/numstates_${K}"
  mkdir -p "${KDIR}"
  java -mx100G -jar "${CHROMHMM_JAR}" LearnModel -b 200 -p 64 -holdcolumnorder \
    "${OUTBIN}" "${KDIR}" "${K}" mm10
done
