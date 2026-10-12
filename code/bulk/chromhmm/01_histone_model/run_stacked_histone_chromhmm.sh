#!/usr/bin/env bash
set -euo pipefail
: "${CHROMHMM_JAR:?Set CHROMHMM_JAR}"
: "${MM10_CHROMSIZES:?Set MM10_CHROMSIZES}"
: "${HISTONE_BAM_DIR:?Set HISTONE_BAM_DIR}"
: "${HISTONE_CHROMHMM_OUT:?Set HISTONE_CHROMHMM_OUT}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TABLE="${SCRIPT_DIR}/../configs/histone_cellmarkfiletable.txt"
OUTBIN="${HISTONE_CHROMHMM_OUT}/bin_mm10_200"
mkdir -p "${HISTONE_CHROMHMM_OUT}" "${OUTBIN}"

java -mx100G -jar "${CHROMHMM_JAR}" BinarizeBam -b 200 -mixed \
  "${MM10_CHROMSIZES}" "${HISTONE_BAM_DIR}" "${TABLE}" "${OUTBIN}"

# Original final analysis script evaluated K=10..40; K=22 was selected.
for K in $(seq 10 40); do
  KDIR="${HISTONE_CHROMHMM_OUT}/numstates_${K}"
  mkdir -p "${KDIR}"
  java -mx100G -jar "${CHROMHMM_JAR}" LearnModel -b 200 -p 64 -holdcolumnorder \
    "${OUTBIN}" "${KDIR}" "${K}" mm10
done
