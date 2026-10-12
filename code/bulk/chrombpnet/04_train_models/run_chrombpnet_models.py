#!/usr/bin/env python3
"""Self-contained Parsl launcher for one bulk ChromBPNet task.

This file keeps the task definition and the orchestration used to submit that
specific task together.  The only shared runtime dependency is the Aurora
Parsl executor configuration in ../aurora/parsl_config.py.

Adapt PROJECT_SUBDIR / paths before running outside the original Aurora layout.
"""
import os
import sys
from pathlib import Path
from datetime import datetime
import parsl
from parsl import bash_app

# Load the shared Aurora executor configuration.
HERE = Path(__file__).resolve().parent
AURORA_DIR = (HERE / "../aurora").resolve()
sys.path.insert(0, str(AURORA_DIR))
from parsl_config import aurora_single_tile_config

PROJECT_SUBDIR = "project_folder/CEBPB_seq/chrombpnet"
TMPDIR = f"/tmp/{os.environ.get('USER', 'parsl_user')}"
TODAY = datetime.today().strftime("%Y-%m-%d")
SELECTED_BIAS = {'Bcell_Graf_Cebpa18h_merged_mm10': ('Bcell_Graf_OKSMd1_merged_mm10', 0.7), 'Bcell_Graf_merged_mm10': ('Bcell_Graf_merged_mm10', 0.5), 'Bcell_Graf_OKSMd1_merged_mm10': ('Bcell_Graf_OKSMd1_merged_mm10', 0.7), 'Bcell_Graf_OKSMd2_merged_mm10': ('Bcell_Graf_OKSMd1_merged_mm10', 0.7), 'ESC_ATAC-Seq_mm10': ('ESC_ATAC-Seq_mm10', 0.7), 'ESC_d1_hCEBPA_ATAC_mm10': ('ESC_d1_hCEBPA_ATAC_mm10', 0.6), 'ESC_d2_hCEBPA_ATAC_mm10': ('ESC_d2_hCEBPA_ATAC_mm10', 0.7), 'MEF_ATAC_mm10': ('MEF_ATAC_mm10', 0.5), 'MEF_d2_mCEBPAflag_ATAC_mm10': ('MEF_d2_mCEBPAflag_ATAC_mm10', 0.5), 'MEF_d2_mCEBPBflag_ATAC_mm10': ('MEF_d2_mCEBPBflag_ATAC_mm10', 0.5), 'MEF_d2_OSKM_ATAC_L002_mm10': ('MEF_d2_OSKM_ATAC_L002_mm10', 0.5), 'MEF_d2_OSKM_hCEBPDflag_ATAC_mm10': ('MEF_d2_OSKM_hCEBPDflag_ATAC_mm10', 0.5), 'MEF_d2_OSKM_mCEBPAflag_ATAC_mm10': ('MEF_d2_OSKM_mCEBPAflag_ATAC_mm10', 0.5), 'MEF_d2_OSKM_mCEBPBflag_ATAC_mm10': ('MEF_d2_OSKM_mCEBPBflag_ATAC_mm10', 0.5), 'NSC_ATAC_S23_L005_mm10': ('NSC_ATAC_S23_L005_mm10', 0.6), 'NSC_d2_mCEBPAflag_ATAC_S15_mm10': ('NSC_d2_mCEBPAflag_ATAC_S15_mm10', 0.6), 'NSC_d2_mCEBPBflag_ATAC_mm10': ('NSC_d2_mCEBPBflag_ATAC_mm10', 0.6), 'NSC_d2_OSKM_ATAC_mm10': ('NSC_d2_OSKM_ATAC_mm10', 0.6), 'NSC_d2_OSKM_mCEBPAflag_ATAC_mm10': ('NSC_d2_OSKM_mCEBPAflag_ATAC_mm10', 0.6), 'NSC_d2_OSKM_mCEBPBflag_ATAC_mm10': ('NSC_d2_OSKM_mCEBPBflag_ATAC_mm10', 0.6)}

FOLDS = range(5)

@bash_app
def train_chrombpnet_model(chrombpnet_dir, sample_name, fold, bias_sample, bias_threshold, tmp_dir, stdout=None, stderr=None):
    return f"""
    mkdir -p {tmp_dir}
    export TMPDIR={tmp_dir}
    chrombpnet pipeline \
      -ibam $HOME/{chrombpnet_dir}/bam_filtered_final/{sample_name}.filtered.bam \
      -d ATAC \
      -g $HOME/{chrombpnet_dir}/reference_data/mm10.fa \
      -c $HOME/{chrombpnet_dir}/reference_data/mm10.chrom.sizes \
      -p $HOME/{chrombpnet_dir}/macs2_peaks_no_blacklist/{sample_name}_peaks.top200000.no_blacklist.bed \
      -n $HOME/{chrombpnet_dir}/chrombpnet_nonpeak_background/{sample_name}_peaks.top200000_fold_{fold}_nonpeaks_negatives.bed \
      -fl $HOME/{chrombpnet_dir}/splits/fold_{fold}.json \
      -b $HOME/{chrombpnet_dir}/6_hrs_bias_models_{bias_threshold}/{bias_sample}/fold_{fold}/models/bias.h5 \
      -o $HOME/{chrombpnet_dir}/12_hrs_chrombpnet_model_selected_bias/{sample_name}_Train_{bias_threshold}/fold_{fold}/
    """

if __name__ == "__main__":
    log_dir = HERE / "logs" / TODAY
    log_dir.mkdir(parents=True, exist_ok=True)
    futures = []
    with parsl.load(aurora_single_tile_config):
        for sample, (bias_sample, threshold) in SELECTED_BIAS.items():
            for fold in FOLDS:
                tag = f"{sample}_fold{fold}"
                futures.append(train_chrombpnet_model(PROJECT_SUBDIR, sample, fold, bias_sample, threshold, TMPDIR,
                    stdout=str(log_dir/f"stdout_{tag}.log"), stderr=str(log_dir/f"stderr_{tag}.log")))
        for i, fut in enumerate(futures):
            try:
                result = fut.result()
                print(f"✅ Task {i} completed successfully.")
            except Exception as e:
                print(f"❌ Task {i} failed: {type(e).__name__}: {e}")

        print("All chrombpnet tasks completed.")
