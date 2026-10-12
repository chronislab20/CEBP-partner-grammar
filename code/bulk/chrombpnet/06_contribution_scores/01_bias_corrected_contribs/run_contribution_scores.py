#!/usr/bin/env python3
"""Self-contained Parsl launcher for one bulk ChromBPNet task.

This file keeps the task definition and the orchestration used to submit that
specific task together.  The only shared runtime dependency is the Aurora
Parsl executor configuration in ../../aurora/parsl_config.py.

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
AURORA_DIR = (HERE / "../../aurora").resolve()
sys.path.insert(0, str(AURORA_DIR))
from parsl_config import aurora_single_tile_config

PROJECT_SUBDIR = "project_folder/CEBPB_seq/chrombpnet"
TMPDIR = f"/tmp/{os.environ.get('USER', 'parsl_user')}"
TODAY = datetime.today().strftime("%Y-%m-%d")
SELECTED_BIAS = {'Bcell_Graf_Cebpa18h_merged_mm10': ('Bcell_Graf_OKSMd1_merged_mm10', 0.7), 'Bcell_Graf_merged_mm10': ('Bcell_Graf_merged_mm10', 0.5), 'Bcell_Graf_OKSMd1_merged_mm10': ('Bcell_Graf_OKSMd1_merged_mm10', 0.7), 'Bcell_Graf_OKSMd2_merged_mm10': ('Bcell_Graf_OKSMd1_merged_mm10', 0.7), 'ESC_ATAC-Seq_mm10': ('ESC_ATAC-Seq_mm10', 0.7), 'ESC_d1_hCEBPA_ATAC_mm10': ('ESC_d1_hCEBPA_ATAC_mm10', 0.6), 'ESC_d2_hCEBPA_ATAC_mm10': ('ESC_d2_hCEBPA_ATAC_mm10', 0.7), 'MEF_ATAC_mm10': ('MEF_ATAC_mm10', 0.5), 'MEF_d2_mCEBPAflag_ATAC_mm10': ('MEF_d2_mCEBPAflag_ATAC_mm10', 0.5), 'MEF_d2_mCEBPBflag_ATAC_mm10': ('MEF_d2_mCEBPBflag_ATAC_mm10', 0.5), 'MEF_d2_OSKM_ATAC_L002_mm10': ('MEF_d2_OSKM_ATAC_L002_mm10', 0.5), 'MEF_d2_OSKM_hCEBPDflag_ATAC_mm10': ('MEF_d2_OSKM_hCEBPDflag_ATAC_mm10', 0.5), 'MEF_d2_OSKM_mCEBPAflag_ATAC_mm10': ('MEF_d2_OSKM_mCEBPAflag_ATAC_mm10', 0.5), 'MEF_d2_OSKM_mCEBPBflag_ATAC_mm10': ('MEF_d2_OSKM_mCEBPBflag_ATAC_mm10', 0.5), 'NSC_ATAC_S23_L005_mm10': ('NSC_ATAC_S23_L005_mm10', 0.6), 'NSC_d2_mCEBPAflag_ATAC_S15_mm10': ('NSC_d2_mCEBPAflag_ATAC_S15_mm10', 0.6), 'NSC_d2_mCEBPBflag_ATAC_mm10': ('NSC_d2_mCEBPBflag_ATAC_mm10', 0.6), 'NSC_d2_OSKM_ATAC_mm10': ('NSC_d2_OSKM_ATAC_mm10', 0.6), 'NSC_d2_OSKM_mCEBPAflag_ATAC_mm10': ('NSC_d2_OSKM_mCEBPAflag_ATAC_mm10', 0.6), 'NSC_d2_OSKM_mCEBPBflag_ATAC_mm10': ('NSC_d2_OSKM_mCEBPBflag_ATAC_mm10', 0.6)}

FOLDS = range(5)
PROFILE_OR_COUNTS = ["counts"]

@bash_app
def contribution_scores(chrombpnet_dir, sample_name, fold, profile_or_counts, bias_threshold, tmp_dir, stdout=None, stderr=None):
    return f"""
    mkdir -p {tmp_dir}
    export TMPDIR={tmp_dir}
    mkdir -p $HOME/{chrombpnet_dir}/12_hrs_chrombpnet_model_selected_bias/{sample_name}_Train_{bias_threshold}/fold_{fold}/Pr_peaks_contribs_bw
    chrombpnet contribs_bw \
      -g $HOME/{chrombpnet_dir}/reference_data/mm10.fa \
      --model-h5 $HOME/{chrombpnet_dir}/12_hrs_chrombpnet_model_selected_bias/{sample_name}_Train_{bias_threshold}/fold_{fold}/models/chrombpnet_nobias.h5 \
      --output-prefix $HOME/{chrombpnet_dir}/12_hrs_chrombpnet_model_selected_bias/{sample_name}_Train_{bias_threshold}/fold_{fold}/Pr_peaks_contribs_bw/peaks_shape \
      --regions $HOME/{chrombpnet_dir}/macs2_peaks_no_blacklist/{sample_name}_peaks.top200000.no_blacklist.bed \
      -c $HOME/{chrombpnet_dir}/reference_data/mm10.chrom.sizes \
      --profile-or-counts {profile_or_counts}
    """

if __name__ == "__main__":
    log_dir = HERE / "logs" / TODAY
    log_dir.mkdir(parents=True, exist_ok=True)
    futures = []
    with parsl.load(aurora_single_tile_config):
        for sample, (_, threshold) in SELECTED_BIAS.items():
            for fold in FOLDS:
                for mode in PROFILE_OR_COUNTS:
                    tag = f"{sample}_fold{fold}_{mode}"
                    futures.append(contribution_scores(PROJECT_SUBDIR, sample, fold, mode, threshold, TMPDIR,
                        stdout=str(log_dir/f"stdout_{tag}.log"), stderr=str(log_dir/f"stderr_{tag}.log")))
        for i, fut in enumerate(futures):
            try:
                result = fut.result()
                print(f"✅ Task {i} completed successfully.")
            except Exception as e:
                print(f"❌ Task {i} failed: {type(e).__name__}: {e}")

        print("All chrombpnet tasks completed.")
