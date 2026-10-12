"""Call FiNeMo motif instances from the merged manuscript motif lexicon.

Public/repository version: local absolute paths are supplied through environment
variables rather than embedded in the script.
"""

import os
import subprocess
from pathlib import Path

pattern_group = "pos_patterns"
project_path = Path(os.environ.get("MERGED_MOTIF_DIR", "/path/to/merged_motif_lexicon"))
pattern_dir = project_path / pattern_group
pattern_dir.mkdir(parents=True, exist_ok=True)

merged_npz_file = pattern_dir / "merged_allchrs.countsshap_scores_inputs.npz"
merged_modisco_h5_path = pattern_dir / "merged_modisco.h5"
output_prefix = pattern_dir / "finemo/merged"
output_prefix.parent.mkdir(parents=True, exist_ok=True)

env = os.environ.copy()
# Override externally if another device should be used.
env.setdefault("CUDA_VISIBLE_DEVICES", "0")

subprocess.run([
    "finemo", "call-hits",
    "-r", str(merged_npz_file),
    "-m", str(merged_modisco_h5_path),
    "--global-lambda", "0.7",
    "-b", "200",
    "--compile",
    "-o", str(output_prefix),
], check=True, env=env)

subprocess.run([
    "finemo", "report",
    "-r", str(merged_npz_file),
    "-m", str(merged_modisco_h5_path),
    "--no-recall",
    "--no-seqlets",
    "-o", str(output_prefix),
    "-H", str(output_prefix),
], check=True, env=env)
