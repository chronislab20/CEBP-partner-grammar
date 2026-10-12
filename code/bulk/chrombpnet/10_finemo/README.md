# 10 — FiNeMo motif-instance calling

`CEBPB_All_Merge_finemo.py` calls predictive motif instances from the merged motif lexicon using the manuscript parameters:

```text
finemo call-hits --global-lambda 0.7 -b 200 --compile
```

The analysis was performed on a local GPU server using NVIDIA RTX PRO 6000 Blackwell GPUs. The public script does not contain the original server path; set:

```bash
export MERGED_MOTIF_DIR=/path/to/merged_motif_lexicon
```

`CUDA_VISIBLE_DEVICES` may be set externally; if unset, the script defaults to device `0`.
