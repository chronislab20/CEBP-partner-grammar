# 08 — TF-MoDISco

This directory contains the TF-MoDISco workflow corresponding to the final manuscript analysis: **counts contribution scores**, a **1,000-bp search window**, and up to **5,000,000 seqlets**.

`run_tfmodisco_1000bp_5m.sh` operates on the 20 manuscript model directories. Set `MODEL_ROOT` and `JASPAR_MEME` before running. Machine-specific GPU/CUDA settings are intentionally not embedded in the public script.

The older exploratory 400-bp TF-MoDISco workflow and historical combined TF-MoDISco/FiNeMo task are intentionally excluded.
