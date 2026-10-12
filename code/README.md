# Analysis code

The public code is organized by analysis branch rather than by computing server.

- `bulk/chrombpnet/` — ATAC preprocessing for ChromBPNet, bias/model training, model evaluation, contribution scores, TF-MoDISco, merged motif lexicon, and FiNeMo.
- `bulk/chromhmm/` — stacked ATAC and histone ChromHMM state models plus state refinement.
- `bulk/integration/chromhmm_chrombpnet/` — integration of FiNeMo motif instances with ChromHMM states and RNA-defined promoter/TSS groups for Fig. 5c and Fig. S5e.
- `single_cell/` — the authors' existing CEBPA/SPI1 CITE-seq and single-cell multiome analysis archive supporting Figs. 7–8 and Figs. S7–S8. The supplied archive is retained unchanged.

See the README inside each directory for workflow-specific inputs, outputs, and scope.
