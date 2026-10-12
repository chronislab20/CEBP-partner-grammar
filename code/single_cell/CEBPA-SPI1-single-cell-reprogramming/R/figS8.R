source(here::here("R", "setup.R"))
source(here::here("R", "atac_helpers.R"))

##############FigS8a####################

atac$sample_plot <- factor(
  sub("_multiome$", "", as.character(atac$sample)),
  levels = sample_order
)

cluster_levels <- as.character(0:14)
atac$ATAC_cluster <- factor(as.character(atac$ATAC_cluster), levels = cluster_levels)
cluster_cols <- setNames(viridisLite::turbo(15), cluster_levels)

figS8a_sample <- DimPlot(
  atac,
  reduction = "umap.peaks",
  group.by = "sample_plot",
  cols = sample_cols,
  pt.size = 0.001,
  raster = FALSE
) +
  ggtitle("Samples") +
  labs(color = NULL) +
  centered_title

figS8a_cluster <- DimPlot(
  atac,
  reduction = "umap.peaks",
  group.by = "ATAC_cluster",
  cols = cluster_cols,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.01,
  raster = FALSE
) +
  ggtitle("ATAC clusters") +
  labs(color = NULL) +
  centered_title

ggsave(
  file.path(figS8_dir, "FigS8a_samples_and_clusters.pdf"),
  figS8a_sample / figS8a_cluster,
  width = 9,
  height = 12,
  useDingbats = FALSE
)

##############FigS8b####################

atac <- add_chromhmm_scores(atac)

save_state_pdf(
  c("St1", "St2", "St3", "St4", "St7", "St9", "St12", "St14", "St17"),
  file.path(figS8_dir, "FigS8b_remaining_ChromHMM_states.pdf")
)

message("Fig S8 panels written to ", figS8_dir)
