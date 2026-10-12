source(here::here("R", "setup.R"))
source(here::here("R", "atac_helpers.R"))

##############Fig8a####################

draw_order <- c(
  "Skeletal/mesenchyme",
  "Adipogenic transition (early)",
  "Sca1+Cd34+CD140+ Fibro-adipogenic Progenitor",
  "Fap+ activated adipogenic progenitor",
  "Residual Fibroblast",
  "Adhesion-remodel mesenchyme"
)

fig8a <- DimPlot(
  atac,
  reduction = "umap.peaks",
  group.by = "ATAC_final",
  cols = identity_cols[atac_order],
  order = draw_order,
  pt.size = 0.1,
  raster = FALSE
) +
  ggtitle("ATAC identities") +
  labs(color = NULL) +
  centered_title +
  theme(legend.title = element_blank())

ggsave(
  file.path(fig8_dir, "Fig8a_ATAC_identities.pdf"),
  fig8a,
  width = 10,
  height = 6,
  useDingbats = FALSE
)

##############Fig8c####################

atac <- add_chromhmm_scores(atac)

save_state_pdf(
  c("St5", "St6", "St8", "St13", "St15", "St16", "St10", "St11"),
  file.path(fig8_dir, "Fig8c_ChromHMM_states.pdf")
)

message("Fig 8 panels written to ", fig8_dir)
