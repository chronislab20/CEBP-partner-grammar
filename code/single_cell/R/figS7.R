source(here::here("R", "setup.R"))
source(here::here("R", "citeseq_helpers.R"))

##############FigS7a####################

check_inputs(c("human_mouse", "zsgreen", "flag"))

human_mouse <- readRDS(paths$human_mouse)
zsgreen_obj <- readRDS(paths$zsgreen)
flag_obj <- readRDS(paths$flag)

human_mouse_umap <- as.data.frame(Embeddings(human_mouse, "rna.adt.wnn.umap")[, 1:2])
colnames(human_mouse_umap) <- c("WNNUMAP_1", "WNNUMAP_2")

mouse_data <- LayerData(human_mouse, assay = "RNA", layer = "data")
human_data <- LayerData(human_mouse, assay = "human_RNA", layer = "data")

plot_species_gene <- function(gene, species, title) {
  mat <- if (species == "human") human_data else mouse_data
  check_features(gene, mat, paste(species, "genes"))

  df <- human_mouse_umap
  df$expr <- as.numeric(mat[gene, rownames(df)])
  df <- arrange(df, expr)

  ggplot(df, aes(WNNUMAP_1, WNNUMAP_2, color = expr)) +
    geom_point(size = 0.001) +
    scale_color_gradientn(colours = rocket_cols) +
    coord_equal() +
    labs(title = title, color = NULL) +
    umap_theme +
    centered_title
}

reporters <- umap
reporters$ZsGreen <- 0
reporters$Flag <- 0

zsgreen_idx <- match(rownames(reporters), rownames(zsgreen_obj@meta.data))
flag_idx <- match(rownames(reporters), rownames(flag_obj@meta.data))

reporters$ZsGreen[!is.na(zsgreen_idx)] <- zsgreen_obj$ZsGreen_count[zsgreen_idx[!is.na(zsgreen_idx)]]
reporters$Flag[!is.na(flag_idx)] <- flag_obj$FLAG_hCEBPA_count[flag_idx[!is.na(flag_idx)]]

plot_reporter <- function(column, colour, title) {
  ggplot(reporters, aes(WNNUMAP_1, WNNUMAP_2)) +
    geom_point(color = "grey95", size = 0.001) +
    geom_point(
      data = filter(reporters, .data[[column]] > 0),
      color = colour,
      size = 0.001
    ) +
    coord_equal() +
    ggtitle(title) +
    umap_theme +
    centered_title
}

figS7a <- (
  plot_species_gene("CEBPA", "human", "hCEBPα") +
    plot_species_gene("SPI1", "human", "hSPI1") +
    plot_reporter("ZsGreen", "#35A853", "ZsGreen")
) / (
  plot_species_gene("Cebpa", "mouse", "mCEBPα") +
    plot_species_gene("Spi1", "mouse", "mSPI1") +
    plot_reporter("Flag", "#CF3046", "Flag")
)

ggsave(
  file.path(figS7_dir, "FigS7a_transgenes_reporters.pdf"),
  figS7a,
  width = 12,
  height = 8,
  device = cairo_pdf
)

##############FigS7b####################

figS7b_umap <- DimPlot(
  seu,
  reduction = "rna.adt.wnn.umap",
  group.by = "wnn_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.001,
  raster = FALSE
) +
  NoLegend() +
  umap_theme

cluster_order <- c(
  "28", "10", "31", "26", "25", "0", "30", "3", "2", "12",
  "21", "6", "32", "19", "22", "1", "37", "34", "4", "36",
  "33", "27", "18", "16", "15", "24", "29", "8", "9", "7",
  "23", "5", "35", "20", "17", "14", "13", "11"
)

figS7b_df <- data.frame(
  cluster = as.character(seu$wnn_clusters),
  sample = as.character(seu$sample)
) %>%
  count(cluster, sample, name = "n_cells") %>%
  complete(cluster = cluster_order, sample = sample_order, fill = list(n_cells = 0)) %>%
  group_by(cluster) %>%
  mutate(pct = 100 * n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    cluster = factor(cluster, levels = rev(cluster_order)),
    sample = factor(sample, levels = sample_order)
  )

figS7b_bar <- ggplot(figS7b_df, aes(cluster, pct, fill = sample)) +
  geom_col(width = 0.9) +
  geom_text(
    aes(label = ifelse(pct >= 5, paste0(round(pct), "%"), "")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 2.3
  ) +
  scale_fill_manual(values = sample_cols, labels = sample_labels) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL, y = "Cells (%)", fill = NULL) +
  theme_classic() +
  theme(axis.ticks.y = element_blank())

ggsave(file.path(figS7_dir, "FigS7b_WNN_clusters.pdf"), figS7b_umap, width = 6, height = 6)
ggsave(file.path(figS7_dir, "FigS7b_cluster_sample_composition.pdf"), figS7b_bar, width = 6, height = 8)

##############FigS7c-f####################

save_gene_panel <- function(genes, file, width, height, ...) {
  ggsave(
    file.path(figS7_dir, file),
    wrap_plots(lapply(genes, plot_gene), ...),
    width = width,
    height = height
  )
}

save_gene_panel(
  c("Snai1", "Meox1", "Twist1", "Scx", "Hmga2", "Fosb"),
  "FigS7c_fibroblast_mesenchymal_markers.pdf",
  width = 15, height = 3, nrow = 1
)

save_gene_panel(
  c("Mki67", "Top2a"),
  "FigS7d_proliferation_markers.pdf",
  width = 4, height = 7, ncol = 1
)

save_gene_panel(
  c("Ucp1", "Dio2", "Sox30", "Pparg"),
  "FigS7e_adipogenic_markers.pdf",
  width = 10, height = 3, nrow = 1
)

save_gene_panel(
  c("Cd74", "H2-Eb1", "Mrc1", "Fcer1g"),
  "FigS7f_myeloid_markers.pdf",
  width = 10, height = 3, nrow = 1
)

##############FigS7g####################

rna_pairs <- c("Ptprc", "Itgam", "Clec12a", "Cd38", "Fcgr3", "Cd86", "Kit", "Cd34", "Pdgfra", "Cd48")
adt_pairs <- c("CD45", "Mac1", "CD371", "CD38", "CD16-CD32", "CD86", "cKit", "CD34", "CD140a", "CD48")

adt_data <- LayerData(seu, assay = "cellbender_adt_CLR", layer = "data")
check_features(adt_pairs, adt_data, "ADT features")

adt_cols <- c("#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#084594")

plot_adt <- function(feature) {
  df <- umap
  df$expr <- as.numeric(adt_data[feature, rownames(df)])
  df$expr[df$expr < 1] <- 1
  df <- arrange(df, expr)

  ggplot(df, aes(WNNUMAP_1, WNNUMAP_2, color = expr)) +
    geom_point(size = 0.001) +
    scale_color_gradientn(colours = adt_cols) +
    coord_equal() +
    labs(title = feature, color = NULL) +
    umap_theme +
    centered_title
}

figS7g <- wrap_plots(
  c(lapply(rna_pairs, plot_gene), lapply(adt_pairs, plot_adt)),
  ncol = 10
)

ggsave(file.path(figS7_dir, "FigS7g_RNA_ADT_pairs.pdf"), figS7g, width = 30, height = 10)

##############FigS7h####################

figS7h_populations <- c(
  "Dual-prog transition",
  "Early adipocyte",
  "Pparg+ Adipo-mac hybrid",
  "Cycling myeloid-like",
  "Inflammatory Mac-like"
)

day_levels <- c("D3", "D6", "D10")
day_cols <- c("D3" = "#C0B1D8", "D6" = "#8E75B8", "D10" = "#4B376C")

figS7h_df <- seu@meta.data %>%
  filter(sample == "MEF_CEBPA_SPI1", RNA_final %in% figS7h_populations) %>%
  count(RNA_final, timepoint, name = "n_cells") %>%
  complete(RNA_final = figS7h_populations, timepoint = day_levels, fill = list(n_cells = 0)) %>%
  group_by(RNA_final) %>%
  mutate(pct = 100 * n_cells / sum(n_cells)) %>%
  ungroup() %>%
  mutate(
    RNA_final = factor(RNA_final, levels = figS7h_populations),
    timepoint = factor(timepoint, levels = day_levels)
  )

figS7h <- ggplot(figS7h_df, aes(RNA_final, pct, fill = timepoint)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = ifelse(pct >= 5, paste0(round(pct), "%"), "")),
    position = position_stack(vjust = 0.5),
    color = "white",
    size = 3
  ) +
  scale_fill_manual(values = day_cols) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  labs(x = NULL, y = "% of population", fill = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks.x = element_blank()
  )

ggsave(
  file.path(figS7_dir, "FigS7h_CEBPA_SPI1_population_timepoints.pdf"),
  figS7h,
  width = 7,
  height = 6
)

##############FigS7i####################

check_inputs("metascape_figS7i")

figS7i_df <- read_excel(paths$metascape_figS7i, sheet = "Plot_data")
figS7i_terms <- read_excel(paths$metascape_figS7i, sheet = "Term_plot_order")

figS7i_order <- c(
  "Residual Fibroblast",
  "Adhesion-remodel mesenchyme",
  "Cycling",
  "Adipogenic transition (early)",
  "Chromatin remodeling-high adipogenic",
  "Cycling myeloid-like",
  "CMP-like progenitor",
  "Myeloid-like",
  "Skeletal/mesenchyme"
)

warn_unmatched(figS7i_df$New_identity, figS7i_order, "Fig S7i population")
warn_unmatched(figS7i_df$Plot_term, figS7i_terms$Plot_term, "Fig S7i GO term")

figS7i_df$New_identity <- factor(figS7i_df$New_identity, levels = figS7i_order)
figS7i_df$Plot_term <- factor(figS7i_df$Plot_term, levels = rev(figS7i_terms$Plot_term))

ggsave(
  file.path(figS7_dir, "FigS7i_remaining_population_GO.pdf"),
  plot_go_dots(figS7i_df),
  width = 10,
  height = 10
)

##############FigS7j####################

ggsave(
  file.path(figS7_dir, "FigS7j_Il12a_UMAP.pdf"),
  plot_gene("Il12a"),
  width = 6,
  height = 6
)

message("Fig S7 panels written to ", figS7_dir)
