source(here::here("R", "setup.R"))
source(here::here("R", "citeseq_helpers.R"))

##############Fig7a####################

timepoint_cols <- list(
  MEF_CEBPA = c("D3" = "#ABDED8", "D6" = "#69C3B9", "D10" = "#21544E"),
  MEF_CEBPA_SPI1 = c("D3" = "#C0B1D8", "D6" = "#8E75B8", "D10" = "#4B376C"),
  MEF_SPI1 = c("D3" = "#FFD3A1", "D6" = "#D99F5D", "D10" = "#9E5500")
)

plot_timecourse <- function(sample_name) {
  ggplot() +
    geom_point(
      data = filter(umap, sample != sample_name),
      aes(WNNUMAP_1, WNNUMAP_2),
      color = "grey95",
      size = 0.001
    ) +
    geom_point(
      data = filter(umap, sample == sample_name),
      aes(WNNUMAP_1, WNNUMAP_2, color = timepoint),
      size = 0.001
    ) +
    scale_color_manual(values = timepoint_cols[[sample_name]]) +
    coord_equal() +
    umap_theme
}

ggsave(
  file.path(fig7_dir, "Fig7a_timecourses.pdf"),
  plot_timecourse("MEF_CEBPA") +
    plot_timecourse("MEF_CEBPA_SPI1") +
    plot_timecourse("MEF_SPI1"),
  width = 12,
  height = 4
)

average_color <- function(cols) {
  x <- rowMeans(grDevices::col2rgb(cols))
  grDevices::rgb(x[1], x[2], x[3], maxColorValue = 255)
}

cluster_combo <- umap %>%
  count(WNN_cluster, sample, name = "n_cells") %>%
  complete(WNN_cluster, sample = sample_order, fill = list(n_cells = 0)) %>%
  group_by(WNN_cluster) %>%
  mutate(pct = 100 * n_cells / sum(n_cells), present = pct >= 5) %>%
  summarise(samples = list(sample[present]), .groups = "drop") %>%
  mutate(
    combo = vapply(samples, function(x) paste(x, collapse = "__"), character(1)),
    label = vapply(samples, function(x) paste(sample_labels[x], collapse = " + "), character(1)),
    color = vapply(samples, function(x) average_color(sample_cols[x]), character(1))
  )

combo_map <- distinct(cluster_combo, combo, label, color)

umap$sample_combo <- setNames(cluster_combo$combo, cluster_combo$WNN_cluster)[umap$WNN_cluster]

fig7a_sample <- ggplot(umap, aes(WNNUMAP_1, WNNUMAP_2, color = sample_combo)) +
  geom_point(size = 0.001, alpha = 0.5) +
  scale_color_manual(
    values = setNames(combo_map$color, combo_map$combo),
    labels = setNames(combo_map$label, combo_map$combo)
  ) +
  coord_equal() +
  umap_theme

ggsave(
  file.path(fig7_dir, "Fig7a_sample_composition.pdf"),
  fig7a_sample,
  width = 7,
  height = 6
)

##############Fig7b####################

umap$MEF_to_Adipocyte_score <- seu$MEF_to_Adipocyte_score
umap$MEF_to_Myeloid_score <- seu$MEF_to_Myeloid_score

plot_score <- function(score, title) {
  df <- arrange(umap, .data[[score]])

  ggplot(df, aes(WNNUMAP_1, WNNUMAP_2, color = .data[[score]])) +
    geom_point(size = 0.001, alpha = 0.9) +
    scale_color_gradient2(
      low = "blue3",
      mid = "grey85",
      high = "red3",
      midpoint = 0.5,
      limits = c(0, 1)
    ) +
    coord_equal() +
    labs(title = title, color = NULL) +
    umap_theme +
    centered_title
}

fig7b <- (
  plot_score("MEF_to_Adipocyte_score", "MEF to Adipocyte-like Trajectory score") +
    plot_gene("Ppargc1a") +
    plot_gene("Cd36")
) / (
  plot_score("MEF_to_Myeloid_score", "MEF to Myeloid Trajectory score") +
    plot_gene("Irf8") +
    plot_gene("Ciita")
)

ggsave(
  file.path(fig7_dir, "Fig7b_trajectory_scores_markers.pdf"),
  fig7b,
  width = 12,
  height = 8
)

##############Fig7c####################

check_inputs(c("metascape_adipo", "metascape_myeloid"))

read_common_go <- function(file, signature) {
  read_excel(file, sheet = "Enrichment") %>%
    mutate(
      Signature = signature,
      NegLog10P = -as.numeric(LogP),
      NegLog10Q = -as.numeric(`Log(q-value)`)
    ) %>%
    filter(str_detect(GroupID, "_Summary$")) %>%
    group_by(Description) %>%
    slice_max(NegLog10P, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(NegLog10P)) %>%
    slice_head(n = 10)
}

go_replacements <- c(
  "triglyceride biosynthetic process" = "triglyceride biosynthesis",
  "lipid biosynthetic process" = "lipid biosynthesis",
  "sulfur compound metabolic process" = "sulfur compound met. process",
  "deoxyguanosine catabolic process" = "deoxyguanosine cat. process",
  "brown fat cell differentiation" = "brown fat cell diff",
  "inflammatory response" = "inflam. response",
  "positive regulation of" = "pos reg of",
  "negative regulation of" = "neg reg of",
  "regulation of" = "reg of",
  "external stimulus" = "ext. stim.",
  "B cell receptor signaling pathway" = "B cell receptor sig.",
  "cytokine-mediated signaling pathway" = "cytokine-mediated sig.",
  "antigen processing and presentation of exogenous peptide antigen via MHC class II" =
    "antigen process & presentation of\nexo. peptide antigen via MHC class II"
)

fig7c_df <- bind_rows(
  read_common_go(paths$metascape_adipo, "Adipo signature"),
  read_common_go(paths$metascape_myeloid, "Myeloid signature")
) %>%
  mutate(Short_description = str_replace_all(Description, go_replacements))

x_max <- ceiling(max(fig7c_df$NegLog10P))
q_limits <- range(fig7c_df$NegLog10Q)

plot_common_go <- function(signature) {
  df <- fig7c_df %>%
    filter(Signature == signature) %>%
    arrange(desc(NegLog10P))

  df$Short_description <- factor(df$Short_description, levels = rev(df$Short_description))

  ggplot(df, aes(NegLog10P, Short_description, fill = NegLog10Q)) +
    geom_col(width = 0.7) +
    scale_fill_gradient(
      low = "blue3",
      high = "red3",
      limits = q_limits,
      name = expression(-log[10](q))
    ) +
    scale_x_continuous(limits = c(0, x_max)) +
    labs(x = expression(-log[10](P)), y = NULL, title = signature) +
    theme_classic() +
    theme(axis.ticks.y = element_blank()) +
    centered_title
}

fig7c <- plot_common_go("Adipo signature") / plot_common_go("Myeloid signature") +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(
  file.path(fig7_dir, "Fig7c_common_signature_GO.pdf"),
  fig7c,
  width = 6.6,
  height = 8.2
)

##############Fig7d####################

fig7d <- DimPlot(
  seu,
  reduction = "rna.adt.wnn.umap",
  group.by = "RNA_final",
  cols = identity_cols,
  label = FALSE,
  pt.size = 0.001,
  raster = FALSE
) +
  NoLegend() +
  umap_theme

ggsave(
  file.path(fig7_dir, "Fig7d_final_identities_UMAP.pdf"),
  fig7d,
  width = 8,
  height = 6
)

##############Fig7e####################

check_inputs("heatmap")
heatmap_data <- readRDS(paths$heatmap)

rna_z <- heatmap_data$RNA_heatmap_Z
adt_z <- heatmap_data$ADT_heatmap_Z
gene_order <- heatmap_data$Gene_order
cell_order <- heatmap_data$Cell_order

heatmap_identity <- as.character(cell_order$RNA_final)
heatmap_sample <- as.character(cell_order$sample)
identity_order <- unique(heatmap_identity)
heatmap_pop <- factor(heatmap_identity, levels = identity_order)

rna_col_fun <- colorRamp2(c(-3, 0, 2), c("#24D6E0", "black", "#E3D400"))
adt_col_fun <- colorRamp2(c(-3, 0, 2), c("#8FDDE6", "black", "#F3B167"))

row_split <- factor(
  case_when(
    gene_order$Heatmap_signature_old == "Common_adipocyte_signature" ~ "Adipo",
    gene_order$Heatmap_signature_old == "Common_myeloid_signature" ~ "Myeloid",
    TRUE ~ " "
  ),
  levels = c(" ", "Adipo", "Myeloid")
)

sample_count <- table(
  Identity = factor(heatmap_identity, levels = identity_order),
  Sample = factor(heatmap_sample, levels = sample_order)
)
sample_prop <- sample_count / rowSums(sample_count)
sample_prop[sample_prop < 0.05] <- 0
sample_prop <- sample_prop / rowSums(sample_prop)
sample_prop_cells <- sample_prop[heatmap_identity, sample_order]

top_anno <- HeatmapAnnotation(
  population = heatmap_pop,
  sample = anno_barplot(
    sample_prop_cells,
    gp = gpar(fill = sample_cols[sample_order], col = NA),
    bar_width = 1,
    border = FALSE,
    axis = FALSE,
    height = unit(3, "mm")
  ),
  col = list(population = identity_cols),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

gene_labels <- c(
  "Acta2", "Col11a1", "Top2a", "Mmp3", "Gli2", "Esf1", "Fap", "Robo1", "Pparg",
  "Cybb", "Cenpa", "Saa3", "Aldh1a2", "Cd74", "C1qa", "Penk", "Mafb"
)
check_features(gene_labels, rna_z, "heatmap label genes")
gene_rows <- match(gene_labels, rownames(rna_z))

gene_anno <- rowAnnotation(
  genes = anno_mark(
    at = gene_rows,
    labels = gene_labels,
    side = "right",
    labels_gp = gpar(fontsize = 7),
    link_gp = gpar(col = "grey35", lwd = 0.7)
  ),
  show_annotation_name = FALSE
)

rna_heatmap <- Heatmap(
  rna_z,
  name = "RNA row Z score",
  col = rna_col_fun,
  top_annotation = top_anno,
  right_annotation = gene_anno,
  column_split = heatmap_pop,
  row_split = row_split,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_rows = FALSE,
  row_order = seq_len(nrow(rna_z)),
  column_gap = unit(0.45, "mm"),
  row_gap = unit(1, "mm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_title_rot = 90,
  column_title_gp = gpar(fontsize = 7),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 8),
  use_raster = TRUE,
  raster_quality = 3,
  height = unit(5, "in")
)

adt_heatmap <- Heatmap(
  adt_z,
  name = "ADT expression",
  col = adt_col_fun,
  column_split = heatmap_pop,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_rows = FALSE,
  row_order = seq_len(nrow(adt_z)),
  column_gap = unit(0.45, "mm"),
  show_row_names = TRUE,
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 7),
  show_column_names = FALSE,
  use_raster = TRUE,
  raster_quality = 3,
  height = unit(2, "in")
)

pdf(
  file.path(fig7_dir, "Fig7e_RNA_ADT_heatmap.pdf"),
  width = 20,
  height = 9,
  useDingbats = FALSE
)
draw(
  rna_heatmap %v% adt_heatmap,
  heatmap_legend_side = "left",
  annotation_legend_side = "left",
  ht_gap = unit(1, "mm"),
  padding = unit(c(12, 8, 4, 4), "mm")
)
invisible(dev.off())

##############Fig7f####################

check_inputs("metascape_fig7f")

fig7f_df <- read_excel(paths$metascape_fig7f, sheet = "Plot_data")
fig7f_terms <- read_excel(paths$metascape_fig7f, sheet = "Term_plot_order")

fig7f_order <- c(
  "FAP-like",
  "Fap+ activated adipogenic progenitor",
  "Early adipocyte",
  "Pparg+ Adipo-mac hybrid",
  "Dual-prog transition",
  "APC Mac-like",
  "Inflammatory Mac-like"
)

warn_unmatched(fig7f_df$New_identity, fig7f_order, "Fig 7f population")
warn_unmatched(fig7f_df$Plot_term, fig7f_terms$Description, "Fig 7f GO term")

fig7f_df$New_identity <- factor(fig7f_df$New_identity, levels = fig7f_order)
fig7f_df$Plot_term <- factor(fig7f_df$Plot_term, levels = rev(fig7f_terms$Description))

ggsave(
  file.path(fig7_dir, "Fig7f_selected_population_GO.pdf"),
  plot_go_dots(fig7f_df),
  width = 10,
  height = 10
)

##############Fig7g####################

ggsave(
  file.path(fig7_dir, "Fig7g_Il12b_UMAP.pdf"),
  plot_gene("Il12b"),
  width = 6,
  height = 6
)

message("Fig 7 panels written to ", fig7_dir)
