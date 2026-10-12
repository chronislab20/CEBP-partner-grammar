check_inputs("seurat")
seu <- readRDS(paths$seurat)

rocket_cols <- viridisLite::rocket(n = 9, direction = -1)

umap <- as.data.frame(Embeddings(seu, "rna.adt.wnn.umap")[, 1:2])
colnames(umap) <- c("WNNUMAP_1", "WNNUMAP_2")
umap$sample <- as.character(seu$sample)
umap$timepoint <- as.character(seu$timepoint)
umap$WNN_cluster <- as.character(seu$wnn_clusters)

rna_data <- LayerData(seu, assay = "cellbender_rna", layer = "data")

plot_gene <- function(gene) {
  check_features(gene, rna_data, "genes")
  df <- umap
  df$expr <- as.numeric(rna_data[gene, rownames(df)])
  df <- arrange(df, expr)

  ggplot(df, aes(WNNUMAP_1, WNNUMAP_2, color = expr)) +
    geom_point(size = 0.001) +
    scale_color_gradientn(colours = rocket_cols) +
    coord_equal() +
    labs(title = gene, color = NULL) +
    umap_theme +
    centered_title
}

plot_go_dots <- function(df) {
  ggplot(df, aes(New_identity, Plot_term, size = NegLog10Q, color = NegLog10P)) +
    geom_point() +
    scale_color_gradient(low = "blue3", high = "red3") +
    scale_size_continuous(range = c(2.5, 8)) +
    labs(
      x = NULL,
      y = NULL,
      size = expression(-log[10](q)),
      color = expression(-log[10](P))
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
      axis.text.y = element_text(colour = "black"),
      axis.ticks = element_blank()
    )
}
