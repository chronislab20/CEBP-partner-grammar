check_inputs("atac")
atac <- readRDS(paths$atac)

atac_order <- c(
  "Adhesion-remodel mesenchyme",
  "Sca1+Cd34+CD140+ Fibro-adipogenic Progenitor",
  "Adipogenic transition (early)",
  "Skeletal/mesenchyme",
  "Residual Fibroblast",
  "Fap+ activated adipogenic progenitor",
  "Early adipocyte",
  "Pparg+ Adipo-mac hybrid",
  "Inflammatory Mac-like",
  "CMP-like progenitor",
  "APC Mac-like",
  "Low-signature-SPI1"
)

atac$ATAC_final <- factor(as.character(atac$ATAC_final), levels = atac_order)

compute_chromhmm_rpkm <- function() {
  check_inputs("chromhmm_bed")

  chrom <- readr::read_tsv(
    paths$chromhmm_bed,
    col_names = c("chr", "start", "end", "state"),
    col_types = "ciic",
    show_col_types = FALSE
  )

  chrom_gr <- makeGRangesFromDataFrame(
    chrom,
    seqnames.field = "chr",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE,
    starts.in.df.are.0based = TRUE
  )

  region_ids <- paste0(seqnames(chrom_gr), "-", start(chrom_gr), "-", end(chrom_gr))
  names(chrom_gr) <- region_ids

  DefaultAssay(atac) <- "peaks"

  feature_matrix <- FeatureMatrix(
    fragments = Fragments(atac),
    features = chrom_gr,
    cells = colnames(atac),
    process_n = 2000
  )

  feature_match <- match(rownames(feature_matrix), region_ids)
  region_state <- as.character(mcols(chrom_gr)$state[feature_match])
  region_length_kb <- width(chrom_gr)[feature_match] / 1000
  total_fragments_million <- Matrix::colSums(feature_matrix) / 1e6

  rpkm <- do.call(cbind, lapply(paste0("U", 1:20), function(state_now) {
    ii <- which(region_state == state_now)
    log2(1 + Matrix::colSums(feature_matrix[ii, , drop = FALSE]) /
           (total_fragments_million * sum(region_length_kb[ii])))
  }))

  rownames(rpkm) <- colnames(feature_matrix)
  colnames(rpkm) <- paste0("St", 1:20)

  dir.create(dirname(paths$chromhmm_rpkm), recursive = TRUE, showWarnings = FALSE)
  saveRDS(rpkm, paths$chromhmm_rpkm)
  rpkm
}

add_chromhmm_scores <- function(obj) {
  rpkm <- if (file.exists(paths$chromhmm_rpkm)) readRDS(paths$chromhmm_rpkm) else compute_chromhmm_rpkm()
  for (state in paste0("St", 1:20)) {
    obj[[paste0("ChromHMM_", state)]] <- rpkm[colnames(obj), state]
  }
  obj
}

state_cols <- rev(RColorBrewer::brewer.pal(11, "PuOr"))

plot_state <- function(state) {
  FeaturePlot(
    atac,
    features = paste0("ChromHMM_", state),
    reduction = "umap.peaks",
    order = TRUE,
    pt.size = 0.01,
    raster = FALSE
  ) +
    scale_color_gradientn(colors = state_cols) +
    ggtitle(state) +
    centered_title
}

save_state_pdf <- function(states, file) {
  pdf(file, width = 6, height = 6, useDingbats = FALSE)
  for (state in states) print(plot_state(state))
  invisible(dev.off())
}
