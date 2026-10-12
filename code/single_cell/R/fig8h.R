source(here::here("R", "setup.R"))
library(rvest); library(ggseqlogo); library(universalmotif); library(scales)

##############Fig8h####################

check_inputs(c("fig8h_values", "modisco_html", "modisco_cwm", "motifs_meme"))

pops <- c("Residual Fibroblast", "Adhesion-remodel mesenchyme", "Skeletal/mesenchyme", "Adipogenic transition (early)", "Sca1+Cd34+CD140+ Fibro-adipogenic Progenitor", "Fap+ activated adipogenic progenitor", "Early adipocyte", "Pparg+ Adipo-mac hybrid", "Inflammatory Mac-like", "APC Mac-like", "CMP-like progenitor", "Low-signature-SPI1")
pop_short <- setNames(c("Res Fibro", "Adh Mes", "Skeletal mes", "Adip trans", "Sca+ FAP", "Fap+ Adip", "Early Adip", "Adipo-Mac", "Infl Mac", "APC Mac", "CMP prog", "Low-SPI"), pops)
states <- paste0("E", 5:17)
access_cols <- c("Fibro access" = "#C3D0C0", "Adipo access" = "#D9A996", "Myelo access" = "#9AACCD")
state_access <- setNames(rep(names(access_cols), c(4, 5, 4)), states)

fig8 <- read.csv(paths$fig8h_values, check.names = FALSE, stringsAsFactors = FALSE)
get_mat <- function(prefix, cols) { m <- as.matrix(fig8[, paste0(prefix, cols), drop = FALSE]); storage.mode(m) <- "double"; dimnames(m) <- list(fig8$pattern, cols); m }
peak_pct <- get_mat("Pct_peaks_containing_motif_", pops); pop_z <- get_mat("Population_importance_rowZ_", pops); state_z <- get_mat("ChromHMM_importance_rowZ_", states)
peak_pct_z <- t(scale(t(peak_pct))); peak_pct_z[!is.finite(peak_pct_z)] <- 0
n_inst <- as.numeric(fig8$Total_FiNeMo_instances)

motif_tbl <- read_html(paths$modisco_html) %>% html_table(fill = TRUE) %>% .[[1]] %>% as.data.frame() %>%
  dplyr::filter(grepl("^pos_patterns\\.pattern_[0-9]+$", pattern)) %>% dplyr::select(pattern, match0) %>% dplyr::distinct(pattern, .keep_all = TRUE)
fig8 <- dplyr::left_join(fig8, motif_tbl, by = "pattern")

fix_pwm <- function(x) {
  if (is.null(x) || length(x) == 0) return(NULL); x <- as.matrix(x)
  if (nrow(x) != 4 && ncol(x) == 4) x <- t(x); if (nrow(x) != 4 || ncol(x) == 0) return(NULL)
  storage.mode(x) <- "double"; x[!is.finite(x)] <- 0; dimnames(x) <- list(c("A", "C", "G", "T"), seq_len(ncol(x))); x
}
norm_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", as.character(x)))
meme <- read_meme(paths$motifs_meme); meme_names <- vapply(meme, function(x) x@name, character(1))
find_known <- function(tf) {
  if (is.na(tf) || tf %in% c("", "NaN")) return(NULL)
  ii <- match(tf, meme_names); if (is.na(ii)) ii <- match(norm_name(tf), norm_name(meme_names)); if (is.na(ii)) ii <- match(norm_name(sub("_.*$", "", tf)), norm_name(meme_names))
  if (is.na(ii)) NULL else meme[[ii]]
}
denovo_mat <- function(p) {
  ff <- file.path(paths$modisco_cwm, p, "modisco_fc.txt"); if (!file.exists(ff)) return(NULL)
  m <- fix_pwm(read.table(ff, header = FALSE)); if (is.null(m)) return(NULL); m[m < 0] <- 0; if (max(m) <= 0) NULL else m
}
known_mat <- function(tf, ori) { m <- find_known(tf); if (is.null(m)) return(NULL); if (!is.na(ori) && ori == "-") m <- motif_rc(m); fix_pwm(m@motif) }
logo_grob <- function(m, ...) if (is.null(m)) NULL else ggplotGrob(ggseqlogo(m, ...) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0)))
denovo_grobs <- lapply(lapply(fig8$pattern, denovo_mat), logo_grob, method = "custom")
known_grobs <- lapply(Map(known_mat, fig8$match0, fig8$Known_orientation), logo_grob, seq_type = "dna")

heat_h <- unit(nrow(fig8) * 6.3, "mm"); title_gp <- gpar(fontsize = 7, fontface = "plain"); cell_border <- gpar(col = "#AAAAAA", lwd = .25)
leg <- function(title, w) list(title = title, direction = "horizontal", legend_width = unit(w, "mm"), title_position = "topcenter")
three_breaks <- function(x) { b <- as.numeric(quantile(x, c(0, .5, 1), na.rm = TRUE)); if (length(unique(b)) < 3) c(min(x), mean(range(x)), max(x)) else b }
text_col <- function(fill) { x <- as.numeric(col2rgb(fill)) / 255; if (sum(c(.2126, .7152, .0722) * x) < .52) "white" else "black" }

side_panel <- function(n, title, w, cell_fun) Heatmap(matrix(0, n, 1), col = colorRamp2(c(0, 1), c("white", "white")), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE, column_title = title, column_title_gp = title_gp, width = unit(w, "mm"), height = heat_h, border = FALSE, rect_gp = gpar(col = NA), cell_fun = cell_fun)
text_panel <- function(lab, title, w, left = FALSE) side_panel(length(lab), title, w, function(j, i, x, y, width, height, fill) grid.text(lab[i], x = if (left) x - width / 2 + unit(1, "mm") else x, y = y, just = if (left) "left" else "centre", gp = gpar(fontsize = 5.8)))
logo_panel <- function(g, title, w) side_panel(length(g), title, w, function(j, i, x, y, width, height, fill) if (!is.null(g[[i]])) { pushViewport(viewport(x = x, y = y, width = width * .97, height = height * .9, clip = "on")); grid.draw(g[[i]]); popViewport() })

count_log <- matrix(log10(n_inst + 1), ncol = 1)
h_count <- Heatmap(count_log, name = "Total_instances", col = colorRamp2(three_breaks(count_log), c("#F7FBFF", "#6BAED6", "#08306B")), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE, show_heatmap_legend = FALSE,
  column_title = "Total FiNeMo\ninstances", column_title_gp = title_gp, width = unit(17, "mm"), height = heat_h, border = TRUE, rect_gp = gpar(col = "#AAAAAA", lwd = .3),
  cell_fun = function(j, i, x, y, width, height, fill) grid.text(comma(n_inst[i], accuracy = 1), x, y, gp = gpar(fontsize = 4.8, col = text_col(fill))))

lim <- max(abs(peak_pct_z), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1
h_peak <- Heatmap(peak_pct_z, name = "Peak_prevalence", col = colorRamp2(c(0, lim), c("grey95", "firebrick3")), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  column_labels = unname(pop_short[pops]), column_names_rot = 90, column_names_gp = gpar(fontsize = 5.2), column_title = "% peaks containing motif", column_title_gp = title_gp,
  width = unit(length(pops) * 5, "mm"), height = heat_h, border = TRUE, rect_gp = cell_border, heatmap_legend_param = leg("Within-motif row Z", 30),
  cell_fun = function(j, i, x, y, width, height, fill) grid.text(sprintf("%.0f", peak_pct[i, j]), x, y, gp = gpar(fontsize = 4, col = "white")))

z_heat <- function(m, name, title, high, anno, labs = colnames(m)) Heatmap(m, name = name, col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#FFFFFF", high)), cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  column_labels = labs, column_names_rot = 90, column_names_gp = gpar(fontsize = 5.2), column_title = title, column_title_gp = title_gp, border = TRUE, rect_gp = cell_border,
  width = unit(ncol(m) * 5, "mm"), height = heat_h, bottom_annotation = anno, heatmap_legend_param = leg("Row z-score", 32))
h_pop <- z_heat(pop_z, "Population_zscore", "ATAC populations\nImportance", "#B2182B",
  HeatmapAnnotation(Population = factor(pops, levels = pops), col = list(Population = identity_cols[pops]), show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(4, "mm")), unname(pop_short[pops]))
h_state <- z_heat(state_z, "ChromHMM_zscore", "ChromHMM E5-E17\nImportance", "#008B8B",
  HeatmapAnnotation(Access = factor(state_access, levels = names(access_cols)), col = list(Access = access_cols), show_annotation_name = FALSE, show_legend = FALSE, simple_anno_size = unit(4, "mm")))

pdf(file.path(fig8_dir, "Fig8h_ChromBPNet_motifs.pdf"), width = 20, height = max(14, nrow(fig8) * 6.3 / 25.4 + 4.8), useDingbats = FALSE)
draw(text_panel(fig8$Name, "Name", 54, TRUE) + text_panel(as.character(fig8$Pattern), "Pattern", 11) + logo_panel(denovo_grobs, "De novo motif", 38) + logo_panel(known_grobs, "Most similar known motif", 38) + h_count + h_peak + h_pop + h_state,
  heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legends = FALSE, legend_grouping = "original", gap = unit(1.2, "mm"), padding = unit(c(18, 10, 22, 10), "mm"))
invisible(dev.off())

message("Fig 8h written to ", fig8_dir)
