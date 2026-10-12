# Run ChromBPNet from the Greenleaf lab
# Liu, B. B., Jessa, S., Kim, S. H., Ng, Y. T., Higashino, S. I., Marinov, G. K., ... & Greenleaf, W. J. (2025). Dissecting regulatory syntax in human development with scalable multiomics and deep learning. bioRxiv, 2025-04.
# https://github.com/GreenleafLab/HDMA/tree/main
# chrombpnet: https://github.com/GreenleafLab/HDMA/blob/main/code/03-chrombpnet/


# Load libraries
library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(BiocParallel)
library(rjson)
library(rvest)
library(reshape2)
library(ComplexHeatmap)

register(MulticoreParam(100))
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

# Set your project paths
# project_path <- "/path/to/bias_models"
# report_path <- "/path/to/bias_models/ChromBPNet_report/"

if (!dir.exists(report_path)) {
  dir.create(report_path, recursive = TRUE)
}


# Use sample names as samples_order
# samples_order <- list.files(project_path)
samples_order <- list.files(project_path, full.names = TRUE)
samples_order <- basename(samples_order[file.info(samples_order)$isdir])
folds <- c("0", "1", "2", "3", "4")

qc_all <- data.frame()
for (prog_sel in samples_order) {
  for (fold_sel in folds) {
    json_path <- file.path(project_path, prog_sel, paste0("fold_", fold_sel), "evaluation", "bias_metrics.json")
    if (!file.exists(json_path)) {
      print(paste("Missing:", json_path))
    } else {
      # qc_json <- rjson::fromJSON(file = json_path)
      if (file.exists(json_path) && file.info(json_path)$size > 0) {
        qc_json <- jsonlite::fromJSON(json_path)
      } else {
        print(paste("Skipping: empty or invalid file:", json_path))
        next
      }


      qc_all <- rbind(qc_all, data.frame(
        Program = prog_sel,
        Fold = fold_sel,
        counts_peaks_pearsonr = qc_json$counts_metrics$peaks$pearsonr,
        counts_non_peaks_pearsonr=qc_json$counts_metrics$nonpeaks$pearsonr,
        profile_peaks_median_jsd = qc_json$profile_metrics$peaks$median_jsd,
        profile_non_peaks_median_jsd = qc_json$profile_metrics$nonpeaks$median_jsd

      ))
    }
  }
}

pdf(file.path(report_path, "evaluation_Bias_metrics.pdf"), width = 15, height = 7, useDingbats = FALSE)
qc_all$counts_peaks_pearsonr_flag <- ifelse(qc_all$counts_peaks_pearsonr < -0.3, "Flagged (less than -0.3)", "Not flagged")
ggplot(qc_all, aes(Fold, counts_peaks_pearsonr, fill = Fold, color = counts_peaks_pearsonr_flag)) +
  geom_bar(stat = "identity", size = 1.2) +
  facet_wrap(~Program) +
  scale_color_manual(
    values = c("Flagged (less than -0.3)" = "black", "Not flagged" = "white"),
    name = "Bar outline"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal()


qc_all$counts_non_peaks_pearsonr_flag <- ifelse(qc_all$counts_non_peaks_pearsonr < 0, "Flagged (less than 0)", "Not flagged")
ggplot(qc_all, aes(Fold, counts_non_peaks_pearsonr, fill = Fold, color = counts_non_peaks_pearsonr_flag)) +
  geom_bar(stat = "identity", size = 1.2) +
  facet_wrap(~Program) +
  scale_color_manual(
    values = c("Flagged (less than 0)" = "black", "Not flagged" = "white"),
    name = "Bar outline"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal()

ggplot(qc_all, aes(Fold, profile_peaks_median_jsd, fill = Fold)) + geom_bar(stat = "identity") + facet_wrap(~Program) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
ggplot(qc_all, aes(Fold, profile_non_peaks_median_jsd, fill = Fold)) + geom_bar(stat = "identity") + facet_wrap(~Program) + scale_y_continuous(expand = expansion(mult = c(0, 0.1)))


dev.off()






# -------- settings --------
need_cols      <- c("pattern", "match0", "qval0")   # columns to show in the table
top_n_perfold  <- 10                                 # how many rows per fold to display
wrap_width     <- 40                                 # wrap long text (pattern) to this many chars
page_size      <- c(width = 11, height = 8.5)        # landscape page helps fit more text
pdf_file       <- file.path(report_path, "Motif_report_folds_per_program.pdf")

if (!dir.exists(report_path)) dir.create(report_path, recursive = TRUE, showWarnings = FALSE)

# -------- open PDF (one page per Program) --------
pdf(pdf_file, width = page_size["width"], height = page_size["height"], onefile = TRUE)



for (prog_sel in samples_order) {
  # gather per-fold tables for this Program
  per_fold <- list()

  for (fold_sel in folds) {
    html_path <- file.path(
      project_path, prog_sel, paste0("fold_", fold_sel),
      "evaluation", "overall_report.html"
    )
    if (!file.exists(html_path)) {
      message("Missing: ", html_path)
      next
    }

    html_tables <- read_html(html_path) %>% html_table(fill = TRUE)
    if (length(html_tables) < 3) {
      message("No 3rd table in: ", html_path)
      next
    }

    df <- html_tables[[3]]
    # check required columns exist
    if (!all(need_cols %in% names(df))) {
      message("Skipping (missing cols) in ", html_path, ": ",
              paste(setdiff(need_cols, names(df)), collapse = ", "))
      next
    }

    # keep only desired columns + Fold for splitting
    subdf <- df %>%
      select(all_of(need_cols)) %>%
      mutate(Fold = fold_sel, .before = 1L) %>%   # put Fold as first column (not printed)
      head(top_n_perfold)

    per_fold[[as.character(fold_sel)]] <- subdf
  }

  if (!length(per_fold)) {
    message("No data for Program ", prog_sel, " — skipping page.")
    next
  }

  # combine all folds (Fold column will be used for row_split, not printed)
  prn <- bind_rows(per_fold)

  # format numbers, wrap long text
  prn <- prn %>%
    mutate(
      across(where(is.numeric), ~ signif(., 3))
    )

  # wrap long strings (e.g., pattern)
  wrap_cols <- intersect(c("pattern"), names(prn))
  for (cc in wrap_cols) {
    prn[[cc]] <- vapply(
      prn[[cc]],
      function(z) paste(strwrap(as.character(z), width = wrap_width), collapse = "\n"),
      character(1)
    )
  }

  # columns to actually show in the grid (Fold used only for splitting)
  show_cols <- need_cols

  # build the character matrix for printing
  mat_chr <- as.matrix(prn[, show_cols, drop = FALSE] %>% mutate(across(everything(), as.character)))
  rownames(mat_chr) <- paste0("Row ", sequence(rle(prn$Fold)$lengths))  # 1..10 within each fold

  # dummy numeric matrix just to define the grid
  m <- matrix(0, nrow = nrow(mat_chr), ncol = ncol(mat_chr))
  colnames(m) <- colnames(mat_chr)
  rownames(m) <- rownames(mat_chr)

  # font size heuristic based on total rows
  nR <- nrow(mat_chr)
  fs <- max(6, 11 - 0.07 * pmax(0, nR - 20))

  # page title
  page_title <- paste0("Program: ", prog_sel, "   |   Top ", top_n_perfold, " per Fold")

  # one heatmap with row splits = folds -> all folds on ONE page
  ht <- Heatmap(
    m,
    name = NULL,
    col = "white",
    show_heatmap_legend = FALSE,
    cluster_rows = FALSE, cluster_columns = FALSE,
    rect_gp = gpar(col = "grey85", fill = NA),
    row_names_side = "left",
    column_names_side = "top",
    column_title = page_title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_split = factor(prn$Fold, levels = folds),
    row_title_gp = gpar(fontsize = 11, fontface = "bold"),
    gap = unit(2, "mm"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.text(mat_chr[i, j], x, y, gp = gpar(fontsize = fs))
    }
  )

  # draw() starts a fresh page for each call when a PDF device is open
  draw(ht)
}

dev.off()
message("Saved PDF: ", pdf_file)






















# # where to save the multi-page PDF
# if (!dir.exists(report_path)) dir.create(report_path, recursive = TRUE, showWarnings = FALSE)
# pdf_file <- file.path(report_path, "overall_report_per_fold.pdf")

# # open a single PDF device (each draw() will be a new page)
# pdf(pdf_file, width = 8.5, height = 11, onefile = TRUE)
# on.exit(dev.off(), add = TRUE)   # make sure the device closes even if errors occur

# html_overall_report_all <- list()
# page_opened <- FALSE
# idx <- 1

# need <- c("Program", "Fold", "pattern", "match0", "qval0")

# for (prog_sel in samples_order) {
#   for (fold_sel in folds) {
#     html_path <- file.path(project_path, prog_sel, paste0("fold_", fold_sel),
#                            "evaluation", "overall_report.html")

#     if (!file.exists(html_path)) {
#       message("Missing: ", html_path)
#       next
#     }

#     html_tables <- read_html(html_path) %>% html_table(fill = TRUE)
#     if (length(html_tables) < 3) {
#       message("No 3rd table in: ", html_path)
#       next
#     }

#     df <- html_tables[[3]]
#     df=head(df, 10)
#     df$Program <- prog_sel
#     df$Fold    <- fold_sel

#     # ensure required columns exist
#     if (!all(need %in% names(df))) {
#       message("Skipping (missing cols) in ", html_path, ": ",
#               paste(setdiff(need, names(df)), collapse = ", "))
#       next
#     }

#     subdf <- df[need]

#     # ---- accumulate for later if you want the combined table too ----
#     html_overall_report_all[[idx]] <- subdf
#     idx <- idx + 1

#     # ---- format values for printing ----
#     prn <- subdf %>%
#       mutate(
#         across(where(is.numeric), ~ signif(., 3))
#       )

#     # wrap long 'pattern' (and others) so they fit in a cell
#     wrap_cols <- c("pattern")  # add more names if needed
#     for (cc in intersect(wrap_cols, names(prn))) {
#       prn[[cc]] <- vapply(
#         prn[[cc]],
#         function(z) paste(strwrap(as.character(z), width = 40), collapse = "\n"),
#         character(1)
#       )
#     }

#     # matrix of strings for printing
#     mat_chr <- as.matrix(prn %>% mutate(across(everything(), as.character)))
#     rownames(mat_chr) <- paste0("Row ", seq_len(nrow(mat_chr)))  # optional row labels

#     # dummy numeric matrix just to define the grid
#     m <- matrix(0, nrow = nrow(mat_chr), ncol = ncol(mat_chr))
#     colnames(m) <- colnames(mat_chr)
#     rownames(m) <- rownames(mat_chr)

#     # pick a readable font size depending on number of rows
#     nR <- nrow(mat_chr)
#     fs <- max(6, 10 - 0.08 * pmax(0, nR - 10))  # tweak as needed

#     # title for the page
#     page_title <- paste0("Program: ", prog_sel, "   |   Fold: ", fold_sel)

#     # build heatmap-like table (white cells, borders, text in each cell)
#     ht <- Heatmap(
#       m,
#       name = NULL,
#       col = "white",
#       show_heatmap_legend = FALSE,
#       cluster_rows = FALSE, cluster_columns = FALSE,
#       rect_gp = gpar(col = "grey85", fill = NA),
#       row_names_side = "left",
#       column_names_side = "top",
#       column_title = page_title,
#       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
#       cell_fun = function(j, i, x, y, w, h, fill) {
#         grid.text(mat_chr[i, j], x, y, gp = gpar(fontsize = fs), default.units = "npc")
#       }
#     )

#     # draw one page
#     if (page_opened) grid.newpage()
#     draw(ht)
#     page_opened <- TRUE
#   }
# }
# dev.off()
# # bind all rows into one data.frame (optional, if you still want it later)
# html_overall_report_all <- dplyr::bind_rows(html_overall_report_all)

# message("Saved multi-page PDF to: ", pdf_file)
















# html_overall_report_all <- data.frame()
# for (prog_sel in samples_order) {
#   for (fold_sel in folds) {
#     html_path <- file.path(project_path, prog_sel, paste0("fold_", fold_sel), "evaluation", "overall_report.html")
#     if (!file.exists(html_path)) {
#       print(paste("Missing:", html_path))
#     } else {
#       html_tables <- read_html(html_path) %>% html_table(fill = TRUE)
#       if (length(html_tables) >= 3) {
#         df <- html_tables[[3]]
#         df$Program <- prog_sel
#         df$Fold <- fold_sel
#         subdf <- df[c("Program", "Fold",'pattern', "match0", "qval0")]

#         html_overall_report_all <- rbind(html_overall_report_all, subdf)
#       }
#     }
#   }
# }





# ## 1) Select and lightly format
# html_overall_report_all <- html_overall_report_all %>%
#   select(Program, Fold, match0, qval0) %>%
#   mutate(
#     across(where(is.numeric), ~ signif(., 3)),   # adjust digits as you like
#     across(everything(), as.character)
#   ) %>%
#   mutate(Row = row_number())                     # keep a row index for plotting

# ## 2) Long format for ggplot
# long_tbl <- melt(html_overall_report_all, id.vars = "Row",
#                  variable.name = "Column", value.name = "Value")

# ## Optional: wrap long strings so they don’t overflow cells
# wrap_width <- 30  # characters per line; tune for your page width
# long_tbl$Value <- vapply(
#   long_tbl$Value,
#   function(z) paste(strwrap(z, width = wrap_width), collapse = "\n"),
#   character(1)
# )

# ## 3) Build the plot
# p <- ggplot(long_tbl, aes(x = Column, y = Row)) +
#   geom_tile(fill = NA, color = "grey80") +
#   geom_text(aes(label = Value), size = 3, lineheight = 0.9) +
#   scale_y_reverse(expand = c(0, 0)) +                 # top row at the top
#   scale_x_discrete(position = "top") +
#   coord_cartesian(clip = "off") +
#   theme_minimal(base_size = 11) +
#   theme(
#     panel.grid = element_blank(),
#     axis.title = element_blank(),
#     axis.text.y = element_blank(),                    # hide numeric row labels
#     axis.ticks = element_blank(),
#     plot.margin = margin(10, 10, 10, 10)
#   )

# ## 4) Save to PDF (height adapts to number of rows)
# pdf_height <- max(2.5, min(11, 1 + 0.30 * nrow(html_overall_report_all)))  # inches
# ggsave(file.path(report_path, "profile_motifs.pdf"), p, width = 8.5, height = pdf_height, units = "in")
