###############################################################################
# ChromBPNet / FiNeMo post processing script
# Motif enrichment, CEBPB 7 group overlaps, clustering, and co occurrence
#
# Source methods from:
# Liu et al, 2025, "Dissecting regulatory syntax in human development with
# scalable multiomics and deep learning"
#
# Greenleaf Lab HDMA code:
# https://github.com/GreenleafLab/HDMA/tree/main/code/03-chrombpnet/
###############################################################################




suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(data.table)
  library(GenomicRanges)
  library(BiocParallel)
  library(ComplexHeatmap)
  library(patchwork)
  library(reticulate)
  library(rtracklayer)
  library(tidyr)
})

###############################################################################
# Parallel and thread settings
###############################################################################

# Use 10 cores for BiocParallel
register(MulticoreParam(10))

# Force single threaded BLAS and OpenMP to avoid over subscription
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)

Sys.setenv(
  OMP_NUM_THREADS          = "1",
  OPENBLAS_NUM_THREADS     = "1",
  MKL_NUM_THREADS          = "1",
  VECLIB_MAXIMUM_THREADS   = "1",
  NUMEXPR_NUM_THREADS      = "1"
)

###############################################################################
# Paths, global options, and program definitions
###############################################################################

pattern_group      <- "pos_patterns"    # "pos_patterns" or "neg_patterns"
counts_or_profile  <- "counts"          # "counts" or "profile"

finemo_merged_dir <- Sys.getenv("FINEMO_MERGED_DIR")
project_path <- Sys.getenv("FINEMO_MERGED_DIR")

# Programs to keep and their colors
program_order <- c(
  "MEF",
  # "MEF_d2_OSKM",
  # "MEF_d2_OSKM_mCEBPA",
  # "MEF_d2_OSKM_mCEBPB",
  # "MEF_d2_mCEBPA",
  # 'MEF_d2_mCEBPB' ,
  'NSC' ,
# 'NSC_d2_mCEBPA' ,
# 'NSC_d2_mCEBPB' ,
# 'NSC_d2_OSKM' ,
# 'NSC_d2_OSKM_mCEBPA' ,
# 'NSC_d2_OSKM_mCEBPB' ,
  'ESC' ,
# 'ESC_d1_hCEBPA' ,
# 'ESC_d2_hCEBPA'  ,
'Bcell_Graf' 
# "Bcell_Graf_Cebpa18h",
# 'Bcell_Graf_OKSMd1' ,
# 'Bcell_Graf_OKSMd2'
)



program_colors <- c(
  "MEF"                = "#1f8fdb",
  # "MEF_d2_OSKM"        = "#903498",
  # "MEF_d2_OSKM_mCEBPA" = "#ffb703",
  # "MEF_d2_OSKM_mCEBPB" = "#8ecae6",
  #   "MEF_d2_mCEBPA"       = "#e63946",
  # "MEF_d2_mCEBPB"       = "#f4a261",
  # "ESC"                = "#09989b",
  "NSC"            = "#118ab2",
  # "NSC_d2_mCEBPA"       = "#ef476f",
  # "NSC_d2_mCEBPB"       = "#ffd166",
  # "NSC_d2_OSKM"         = "#06d6a0",
  # "NSC_d2_OSKM_mCEBPA"  = "#073b4c",
  # "NSC_d2_OSKM_mCEBPB"  = "#8338ec",
  # "Bcell_Graf_Cebpa18h"="#3a86ff",
  "Bcell_Graf"="#ffbe0b",
  # "Bcell_Graf_OKSMd1"="#9b5de5",
  # "Bcell_Graf_OKSMd2"="#00f5d4",
    "ESC"            = "#09989b"
  # "ESC_d1_hCEBPA"       = "#ff006e",
  # "ESC_d2_hCEBPA"       = "#fb5607"
)



# Output directories for figures and tables
figure_dir <- file.path(project_path, pattern_group, "figures_Starting_cell_types_clusters")
tables_dir <- file.path(project_path, pattern_group, "figures_Starting_cell_types_clusters")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
# 1. Parse MoDISco motif table (motifs.html) to get all_matches labels
###############################################################################

modisco_html <- rvest::read_html(
  file.path(finemo_merged_dir, "pos_patterns/merged_modisco.html/motifs.html")
)

modisco_html_tables <- modisco_html %>%
  rvest::html_table(fill = TRUE)

# Use the first table which contains the motif summary
modisco_table <- modisco_html_tables[[1]]

# Tag patterns with q value > 0.05
modisco_table$pval_tag <- ""
modisco_table$pval_tag[modisco_table$qval0 > 0.05] <- ";*"

# Build a compact motif label that summarizes the top three TOMTOM matches
modisco_table$all_matches <- paste(
  gsub("_.*", "", modisco_table$match0),
  gsub("_.*", "", modisco_table$match1),
  gsub("_.*", "", modisco_table$match2),
  gsub(".*_", "", modisco_table$pattern),
  sep = ";"
)

# Append significance tag
modisco_table$all_matches <- paste0(
  modisco_table$all_matches,
  modisco_table$pval_tag
)

###############################################################################
# 2. Load FiNeMo hits for merged run and count hits per program and motif
###############################################################################

# In this section we work only with the "merged" FiNeMo run
selected_finemodel <- "merged"

finemo_hits <- fread(
  file.path(finemo_merged_dir, "pos_patterns/finemo", selected_finemodel, "hits.tsv")
)

# Clean program label (Pr) from peak_name
# Example peak_name: "MEF_mm10_chr1_..." or "MEF_ATAC_mm10_..."
finemo_hits$Pr <- finemo_hits$peak_name
finemo_hits$Pr <- sub("(_mm10.*)$", "", finemo_hits$Pr)
finemo_hits$Pr <- sub("[-_]?ATAC.*$", "", finemo_hits$Pr, ignore.case = TRUE)
finemo_hits$Pr <- sub("[-_]?merged*$", "", finemo_hits$Pr, ignore.case = TRUE)
finemo_hits$Pr <- sub("[-_]?flag$", "", finemo_hits$Pr, ignore.case = TRUE)

# unique(finemo_hits$Pr)
# length(unique(finemo_hits$Pr))

# Restrict FiNeMo hits to the selected programs
finemo_hits <- finemo_hits %>%
  filter(Pr %in% program_order)

unique(finemo_hits$Pr)
length(unique(finemo_hits$Pr))

# Count hits per program and motif
finemo_hits_counts <- as.data.frame(
  table(finemo_hits[, c("Pr", "motif_name")])
)

colnames(finemo_hits_counts) <- c("sample", "pattern", "Freq")

# Attach MoDISco labels
finemo_hits_counts <- merge(
  finemo_hits_counts,
  modisco_table[, c("pattern", "all_matches")],
  by.x  = "pattern",
  by.y  = "pattern",
  all.x = TRUE
)

# Turn into matrix: rows are samples, columns are all_matches
finemo_hits_counts_mat <- reshape2::acast(
  finemo_hits_counts,
  sample ~ all_matches,
  value.var     = "Freq",
  fun.aggregate = sum
)

finemo_hits_counts_mat[is.na(finemo_hits_counts_mat)] <- 0

###############################################################################
# 3. Heatmaps of FiNeMo hit counts per program and motif
###############################################################################

figure_file <- file.path(
  figure_dir,
  paste0("Finemo_merged_counts_", counts_or_profile, ".pdf")
)

pdf(figure_file, width = 30, height = 6, useDingbats = FALSE)

# Use only motifs with at least 5000 hits across all programs
selected_cols <- colSums(finemo_hits_counts_mat) > 5000
counts_mat    <- finemo_hits_counts_mat[, selected_cols, drop = FALSE]

ComplexHeatmap::Heatmap(
  counts_mat,
  cluster_columns        = TRUE,
  cluster_rows           = TRUE,
  heatmap_legend_param   = list(title = "Freq")
)

ComplexHeatmap::Heatmap(
  log10(counts_mat + 1),
  cluster_columns        = TRUE,
  cluster_rows           = TRUE,
  heatmap_legend_param   = list(title = "log10(Freq)")
)

ComplexHeatmap::Heatmap(
  scale(counts_mat),
  cluster_columns        = TRUE,
  cluster_rows           = TRUE,
  heatmap_legend_param   = list(title = "scaled(Freq)")
)

dev.off()

###############################################################################
# 4. Summed FiNeMo importance per program and motif
###############################################################################

# Sum hit_importance per (Pr, motif_name)
sum_importance_dt <- as.data.table(finemo_hits)[
  ,
  .(sum_importance = sum(hit_importance, na.rm = TRUE)),
  by = .(Pr, motif_name)
][order(motif_name, Pr)]

# Order programs on the x axis using the predefined program_order
sum_importance_dt$Pr <- factor(sum_importance_dt$Pr, levels = program_order)

# Output file for the summed importance barplot
figure_file <- file.path(
  figure_dir,
  paste0("Finemo_merged_scores_", counts_or_profile, ".pdf")
)

pdf(figure_file, width = 20, height = 20, useDingbats = FALSE)

# Bar plot of summed importance per program, faceted by motif
ggplot(
  sum_importance_dt,
  aes(x = Pr, y = sum_importance, fill = Pr)
) +
  geom_bar(stat = "identity") +
  facet_wrap(
    ~ sub("^pos_patterns\\.pattern_", "", motif_name),
    scale = "free_y"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

dev.off()


###############################################################################
# 5. CEBPB 7 groups: annotation of overlaps between CEBPB groups and ATAC peaks
###############################################################################

# Group codes to use
selected_groups <- c("RNA_cluster_1", "RNA_cluster_2", "RNA_cluster_3", "RNA_cluster_4")

# Directory with one BED file per CEBPB group
cebpb_groups_bed_dir <- Sys.getenv("TSS_GROUP_BED_DIR")

# Map CEBPB group code to Pr programs
group_to_programs <- list(
  "RNA_cluster_1" = "MEF",
  "RNA_cluster_2" = "NSC",
  "RNA_cluster_3" = "ESC",
  "RNA_cluster_4" = "Bcell_Graf"
)



# Map Pr labels to ATAC sample names in the interpreted regions file
program_to_sample_map <- c(
  "MEF"                = "MEF",
  "NSC"        = "NSC",
  "ESC" = "ESC",
  "Bcell_Graf" = "Bcell_Graf"

)


# Read ATAC interpreted regions (one row per peak)
# Expected columns: chr, peak_start, peak_end, peak_name, sample, Pr, ...
atac_peaks_df <- fread(
  file.path(finemo_merged_dir, "pos_patterns/merged_allchrs.interpreted_regions.bed")
)

# Store total number of CEBPB peaks per group
group_sizes <- data.frame(
  group      = character(),
  group_size = integer(),
  stringsAsFactors = FALSE
)

# Store all overlaps between ATAC peaks and CEBPB groups
atac_groups_overlap_df <- data.frame(
  peak_name = character(),
  group     = character(),
  Pr        = character(),
  cebpb_id  = character(),
  stringsAsFactors = FALSE
)

# For each CEBPB group: import BED, expand to all relevant programs, and intersect
for (group_code in selected_groups) {
  message("Processing CEBPB group: ", group_code)
  
  # Import CEBPB group peaks
  cebpb_gr <- rtracklayer::import(
    file.path(cebpb_groups_bed_dir, paste0(group_code, ".bed"))
  )
  
  # Add integer id and group label
  cebpb_gr$id    <- seq_along(cebpb_gr)
  cebpb_gr$group <- group_code
  
  # Record group size
  group_sizes <- rbind(
    group_sizes,
    data.frame(
      group      = group_code,
      group_size = length(cebpb_gr),
      stringsAsFactors = FALSE
    )
  )
  
  # Determine which programs are active for this group
  pr_string <- group_to_programs[[group_code]]
  pr_values <- strsplit(pr_string, "\\|")[[1]]
  
  # Expand CEBPB peaks per program (Pr)
  cebpb_gr_expanded <- do.call(
    c,
    lapply(pr_values, function(pr_val) {
      g     <- cebpb_gr
      g$Pr  <- pr_val
      g
    })
  )
  
  # Unique id per CEBPB peak and program combination
  cebpb_gr_expanded$uid <- paste0("id", cebpb_gr_expanded$id, "_", cebpb_gr_expanded$Pr)
  
  # Convert programs to ATAC sample names
  sample_values <- unname(program_to_sample_map[pr_values])
  sample_values <- sample_values[!is.na(sample_values)]
  
  # Restrict ATAC peaks to these samples
  atac_peaks_in_pr_df <- atac_peaks_df %>%
    filter(sample %in% sample_values)


  # Invert the map: sample -> program
  sample_to_program <- setNames(
  names(program_to_sample_map),
  program_to_sample_map
  )

  # Create Pr column by mapping sample and removing trailing "flag"
  atac_peaks_in_pr_df$Pr <- sample_to_program[atac_peaks_in_pr_df$sample]
  atac_peaks_in_pr_df$Pr <- sub("flag$", "", atac_peaks_in_pr_df$Pr)

  
  # Convert to GRanges
  atac_peaks_in_pr_gr <- makeGRangesFromDataFrame(
    atac_peaks_in_pr_df,
    keep.extra.columns = TRUE,
    seqnames.field     = "chr",
    start.field        = "peak_start",
    end.field          = "peak_end",
    ignore.strand      = TRUE
  )


  common_pr <- intersect(
  unique(atac_peaks_in_pr_gr$Pr),
  unique(cebpb_gr_expanded$Pr)
  )


  for (Pr in common_pr) {
    print(Pr)

    atac_sub  <- atac_peaks_in_pr_gr[atac_peaks_in_pr_gr$Pr == Pr]
    cebpb_sub <- cebpb_gr_expanded[cebpb_gr_expanded$Pr == Pr]

    # Find overlaps between ATAC peaks and expanded CEBPB peaks
    overlaps <- findOverlaps(
    atac_sub,
    cebpb_sub,
    ignore.strand = TRUE
    )

    # Build overlap table for this group
    atac_group_overlap_this_group <- data.frame(
    peak_name = atac_sub$peak_name[queryHits(overlaps)],
    group     = cebpb_sub$group[subjectHits(overlaps)],
    Pr        = cebpb_sub$Pr[subjectHits(overlaps)],
    cebpb_id  = cebpb_sub$uid[subjectHits(overlaps)],
    stringsAsFactors = FALSE
    )

    # Append to global overlap table
    atac_groups_overlap_df <- rbind(
    atac_groups_overlap_df,
    atac_group_overlap_this_group
    )
  }

}


# Define output file
csv_file <- file.path(
  figure_dir,
  "atac_TSS_4group_overlap.csv"
)

# Ensure output directory exists
dir.create(dirname(csv_file), showWarnings = FALSE, recursive = TRUE)

# Save the data frame
write.csv(
  atac_groups_overlap_df,
  file = csv_file,
  row.names = FALSE
)



###############################################################################
# 6. Merge FiNeMo hits with CEBPB overlaps and annotate motifs
###############################################################################

# Restrict FiNeMo hits to the peaks that overlap a CEBPB group
finemo_hits_cebpb_overlap <- merge(
  finemo_hits,
  atac_groups_overlap_df,
  by = c("peak_name", "Pr")
)

# unique(finemo_hits_cebpb_overlap$Pr)

# Add MoDISco motif annotation (all_matches) based on motif_name
finemo_hits_cebpb_overlap <- merge(
  finemo_hits_cebpb_overlap,
  modisco_table[, c("pattern", "all_matches")],
  by.x  = "motif_name",
  by.y  = "pattern",
  all.x = TRUE
)

# Add a clean motif index and a composite facet label for plotting
finemo_hits_cebpb_overlap <- finemo_hits_cebpb_overlap %>%
  mutate(
    motif_clean = sub("^pos_patterns\\.pattern_", "", motif_name),
    facet_label = paste0("Pattern ", motif_clean, " | ", all_matches)
  )



# Define output file
csv_file <- file.path(
  figure_dir,
  "pattren_names.csv"
)

# Ensure output directory exists
dir.create(dirname(csv_file), showWarnings = FALSE, recursive = TRUE)


# select columns
pattren_names <- finemo_hits_cebpb_overlap[, c("motif_name", "motif_clean", "all_matches")]

# drop duplicate rows
pattren_names <- unique(pattren_names)

# sort by motif_clean as integer
pattren_names <- pattren_names[order(as.integer(pattren_names$motif_clean)), ]

# write to CSV
write.csv(pattren_names, file = csv_file, row.names = FALSE)



###############################################################################
# 7. Fraction of CEBPB peaks with a given motif per group and program
###############################################################################

# Count unique CEBPB peaks with a given motif for each (group, program, motif)
cebpb_overlap_count_motifs <- finemo_hits_cebpb_overlap %>%
  distinct(group, Pr, motif_name, all_matches, cebpb_id) %>%
  count(
    group, Pr, motif_name, all_matches,
    name = "n_cebp_peaks_with_motif"
  ) %>%
  complete(
    group, Pr,
    nesting(motif_name, all_matches),
    fill = list(n_cebp_peaks_with_motif = 0L)
  ) %>%
  left_join(group_sizes, by = "group") %>%
  mutate(
    frac_cebp_with_motif = n_cebp_peaks_with_motif / group_size
  )

# Map motif_name to facet label (used later for plotting)
facet_map <- finemo_hits_cebpb_overlap %>%
  distinct(motif_name, facet_label)

# Attach facet labels and motif indices for ordering
cebpb_overlap_count_motifs_plot <- cebpb_overlap_count_motifs %>%
  left_join(facet_map, by = "motif_name") %>%
  mutate(
    motif_idx   = as.integer(sub(".*_", "", motif_name))
  ) %>%
  arrange(motif_idx) %>%
  mutate(
    facet_label = factor(facet_label, levels = unique(facet_label)),
    group       = factor(group, levels = selected_groups)
  )

###############################################################################
# 8. Heatmap of fraction of CEBPB peaks with motif per group and program
###############################################################################

figure_file <- file.path(
  figure_dir,
  paste0("Frac_TSS_peaks_with_motif_heatmap_", counts_or_profile, ".pdf")
)

dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)

pdf(figure_file, width = 12, height = 50, useDingbats = FALSE)

ggplot(
  cebpb_overlap_count_motifs_plot,
  aes(x = Pr, y = group, fill = frac_cebp_with_motif)
) +
  geom_tile(color = "white", size = 0.2) +
  facet_wrap(~ facet_label, ncol = 6) +
  scale_fill_gradient(
    name = "Fraction of TSS peaks",
    low  = "white",
    high = "steelblue"
  ) +
  theme_bw() +
  theme(
    strip.text       = element_text(size = 5),
    axis.text.x      = element_text(
      angle  = 90,
      vjust  = 0.5,
      hjust  = 1,
      size   = 7
    ),
    axis.text.y      = element_text(size = 7),
    axis.title.x     = element_blank(),
    axis.title.y     = element_blank(),
    panel.grid       = element_blank(),
    legend.key.size  = unit(3, "mm"),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7)
  )

dev.off()

###############################################################################
# 9. Barplot of fraction of CEBPB peaks with motif per group and program
###############################################################################

# Motif annotation for barplot labels
motif_ann <- finemo_hits_cebpb_overlap %>%
  mutate(
    motif_clean = sub("^pos_patterns\\.pattern_", "", motif_name)
  ) %>%
  distinct(motif_name, motif_clean, all_matches)

# Build labels per (motif_name, group)
facet_map_bar <- expand.grid(
  motif_name       = unique(motif_ann$motif_name),
  group            = selected_groups,
  stringsAsFactors = FALSE
) %>%
  left_join(motif_ann, by = "motif_name") %>%
  mutate(
    facet_label = paste0(
      "Group ", group,
      " | Pattern ", motif_clean,
      " | ", all_matches
    )
  )

# Attach these labels to the fraction table
cebpb_frac_for_bar <- cebpb_overlap_count_motifs %>%
  left_join(facet_map_bar, by = c("motif_name", "group")) %>%
  mutate(
    motif_idx   = as.integer(sub(".*_", "", motif_name))
  ) %>%
  arrange(motif_idx, group) %>%
  mutate(
    facet_label = factor(facet_label, levels = unique(facet_label)),
    group       = factor(group, levels = selected_groups)
  )

figure_file <- file.path(
  figure_dir,
  paste0("Frac_TSS_peaks_with_motif_bar_", counts_or_profile, ".pdf")
)

pdf(figure_file, width = 25, height = 500, useDingbats = FALSE)

ggplot(
  cebpb_frac_for_bar,
  aes(x = Pr, y = frac_cebp_with_motif, fill = Pr)
) +
  geom_bar(stat = "identity") +
  facet_wrap(
    ~ facet_label,
    ncol = 5
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x   = "Program (Pr)",
    y   = "Fraction of TSS peaks with motif",
    fill = "Program"
  ) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
    axis.text.y      = element_text(size = 7),
    axis.title.x     = element_text(size = 9),
    axis.title.y     = element_text(size = 9),
    strip.text       = element_text(size = 10),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_blank(),
    legend.key.size  = unit(3, "mm"),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 7)
  )

dev.off()

###############################################################################
# 10. Summed FiNeMo importance per group, program, and motif
###############################################################################

# Sum hit_importance for each (Pr, motif, CEBPB group)
sum_scores <- as.data.table(finemo_hits_cebpb_overlap)[
  ,
  .(sum_importance = sum(hit_importance, na.rm = TRUE)),
  by = .(Pr, motif_name, group)
][order(motif_name, group, Pr)]

# Ensure all combinations of Pr, motif_name, and group exist
sum_scores <- sum_scores %>%
  complete(Pr, motif_name, group, fill = list(sum_importance = 0))

# Merge group_sizes so that we can normalize by group_size
sum_scores <- merge(
  sum_scores,
  group_sizes,
  by   = "group",
  all.x = TRUE
)

# Normalized sum importance per CEBPB peak in the group
sum_scores$sum_imp_per_size <- sum_scores$sum_importance / sum_scores$group_size

# Add MoDISco annotation
sum_scores <- merge(
  sum_scores,
  modisco_table[, c("pattern", "all_matches")],
  by.x  = "motif_name",
  by.y  = "pattern",
  all.x = TRUE
)

# Build facet labels for sum importance plots
sum_scores <- sum_scores %>%
  mutate(
    motif_clean = sub("^pos_patterns\\.pattern_", "", motif_name),
    facet_label = paste0(
      "Group ", group,
      " | Pattern ", motif_clean,
      " | ", all_matches
    )
  )

sum_scores$facet_label <- factor(
  sum_scores$facet_label,
  levels = sum_scores %>%
    arrange(as.numeric(motif_clean), group) %>%
    pull(facet_label) %>%
    unique()
)

# Add a column that stores the maximum y value across all groups and Pr for each motif.

sum_scores <- sum_scores %>%
  group_by(motif_clean) %>%
  mutate(row_ymax = max(sum_imp_per_size, na.rm = TRUE)) %>%
  ungroup()

setDT(sum_scores)

# Plot sum_imp_per_size per program, faceted by group and motif
figure_file <- file.path(
  figure_dir,
  paste0("sum_imp_per_size_bar_", counts_or_profile, ".pdf")
)

pdf(figure_file, width = 25, height = 300, useDingbats = FALSE)

ggplot(sum_scores, aes(x = Pr, y = sum_imp_per_size, fill = Pr)) +
  geom_bar(stat = "identity") +
  geom_blank(aes(y = row_ymax)) +
  facet_wrap(
    ~ facet_label,
    ncol   = 5,
    scales = "free_y"
  ) +
  labs(
    x   = "Program (Pr)",
    y   = "(Summed FiNeMo importance) / group size",
    fill = "Program"
  ) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
    axis.text.y  = element_text(size = 7),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    strip.text   = element_text(size = 10)
  )

dev.off()

# ggplot(sum_scores, aes(x = Pr, y = sum_imp_per_size, fill = Pr)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(
#     ~ facet_label,
#     ncol   = 5,
#     scales = "free_y"
#   ) +
#   labs(
#     x   = "Program (Pr)",
#     y   = "(Summed FiNeMo importance) / group size",
#     fill = "Program"
#   ) +
#   theme_bw() +
#   theme(
#     axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1, size = 7),
#     axis.text.y      = element_text(size = 7),
#     axis.title.x     = element_text(size = 9),
#     axis.title.y     = element_text(size = 9),
#     strip.text       = element_text(size = 10),
#     panel.grid.major = element_line(size = 0.2),
#     panel.grid.minor = element_blank(),
#     legend.key.size  = unit(3, "mm"),
#     legend.title     = element_text(size = 8),
#     legend.text      = element_text(size = 7)
#   )

# dev.off()


###############################################################################
# 11. Co occurrence analysis per program and CEBPB group
###############################################################################

# library(data.table)
# library(dplyr)
# library(ComplexHeatmap)

# # Load motif occurrence matrix from FiNeMo
# # Rows: peaks (peak_id), columns: motifs, plus "total" and "peak_id"
# motif_occurrences <- as.data.frame(
#   fread(
#     file.path(
#       finemo_merged_dir,
#       "pos_patterns/finemo",
#       selected_finemodel,
#       "motif_occurrences.tsv"
#     )
#   )
# )

# # Load interpreted regions table to map peak_id to sample and peak_name
# peak_ids <- as.data.frame(
#   fread(
#     file.path(
#       finemo_merged_dir,
#       "pos_patterns/merged_allchrs.interpreted_regions.bed"
#     )
#   )
# )

# # very peak_id in motif_occurrences exists in peak_ids
# stopifnot(all(motif_occurrences$peak_id %in% rownames(peak_ids)))


# # Make rownames of peak_ids correspond to the zero based peak_id used by FiNeMo
# row.names(peak_ids) <- 0:(nrow(peak_ids) - 1)

# # Reorder peak_ids so that rows match motif_occurrences$peak_id
# peak_ids <- peak_ids[as.character(motif_occurrences[, "peak_id"]), ]

# # Use peak_id as rownames for motif_occurrences
# row.names(motif_occurrences) <- motif_occurrences[, "peak_id"]

# # Drop non motif columns
# motif_occurrences <- motif_occurrences[
#   ,
#   !colnames(motif_occurrences) %in% c("total", "peak_id")
# ]

# # Output file for co occurrence heatmaps
# figure_file <- file.path(
#   figure_dir,
#   "CEBPB_7groups_MEF_NSC",
#   paste0("Finemo_merged_counts_cooccurrences_", counts_or_profile, ".pdf")
# )

# dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)

# pdf(figure_file, width = 30, height = 30, useDingbats = FALSE)

# # Helper table to map motif pattern id to all_matches label
# motif_names_df <- as.data.frame(
#   modisco_table[, c("pattern", "all_matches")]
# )
# row.names(motif_names_df) <- motif_names_df$pattern

# # Loop over programs and CEBPB groups
# for (selected_program in c("MEF", "MEF_d2_OSKM", "MEF_d2_OSKM_mCEBPA", "MEF_d2_OSKM_mCEBPB")) {
#   message("Processing program for co occurrence: ", selected_program)
  
#   selected_sample <- program_to_sample_map[[selected_program]]
  
#   for (group_code in selected_groups) {
#     message("  Group: ", group_code)
    
#     # CEBPB peaks for this group and program
#     group_peaks <- atac_groups_overlap_df[
#       atac_groups_overlap_df$Pr    == selected_program &
#       atac_groups_overlap_df$group == group_code,
#       ,
#       drop = FALSE
#     ]
    
#     # Peak ids whose peak_name appears in these overlaps
#     rows_group_pr <- rownames(peak_ids)[peak_ids$peak_name %in% group_peaks$peak_name]
    
#     # Sanity check: restrict to peaks that also have the expected sample label
#     rows_sample <- rownames(peak_ids)[peak_ids$sample == selected_sample]
    
#     selected_rows <- intersect(rows_group_pr, rows_sample)
#     selected_rows <- intersect(selected_rows, rownames(motif_occurrences))
    
#     motif_occurrences_subset <- motif_occurrences[selected_rows, , drop = FALSE]
    
#     message("    n peaks in subset = ", nrow(motif_occurrences_subset))
    
#     # Require a minimum number of peaks for stable co occurrence estimates
#     if (nrow(motif_occurrences_subset) > 1000) {
#       message("    Plotting co occurrence heatmaps for ", selected_program, " / ", group_code)
      
#       # Co occurrence counts: co_counts(i, j) = sum over peaks of (count_i * count_j)
#       mat       <- as.matrix(motif_occurrences_subset)
#       co_counts <- t(mat) %*% mat
      
#       # Diagonal entries and normalization factor
#       d <- diag(co_counts)
#       cov_norm <- 1 / sqrt(d)
#       cov_norm[!is.finite(cov_norm)] <- 0
      
#       # Normalize co occurrence: D^{-1/2} * co_counts * D^{-1/2}
#       Dinv     <- diag(cov_norm)
#       norm_mat <- Dinv %*% co_counts %*% Dinv
      
#       row.names(norm_mat) <- row.names(co_counts)
#       colnames(norm_mat)  <- row.names(co_counts)
      
#       # Map motif ids (pattern) to MoDISco all_matches labels
#       row.names(norm_mat) <- motif_names_df[row.names(co_counts), "all_matches"]
#       colnames(norm_mat)  <- motif_names_df[row.names(co_counts), "all_matches"]
      
#       # Full co occurrence heatmap
#       h_full <- ComplexHeatmap::Heatmap(
#         norm_mat,
#         cluster_columns      = TRUE,
#         cluster_rows         = TRUE,
#         heatmap_legend_param = list(
#           title = paste(selected_program, group_code)
#         )
#       )
#       print(h_full)
      
#       # Filter motifs with strong total co occurrence signal (heuristic threshold rowSums > 7)
#       selected_motifs <- row.names(
#         norm_mat[rowSums(norm_mat) > 7, , drop = FALSE]
#       )
      
#       if (length(selected_motifs) > 1) {
#         norm_mat_filtered <- norm_mat[selected_motifs, selected_motifs, drop = FALSE]
        
#         h_filtered <- ComplexHeatmap::Heatmap(
#           norm_mat_filtered,
#           cluster_columns      = TRUE,
#           cluster_rows         = TRUE,
#           heatmap_legend_param = list(
#             title = paste(selected_program, group_code, "(filtered)")
#           )
#         )
#         print(h_filtered)
#       } else {
#         message(
#           "    Skipping filtered heatmap for ",
#           selected_program,
#           " / ",
#           group_code,
#           " because <= 1 motif passed rowSums > 7 filter"
#         )
#       }
#     } else {
#       message(
#         "    Skipping ",
#         selected_program,
#         " / ",
#         group_code,
#         " because only ",
#         nrow(motif_occurrences_subset),
#         " peaks (< 1000)"
#       )
#     }
#   }
# }

# dev.off()

###############################################################################
# 12. Count based clustering across CEBPB group and program (heatmap)
###############################################################################

# Build matrix: rows are motifs (all_matches), columns are "group_Pr" combinations
cebpb_overlap_counts_table <- as.data.frame(
  table(finemo_hits_cebpb_overlap[, c("Pr", "motif_name", "group")])
)

cebpb_overlap_counts_annot <- merge(
  cebpb_overlap_counts_table,
  modisco_table[, c("pattern", "all_matches")],
  by.x  = "motif_name",
  by.y  = "pattern",
  all.x = TRUE
)

motif_count_mat <- reshape2::acast(
  cebpb_overlap_counts_annot,
  all_matches ~ group + Pr,
  value.var     = "Freq",
  fun.aggregate = sum
)

# Compute max fraction of CEBPB peaks with motif across all group Pr
motif_max_frac <- cebpb_overlap_count_motifs %>%
  dplyr::group_by(all_matches) %>%
  dplyr::summarise(
    max_frac_cebp_with_motif = max(frac_cebp_with_motif, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

# Keep only motifs with max fraction >= 0.05
motifs_to_keep <- motif_max_frac %>%
  dplyr::filter(max_frac_cebp_with_motif >= 0.05) %>%
  dplyr::pull(all_matches)

# motif_count_mat_filtered <- motif_count_mat[
#   rownames(motif_count_mat) %in% motifs_to_keep,
#   ,
#   drop = FALSE
# ]

motif_count_mat_filtered <- motif_count_mat

# Z score across columns (per motif)
finemo_scaled_counts <- t(scale(t(motif_count_mat_filtered)))

# Column annotations
column_names <- colnames(finemo_scaled_counts)

# Extract CEBPB group and program from "group_Pr" column names
# cebpb_group_vec <- sub("_.*$", "", column_names)
cebpb_group_vec <- cluster_names <- sub("^(RNA_cluster_[0-9]+)_.*$", "\\1", column_names)
program_vec     <- sub("^RNA_cluster_[0-9]+_", "", column_names)

cebpb_group_vec <- factor(cebpb_group_vec, levels = selected_groups)
program_vec     <- factor(program_vec,     levels = program_order)

row_annotation_counts <- rowAnnotation(
  foo = anno_mark(
    at     = seq_len(nrow(finemo_scaled_counts)),
    labels = rownames(finemo_scaled_counts)
  )
)

library(circlize)
library(RColorBrewer)

cebpb_group_colors <- structure(
  brewer.pal(length(levels(cebpb_group_vec)), "Set1"),
  names = levels(cebpb_group_vec)
)

# program_colors_top <- c(
#   "MEF"                = "#1f8fdb",
#   "MEF_d2_OSKM"        = "#903498",
#   "MEF_d2_OSKM_mCEBPA" = "#ffb703",
#   "MEF_d2_OSKM_mCEBPB" = "#8ecae6",
#   "ESC"                = "#09989b"
# )

# program_colors_top <- c(
#   "MEF"                = "#1f8fdb",
#   "MEF_d2_OSKM"        = "#903498",
#   "MEF_d2_OSKM_mCEBPA" = "#ffb703",
#   "MEF_d2_OSKM_mCEBPB" = "#8ecae6",
#     "MEF_d2_mCEBPA"       = "#e63946",
#   "MEF_d2_mCEBPB"       = "#f4a261"
# )



program_colors_top <- c(
  "MEF"                = "#1f8fdb",
  # "MEF_d2_OSKM"        = "#903498",
  # "MEF_d2_OSKM_mCEBPA" = "#ffb703",
  # "MEF_d2_OSKM_mCEBPB" = "#8ecae6",
  #   "MEF_d2_mCEBPA"       = "#e63946",
  # "MEF_d2_mCEBPB"       = "#f4a261",
  # "ESC"                = "#09989b",
  "NSC"            = "#118ab2",
  # "NSC_d2_mCEBPA"       = "#ef476f",
  # "NSC_d2_mCEBPB"       = "#ffd166",
  # "NSC_d2_OSKM"         = "#06d6a0",
  # "NSC_d2_OSKM_mCEBPA"  = "#073b4c",
  # "NSC_d2_OSKM_mCEBPB"  = "#8338ec",
  # "Bcell_Graf_Cebpa18h"="#3a86ff",
  "Bcell_Graf"="#ffbe0b",
  # "Bcell_Graf_OKSMd1"="#9b5de5",
  # "Bcell_Graf_OKSMd2"="#00f5d4",
    "ESC"            = "#09989b"
  # "ESC_d1_hCEBPA"       = "#ff006e",
  # "ESC_d2_hCEBPA"       = "#fb5607"
)



top_annotation_counts <- HeatmapAnnotation(
  CEBPB_group = cebpb_group_vec,
  Pr          = program_vec,
  col = list(
    CEBPB_group = cebpb_group_colors,
    Pr          = program_colors_top
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 10)
)

# Column order: sort by group then program
column_order_counts <- order(cebpb_group_vec, program_vec)

cluster_dir_counts <- file.path(
  figure_dir,
  paste0("Finemo_merged_counts_cleaned_", counts_or_profile, "/clusters")
)

dir.create(cluster_dir_counts, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
for (K in 3:20) {
  figure_file <- file.path(
    cluster_dir_counts,
    paste0(
      "Finemo_merged_counts_cleaned_",
      counts_or_profile,
      "K_",
      K,
      ".pdf"
    )
  )
  
  pdf(figure_file, width = 20, height = 30, useDingbats = FALSE)
  
  ht_list <- Heatmap(
    finemo_scaled_counts,
    name              = "Scaled(freq)",
    cluster_columns   = FALSE,
    cluster_rows      = TRUE,
    show_row_names    = FALSE,
    row_names_gp      = gpar(fontsize = 16),
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 16),
    right_annotation  = row_annotation_counts,
    top_annotation    = top_annotation_counts,
    column_split      = cebpb_group_vec,
    column_order      = column_order_counts,
    k                 = K,
    border            = TRUE,
    gap               = unit(1, "mm"),
    heatmap_legend_param = list(
      title    = "Scaled(freq)",
      border   = "black",
      title_gp = gpar(col = "black", fontsize = 15)
    )
  )
  
  ht <- draw(ht_list)
  dev.off()
}

###############################################################################
# 13. Sum importance based clustering across CEBPB group and program
###############################################################################

# Build matrix using summed importance (not normalized)
motif_sum_importance_mat <- reshape2::acast(
  sum_scores,
  all_matches ~ group + Pr,
  value.var     = "sum_importance",
  fun.aggregate = sum
)

motif_sum_importance_mat_filtered <- motif_sum_importance_mat[
  rownames(motif_sum_importance_mat) %in% motifs_to_keep,
  ,
  drop = FALSE
]

finemo_scaled_sum_importance <- t(scale(t(motif_sum_importance_mat_filtered)))

cluster_dir_sum_imp <- file.path(
  figure_dir,
  paste0("Finemo_merged_sum_importance_cleaned_", counts_or_profile, "/clusters")
)

dir.create(cluster_dir_sum_imp, recursive = TRUE, showWarnings = FALSE)


library(circlize)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("blue", "white", "red")
)

set.seed(42)
for (K in 3:20) {
  figure_file <- file.path(
    cluster_dir_sum_imp,
    paste0(
      "Finemo_merged_sum_importance_cleaned_",
      counts_or_profile,
      "K_",
      K,
      ".pdf"
    )
  )
  
  pdf(figure_file, width = 25, height = 25, useDingbats = FALSE)
  
  ht_list <- Heatmap(
    finemo_scaled_sum_importance,
    col = col_fun,
    name              = "Scaled(sum_importance)",
    cluster_columns   = FALSE,
    cluster_rows      = TRUE,
    show_row_names    = FALSE,
    row_names_gp      = gpar(fontsize = 16),
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 16),
    right_annotation  = row_annotation_counts,
    top_annotation    = top_annotation_counts,
    column_split      = cebpb_group_vec,
    column_order      = column_order_counts,
    k                 = K,
    border            = TRUE,
    gap               = unit(1, "mm"),
    heatmap_legend_param = list(
      title    = "Scaled(sum_importance)",
      border   = "black",
      title_gp = gpar(col = "black", fontsize = 15)
    )
  )
  
  ht <- draw(ht_list)
  dev.off()
}

###############################################################################
# 14. Sum importance per group size based clustering across CEBPB group and program
###############################################################################

# Matrix using sum_importance normalized by group_size (sum_imp_per_size)
motif_sum_imp_per_size_mat <- reshape2::acast(
  sum_scores,
  all_matches ~ group + Pr,
  value.var     = "sum_imp_per_size",
  fun.aggregate = sum
)

# motif_sum_imp_per_size_mat_filtered <- motif_sum_imp_per_size_mat[
#   rownames(motif_sum_imp_per_size_mat) %in% motifs_to_keep,
#   ,
#   drop = FALSE
# ]

motif_sum_imp_per_size_mat <- t(scale(t(motif_sum_imp_per_size_mat)))

motif_sum_imp_per_size_mat_filtered <- motif_sum_imp_per_size_mat

library(circlize)

col_fun <- colorRamp2(
  c(-3, 0, 3),
  c("blue", "white", "red")
)

## Clip the matrix values before plotting
Z_plot <- motif_sum_imp_per_size_mat_filtered
Z_plot[Z_plot >  3] <-  3
Z_plot[Z_plot < -3] <- -3


cluster_dir_sum_imp_per_size <- file.path(
  figure_dir,
  paste0("Finemo_merged_sum_imp_per_size_cleaned_", counts_or_profile, "/clusters")
)

dir.create(cluster_dir_sum_imp_per_size, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
for (K in 3:20) {
  figure_file <- file.path(
    cluster_dir_sum_imp_per_size,
    paste0(
      "Finemo_merged_sum_imp_per_size_cleaned_",
      counts_or_profile,
      "K_",
      K,
      ".pdf"
    )
  )
  
  pdf(figure_file, width = 15, height = 25, useDingbats = FALSE)
  
  ht_list <- Heatmap(
    Z_plot,
    col = col_fun,
    name              = "Scaled(sum_importance/group_size)",
    cluster_columns   = FALSE,
    cluster_rows      = TRUE,
    show_row_names    = FALSE,
    row_names_gp      = gpar(fontsize = 16),
    show_column_names = TRUE,
    column_names_gp   = gpar(fontsize = 16),
    right_annotation  = row_annotation_counts,
    top_annotation    = top_annotation_counts,
    column_split      = cebpb_group_vec,
    column_order      = column_order_counts,
    k                 = K,
    border            = TRUE,
    gap               = unit(1, "mm"),
    heatmap_legend_param = list(
      title    = "Scaled(sum_importance/group_size)",
      border   = "black",
      at       = c(-3, -2, -1, 0, 1, 2, 3),
      title_gp = gpar(col = "black", fontsize = 15)
    )
  )
  
  ht <- draw(ht_list)
  dev.off()
}

write.csv(motif_sum_imp_per_size_mat_filtered, 
	paste(figure_dir,"/Finemo_merged_sum_imp_per_size_cleaned_","Starting_cell_types_clusters_",counts_or_profile,".csv",sep=""))

write.csv(motif_max_frac, 
	paste(figure_dir,"/TSS_cluster_motif_max_frac_",counts_or_profile,".csv",sep=""))
  
