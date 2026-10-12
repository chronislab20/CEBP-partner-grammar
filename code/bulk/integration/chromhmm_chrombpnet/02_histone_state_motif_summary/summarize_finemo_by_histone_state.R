# Run ChromBPNet from the Greenleaf lab
# Liu, B. B., Jessa, S., Kim, S. H., Ng, Y. T., Higashino, S. I., Marinov, G. K., ... & Greenleaf, W. J. (2025). Dissecting regulatory syntax in human development with scalable multiomics and deep learning. bioRxiv, 2025-04.
# https://github.com/GreenleafLab/HDMA/tree/main
# chrombpnet: https://github.com/GreenleafLab/HDMA/blob/main/code/03-chrombpnet/

# echo $DISPLAY
# R --no-restore --no-save

library(dplyr)
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(BiocParallel)
library(ComplexHeatmap)
library(patchwork)
library(rtracklayer)
library(reticulate)
library(BSgenome.Hsapiens.UCSC.hg38)


register(MulticoreParam(10)) # Use 100 cores

# Limit the number threads (needs to set to 1)
RhpcBLASctl::blas_set_num_threads(1)
RhpcBLASctl::omp_set_num_threads(1)
Sys.setenv(OMP_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1",
           MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1",
           NUMEXPR_NUM_THREADS="1")


# set params
project_path <- Sys.getenv("FINEMO_MERGED_DIR")
out_project_path <- Sys.getenv("INTEGRATION_OUT")
figure_dir <- paste(out_project_path,"/figures_V6/",sep="")
tables_dir <- paste(out_project_path,"/tables_V6/",sep="")

selected_finemodel <- "merged"
pattern_group      <- "pos_patterns"    # "pos_patterns" or "neg_patterns"
counts_or_profile  <- "counts"          # "counts" or "profile
selected_group <- "Combined_modisco_1000_5000000_all_CEBPB"

# Programs to keep and their colors
program_order <- c(
  "MEF",
  "MEF_d2_OSKM",
  "MEF_d2_OSKM_mCEBPA",
  "MEF_d2_OSKM_mCEBPB",
  "MEF_d2_mCEBPA",
  'MEF_d2_mCEBPB' ,
  'NSC' ,
'NSC_d2_mCEBPA' ,
'NSC_d2_mCEBPB' ,
'NSC_d2_OSKM' ,
'NSC_d2_OSKM_mCEBPA' ,
'NSC_d2_OSKM_mCEBPB' ,
  'ESC' ,
'ESC_d1_hCEBPA' ,
'ESC_d2_hCEBPA'  ,
'Bcell_Graf' ,
"Bcell_Graf_Cebpa18h",
'Bcell_Graf_OKSMd1' ,
'Bcell_Graf_OKSMd2'
)



program_colors <- c(
  "MEF"                = "#1f8fdb",
  "MEF_d2_OSKM"        = "#903498",
  "MEF_d2_OSKM_mCEBPA" = "#ffb703",
  "MEF_d2_OSKM_mCEBPB" = "#8ecae6",
    "MEF_d2_mCEBPA"       = "#e63946",
  "MEF_d2_mCEBPB"       = "#f4a261",
  "ESC"                = "#09989b",
  "NSC"            = "#118ab2",
  "NSC_d2_mCEBPA"       = "#ef476f",
  "NSC_d2_mCEBPB"       = "#ffd166",
  "NSC_d2_OSKM"         = "#06d6a0",
  "NSC_d2_OSKM_mCEBPA"  = "#073b4c",
  "NSC_d2_OSKM_mCEBPB"  = "#8338ec",
  "Bcell_Graf_Cebpa18h"="#3a86ff",
  "Bcell_Graf"="#ffbe0b",
  "Bcell_Graf_OKSMd1"="#9b5de5",
  "Bcell_Graf_OKSMd2"="#00f5d4",
    "ESC"            = "#09989b",
  "ESC_d1_hCEBPA"       = "#ff006e",
  "ESC_d2_hCEBPA"       = "#fb5607"
)

###############################################################################
# 1. Parse MoDISco motif table (motifs.html) to get all_matches labels
###############################################################################

modisco_html <- rvest::read_html(
  file.path(project_path, "pos_patterns/merged_modisco.html/motifs.html")
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
motif_names_df <- as.data.frame(
  modisco_table[, c("pattern", "all_matches")]
)

unique_patterns <- unique(motif_names_df[,c("pattern","all_matches")])
# Get unique motifs
# ordered_patterns <- unique_patterns %>% arrange(desc(all_matches)) %>% pull(pattern)
ordered_patterns <- unique(unique_patterns[,c("pattern")])


#################################################
# 2) Overlap with chromHMM

chr_order <- paste0("chr", c(1:19, "X", "Y"))

# Histone ChromHMM
Histones_chromHMM = makeGRangesFromDataFrame(fread(Sys.getenv("HISTONE_FILTERED_SEGMENTS"), sep='\t', header=F),seqnames.field="V1",start.field="V2",end.field="V3", keep.extra.columns=TRUE, ignore.strand=TRUE)
start(Histones_chromHMM) <- start(Histones_chromHMM) + 1
Histones_chromHMM <- Histones_chromHMM[seqnames(Histones_chromHMM) %in% chr_order]
seqlevels(Histones_chromHMM) <- chr_order
Histones_chromHMM <- sort(Histones_chromHMM, ignore.strand = TRUE)


# assume the state label is V4; change if yours is different
state_col <- "V4"
## 2) total bp per state
Histones_chromHMM_state_bp <- tapply(width(Histones_chromHMM),
                   mcols(Histones_chromHMM)[[state_col]],
                   sum)

Histones_chromHMM_state_bp <- sort(Histones_chromHMM_state_bp, decreasing = TRUE)

## 3) make a nice table (also include fraction)
Histones_chromHMM_state_size_df <- data.frame(
  state = names(Histones_chromHMM_state_bp),
  bp = as.numeric(Histones_chromHMM_state_bp),
  frac = as.numeric(Histones_chromHMM_state_bp) / sum(width(Histones_chromHMM)),
  row.names = NULL
)







# Tile both to 200 bps, and combine
Histones_chromHMM_200 <- unlist(tile(Histones_chromHMM, width = 200))
parent_idx <- rep(seq_along(Histones_chromHMM), elementNROWS(tile(Histones_chromHMM, width = 200)))
Histones_chromHMM_200$Histone_state <- mcols(Histones_chromHMM)[parent_idx, , drop = FALSE]$V4







###############################################################################
# 3. Load FiNeMo hits for merged run and count hits per program and motif
###############################################################################

# In this section we work only with the "merged" FiNeMo run
selected_finemodel <- "merged"

finemo_hits <- fread(
  file.path(project_path, "pos_patterns/finemo", selected_finemodel, "hits.tsv")
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

# Attach MoDISco labels
finemo_hits <- merge(
  finemo_hits,
  motif_names_df[, c("pattern", "all_matches")],
  by.x  = "motif_name",
  by.y  = "pattern",
  all.x = TRUE
)




##############################################
# 4) Motifs in ChromHMM

# send motifs to GR
finemo_hits_GR <- makeGRangesFromDataFrame(finemo_hits, keep.extra.columns=TRUE, ignore.strand=FALSE ,seqnames.field="chr")


# Get overlap between motifs and chromHMM
olaps <- findOverlaps(finemo_hits_GR, Histones_chromHMM_200)
finemo_with_states <- data.frame(pattern=finemo_hits_GR$motif_name[queryHits(olaps)], 
	sample=finemo_hits_GR$Pr[queryHits(olaps)], 
	chr=seqnames(finemo_hits_GR[queryHits(olaps)]), 
	start=start(finemo_hits_GR[queryHits(olaps)]),
	end=end(finemo_hits_GR[queryHits(olaps)]),
	Histone_state=Histones_chromHMM_200$Histone_state[subjectHits(olaps)],
	hit_importance=finemo_hits_GR$hit_importance[queryHits(olaps)])




# make unique (for cases where the motif overlap 2 200bps regions)
finemo_with_states <- unique(finemo_with_states)
finemo_with_states$sample <- factor(finemo_with_states$sample , levels=program_order)




##############################################
# filter finemo_with_states by state_sample_dict





 histone_state_sample_dict <- list(
 'H4'=c('MEF'), 
 'H5'=c('NSC'), 
 'H6'=c('ESC'), 
 'H7'=c('ESC'), 
 
 'H8'=c('Bcell_Graf'), 
 'H9'=c('MEF', 'NSC'), 
 
 
 
 'H10'=c('MEF', 'NSC', 'ESC'), 
 
 'H11'=c('MEF', 'NSC', 'ESC')
 )


# ------------------------------------------------------------
# Purpose:
# Filter finemo_with_states so that for each ATAC_state,
# we keep only rows where the sample is allowed according
# to state_sample_dict.
#
# state_sample_dict structure:
#   names  -> ATAC_state (e.g. "A1", "A2")
#   values -> vector of allowed samples for that state
#
# Strategy:
#   1. Convert the dictionary into a long lookup table
#      with columns (ATAC_state, sample).
#   2. Perform an inner merge to retain only valid
#      ATAC_state–sample combinations.
# ------------------------------------------------------------


# Step 1: Convert dictionary (list) into a lookup data.frame
# This creates one row per valid (ATAC_state, sample) pair
histone_state_sample_lookup <- do.call(
  rbind,
  lapply(names(histone_state_sample_dict), function(st) {
    data.frame(
      Histone_state = st,                     # ATAC state name
      sample = histone_state_sample_dict[[st]],    # allowed samples for that state
      stringsAsFactors = FALSE             # prevent factor conversion
    )
  })
)



# Step 2: Perform inner merge to keep only allowed combinations
# by = c("ATAC_state", "sample") ensures matching on both columns
# all = FALSE means:
#   keep only rows present in both data frames
#   (i.e., valid ATAC_state–sample pairs)


finemo_with_states_filtered <- merge(
  finemo_with_states,
  histone_state_sample_lookup,
  by = c("Histone_state", "sample"),
  all = FALSE
)


# Result:
# finemo_with_states_filtered contains only rows where
# sample is permitted for its corresponding ATAC_state.
# Output remains a standard data.frame.



############################################################
# Find patterns_keep
#
#
# Using the precomputed overlaps (olaps) between FiNeMo motif hits
# and 200 bp ChromHMM bins, this section calculates, for each motif
# pattern and each ChromHMM state:
#
#   Proportion of bins in state j that contain motif i
#   = (number of unique 200 bp bins in state j overlapping motif i)
#     / (total number of bins in state j)
#
# This is computed separately for:
#   1) Histone_state
#   2) ATAC_state
#
# Important details:
# - Each 200 bp bin is counted at most once per (pattern, state),
#   even if multiple motif hits fall in the same bin.
# - The denominator is the total number of 200 bp bins assigned
#   to a given state across the genome.
#
# Finally, motifs are retained (patterns_keep) if they reach
# a minimum proportion threshold (min_prop = 0.05) in at least
# one Histone_state or ATAC_state.
############################################################


## 2) Make a compact overlap table with bin ids
ov <- tibble(
  pattern = finemo_hits_GR$motif_name[queryHits(olaps)],
  bin_id  = subjectHits(olaps),  # index of the 200 bp bin in Histones_chromHMM_200
  Histone_state = Histones_chromHMM_200$Histone_state[subjectHits(olaps)]
)

## 3) Denominators: total bins per state
histone_total <- tibble(Histone_state = Histones_chromHMM_200$Histone_state) %>%
  count(Histone_state, name = "total_bins")

## 4) Numerators: unique bins per (pattern, state)
histone_num <- ov %>%
  distinct(pattern, Histone_state, bin_id) %>%   # ensures bin counted once even if multiple hits in same bin
  count(pattern, Histone_state, name = "bins_with_motif")


## 5) Proportions
histone_prop <- histone_num %>%
  left_join(histone_total, by = "Histone_state") %>%
  mutate(prop = bins_with_motif / total_bins) %>%
  arrange(desc(prop))



min_prop <- 0.05

patterns_keep <- histone_prop %>%
    filter(prop >= min_prop) %>%
    distinct(pattern) %>%
    pull(pattern)


length(patterns_keep)
head(patterns_keep)



############################################################

# # Count hits per program and motif
# finemo_hits_counts <- as.data.frame(
#   table(finemo_hits[, c("Pr", "motif_name")])
# )
# colnames(finemo_hits_counts) <- c("sample", "pattern", "Freq")

############################################################
# Compute Sum importance scores
sum_scores_Pr <- as.data.table(finemo_hits)[, .(sum_importance = sum(hit_importance, na.rm = TRUE)),
             by = .(Pr, motif_name)][order(motif_name, Pr)]
sum_scores_Histone_state <- as.data.table(finemo_with_states_filtered)[, .(sum_importance = sum(hit_importance, na.rm = TRUE)),
             by = .(Histone_state, pattern)][order(pattern, Histone_state)]

sum_scores_Pr_Histone_state <- as.data.table(finemo_with_states_filtered)[, .(sum_importance = sum(hit_importance, na.rm = TRUE)),
             by = .(Histone_state, sample, pattern)][order(pattern, sample, Histone_state)]

############################################################
## # Normalize motif importance by genomic fraction of each ChromHMM state
############################################################
# convert state size table to data.table
Histones_chromHMM_state_size_df_dt <- as.data.table(Histones_chromHMM_state_size_df)

# rename 'state' column to match
setnames(Histones_chromHMM_state_size_df_dt, "state", "Histone_state")

# join and normalize
sum_scores_Histone_state_norm <- merge(
  sum_scores_Histone_state,
  Histones_chromHMM_state_size_df_dt[, .(Histone_state, frac)],
  by = "Histone_state",
  all.x = TRUE
)[
  , sum_importance_over_state_size := sum_importance / frac
][
  order(pattern, Histone_state)
]



## similalry for sum_scores_Pr_Histone_state
# join and normalize
sum_scores_Pr_Histone_state_norm <- merge(
  sum_scores_Pr_Histone_state,
  Histones_chromHMM_state_size_df_dt[, .(Histone_state, frac)],
  by = "Histone_state",
  all.x = TRUE
)[
  , sum_importance_over_state_size := sum_importance / frac
][
  order(pattern, sample, Histone_state)
]






########################################################################################################################



# Change order
# sum_scores_Pr$Pr <- factor(sum_scores_Pr$Pr , levels=program_order)
sum_scores_Histone_state_norm$Histone_state <- factor(sum_scores_Histone_state_norm$Histone_state , levels=paste0("H",1:22))



# make matrices
sum_scores_Pr_heatmap <- reshape2::acast(sum_scores_Pr, motif_name~Pr, value.var="sum_importance")
sum_scores_Pr_heatmap[is.na(sum_scores_Pr_heatmap)] <- 0

sum_scores_Histone_state_heatmap <- reshape2::acast(sum_scores_Histone_state_norm, pattern~Histone_state, value.var="sum_importance_over_state_size")
sum_scores_Histone_state_heatmap[is.na(sum_scores_Histone_state_heatmap)] <- 0

sum_scores_Pr_Histone_state_heatmap <- reshape2::acast(sum_scores_Pr_Histone_state_norm, pattern~Histone_state+sample, value.var="sum_importance_over_state_size")
sum_scores_Pr_Histone_state_heatmap[is.na(sum_scores_Pr_Histone_state_heatmap)] <- 0



# create motifs list
motif_list_finemo <- list()
motif_list_modisco <- list()
prog_sel <- "merged"
for(pattern_sel in ordered_patterns){
	motif_finemo <- read.table(paste0(project_path,"/pos_patterns/finemo/",prog_sel,"/CWMs/",pattern_sel,"/hits_fc.txt")) 
	motif_modisco <- read.table(paste0(project_path,"/pos_patterns/finemo/",prog_sel,"/CWMs/",pattern_sel,"/modisco_fc.txt")) 
	rownames(motif_finemo)  <- c("A","C","G","T")
	rownames(motif_modisco) <- c("A","C","G","T")
	colnames(motif_finemo)  <- 1:ncol(motif_finemo)
	colnames(motif_modisco) <- 1:ncol(motif_modisco)
	if(length(motif_finemo[is.na(motif_finemo)])==0){
		motif_list_finemo[[paste("finemo",pattern_sel)]]  <- as.matrix(motif_finemo)
		motif_list_modisco[[paste("Modisco",pattern_sel)]] <- as.matrix(motif_modisco)
	}
}




rownames(sum_scores_Pr_heatmap)
rownames(sum_scores_Histone_state_heatmap)

sum_scores_Pr_heatmap <-
  sum_scores_Pr_heatmap[
    rownames(sum_scores_Pr_heatmap) %in% 
      rownames(sum_scores_Histone_state_heatmap),
  ]

nrow(sum_scores_Pr_heatmap)

sub_sum_scores_Histone_state_heatmap <- sum_scores_Histone_state_heatmap[,paste0("H",c(
	
	4 , #MEF Enh
    5 , #NSC Enh
    6,7 , #ESC Enh
    8 , # Bcell Enh
	9 , #MEF NSC Enh
	10,11  #MEF NSC ESC Enh
	)
	)]

scaled_sep  <- cbind(t(scale(t((sum_scores_Pr_heatmap[,program_order])))),t(scale(t((sum_scores_Histone_state_heatmap)))) )

scaled_both <-t(scale(t(((sum_scores_Pr_Histone_state_heatmap )))))




sub_histone_prop <- histone_prop[
  histone_prop$Histone_state %in% paste0("H",c(
	
	4 , #MEF Enh
    5 , #NSC Enh
    6,7 , #ESC Enh
    8 , # Bcell Enh
	9 , #MEF NSC Enh
	10,11  #MEF NSC ESC Enh
	)),  ,drop = FALSE]




min_prop <- 0.05

patterns_keep <- sub_histone_prop %>%
    filter(prop >= min_prop) %>%
    distinct(pattern) %>%
    pull(pattern)



length(patterns_keep)
head(patterns_keep)


# Convert patterns to motif names
convert_names <- as.data.frame(motif_names_df[,c("pattern","all_matches")])
row.names(convert_names) <- motif_names_df$pattern
row.names(scaled_sep) <- convert_names[row.names(scaled_sep),"all_matches"]
row.names(scaled_both)<- convert_names[row.names(scaled_both),"all_matches"]



## Filter by patterns_keep

# motifs_to_keep <-  convert_names[patterns_keep,]$all_matches

# scaled_both_filtered <- scaled_both[
#   rownames(scaled_both) %in% motifs_to_keep,
#   ,
#   drop = FALSE
# ]

## No filtering
scaled_both_filtered <- scaled_both




# !!! Make the heatmap and get rows order
# We use this order throughout
set.seed(5)
kmean <- 20
comp_ht <- (ComplexHeatmap::Heatmap(scaled_both_filtered[,matrixStats::colMaxs(scaled_both_filtered)>0], k=kmean))
figure_file<-paste(figure_dir,"/Finemo_merged_Original_map_histones_sum_importance_over_state_size_heatmap_complex","_",counts_or_profile,".pdf",sep="")
pdf(figure_file, width = 15, height = 30, useDingbats=FALSE)
row_order_scaled_both_filtered <- row_order(draw(comp_ht))
dev.off()

saveRDS(row_order_scaled_both_filtered, paste0(tables_dir,"/Pr_Histone_heatmap_sum_importance_over_state_size_",kmean,"_",".rds"))

# matrixStats::rowMaxs(scaled_both_filtered)
# matrixStats::rowMins(scaled_both_filtered)

cluster_order <- 1:kmean

figure_file<-paste(figure_dir,"/Finemo_merged_histones_histones_sum_importance_over_state_size_heatmap_complex","_",counts_or_profile,".pdf",sep="")
pdf(figure_file, width = 60, height = 60, useDingbats=FALSE)
row_order_scaled_both_filtered_mat <- data.frame(
  rowname = unlist(lapply(row_order_scaled_both_filtered, \(i) rownames(scaled_both_filtered)[i])),
  cluster = factor(rep(names(row_order_scaled_both_filtered), lengths(row_order_scaled_both_filtered)), levels=names(row_order_scaled_both_filtered))
)
row_order_scaled_both_filtered_mat$cluster <- factor(row_order_scaled_both_filtered_mat$cluster, levels = cluster_order)



dev.off()


library(circlize)

# col_fun <- colorRamp2(
#   c(-2, 0, 2),
#   # c(-3, 0, 3),
#   c("blue", "white", "red")
# )

# col_fun2 <- colorRamp2(
#   c(-2, 0, 2),
#   c("#5E3C99", "#FFFFFF", "#E66101")
# )

# col_fun2 <- colorRamp2(
#   c(-2, 0, 2),
#   c("#1B9E77", "#FFFFFF", "#D95F02")
# )

# col_fun2 <- colorRamp2(
#   c(-2, 0, 2),
#   c("navy", "grey90", "firebrick3")
# )


col_fun2 <- colorRamp2(
  c(-2, 0, 2),
  c("#3B4CC0", "#FFFFFF", "#E66101")
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#3B4CC0", "#FFFFFF", "#B40426")
)

figure_file<-paste(figure_dir,"/Finemo_merged_Pr_Histone_sum_importance_over_state_size_sep","_",counts_or_profile,".pdf",sep="")
pdf(figure_file, width = 25, height = 60, useDingbats=FALSE)

h1=ComplexHeatmap::Heatmap(scaled_sep[row_order_scaled_both_filtered_mat$rowname, unique(histone_state_sample_lookup$sample)],
	cluster_columns=F, cluster_rows=F,
	#col = circlize::colorRamp2(c(min(mvh),max(mvh)), hcl_palette ="Tropic"),
	# col = circlize::colorRamp2(seq(0,1, length = 3), c("lightgrey", "white", "#3B6FB6")),
  col = col_fun,
	row_split = row_order_scaled_both_filtered_mat$cluster,border = TRUE,gap = unit(1, "mm"),
	heatmap_legend_param = list(title = "Scaled(sum_importance)", border = "black", title_gp = gpar(col = "black", fontsize = 20)),
	width = unit(1, "null"))

h3=ComplexHeatmap::Heatmap(scaled_sep[row_order_scaled_both_filtered_mat$rowname,paste0("H",c(
	4 , #MEF Enh
    5 , #NSC Enh
    6,7 , #ESC Enh
    8 , # Bcell Enh
	9 , #MEF NSC Enh
	10,11  #MEF NSC ESC Enh
	)
	)],
	cluster_columns=F, cluster_rows=F,cluster_column_slices = FALSE,
	# col = circlize::colorRamp2(seq(0,1, length = 3), c("lightgrey", "white", "#D1495B")),
  col = col_fun2,
	row_split = row_order_scaled_both_filtered_mat$cluster,border = TRUE,gap = unit(1, "mm"),
	column_split = factor(c(
		# rep("Prom",3),
		rep("MEF\nEnh",1),
    rep("NSC\nEnh",1),
    rep("ESC\nEnh",2),
		rep("Bcell\nEnh",1),
		rep("MEF+NSC\nEnh",1),
		rep("MEF+NSC+\nESC Enh",2)
		)),
	heatmap_legend_param = list(title = "Scaled((motif_count_over_state_size))", border = "black", title_gp = gpar(col = "black", fontsize = 10)),
	width = unit(2, "null"))

# ht <- h1 + h2 + h3
ht <- h1 + h3
draw(ht, gap = unit(3, "mm"))
dev.off()

write.csv(scaled_sep[row_order_scaled_both_filtered_mat$rowname,], 
	paste(figure_dir,"/Finemo_merged_Pr_Histone_sum_importance_over_state_size_sep","_",counts_or_profile,".csv",sep=""))

write.csv(row_order_scaled_both_filtered_mat, 
	paste(figure_dir,"/Finemo_merged_Pr_Histone_sum_importance_over_state_size_sep","row_order_",counts_or_profile,".csv",sep=""))


write.csv(histone_prop, 
	paste(figure_dir,"/Finemo_merged_Pr_histone_motif_prop_",counts_or_profile,".csv",sep=""))

write.csv(convert_names, 
	paste(figure_dir,"convert_names_pattren_names_",counts_or_profile,".csv",sep=""))

  




#### plot top 5% #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
motifs_to_keep <-  convert_names[patterns_keep,]$all_matches

scaled_both_filtered <- scaled_both[
  rownames(scaled_both) %in% motifs_to_keep,
  ,
  drop = FALSE
]




# !!! Make the heatmap and get rows order
# We use this order throughout
set.seed(5)
kmean <- 20
comp_ht <- (ComplexHeatmap::Heatmap(scaled_both_filtered[,matrixStats::colMaxs(scaled_both_filtered)>0], k=kmean))
row_order_scaled_both_filtered <- row_order(draw(comp_ht))

row_order_scaled_both_filtered_mat <- data.frame(
  rowname = unlist(lapply(row_order_scaled_both_filtered, \(i) rownames(scaled_both_filtered)[i])),
  cluster = factor(rep(names(row_order_scaled_both_filtered), lengths(row_order_scaled_both_filtered)), levels=names(row_order_scaled_both_filtered))
)
row_order_scaled_both_filtered_mat$cluster <- factor(row_order_scaled_both_filtered_mat$cluster, levels = cluster_order)


library(circlize)

col_fun2 <- colorRamp2(
  c(-2, 0, 2),
  c("#3B4CC0", "#FFFFFF", "#E66101")
)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#3B4CC0", "#FFFFFF", "#B40426")
)

figure_file<-paste(figure_dir,"/Finemo_merged_Pr_Histone_sum_importance_over_state_size_sep","_",counts_or_profile,"_",min_prop,".pdf",sep="")
pdf(figure_file, width = 25, height = 30, useDingbats=FALSE)

h1=ComplexHeatmap::Heatmap(scaled_sep[row_order_scaled_both_filtered_mat$rowname, unique(histone_state_sample_lookup$sample)],
	cluster_columns=F, cluster_rows=F,
	#col = circlize::colorRamp2(c(min(mvh),max(mvh)), hcl_palette ="Tropic"),
	# col = circlize::colorRamp2(seq(0,1, length = 3), c("lightgrey", "white", "#3B6FB6")),
  col = col_fun,
	row_split = row_order_scaled_both_filtered_mat$cluster,border = TRUE,gap = unit(1, "mm"),
	heatmap_legend_param = list(title = "Scaled(sum_importance)", border = "black", title_gp = gpar(col = "black", fontsize = 20)),
	width = unit(1, "null"))

h3=ComplexHeatmap::Heatmap(scaled_sep[row_order_scaled_both_filtered_mat$rowname,paste0("H",c(
	4 , #MEF Enh
    5 , #NSC Enh
    6,7 , #ESC Enh
    8 , # Bcell Enh
	9 , #MEF NSC Enh
	10,11  #MEF NSC ESC Enh
	)
	)],
	cluster_columns=F, cluster_rows=F,cluster_column_slices = FALSE,
	# col = circlize::colorRamp2(seq(0,1, length = 3), c("lightgrey", "white", "#D1495B")),
  col = col_fun2,
	row_split = row_order_scaled_both_filtered_mat$cluster,border = TRUE,gap = unit(1, "mm"),
	column_split = factor(c(
		# rep("Prom",3),
		rep("MEF\nEnh",1),
    rep("NSC\nEnh",1),
    rep("ESC\nEnh",2),
		rep("Bcell\nEnh",1),
		rep("MEF+NSC\nEnh",1),
		rep("MEF+NSC+\nESC Enh",2)
		)),
	heatmap_legend_param = list(title = "Scaled((motif_count_over_state_size))", border = "black", title_gp = gpar(col = "black", fontsize = 10)),
	width = unit(2, "null"))

# ht <- h1 + h2 + h3
ht <- h1 + h3
draw(ht, gap = unit(3, "mm"))
dev.off()



