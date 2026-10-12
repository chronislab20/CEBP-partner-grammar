# Adapted from https://github.com/GreenleafLab/HDMA
# code/03-chrombpnet/02-compendium/03-merge_modisco.py
# and https://github.com/jmschrei/tfmodisco-lite/blob/main/examples/MergeMotifsAcrossRuns.ipynb
#
# Public/repository version: server-specific absolute paths have been replaced
# by environment variables. Analysis logic and manuscript parameters are unchanged.

import random
import numpy as np
import h5py as h5
import csv
import subprocess
import modiscolite
import os
import pandas as pd
import polars as pl
from finemo import data_io
from modiscolite import aggregator
from modiscolite.core import TrackSet, Seqlet, SeqletSet
from modiscolite.report import report_motifs, make_logo
from memelite.io import read_meme

import sys
from pathlib import Path
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import Merge_modisco_functions as mmf
import importlib
importlib.reload(mmf)
from Merge_modisco_functions import TrackSet_ID, Seqlet, SimilarPatternsCollapser_ID, save_patterns_modisco_h5
from Merge_modisco_functions import run_tomtom_ID

random.seed(42)
np.random.seed(42)

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
# MODEL_ROOT should contain the 20 trained-model directories, e.g.
#   MEF_ATAC_mm10_Train_0.5/Ave/...
# MERGED_MOTIF_DIR is the output directory for the consolidated motif lexicon.
# JASPAR_MEME is the curated MEME-format known-motif database used by TOMTOM.
# KNOWN_MOTIF_LOGO_DIR stores logos generated from that database.
MODEL_ROOT = Path(os.environ.get("MODEL_ROOT", "/path/to/chrombpnet_model_selected_bias"))
project_path = Path(os.environ.get("MERGED_MOTIF_DIR", "/path/to/merged_motif_lexicon"))
motifs_meme_file = os.environ.get("JASPAR_MEME", "/path/to/curated_JASPAR_motifs.meme")
all_motifs_logo_path = Path(os.environ.get("KNOWN_MOTIF_LOGO_DIR", str(project_path / "known_motif_logos")))

# Exact 20-condition order used in the manuscript motif-merging analysis.
# tuple: (trained-model directory, sample label)
SAMPLES = [
    ("MEF_ATAC_mm10_Train_0.5", "MEF"),
    ("MEF_d2_mCEBPAflag_ATAC_mm10_Train_0.5", "MEF_d2_mCEBPAflag"),
    ("MEF_d2_mCEBPBflag_ATAC_mm10_Train_0.5", "MEF_d2_mCEBPBflag"),
    ("MEF_d2_OSKM_ATAC_L002_mm10_Train_0.5", "MEF_d2_OSKM"),
    ("MEF_d2_OSKM_hCEBPDflag_ATAC_mm10_Train_0.5", "MEF_d2_OSKM_hCEBPDflag"),
    ("MEF_d2_OSKM_mCEBPAflag_ATAC_mm10_Train_0.5", "MEF_d2_OSKM_mCEBPAflag"),
    ("MEF_d2_OSKM_mCEBPBflag_ATAC_mm10_Train_0.5", "MEF_d2_OSKM_mCEBPBflag"),
    ("NSC_ATAC_S23_L005_mm10_Train_0.6", "NSC"),
    ("NSC_d2_mCEBPAflag_ATAC_S15_mm10_Train_0.6", "NSC_d2_mCEBPAflag"),
    ("NSC_d2_mCEBPBflag_ATAC_mm10_Train_0.6", "NSC_d2_mCEBPBflag"),
    ("NSC_d2_OSKM_ATAC_mm10_Train_0.6", "NSC_d2_OSKM"),
    ("NSC_d2_OSKM_mCEBPAflag_ATAC_mm10_Train_0.6", "NSC_d2_OSKM_mCEBPAflag"),
    ("NSC_d2_OSKM_mCEBPBflag_ATAC_mm10_Train_0.6", "NSC_d2_OSKM_mCEBPBflag"),
    ("ESC_ATAC-Seq_mm10_Train_0.7", "ESC"),
    ("ESC_d1_hCEBPA_ATAC_mm10_Train_0.6", "ESC_d1_hCEBPA"),
    ("ESC_d2_hCEBPA_ATAC_mm10_Train_0.7", "ESC_d2_hCEBPA"),
    ("Bcell_Graf_Cebpa18h_merged_mm10_Train_0.7", "Bcell_Graf_Cebpa18h"),
    ("Bcell_Graf_merged_mm10_Train_0.5", "Bcell_Graf"),
    ("Bcell_Graf_OKSMd1_merged_mm10_Train_0.7", "Bcell_Graf_OKSMd1"),
    ("Bcell_Graf_OKSMd2_merged_mm10_Train_0.7", "Bcell_Graf_OKSMd2"),
]

merge_list = []
for sample_dir, label in SAMPLES:
    sample_root = MODEL_ROOT / sample_dir
    merge_list.append((
        str(sample_root / "Ave/modisco_1000_5000000/counts/merged_allchrs.counts_modisco.h5"),
        str(sample_root / "Ave/average_shaps.counts.h5"),
        str(sample_root / "fold_0/Pr_peaks_contribs_bw/peaks_shape.interpreted_regions.bed"),
        label,
    ))

# Define output directories
pattern_group = "pos_patterns"  # pos_patterns or neg_patterns
pattern_dir = project_path / pattern_group
pattern_dir.mkdir(parents=True, exist_ok=True)

log_path = str(pattern_dir / "21.6.log.txt")
merged_modisco_coords = str(pattern_dir / "merged_modisco_coordinates.tsv")
merged_peak_file = str(pattern_dir / "merged_allchrs.interpreted_regions.bed")
merged_npz_file = str(pattern_dir / "merged_allchrs.countsshap_scores_inputs.npz")
merged_modisco_h5_path = str(pattern_dir / "merged_modisco.h5")
merged_modisco_html = str(pattern_dir / "merged_modisco.html")
all_motifs_logo_path = str(all_motifs_logo_path)
project_path = str(project_path)

## Open the log file
#log_file = open(log_path, "a")  # open once, append mode
#print("message", file=log_file, flush=True)

# Define the hyperparameters for merging - these are set to the defaults in tfmodisco-lite
# tfmodisco-lite defaults (as of version 2.2.0 from https://github.com/jmschrei/tfmodisco-lite/blob/main/examples/MergeMotifsAcrossRuns.ipynb)
min_overlap = 0.7
prob_and_pertrack_sim_merge_thresholds = [(0.8,0.8), (0.5, 0.85), (0.2, 0.9)]
prob_and_pertrack_sim_dealbreaker_thresholds = [(0.4, 0.75), (0.2,0.8), (0.1, 0.85), (0.0,0.9)]
min_frac = 0.2 # also called frac_support_to_trim_to
min_num = 30 # also called min_num_to_trim_to
flank_to_add = 5 # also called initial_flank_to_add
window_size = 20 # also called trim_to_window_size
max_seqlets_subsample = 300 # also called merging_max_seqlets_subsample
# input_window = 400 # corresponds to default for tfmodisco-lite args.window
input_window = 1000 # corresponds to default for tfmodisco-lite args.window

def calculate_window_offsets(center: int, window_size: int) -> tuple:
        return (center - window_size // 2, center + window_size // 2)



###############################################
# 1) Load modisco H5 files and merge similar seqlets

# merge_list = []
# with open(merge_inputs, newline="") as f:
#     reader = csv.DictReader(f, delimiter="\t")
#     for row in reader:
#         merge_list.append((row["modisco_path"],
#                            row["shaps_path"],
#                            row["bed_path"],
#                            row["label"]))

# print(merge_list[0])


# for each program, load the modisco h5 files and extract the pattern information
exampleidx_offset = 0 #incremented after each modisco results file
seqlet_unique_index = [-1] 
union_onehot = []
union_hypscores = []
union_projscores = []
union_onehot_full = []
union_hypscores_full = []
union_projscores_full = []
union_bed = []
union_patterncoords = []
new_coordinates_bed =[]
#for (results_file, shaps_file, bed_file, sample_name) in merge_list:
for (results_file, shaps_file, bed_file, sample_name) in merge_list:
  print("On file "+results_file+" "+sample_name)
  results_fh = h5.File(results_file)
  shape_fh = h5.File(shaps_file)
  bed_fh = pd.read_csv(bed_file, sep="\t", header=None, usecols=[0,1,2,3,9], names=["chrom","start","end","peak_name","summit"])
  # adj coords relative to the summit
  motif_width=shape_fh["raw"]["seq"].shape[2]//2
  bed_fh["summit_center_start"] = (bed_fh["start"] + bed_fh["summit"]) - motif_width
  bed_fh["summit_center_end"]   = (bed_fh["start"] + bed_fh["summit"]) + motif_width
  # adj coords relative to the summit using modisco window size
  bed_fh["modisco_start"] = (bed_fh["start"] + bed_fh["summit"]) - (input_window//2)
  bed_fh["modisco_end"]   = (bed_fh["start"] + bed_fh["summit"]) + (input_window//2)
  bed_fh["sample"] = sample_name
  # extract contribution scores
  onehot = shape_fh["raw"]["seq"]
  hypscores = shape_fh["shap"]["seq"]
  projscores = shape_fh["projected_shap"]["seq"]
  allpattern_exampleidxs = []
  # Keep only peaks with seqlets:
  # first, iterate through the patterns and get all the example indices
  # (Note: "example_idx" refers to the index of the sequence that contained the seqlet)
  for pattern_name in results_fh[pattern_group].keys():
    seqlets_grp = results_fh[pattern_group][pattern_name]['seqlets']
    allpattern_exampleidxs.extend(np.array(seqlets_grp['example_idx'])) 
  # figure out the subset of indices that actually have seqlets, sort it.
  surviving_indices = sorted(list(set(allpattern_exampleidxs)))
  print(str(len(surviving_indices))+" indices had seqlets out of "+str(len(onehot)))
  # append all scores
  union_onehot_full.append(onehot)
  union_hypscores_full.append(hypscores)
  union_projscores_full.append(projscores)
  # add the scores for the subset of sequences that have scores to the 'union' list.
  center = onehot.shape[2] // 2
  start, end = calculate_window_offsets(center, input_window)
  for idx in surviving_indices:
    union_onehot.append(onehot[idx,:,start:end])
    union_hypscores.append(hypscores[idx,:,start:end])
    union_projscores.append(projscores[idx,:,start:end])
  # create an index remapping based on the subset of surviving indices
  # (we will add exampleidx_offset later)
  idx_remapping = dict(zip(surviving_indices,np.arange(len(surviving_indices))))
  #union_bed.append(bed_fh.iloc[surviving_indices])
  union_bed.append(bed_fh)
  # Now iterate through the patterns again and prep the seqlet coordinates,
  # remapping the example indices as needed.
  # We also add in exampleidx_offset to account for all the previous seqeuences
  # that have already been added to the 'union' lists
  for pattern_name in results_fh[pattern_group].keys():
    seqlets_grp = results_fh[pattern_group][pattern_name]['seqlets']
    pattern_exampleidxs = np.array(seqlets_grp['example_idx'])
    #remap the example idxs
    pattern_remapped_exampleidxs = np.array([(exampleidx_offset+idx_remapping[idx]) for idx in pattern_exampleidxs])
    pattern_start = np.array(seqlets_grp['start'])
    pattern_end   = np.array(seqlets_grp['end'])
    pattern_isrevcomp = np.array(seqlets_grp['is_revcomp'])
    seqlet_unique_index = range((max(seqlet_unique_index)+1),((max(seqlet_unique_index)+1)+len(pattern_end)))
    # add old and new peak index number
    sub_bed_fh=bed_fh.loc[pattern_exampleidxs]
    sub_bed_fh["old_peak_idx"] = pattern_exampleidxs
    sub_bed_fh["new_peak_idx"] = pattern_remapped_exampleidxs
    sub_bed_fh["sample"] = sample_name
    sub_bed_fh["pattern_name_original"] = pattern_name
    sub_bed_fh["pattern_group"] = pattern_group
    sub_bed_fh["unique_seqlet_index"] = np.array(seqlet_unique_index)
    sub_bed_fh["pattern_start"] = pattern_start
    sub_bed_fh["pattern_end"] = pattern_end
    sub_bed_fh["pattern_isrevcomp"] = pattern_isrevcomp
    new_coordinates_bed.append(sub_bed_fh)
    # set new attr with sample name
    str_dt = h5.string_dtype("utf-8")
    sample_arr = np.full(pattern_end.shape, sample_name, dtype=str_dt)
    union_patterncoords.append((pattern_remapped_exampleidxs, pattern_start, pattern_end, pattern_isrevcomp, pattern_exampleidxs, np.array(seqlet_unique_index), sample_arr))
  # increment exampleidx_offset
  exampleidx_offset = (exampleidx_offset + len(surviving_indices))
  results_fh.close()

# Write new coordinates to file (for debugging)
final_coordinates_bed = pd.concat(new_coordinates_bed, ignore_index=True, copy=False)
final_coordinates_bed.to_csv(merged_modisco_coords, sep="\t", mode="w", index=False)

# Create a merged npz file to run with finemo
# Create merged peak file
peaks_df = pd.concat(union_bed, ignore_index=True)
peaks_df = peaks_df.rename(columns={"chrom": "chr"})
peaks_df = peaks_df.rename(columns={"summit_center_start": "peak_region_start"})
peaks_df = peaks_df.rename(columns={"start": "peak_start"})
peaks_df = peaks_df.rename(columns={"end": "peak_end"})
peaks_df = peaks_df.rename(columns={"summit": "peak_summit"})
# add chr id
chrom_order = []
chrom_order_set = set(chrom_order)
chrom_order_peaks = [i for i in pd.unique(peaks_df["chr"]) if i not in chrom_order_set]
chrom_order.extend(chrom_order_peaks)
chrom_ind_map = {val: ind for ind, val in enumerate(chrom_order)}
peaks_df["chr_id"] = peaks_df["chr"].replace(chrom_ind_map)
peaks_df["peak_id"] = peaks_df.index.astype("uint32")
peaks_df["peak_start"] = peaks_df["peak_start"].astype("int32")
peaks_df["peak_end"] = peaks_df["peak_end"].astype("int32")
peaks_df["peak_score"] = 0
peaks_df["peak_score"] = peaks_df["peak_score"].astype("int32")
peaks_df["peak_strand"] = ""
peaks_df["peak_signal"] = 0
peaks_df["peak_signal"] = peaks_df["peak_signal"].astype("float32")
peaks_df["peak_pval"] = 0
peaks_df["peak_pval"] = peaks_df["peak_pval"].astype("float32")
peaks_df["peak_qval"] = 0
peaks_df["peak_qval"] = peaks_df["peak_qval"].astype("float32")
peaks_df["peak_summit"] = peaks_df["peak_summit"].astype("int32")

peaks_df.to_csv(merged_peak_file, sep="\t", mode="w", index=False)

# Write the npz file 
peaks_pl = pl.DataFrame(peaks_df)
#data_io.write_regions_npz(np.array(union_onehot), np.array(union_hypscores), merged_npz_file, peaks_df=peaks_pl)
data_io.write_regions_npz(np.concatenate(union_onehot_full, axis=0), np.concatenate(union_hypscores_full, axis=0), merged_npz_file, peaks_df=peaks_pl)


# create the trackset object for modisco
union_onehot = np.transpose(np.array(union_onehot), axes=(0,2,1))
union_hypscores = np.transpose(np.array(union_hypscores), axes=(0,2,1))
union_projscores = np.transpose(np.array(union_projscores), axes=(0,2,1))
track_set = TrackSet_ID(one_hot=union_onehot, contrib_scores=union_projscores, hypothetical_contribs=union_hypscores)

# Create pattern objects using the new track_set and modified coordinates
print("Patterns to be merged:")
all_patterns = []
#for (example_idxs, starts, ends, isrevcomps) in union_patterncoords:
for (example_idxs, starts, ends, isrevcomps, old_peak_idxs, seqlet_unique_indexs, sample_names ) in union_patterncoords:
  #tfmlite reuses the same object for representing seqlet coordinates as well as seqlets
  seqlet_coords = [Seqlet(example_idx, start, end, isrevcomp, old_peak_idx, seqlet_unique_index, sample_name) for
                   (example_idx, start, end, isrevcomp, old_peak_idx, seqlet_unique_index, sample_name) in zip(
                   example_idxs, starts, ends, isrevcomps, old_peak_idxs, seqlet_unique_indexs, sample_names)]
  seqlets = track_set.create_seqlets_ID(seqlet_coords)
  pattern = SeqletSet(seqlets) #SeqletSet in tfm lite = AggregatedSeqlet in tfm
  all_patterns.append(pattern)
  print("numseqlets:",len(pattern.seqlets))

#bg_freq is used for identifying the region within a merged motif of
# highest information content and then expanding around that region.
bg_freq = np.mean(union_onehot, axis=(0, 1))

# Do the merging
merged_patterns, pattern_merge_hierarchy = SimilarPatternsCollapser_ID(
    patterns=all_patterns,
    #patterns=all_patterns[0:2],
    track_set=track_set,
    min_overlap=min_overlap,
    prob_and_pertrack_sim_merge_thresholds=prob_and_pertrack_sim_merge_thresholds,
    prob_and_pertrack_sim_dealbreaker_thresholds=prob_and_pertrack_sim_dealbreaker_thresholds,
    min_frac=min_frac,
    min_num=min_num,
    flank_to_add=flank_to_add,
    window_size=window_size,
    bg_freq=bg_freq,
    max_seqlets_subsample=max_seqlets_subsample)

# Save all data into one merged pattern file
save_patterns_modisco_h5(merged_patterns, merged_modisco_h5_path, sign=pattern_group, sample_name="")

# Run modisco report
cmd = [
    "modisco", "report",
    "-m", motifs_meme_file,              
    "-s", merged_modisco_html,         
    "-i", merged_modisco_h5_path,          
    "-o", merged_modisco_html             
]
res = subprocess.run(cmd, check=True, text=True, capture_output=True)
print("STDOUT:\n", res.stdout)
print("STDERR:\n", res.stderr)

# Optional (if one wants to keep tomtom output)
tomtom_output = f"{project_path}/{pattern_group}/tomtom_cwm"
os.makedirs(tomtom_output, exist_ok=True)
tomtom_df = run_tomtom_ID(merged_modisco_h5_path, motifs_meme_file, tomtom_output)


# Optional, create a folder with all meme logos
os.makedirs(f"{all_motifs_logo_path}", exist_ok=True)
p = Path(all_motifs_logo_path)
if p.exists() and p.is_dir() and not any(p.iterdir()):
    os.makedirs(f"{all_motifs_logo_path}/fwd", exist_ok=True)
    os.makedirs(f"{all_motifs_logo_path}/rev", exist_ok=True)
    targets = read_meme(motifs_meme_file)
    for motif_id, P in targets.items():
        targets[motif_id]=targets[motif_id].T
        make_logo(motif_id, f"{all_motifs_logo_path}/fwd", targets)
        targets[motif_id]=targets[motif_id][::-1, ::-1]
        make_logo(motif_id, f"{all_motifs_logo_path}/rev", targets)


