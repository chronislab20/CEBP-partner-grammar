
from modiscolite.aggregator import *
from modiscolite.aggregator import _align_patterns, PatternMergeHierarchyNode
import os
import h5py as h5
import pandas
import tempfile
from modiscolite.report import write_meme_file

class Seqlet(object):
	def __init__(self, example_idx, start, end, is_revcomp, old_peak_idx, seqlet_unique_index, sample_name):
		self.example_idx = example_idx
		self.start = start
		self.end = end
		self.is_revcomp = is_revcomp
		self.old_peak_idx = old_peak_idx
		self.seqlet_unique_index = seqlet_unique_index
		self.sample_name = sample_name

		self.sequence = None
		self.contrib_scores = None
		self.hypothetical_contribs = None

		super(Seqlet, self).__init__()

	def __str__(self):
		return ("example:"+str(self.example_idx)
				+",start:"+str(self.start)+",end:"+str(self.end)
				+",rc:"+str(self.is_revcomp)
				+",old_peak_idx:"+str(self.old_peak_idx)+",seqlet_unique_index:"+str(self.seqlet_unique_index)+",sample_name:"+str(self.sample_name))

	def __len__(self):
		return self.end - self.start

	@property
	def string(self):
		return str(self.example_idx)+"_"+str(self.start)+"_"+str(self.end)

	def revcomp(self):
		new_seqlet = Seqlet(
				example_idx=self.example_idx,
				start=self.start, end=self.end,
				is_revcomp=(self.is_revcomp==False),
				old_peak_idx=self.old_peak_idx,seqlet_unique_index=self.seqlet_unique_index,sample_name=self.sample_name)

		new_seqlet.sequence = self.sequence[::-1, ::-1]
		new_seqlet.contrib_scores = self.contrib_scores[::-1, ::-1]
		new_seqlet.hypothetical_contribs = self.hypothetical_contribs[::-1, ::-1]
		return new_seqlet

	def shift(self, shift_amt):
		return Seqlet(
				example_idx=self.example_idx,
				start=self.start+shift_amt, end=self.end+shift_amt,
				is_revcomp=self.is_revcomp, 
				old_peak_idx=self.old_peak_idx,seqlet_unique_index=self.seqlet_unique_index,sample_name=self.sample_name)

	def trim(self, start_idx, end_idx):
		if self.is_revcomp == False:
			new_start = self.start + start_idx 
			new_end = self.start + end_idx
		else:
			new_start = self.end - end_idx
			new_end = self.end - start_idx

		new_seqlet = Seqlet(example_idx=self.example_idx,
			start=new_start, end=new_end, is_revcomp=self.is_revcomp,
			old_peak_idx=self.old_peak_idx,seqlet_unique_index=self.seqlet_unique_index,sample_name=self.sample_name)

		s, e = start_idx, end_idx
		new_seqlet.sequence = self.sequence[s:e]
		new_seqlet.contrib_scores = self.contrib_scores[s:e]
		new_seqlet.hypothetical_contribs = self.hypothetical_contribs[s:e]
		return new_seqlet


class TrackSet_ID(object):
	def __init__(self, one_hot, contrib_scores, hypothetical_contribs):
		self.one_hot = one_hot
		self.contrib_scores = contrib_scores
		self.hypothetical_contribs = hypothetical_contribs
		self.length = len(one_hot[0])

	def create_seqlets_ID(self, seqlets):
		for seqlet in seqlets:
			idx = seqlet.example_idx
			s, e = seqlet.start, seqlet.end
			old_peak_idx=seqlet.old_peak_idx
			seqlet_unique_index=seqlet.seqlet_unique_index
			sample_name=seqlet.sample_name

			if seqlet.is_revcomp:
				seqlet.sequence = self.one_hot[idx][s:e][::-1, ::-1]
				seqlet.contrib_scores = self.contrib_scores[idx][s:e][::-1, ::-1]
				seqlet.hypothetical_contribs = self.hypothetical_contribs[idx][s:e][::-1, ::-1]
			else:
				seqlet.sequence = self.one_hot[idx][s:e]
				seqlet.contrib_scores = self.contrib_scores[idx][s:e]
				seqlet.hypothetical_contribs = self.hypothetical_contribs[idx][s:e]				

		return seqlets


def polish_pattern(pattern, min_frac, min_num, track_set, flank, window_size, bg_freq):
	pattern = pattern.trim_to_support(min_frac=min_frac, min_num=min_num)

	pattern = _expand_seqlets_to_fill_pattern(pattern, 
		track_set=track_set, left_flank_to_add=flank,
		right_flank_to_add=flank)

	if pattern is None:
		return None

	# Trim by IC
	ppm = pattern.sequence
	per_pos_ic = util.compute_per_position_ic(
		ppm=ppm, background=bg_freq, pseudocount=0.001)

	start_idx = np.argmax(util.cpu_sliding_window_sum(
		arr=per_pos_ic, window_size=window_size))

	end_idx = start_idx + window_size
	pattern = pattern.trim_to_idx(start_idx=start_idx, end_idx=end_idx)

	pattern = _expand_seqlets_to_fill_pattern(pattern, 
		track_set=track_set, left_flank_to_add=flank,
		right_flank_to_add=flank)

	return pattern

def _expand_seqlets_to_fill_pattern(pattern, track_set, left_flank_to_add,
	right_flank_to_add):

	new_seqlets = []
	for seqlet in pattern.seqlets:
		left_expansion = left_flank_to_add 
		right_expansion = len(pattern) - len(seqlet) + right_flank_to_add

		if seqlet.is_revcomp == False:
			start = seqlet.start - left_expansion
			end = seqlet.end + right_expansion
		else:
			start = seqlet.start - right_expansion
			end = seqlet.end + left_expansion
		
		if start >= 0 and end <= track_set.length:
			seqlet = track_set.create_seqlets_ID(
				seqlets=[Seqlet(example_idx=seqlet.example_idx,
					start=start, end=end, is_revcomp=seqlet.is_revcomp, old_peak_idx=seqlet.old_peak_idx, seqlet_unique_index=seqlet.seqlet_unique_index, sample_name=seqlet.sample_name)])[0]

			new_seqlets.append(seqlet)

	if len(new_seqlets) > 0:
		return core.SeqletSet(seqlets=new_seqlets)
	else:
		return None



def merge_in_seqlets_filledges_ID(parent_pattern, seqlets_to_merge,
	track_set, metric, min_overlap, transformer='l1', 
	include_hypothetical=True):

	parent_pattern = parent_pattern.copy()

	for seqlet in seqlets_to_merge:
		alnmt, revcomp_match, alnmt_score = _align_patterns(parent_pattern, 
			seqlet, metric, min_overlap, transformer, include_hypothetical)
		
		if revcomp_match:
			seqlet = seqlet.revcomp()

		preexpansion_seqletlen = len(seqlet)

		left_expansion = max(alnmt,0)
		right_expansion = max((len(parent_pattern) - (alnmt+len(seqlet))), 0)

		if seqlet.is_revcomp == False:
			start = seqlet.start - left_expansion
			end = seqlet.end + right_expansion
		else:
			start = seqlet.start - right_expansion
			end = seqlet.end + left_expansion

		if start >= 0 and end <= track_set.length:
			seqlet = track_set.create_seqlets_ID(
				seqlets=[Seqlet(example_idx=seqlet.example_idx,
					start=start, end=end, is_revcomp=seqlet.is_revcomp, old_peak_idx=seqlet.old_peak_idx, seqlet_unique_index=seqlet.seqlet_unique_index, sample_name=seqlet.sample_name)])[0]
		else:
			continue #don't try adding this seqlet

		#also expand the pattern (if needed) so that the seqlet
		# doesn't go over the edge
		parent_left_expansion = max(0, -alnmt)
		parent_right_expansion = max(0, (alnmt+preexpansion_seqletlen)
										 - len(parent_pattern))

		if (parent_left_expansion > 0) or (parent_right_expansion > 0):
			candidate_parent_pattern = _expand_seqlets_to_fill_pattern(
				parent_pattern, track_set=track_set, 
				left_flank_to_add=parent_left_expansion,
				right_flank_to_add=parent_right_expansion)

			if candidate_parent_pattern is not None:
				parent_pattern = candidate_parent_pattern
			else:
				continue

		#add the seqlet in at alignment 0, assuming it's not already
		# part of the pattern
		if seqlet.string not in parent_pattern.unique_seqlets:
			parent_pattern._add_seqlet(seqlet=seqlet)

	return parent_pattern

def SimilarPatternsCollapser_ID(patterns, track_set,
	min_overlap, prob_and_pertrack_sim_merge_thresholds,
	prob_and_pertrack_sim_dealbreaker_thresholds,
	min_frac, min_num, flank_to_add, window_size, bg_freq,
	max_seqlets_subsample=1000):
	patterns = [x.copy() for x in patterns]

	merge_hierarchy_levels = []        
	current_level_nodes = [
		PatternMergeHierarchyNode(pattern=x) for x in patterns]
	merge_hierarchy_levels.append(current_level_nodes)

	merge_occurred_last_iteration = True

	#negative numbers to indicate which
	# entries need to be filled (versus entries we can infer
	# from the previous iteration of the while loop)
	pairwise_aurocs = -np.ones((len(patterns), len(patterns)))
	pairwise_sims = np.zeros((len(patterns), len(patterns)))

	#loop until no more patterns get merged
	while merge_occurred_last_iteration:
		merge_occurred_last_iteration = False

		#Let's subsample 'patterns' to prevent runtime from being too
		# large in calculating pairwise sims. 
		subsample_patterns = []
		for pattern in patterns:
			if len(pattern.seqlets) > max_seqlets_subsample:
				subsample = np.random.RandomState(1234).choice(
					a=pattern.seqlets, replace=False, size=max_seqlets_subsample)
				pattern = core.SeqletSet(seqlets=subsample)

			subsample_patterns.append(pattern)

		n = len(patterns)
		for i in range(n):
			for j in range(n):
				#Note: I compute both i,j AND j,i because although
				# the result is the same for the sim, it can be different
				# for the auroc because a different motif is getting
				# shifted over.
				if j == i:
					pairwise_aurocs[i, j] = 0.5
					pairwise_sims[i, j] = 1.0
					continue

				if pairwise_aurocs[i, j] >= 0: #filled in from previous iter
					continue 

				#Compute best alignment between pattern pair
				alnmt, rc, aligner_sim =\
					_align_patterns(parent_pattern=patterns[i], 
						child_pattern=patterns[j], 
						metric=affinitymat.pearson_correlation, 
						min_overlap=min_overlap, 
						include_hypothetical=False,
						transformer='magnitude') 

				pairwise_sims[i, j] = aligner_sim

				#get realigned pattern2
				pattern2_coords = subsample_patterns[j].seqlets

				if rc: #flip strand if needed to align
					pattern2_coords  = [x.revcomp() for x in pattern2_coords]

				#now apply the alignment
				pattern2_coords = [x.shift((1 if x.is_revcomp else -1)*alnmt)
					for x in pattern2_coords] 

				# Filter out bad seqlets
				pattern2_coords = [seqlet for seqlet in pattern2_coords 
					if seqlet.start >= 0 and seqlet.end < track_set.length]

				if len(pattern2_coords) == 0:
					pairwise_sims[i, j] = 0.0
					pairwise_aurocs[i, j] = 0.5
					continue

				pattern2_shifted_seqlets = track_set.create_seqlets_ID(
					seqlets=pattern2_coords)

				pattern1_fwdseqdata, _ =\
				  util.get_2d_data_from_patterns(subsample_patterns[i].seqlets)

				pattern2_fwdseqdata, _ =\
				  util.get_2d_data_from_patterns(pattern2_shifted_seqlets)

				#Flatten, compute continjacc sim at this alignment
				flat_pattern1_fwdseqdata = pattern1_fwdseqdata.reshape(
					(len(pattern1_fwdseqdata), -1))
				flat_pattern2_fwdseqdata = pattern2_fwdseqdata.reshape(
					(len(pattern2_fwdseqdata), -1))

				between_pattern_sims = affinitymat.jaccard(
					flat_pattern1_fwdseqdata[:, :, None], 
					flat_pattern2_fwdseqdata[:, :, None])[:, :, 0].flatten()

				within_pattern1_sims = affinitymat.jaccard(
					flat_pattern1_fwdseqdata[:, :, None], 
					flat_pattern1_fwdseqdata[:, :, None])[:, :, 0].flatten()

				auroc = roc_auc_score(
					y_true=[0 for x in between_pattern_sims]
						   +[1 for x in within_pattern1_sims],
					y_score=list(between_pattern_sims)
							+list(within_pattern1_sims))

				#The symmetrization over i,j and j,i is done later
				pairwise_aurocs[i,j] = auroc


		#pairwise_sims is not symmetric; differ based on which pattern is
		# padded with zeros.
		patterns_to_patterns_aligner_sim = 0.5*(pairwise_sims+pairwise_sims.T)
		cross_contamination = 2*(1-np.maximum(pairwise_aurocs,0.5))

		indices_to_merge = []
		merge_partners_so_far = dict([(i, set([i])) for i in
									  range(len(patterns))])

		#merge patterns with highest similarity first
		sorted_pairs = sorted([(i,j,patterns_to_patterns_aligner_sim[i,j])
						for i in range(len(patterns))
						for j in range(len(patterns)) if (i < j)],
						key=lambda x: -x[2])

		#iterate over pairs
		for i, j, aligner_sim in sorted_pairs:
			#symmetrize asymmetric crosscontam
			# take min rather than avg to avoid aggressive merging
			cross_contam = min(cross_contamination[i,j],
								cross_contamination[j,i])

			collapse = any([(cross_contam >= x[0] and aligner_sim >= x[1])
				for x in prob_and_pertrack_sim_merge_thresholds])

			if collapse:
				collapse_passed = True
				#check compatibility for all indices that are
				#about to be merged
				merge_under_consideration = set(
					list(merge_partners_so_far[i])
					+list(merge_partners_so_far[j]))

				for m1 in merge_under_consideration:
					for m2 in merge_under_consideration:
						if (m1 < m2):
							cross_contam_here =\
								0.5*(cross_contamination[m1, m2]+
									 cross_contamination[m2, m1])
							aligner_sim_here =\
								patterns_to_patterns_aligner_sim[
									m1, m2]

							dealbreaker = any([(cross_contam_here <= x[0] and aligner_sim_here <= x[1])              
								for x in prob_and_pertrack_sim_dealbreaker_thresholds])

							if dealbreaker:
								collapse_passed = False
								break

				if collapse_passed:
					indices_to_merge.append((i,j))
					for an_idx in merge_under_consideration:
						merge_partners_so_far[an_idx]=\
							merge_under_consideration 

		for i, j in indices_to_merge:
			pattern1 = patterns[i]
			pattern2 = patterns[j]

			if pattern1 != pattern2: #if not the same object
				if len(pattern1.seqlets) < len(pattern2.seqlets):
					parent_pattern, child_pattern = pattern2, pattern1
				else:
					parent_pattern, child_pattern = pattern1, pattern2

				new_pattern = merge_in_seqlets_filledges_ID(
					parent_pattern=parent_pattern,
					seqlets_to_merge=child_pattern.seqlets,
					include_hypothetical=False,
					metric=affinitymat.pearson_correlation,
					min_overlap=min_overlap,
					transformer='magnitude',
					track_set=track_set)

				new_pattern = polish_pattern(new_pattern, min_frac=min_frac, 
					min_num=min_num, track_set=track_set, flank=flank_to_add, 
					window_size=window_size, bg_freq=bg_freq)

				if new_pattern is not None:
					for k in range(len(patterns)):
						#Replace EVERY case where the parent or child
						# pattern is present with the new pattern. This
						# effectively does single-linkage.
						if (patterns[k]==parent_pattern or
							patterns[k]==child_pattern):
							patterns[k]=new_pattern

		merge_occurred_last_iteration = (len(indices_to_merge) > 0)

		if merge_occurred_last_iteration:
			#Once we are here, each element of 'patterns'
			#will have the new parent of the corresponding element
			#of 'old_patterns'
			old_to_new_pattern_mapping = patterns

			#sort by size and remove redundant patterns
			patterns = sorted(patterns, key=lambda x: -len(x.seqlets))
			patterns = list(OrderedDict([(x,1) for x in patterns]).keys())

			#let's figure out which indices don't require recomputation
			# and use it to repopulate pairwise_sims and pairwise_aurocs
			old_to_new_index_mappings = OrderedDict()
			for old_pattern_idx,(old_pattern_node, corresp_new_pattern)\
				in enumerate(zip(current_level_nodes,
								 old_to_new_pattern_mapping)):
				#if the old pattern was NOT changed in this iteration
				if old_pattern_node.pattern == corresp_new_pattern:
					new_idx = patterns.index(corresp_new_pattern) 
					old_to_new_index_mappings[old_pattern_idx] = new_idx

			new_pairwise_aurocs = -np.ones((len(patterns), len(patterns)))
			new_pairwise_sims = np.zeros((len(patterns), len(patterns)))
			for old_idx_i, new_idx_i in\
				old_to_new_index_mappings.items():
				for old_idx_j, new_idx_j in\
					old_to_new_index_mappings.items():
					new_pairwise_aurocs[new_idx_i, new_idx_j] =\
						pairwise_aurocs[old_idx_i, old_idx_j]
					new_pairwise_sims[new_idx_i, new_idx_j] =\
						pairwise_sims[old_idx_i, old_idx_j]
			pairwise_aurocs = new_pairwise_aurocs 
			pairwise_sims = new_pairwise_sims
				 

			#update the hierarchy
			#the current 'top level' will consist of all the current
			# nodes that didn't get a new parent, plus any new parents
			# created                
			next_level_nodes = []
			for frontier_pattern in patterns:
				#either this pattern is in old_pattern_nodes, in which
				# case take the old_pattern_node entry, or it's a completely
				# new pattern in which case make a node for it
				old_pattern_node_found = False
				for old_pattern_node in current_level_nodes:
					if (old_pattern_node.pattern==frontier_pattern):
						#sanity check..there should be only one node
						# per pattern
						next_level_nodes.append(old_pattern_node)
						old_pattern_node_found = True 
				if (old_pattern_node_found==False):
				   next_level_nodes.append(
					PatternMergeHierarchyNode(frontier_pattern)) 

			for next_level_node in next_level_nodes:
				#iterate over all the old patterns and their new parent
				# in order to set up the child nodes correctly
				for old_pattern_idx,(old_pattern_node, corresp_new_pattern)\
					in enumerate(zip(current_level_nodes,
									 old_to_new_pattern_mapping)):

					#if the node has a new parent
					if (old_pattern_node.pattern != corresp_new_pattern):
						if (next_level_node.pattern==corresp_new_pattern):

							
							#corresp_new_pattern should be comprised of a 
							# merging of all the old patterns at
							# indices_merged_with
							indices_merged = tuple(sorted(
								merge_partners_so_far[old_pattern_idx])) 
							#get the relevant slice         
							submat_crosscontam =\
							 cross_contamination[indices_merged,:][:,
												 indices_merged]
							submat_alignersim =\
							 patterns_to_patterns_aligner_sim[
								indices_merged, :][:,indices_merged]

							if (next_level_node.indices_merged is not None):
								assert (next_level_node.indices_merged
										==indices_merged),\
								 (next_level_node.indices_merged,
								  indices_merged)
							else:
								next_level_node.indices_merged =\
									indices_merged
								next_level_node.submat_crosscontam =\
									submat_crosscontam
								next_level_node.submat_alignersim =\
									submat_alignersim

							next_level_node.child_nodes.append(
											old_pattern_node) 

							old_pattern_node.parent_node = next_level_node

			current_level_nodes = next_level_nodes

	return patterns, PatternMergeHierarchy(root_nodes=current_level_nodes)


def _as_float32(x):
    return np.asarray(x, dtype=np.float32)

def write_seqlets_old(pg, pat, sample_name):
    #Create seqlets/ with at least n_seqlets (len-1 int array)
    seqlets = getattr(pat, "seqlets", None)
    if seqlets is None:
        # still create the group with n_seqlets=[0]
        sg = pg.create_group("seqlets")
        sg.create_dataset("n_seqlets", data=np.array([0], dtype=np.int64))
        return
    sg = pg.create_group("seqlets")
    n = len(seqlets)
    # Try to collect common fields if they exist
    def collect(field, cast=None):
        vals = []
        for s in seqlets:
            v = getattr(s, field, None)
            if v is None and isinstance(s, dict):
                v = s.get(field, None)
            if v is not None:
                vals.append(v)
        if vals:
            arr = np.asarray(vals, dtype=cast) if cast is not None else np.asarray(vals)
            sg.create_dataset(field, data=arr, compression="gzip")
    def collect_str(field):
        vals = []
        for s in seqlets:
            v = getattr(s, field, None)
            if v is None and isinstance(s, dict):
                v = s.get(field, None)
            if v is not None:
                vals.append(str(v))
        if vals:
            sg.create_dataset(field, data=np.asarray(vals, dtype=str_dt), compression="gzip")        
    collect("start", np.int64)
    collect("end", np.int64)
    collect("example_idx", np.int64)
    #collect_str("sample_name")
    collect("old_peak_idx", np.int64)
    collect("sample_name", str_dt)
    collect("seqlet_unique_index", np.int64)
    # store booleans as 0/1 small ints
    rc_vals = []
    for s in seqlets:
        v = getattr(s, "is_revcomp", None)
        if v is None and isinstance(s, dict):
            v = s.get("is_revcomp", None)
        if v is not None:
            rc_vals.append(1 if bool(v) else 0)
    if rc_vals:
        sg.create_dataset("is_revcomp", data=np.asarray(rc_vals, dtype=np.int8), compression="gzip")
    sg.create_dataset("n_seqlets", data=np.array([n], dtype=np.int64))

def write_seqlets(pg, pat, sample_name, str_dt=None):
    if str_dt is None:
        str_dt = np.dtype("S64")
    seqlets_all = getattr(pat, "seqlets", None)
    # always create the group
    sg = pg.create_group("seqlets")
    # handle no seqlets on the pattern at all
    if seqlets_all is None:
        sg.create_dataset("n_seqlets", data=np.array([0], dtype=np.int64))
        return 0  # wrote zero
    # optional filtering by sample name
    if sample_name:
        seqlets = []
        for s in seqlets_all:
            sn = getattr(s, "sample_name", None)
            if sn is None and isinstance(s, dict):
                sn = s.get("sample_name", None)
            if sn is not None and str(sn) == str(sample_name):
                seqlets.append(s)
    else:
        seqlets = list(seqlets_all)
    n = len(seqlets)
    # helper collectors
    def collect(field, cast=None, override=None):
        vals = []
        for s in seqlets:
            if override is not None:
                v = override(s)
            else:
                v = getattr(s, field, None)
                if v is None and isinstance(s, dict):
                    v = s.get(field, None)
            if v is not None:
                vals.append(v)
        if vals:
            arr = np.asarray(vals, dtype=cast) if cast is not None else np.asarray(vals)
            sg.create_dataset(field, data=arr, compression="gzip")
    def collect_str(field):
        vals = []
        for s in seqlets:
            v = getattr(s, field, None)
            if v is None and isinstance(s, dict):
                v = s.get(field, None)
            if v is not None:
                vals.append(str(v))
        if vals:
            sg.create_dataset(field, data=np.asarray(vals, dtype=str_dt), compression="gzip")
    # write fields
    collect("start", np.int64)
    collect("end", np.int64)
    # if filtering by sample_name, override example_idx with old_peak_idx
    if sample_name:
        collect("example_idx", np.int64, override=lambda s: getattr(s, "old_peak_idx", None) if not isinstance(s, dict) else s.get("old_peak_idx"))
    else:
        collect("example_idx", np.int64)
    collect("old_peak_idx", np.int64)
    collect_str("sample_name")
    collect("seqlet_unique_index", np.int64)
    # store booleans as tiny ints
    rc_vals = []
    for s in seqlets:
        v = getattr(s, "is_revcomp", None)
        if v is None and isinstance(s, dict):
            v = s.get("is_revcomp", None)
        if v is not None:
            rc_vals.append(1 if bool(v) else 0)
    if rc_vals:
        sg.create_dataset("is_revcomp", data=np.asarray(rc_vals, dtype=np.int8), compression="gzip")
    # always record count
    sg.create_dataset("n_seqlets", data=np.array([n], dtype=np.int64))
    # add a hint for traceability
    if sample_name:
        sg.attrs["filtered_by_sample_name"] = str(sample_name)
    return n


def save_patterns_modisco_h5(merged_patterns, out_h5, sign, sample_name):
    if "neg" in sign:
        sign_short = "neg"
    else:
        sign_short = "pos"

	# ensure iterable
    if not hasattr(merged_patterns, "__iter__") or isinstance(merged_patterns, (str, bytes, dict)):
        merged_patterns = [merged_patterns]
    with h5.File(out_h5, "w") as f:
        root = f.create_group(f"{sign_short}_patterns")
        n_written = 0
        for i, pat in enumerate(merged_patterns):
            pg_name = f"pattern_{i}"
            pg = root.create_group(pg_name)
            # pattern-level arrays
            if getattr(pat, "sequence", None) is not None:
                pg.create_dataset("sequence", data=_as_float32(pat.sequence), compression="gzip")
            if getattr(pat, "contrib_scores", None) is not None:
                pg.create_dataset("contrib_scores", data=_as_float32(pat.contrib_scores), compression="gzip")
            if getattr(pat, "hypothetical_contribs", None) is not None:
                pg.create_dataset("hypothetical_contribs", data=_as_float32(pat.hypothetical_contribs), compression="gzip")
            # Optional name
            name = getattr(pat, "name", None)
            if name is not None:
                pg.attrs["name"] = str(name)
            # Write seqlets (may be filtered by sample_name); returns how many written
            n_seq = write_seqlets(pg, pat, sample_name=sample_name)
            # Skip / remove empty patterns
            if n_seq == 0:
                del root[pg_name]   # remove the group we just created
                continue
            n_written += 1
        root.attrs["n_patterns"] = n_written
    return n_written

import numpy as np

def cwm_to_ppm(cwm, method="relu", alpha=1.0, pseudocount=1e-6, temperature=1.0,
                    return_shape="Lx4"):
    """
    Convert a CWM to a PPM with columns (positions) summing to 1.
    Accepts either (4, L) or (L, 4) input and auto-orients internally.

    return_shape: "4xL" (default) or "Lx4" for MEME-style PSPM rows.
    """
    x = np.asarray(cwm, dtype=float)

    # Auto-orient so rows are bases (A,C,G,T) and cols are positions
    if x.shape[0] == 4:
        x4L = x
    elif x.shape[1] == 4:
        x4L = x.T
    else:
        raise ValueError(f"Need one dimension to be 4 (A,C,G,T). Got {x.shape}.")

    # Turn contributions into non-negative weights
    if method == "relu":
        w = np.maximum(x4L, 0.0)**alpha + float(pseudocount)
    elif method == "abs":
        w = np.abs(x4L)**alpha + float(pseudocount)
    elif method == "softmax":
        m = x4L.max(axis=0, keepdims=True)
        e = np.exp((x4L - m) / float(temperature))
        denom = np.clip(e.sum(axis=0, keepdims=True), 1e-300, None)
        ppm = e / denom
    else:
        raise ValueError("method must be one of {'relu','abs','softmax'}")

    if method in {"relu","abs"}:
        denom = np.clip(w.sum(axis=0, keepdims=True), 1e-300, None)
        ppm = w / denom

    # Return in requested shape
    if return_shape == "4xL":
        return ppm
    elif return_shape == "Lx4":
        return ppm.T
    else:
        raise ValueError("return_shape must be '4xL' or 'Lx4'")


# motifs_db=meme_motif_db, tomtom_dir_temp=tomtom_dir,
def fetch_tomtom_matches_ID(ppm, cwm, motifs_db, tomtom_dir_temp, pattern_name,
        background=[0.25, 0.25, 0.25, 0.25], tomtom_exec_path='tomtom',
        trim_threshold=0.3, trim_min_length=3):

        """Fetches top matches from a motifs database using TomTom.
        Args:
                ppm: position probability matrix- numpy matrix of dimension (N,4)
                background: list with ACGT background probabilities
                tomtom_exec_path: path to TomTom executable
                motifs_db: path to motifs database in meme format
                n: number of top matches to return, ordered by p-value
                temp_dir: directory for storing temp files
                trim_threshold: the ppm is trimmed from left till first position for which
                        probability for any base pair >= trim_threshold. Similarly from right.
        Returns:
                list: a list of up to n results returned by tomtom, each entry is a
                        dictionary with keys 'Target ID', 'p-value', 'E-value', 'q-value'
        """

        #_, fname = tempfile.mkstemp(dir=tomtom_dir_temp)
        #_, tomtom_fname = tempfile.mkstemp(dir=tomtom_dir_temp)
        fname = f"{tomtom_dir_temp}/{pattern_name}_cwm.txt"
        tomtom_fname = f"{tomtom_dir_temp}/{pattern_name}_tomtom.csv"

        score = np.sum(np.abs(cwm), axis=1)
        trim_thresh = np.max(score) * trim_threshold  # Cut off anything less than 30% of max score
        pass_inds = np.where(score >= trim_thresh)[0]
        trimmed = ppm[np.min(pass_inds): np.max(pass_inds) + 1]

        # can be None of no base has prob>t
        if trimmed is None:
                return []

        # convert cwm to ppm
        trimmed = cwm_to_ppm(cwm[np.min(pass_inds): np.max(pass_inds) + 1])

        # trim and prepare meme file
        write_meme_file(trimmed, background, fname)

        # run tomtom
        #cmd = '%s -no-ssc -oc . --verbosity 1 -text -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 %s %s > %s' % (tomtom_exec_path, fname, motifs_db, tomtom_fname)
        cmd = '%s -no-ssc -oc . --verbosity 1 -text -min-overlap 5 -mi 1 -dist kullback -evalue -thresh 10.0 %s %s > %s' % (tomtom_exec_path, fname, motifs_db, tomtom_fname)

        os.system(cmd)
        tomtom_results = pandas.read_csv(tomtom_fname, sep="\t", usecols=(1, 5))
        #os.system('rm ' + tomtom_fname)
        #os.system('rm ' + fname)
        return tomtom_results

def run_tomtom_ID(modisco_h5py, meme_motif_db, tomtom_dir, top_n_matches=3,
        tomtom_exec="tomtom", trim_threshold=0.3, trim_min_length=3):
        modisco_results = h5.File(modisco_h5py, 'r')

        tomtom_results = {'pattern': [], 'num_seqlets': []}
        for i in range(top_n_matches):
                tomtom_results['match{}'.format(i)] = []
                tomtom_results['qval{}'.format(i)] = []

        for name in ['pos_patterns', 'neg_patterns']:
                if name not in modisco_results.keys():
                        continue

                metacluster = modisco_results[name]
                key = lambda x: int(x[0].split("_")[-1])
                for pattern_name, pattern in sorted(metacluster.items(), key=key):
                        ppm = np.array(pattern['sequence'][:])
                        cwm = np.array(pattern["contrib_scores"][:])

                        num_seqlets = pattern['seqlets']['n_seqlets'][:][0]
                        tag = '{}.{}'.format(name, pattern_name)

                        r = fetch_tomtom_matches_ID(ppm, cwm, motifs_db=meme_motif_db, tomtom_dir_temp=tomtom_dir,pattern_name=tag,
                                tomtom_exec_path=tomtom_exec, trim_threshold=trim_threshold,
                                trim_min_length=trim_min_length)


        modisco_results.close()
        return pandas.DataFrame(tomtom_results)



###################
# Testing motif matches instead of tomtom (Can ignore for now)

import numpy as np
import pandas as pd
from math import sqrt
from functools import lru_cache
import math

# ===================== basic helpers =====================

ALPH = np.array(list("ACGT"))

def to_4xL(M):
    M = np.asarray(M)
    if M.ndim != 2 or 4 not in M.shape:
        raise ValueError(f"Expected 2D with a 4, got {M.shape}")
    return M if M.shape[0] == 4 else M.T

def revcomp(M):  # A,C,G,T -> T,G,C,A and reverse columns
    return M[::-1, ::-1]

#def colweights_from_cwm(C, mode="abs", eps=1e-12):
#    w = np.sum(np.abs(C), axis=0) if mode == "abs" else np.sum(C**2, axis=0)
#    s = w.sum()
#    return w / (s + eps) if s > 0 else np.full(C.shape[1], 1.0 / C.shape[1])

def colweights_from_cwm(C, mode="abs", eps=1e-12):
    """
    Compute per-column weights from a CWM (4×L) for use in weighted similarity.
    Parameters
    ----------
    C : array-like, shape (4, L) or (L, 4)
        Contribution Weight Matrix (signed importance scores).
    mode : str
        How to score each column before normalizing to sum=1.
        - "abs"      : sum of absolute contributions per column (default).
        - "sq"       : sum of squared contributions per column (emphasizes large magnitudes).
        - "pos_patterns"      : sum of positive contributions only.
        - "neg_patterns"      : sum of negative contributions (as positive magnitudes).
    eps : float
        Small constant to avoid division by zero during normalization.
    Returns
    -------
    w : np.ndarray, shape (L,)
        Non-negative column weights normalized to sum to 1. If all columns score 0,
        returns a uniform distribution over columns.
    """
    # Split into positive and negative (as positive magnitudes)
    P = np.maximum(C, 0.0)   # positive contributions
    N = np.maximum(-C, 0.0)  # negative contributions
    if mode == "abs":
        w = np.sum(np.abs(C), axis=0)
    elif mode == "sq":
        w = np.sum(C*C, axis=0)
    elif mode == "pos_patterns":
        w = np.sum(P, axis=0)
    elif mode == "neg_patterns":
        w = np.sum(N, axis=0)
    else:
        raise ValueError(f"Unknown mode: {mode}")
    s = float(w.sum())
    return w/(s + eps) if s > 0 else np.full(C.shape[1], 1.0 / C.shape[1])

def pwm_to_signed(P, mode="logodds"):
    if mode == "center":   return P - 0.25
    if mode == "logodds":  return np.log2((P + 1e-9)/0.25)
    return P  # assume already signed

def consensus_from_cwm(C):
    idx = np.argmax(np.abs(C), axis=0)
    return "".join(ALPH[idx])

def consensus_from_ppm(P):
    idx = np.argmax(P, axis=0)
    return "".join(ALPH[idx])

def count_valid_offsets(L1, L2, min_overlap=6):
    # exact count of ungapped slides whose overlap >= min_overlap
    K = 0
    #for off in range(-(L2-1), L1):
    for off in range(-(math.ceil( (L2/1.5))-1), L1):
        s1, s2 = max(0, off), max(0, -off)
        ov = min(L1 - s1, L2 - s2)
        if ov >= min_overlap:
        #if ov >= max(min_overlap,min((L2/1.5),(L1/1.5))):
            K += 1
    return K

# ===================== alignment (fast) =====================

def _weighted_cosine_symmetric(A, B, w, eps=1e-12):
    ws = np.sqrt(w)[None, :]
    vA = (A * ws).ravel(order="F")
    vB = (B * ws).ravel(order="F")
    nA = np.linalg.norm(vA); nB = np.linalg.norm(vB)
    return float(vA @ vB) / (nA*nB + eps)

def best_alignment_fast(A, B, w, min_overlap=6):
    """
    A,B: (4,L); w: (L,) normalized column-weights (sum=1)
    Returns dict: score_dot, score_cos, offset, strand, overlap
    """
    L1, L2 = A.shape[1], B.shape[1]
    best = {"score_dot": -np.inf, "score_cos": -2.0, "offset": None, "strand": "+", "overlap": 0}
    for strand, M2 in [("+", B), ("-", revcomp(B))]:
        #for off in range(-L2+1, L1):
        for off in range(-(math.ceil( (L2/1.5))-1), L1):
            s1, s2 = max(0, off), max(0, -off)
            ov = min(L1 - s1, L2 - s2)
            if ov < min_overlap:
            #if ov < max(min_overlap,min((L2/1.5),(L1/1.5))):
                continue
            Awin, Bwin = A[:, s1:s1+ov], M2[:, s2:s2+ov]
            #wwin = w[s1:s1+ov]; wwin = wwin / (wwin.sum() + 1e-12)
            wwin = w[s1:s1+ov]
            # dot (weights apply once; A is CWM-signed)
            dot = float(np.sum(Awin * Bwin * wwin[None, :]))
            if dot > best["score_dot"]:
                cos = _weighted_cosine_symmetric(Awin, Bwin, wwin)
                best.update(score_dot=dot, score_cos=cos, offset=off, strand=strand, overlap=ov)
    return best

# ===================== Tomtom-like p for weighted_dot =====================

def _sample_dirichlet_columns(L, bg=(0.25,0.25,0.25,0.25), alpha=20.0, rng=None):
    rng = np.random.default_rng() if rng is None else rng
    a = np.asarray(bg, float) * float(alpha)
    return rng.dirichlet(a, size=L).T  # (4,L) probabilities

def _build_single_col_pmf(A, w, *, samples=200_000, bg=(0.25,)*4, alpha=20.0,
                          mode_signed="logodds", bin_width=1e-3, clip=0.5, seed=None):
    """
    PMF of single-column contribution: S = w_i * <A[:,i], B_rand_signed>, with i ~ w.
    Returns (edges, pmf1, delta, clip)
    """
    rng = np.random.default_rng(seed)
    A, w = to_4xL(A), np.asarray(w, float)
    w = w / (w.sum() + 1e-12)
    L1 = A.shape[1]
    idx = rng.choice(L1, size=samples, p=w)
    Brand = _sample_dirichlet_columns(samples, bg=bg, alpha=alpha, rng=rng)
    if mode_signed == "center":
        Brand = Brand - 0.25
    elif mode_signed == "logodds":
        Brand = np.log2((Brand + 1e-9)/0.25)
    dots = np.sum(A[:, idx] * Brand, axis=0) * w[idx]
    lo, hi = -clip, clip
    nbins = int(np.ceil((hi - lo) / bin_width)) + 1
    edges = np.linspace(lo, hi, nbins)
    hist, _ = np.histogram(np.clip(dots, lo, hi), bins=edges, density=True)
    pmf = hist / (hist.sum() + 1e-12)
    return edges, pmf, float(edges[1]-edges[0]), clip

def _pmf_convolve_power_fft(pmf, L):
    """
    L-fold convolution of a discrete PMF using FFT.
    Returns the FULL support of length (m-1)*L + 1 with NO trimming.
    """
    m = len(pmf)                 # single-column pmf length (bins)
    if L <= 1:
        return pmf.copy()
    need = (m - 1) * L + 1       # exact length of L-fold convolution
    # FFT size: next power of two >= need
    n = 1
    while n < need:
        n <<= 1
    spec = np.fft.rfft(pmf, n)
    out  = np.fft.irfft(spec ** L, n)[:need]   # keep EXACT needed length
    out[out < 0] = 0
    s = out.sum()
    if s > 0:
        out /= s
    return out

@lru_cache(maxsize=None)
def _cdf_for_L_cached(pmf_bytes, L):
    pmf = np.frombuffer(pmf_bytes, dtype=np.float64)
    pmfL = _pmf_convolve_power_fft(pmf, L)
    cdf = np.minimum(np.cumsum(pmfL), 1.0)
    start_bin = 0
    return pmfL, cdf, start_bin

def tomtom_like_p_for_target(A, w, B_prob, *,
                             mode_signed="logodds", min_overlap=6,
                             bg=(0.25,)*4, alpha=20.0,
                             pmf_pack=None,
                             seed=None):
    A = to_4xL(A); w = np.asarray(w, float); w = w / (w.sum() + 1e-12)
    A_w = A
    Bp = to_4xL(np.asarray(B_prob))
    if mode_signed == "center":
        Bf = Bp - 0.25
    elif mode_signed == "logodds":
        Bf = np.log2((Bp + 1e-9)/0.25)
    else:
        Bf = Bp
    best_fwd = best_alignment_fast(A_w, Bf, w, min_overlap=min_overlap)
    best_rc  = best_alignment_fast(A_w, revcomp(Bf), w, min_overlap=min_overlap)
    edges, pmf1, delta, clip = pmf_pack
    pmf1 = pmf1.astype(np.float64, copy=False)
    def p_for(best):
        L = int(best["overlap"])
        if L < min_overlap:
        #if L < max(min_overlap,min((L2/1.5),(L1/1.5))):
            return 1.0
        # one-sided (upper-tail): if the observed dot is <= 0, it's not significant
        if best["score_dot"] <= 0:
            return 1.0
        # UNPACK THREE VALUES HERE:
        pmfL, cdfL, start_bin = _cdf_for_L_cached(pmf1.astype(np.float64).tobytes(), L)
        # map score -> index on sum support; note the start_bin correction
        idx_raw = int(np.floor((best["score_dot"] + L*clip) / delta + 1e-9))
        idx = int(np.clip(idx_raw - start_bin, 0, len(cdfL) - 1))
        Fsum = float(cdfL[idx])       # CDF at the sum
        p_sum = 1.0 - Fsum            # upper tail
        K = count_valid_offsets(A.shape[1], Bp.shape[1], min_overlap=min_overlap)
        p_max = 1.0 - (1.0 - p_sum)**K   # Šidák for max over shifts on this strand
        return float(np.clip(p_max, 0.0, 1.0))
    p_fwd = p_for(best_fwd)
    p_rc  = p_for(best_rc)
    p_both = 1.0 - (1.0 - min(p_fwd, p_rc))**2
    best = best_fwd if best_fwd["score_dot"] >= best_rc["score_dot"] else best_rc
    strand = "+" if best is best_fwd else "-"
    return {
        "p_dot_tomtom": p_both,
        "best": {**best, "strand": strand}
    }

# ===================== main driver =====================

def run_cwm_vs_targets_tomtom_like(
    cwm, targets, *,
    pwm_mode="logodds",
    min_overlap=6, wmode="abs",
    # Tomtom-like null options:
    bg=(0.25,0.25,0.25,0.25), alpha=20.0,
    pmf_samples=200_000, bin_width=1e-3, clip=0.5, seed=1
):
    """
    targets: dict {motif_id: PWM_prob_matrix (4xL or Lx4)}
    Returns pandas DataFrame with required columns including p_dot_tomtom and qvalue (BH-FDR).
    """
    A = to_4xL(cwm)
    w = colweights_from_cwm(A, mode=wmode)
    # Build single-column PMF once for the query (Tomtom-like null)
    pmf_pack = _build_single_col_pmf(
        A, w, samples=pmf_samples, bg=bg, alpha=alpha,
        mode_signed=pwm_mode, bin_width=bin_width, clip=clip, seed=seed
    )  # (edges, pmf1, delta, clip)
    rows = []
    for motif_id, P in targets.items():
        Pprob = to_4xL(np.asarray(P))
        # observed best alignment (scores + offset/strand/overlap)
        Bsigned = pwm_to_signed(Pprob, pwm_mode)
        best = best_alignment_fast(A, Bsigned, w, min_overlap=min_overlap)
        # consensus over the aligned window, oriented to reported strand
        s1, s2 = max(0, best["offset"]), max(0, -best["offset"])
        ov = best["overlap"]
        Awin = A[:, s1:s1+ov]
        Poriented_prob = Pprob if best["strand"] == "+" else revcomp(Pprob)
        Pwin_prob = Poriented_prob[:, s2:s2+ov]
        cons_query  = consensus_from_cwm(Awin)
        cons_target = consensus_from_ppm(Pwin_prob)
        # Tomtom-like p for weighted_dot (max over shifts & both strands)
        tt = tomtom_like_p_for_target(
            A, w, Pprob, mode_signed=pwm_mode, min_overlap=min_overlap,
            bg=bg, alpha=alpha, pmf_pack=pmf_pack, seed=None
        )
        rows.append({
            "motif_id": motif_id,
            "offset": best["offset"],
            "strand": best["strand"],
            "overlap": best["overlap"],
            "weighted_dot": best["score_dot"],
            "weighted_cosine": best["score_cos"],
            "consensus_query": cons_query,
            "consensus_target": cons_target,
            "p_dot_tomtom": tt["p_dot_tomtom"]
        })
    df = pd.DataFrame(rows).sort_values(
        by=["p_dot_tomtom", "weighted_dot", "weighted_cosine"],
        ascending=[True, False, False], ignore_index=True
    )
    # Benjamini-Hochberg q-value over p_dot_tomtom
    p = df["p_dot_tomtom"].values
    m = len(p)
    order = np.argsort(p)
    ranks = np.empty_like(order); ranks[order] = np.arange(1, m+1)
    q = p * m / ranks
    q = np.minimum.accumulate(q[order[::-1]])[::-1]
    df["qvalue"] = np.clip(q, 0, 1)
    return df


##################
# html report
import os, base64, mimetypes

def image_to_data_uri(path):
    # infer MIME type from extension; default to png if unknown
    mime, _ = mimetypes.guess_type(path)
    if mime is None:
        mime = "image/png"
    with open(path, "rb") as f:
        b64 = base64.b64encode(f.read()).decode("ascii")
    return f"data:{mime};base64,{b64}"

def path_to_image_html(path, max_width_px=160):
    """
    Formatter for DataFrame.to_html that embeds the image file into the HTML.
    If the file is missing, returns the raw path as text (so you can see what broke).
    """
    if isinstance(path, str) and os.path.isfile(path):
        src = image_to_data_uri(path)
        return f'<img src="{src}" style="max-width:{max_width_px}px; height:auto;"/>'
    if isinstance(path, str) and path.startswith("data:"):
        return f'<img src="{path}" style="max-width:{max_width_px}px; height:auto;"/>'
    return f'<span style="color:#b91c1c;">(missing)</span> <code>{path}</code>'