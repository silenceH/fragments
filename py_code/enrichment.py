from Bioisosteres import *

class Pair:
	def __init__(self, frags, twoDim=None, threeDim=None):
		self.frags=frags
		if twoDim is None:
			self.twoDim=twoDim
		self.twoDim=twoDim
		if threeDim is None:
			self.threeDim=threeDim
		self.threeDim=threeDim

		self.active = 0

	def set_score(self,score,dim):
		if dim == 2:
			self.twoDim = score
		if dim == 3: 
			self.threeDim = score
		if dim != 2 or dim != 3:
			print str(dim) + " not a valid dimension"

## define return list of frags that have the query as a bioisosteric pair
def pairs_with_query(list_of_pairs, query):
	## test for group objects
	try: 
		## take pairs from Group data type 
		list_of_pairs = [x.group for x in list_of_pairs]
	except AttributeError:
		list_of_pairs = list_of_pairs
	## get the first pair if the second is equal to the query
	final_list_a = [x[0] for x in list_of_pairs if x[1].are_similar(query,1.0)]	
	final_list_b = [x[1] for x in list_of_pairs if x[0].are_similar(query,1.0)]	

	## combine two lists
	final_list = final_list_a + final_list_b
	return final_list

def score_similarity(query,list_of_fragments,threeDim=False):
	## returna  list of all the fragments in a set of overlays
	similarity_scores = [x.tanimoto_score(query) for x in list_of_fragments]
	if threeDim:
		threeDim_scores = [get_3D_score(query.frag,x.frag) for x in list_of_fragments]
		return similarity_scores,threeDim_scores
	return similarity_scores

def get_3D_score(query, ref):
	#write to tmp sdf file
	f = Chem.SDWriter("temp_ref.sdf")
	f.write(query)
	f.flush()

	g = Chem.SDWriter("temp_dat.sdf")
	g.write(ref)
	g.flush()

	# get score
	os.system("shape-it -r temp_ref.sdf -d temp_dat.sdf -s output.dat")
	f = open('output.dat','r')
	results = f.readlines()[1]
	first_tab = results.find('\t')
	second_tab = results.find('\t',first_tab+1)
	third_tab = results.find('\t',second_tab+1)
	score = float(results[second_tab:third_tab])

	return score


def rank_list_of_fragments(query,list_of_fragments,return_scores = False):
	## rank a list of fragments by the similarity to a query
	scores = score_similarity(query,list_of_fragments)
	changed = True
	while changed:
		changed = False
		for i in range(len(list_of_fragments)-1):
			if scores[i] < scores[i+1]:
				scores[i], scores[i+1] = scores[i+1], scores[i]
				list_of_fragments[i], list_of_fragments[i+1] = list_of_fragments[i+1], list_of_fragments[i]
				changed = True
	if return_scores:
		return list_of_fragments,scores
	return list_of_fragments

def score_list_frag_pairs(query, list_of_fragments,threeDim=False):
	## create list of objects
	pair_list = []
	if threeDim:
		scores,threeDim_scores = score_similarity(query,list_of_fragments,True)
	else:
		scores = score_similarity(query,list_of_fragments)
	for i in range(len(list_of_fragments)):
		if threeDim:
			new_pair = Pair((query,list_of_fragments[i]),twoDim=scores[i],threeDim=threeDim_scores[i])
		else:
			new_pair = Pair((query,list_of_fragments[i]),twoDim=scores[i])
		pair_list.append(new_pair)
	return pair_list

def mergeSorted(scores_1,fragments_1, scores_2, fragments_2):
	final_score= []
	final_frag= []

	assert len(scores_2) == len(fragments_2)
	assert len(scores_1) == len(fragments_1)

	# check to see whether either scores is empty
	if len(scores_1) == 0:
		return scores_2,fragments_2
	if len(scores_2) == 0:
		return scores_1,fragments_1

	score_1 = scores_1.pop(0)	# gets the biggest element
	score_2 = scores_2.pop(0)

	frag_pair_1 = fragments_1.pop(0)
	frag_pair_2 = fragments_2.pop(0)
	
	assert len(scores_2) == len(fragments_2)
	assert len(scores_1) == len(fragments_1)
	## note:: this needs fixing for scoress of one element 
	while len(scores_1)!=0 and len(scores_2)!=0 :
		if score_1>=score_2:
			final_score.append(score_1)
			final_frag.append(frag_pair_1)
			score_1 = scores_1.pop(0)	
			frag_pair_1 = fragments_1.pop(0)

		else:
			final_score.append(score_2)
			final_frag.append(frag_pair_2)
			score_2 = scores_2.pop(0)	
			frag_pair_2 = fragments_2.pop(0)

	if len(scores_1)==0 or len(scores_2)==0:

		if len(scores_1) == 0:
			empty_scores = scores_1
			empty_pairs = fragments_1
			last_score = score_1
			last_pair = frag_pair_1

			non_empty_scores = scores_2
			non_empty_pairs = fragments_2
			non_last_score = score_2
			non_last_pair = frag_pair_2

		if len(scores_2) == 0:
			empty_scores = scores_2
			empty_pairs = fragments_2
			last_score = score_2
			last_pair = frag_pair_2

			non_empty_scores = scores_1
			non_empty_pairs = fragments_1
			non_last_score = score_1
			non_last_pair = frag_pair_1
		
		finished = False

		while not finished:
			if last_score >= non_last_score:
				final_score.extend([last_score,non_last_score])
				final_frag.extend([last_pair,non_last_pair])
				final_score.extend(non_empty_scores)
				final_frag.extend(non_empty_pairs)
				finished = True

			else: 
				final_score.append(non_last_score)
				final_frag.append(non_last_pair)
				non_last_score = non_empty_scores.pop(0)	
				non_last_pair = non_empty_pairs.pop(0)

			if len(non_empty_scores) == 0:
				if last_score >= non_last_score:
					final_score.extend([last_score,non_last_score])
					final_frag.extend([last_pair,non_last_pair])
					finished = True
				else:
					final_score.extend([non_last_score,last_score])
					final_frag.extend([non_last_pair,last_pair])
					finished = True

	return final_frag, final_score

def rank_all_frags(list_of_fragments,pair_obj=False):
	if pair_obj:
		final_list = []
		for frag_1 in list_of_fragments:
			query_list = [x for x in list_of_fragments if x is not frag_1]
			list_of_pairs = score_list_frag_pairs(frag_1,query_list)
			final_list.extend(list_of_pairs)

		final_list.sort(key=lambda x: x.twoDim, reverse=True)
		return final_list

	final_list = []
	final_scores = []


	for frag_1 in list_of_fragments:
		query_list = [x for x in list_of_fragments if x is not frag_1]
		ranked_list, scores = rank_list_of_fragments(frag_1,query_list,True)
		ranked_list_pairs = [(frag_1,x) for x in ranked_list]
		if len(final_list) == 0:
			final_list = ranked_list_pairs
			final_scores = scores
		else:
			final_list,final_scores = mergeSorted(final_scores,final_list,scores,ranked_list_pairs)
	return final_list,final_scores

def rank_all_frags_threeDim(list_of_fragments):
	final_list = []
	for frag_1 in list_of_fragments:
		query_list = [x for x in list_of_fragments if x is not frag_1]
		list_of_pairs = score_list_frag_pairs(frag_1,query_list,True)
		final_list.extend(list_of_pairs)

	return final_list

		

def ranked_list_of_all_fragments(list_of_fragments):
	## define final list of fragments and scores
	final_list = []
	final_scores = []

	# iterate for fragments to get a ranked list of pairs for each fragment
	for frag_1 in list_of_fragments:
		query_list = [x for x in list_of_fragments if x is not frag_1]
		ranked_list, scores = rank_list_of_fragments(frag_1,query_list,True)
		if len(final_list) == 0:
			final_list = [(frag_1,frag_2) for frag_2 in ranked_list]
			final_scores = scores
		else: 
			while(len(scores)>0):
				test_score = scores.pop(0)
				test_frag = ranked_list.pop(0)
				for i in range(len(final_scores)):
					if test_score >= final_scores[i]:
						final_list.insert(i,(frag_1,test_frag))
						final_scores.insert(i,test_score)
						break
	return final_list,final_scores

	
def rank_targets(*args):
	all_fragments = get_unique_fragments_from_files(*args)
	return ranked_list_of_all_fragments(all_fragments)

def rank_targets_test(threeDim,pair_obj,*args):
	all_fragments = get_unique_fragments_from_files(*args)
	if pair_obj:
		return rank_all_frags(all_fragments,True)
	if threeDim:
		return rank_all_frags_threeDim(all_fragments)
	return rank_all_frags(all_fragments)

def pairs_are_the_same(pair1,pair2):
	return True

def find_rank(fragment,list_of_fragments):
	for i in range(len(list_of_fragments)):
		if fragment.are_similar(list_of_fragments[i],1.0):
			return i
	return None

def pair_enrichment(num,target,*args):
	# give the ranks of the bioisosteric pairs 
	
	# first get the pairs from a 
	target_pairs = get_bioisosteres(target,return_pairs=True)

	# get all the unique fragments from the files
	all_fragments = get_unique_fragments_from_files(*args)

	if num>0:
		# select num pairs 
		import random
		pairs = []
		for n in range(num):
			choice = random.randint(0,len(all_fragments)) - 1
			pairs.append(choice)
	else:
		pairs = range(len(target_pairs))

	for pair in pairs:
		query_pair = target_pairs[pair].group
		## rank by similarity to first fragment in the pair
		## then find the rank of the second fragment
		ranked_list = rank_list_of_fragments(query_pair[0],all_fragments)
		rank_0 = find_rank(query_pair[1],ranked_list)
		print "for pair " + str(pair)
		print "the rank of the pair wrt first fragment: " + str(rank_0)

		## rank by similarity to the second fragment in the pair
		## then find the rank of the second fragment
		ranked_list = rank_list_of_fragments(query_pair[1],all_fragments)
		rank_1 = find_rank(query_pair[0],ranked_list)
		print "the rank of the pair wrt second fragment: " + str(rank_1)

	
def group_enrichment(target,*args):
	## find the enrichment for molecules in the groups of an overlay

	## get all the groups
	target_groups = get_bioisosteres(target)
	num_groups = len(target_groups) 	# number of groups

	## get all the fragments in the test set
	all_fragments = get_unique_fragments_from_files(*args)
	
	total_fragments = len(all_fragments)
	print "total number of fragments: " + str(total_fragments)

	## define statistics
	av_ranks = [] # average ranks
	pair_ranks = [] # pairwise ranks
	grp_size = [g.size() for g in target_groups] # group sizes
	grp_num_pairs = [2*x for x in grp_size]
	
	## find the 2D rank for each group in the target groups
	for grp in target_groups:
		pairs = grp.get_pairs_from_group()
		total_rank = 0
		grp_rank = [] 	# all ranks of the pairs
		for (frag_1,frag_2) in pairs: 
			ranked_list = rank_list_of_fragments(frag_1,all_fragments)
			rank_0 = find_rank(frag_2,ranked_list)
			grp_rank.append(rank_0)
			total_rank += rank_0
			ranked_list = rank_list_of_fragments(frag_2,all_fragments)
			rank_1 = find_rank(frag_1,ranked_list)
			grp_rank.append(rank_1)
			total_rank += rank_1
		av_rank = total_rank / (2*len(pairs))
		av_ranks.append(av_rank)
		pair_ranks.append(grp_rank)
	
	## test to see whether the lists are the same size as the number of groups
	assert num_groups == len(av_ranks)
	
	## write stats to file
	directory = '../test_output/enrichment/'
	try: 
		os.makedirs(directory)
	except OSError:
		print "enrichment dir already exist"
	
	f = open(directory+"group_stats.csv",'w')

	# print summary statistics
	f.write("group ranks by pairs\n")
	for i in xrange(max(grp_num_pairs)):
		for j in xrange(num_groups):
			if i >= grp_num_pairs[j]:
				f.write('\t')
			else:
				f.write(str(pair_ranks[j][i])+'\t')
		f.write('\n')

	for i in xrange(num_groups):
		f.write(str(av_ranks[i]) + "\t")
	f.write("\n\n")
	f.write(str(total_fragments)+ '\n')
	f.write("average group rating\t" + str(sum(av_ranks)/num_groups)+'\n')
	f.close()
	print "statistics written to " + directory
	return
