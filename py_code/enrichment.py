from Bioisosteres import *

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

def score_similarity(query,list_of_fragments):
	## returna  list of all the fragments in a set of overlays
	similarity_scores = [x.tanimoto_score(query) for x in list_of_fragments]
	return similarity_scores

def rank_list_of_fragments(query,list_of_fragments):
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
	return list_of_fragments

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
	
	print "total number of fragments: " + str(len(all_fragments))

	## define statistics
	av_ranks = [] # average ranks
	pair_ranks = [] # pairwise ranks
	grp_size = [g.size() for g in target_groups] # group sizes
	grp_num_pairs = [2*x for x in grp_size]
	
	## find the enrichment for each group in the target groups
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
	f.write("average group rating\t" + str(sum(av_ranks)/num_groups)+'\n')
	f.close()
	print "statistics written to " + directory
	print grp_size
