import enrichment
import Bioisosteres

# get enrichment over all targets

# TODO: add 3D enrichment

def get_enrichment(args):
	# get a list of all the ranked 2D pairs
	all_fragments_2D = []
	all_bioisosteres = []

	for target in args: 
		# get ranked list of 2D targets
		target_2D_ranks = enrichment.rank_targets_test(False,True, target) # True for pair objects

		# identify bioisosteres from the target
		bioisostere_pairs = [x.group for x in  Bioisosteres.get_bioisosteres(target,test=True,return_pairs=True)]

		# add to overall lists
		all_fragments_2D.extend(target_2D_ranks)
		all_bioisosteres.extend(bioisostere_pairs)
		
		# include 3D shape comparison
		# threeDim_frags = enrichment.rank_targets_test(True,False,target)

	# sort list of fragment pairs into order
	all_fragments_2D.sort(key=lambda x: x.twoDim, reverse=True)
	
	ranked_list = [0 for x in range(len(all_fragments_2D))]

	for i in range(len(all_fragments_2D)):
		test_pair = all_fragments_2D[i].frags
		for ranked_pair in all_bioisosteres:
			if (test_pair[0].are_similar(ranked_pair[0],1) and test_pair[1].are_similar(ranked_pair[1],1)) or \
					(test_pair[0].are_similar(ranked_pair[0],1) and test_pair[1].are_similar(ranked_pair[1],1)):
				ranked_list[i] = 1
	
	actives = len(all_bioisosteres)
	total_pairs = len(all_fragments_2D)

	print "total actives = " + str(actives)
	print "total inactives = " + str(total_pairs-actives)
	print "enrichment factor threshold" 

	efs = []

	for threshold in [0.01,0.05,0.10]:
		try:
			n = int(round(total_pairs * threshold))
			pick_n = ranked_list[:n]
			t_p = sum(pick_n)
			ef = (float(t_p)/float(n))/(float(actives)/float(total_pairs))
			print str(threshold) +'\t' + str(ef) 
			efs.append(ef)
		except ZeroDivisionError:
			print "zero division"
			print total_pairs * threshold
			efs.append('nan')
	return [x.twoDim for x in all_fragments_2D], ranked_list 

	




import os
import fnmatch
try:
	data = os.environ['DATA']		## get data env
except KeyError:
	print "cannot find data environment variable"
data_dir = data + '/validation_overlays/'
files = []
for file in os.listdir(data_dir):
	if fnmatch.fnmatch(file, '*.sdf'):
		files.append(file)

targets = [target[:-4] for target in files]

stats_file = open('ranked_list.csv','w')
scores,ranks = get_enrichment(targets[:10])


for i in range(len(ranks)):
	stats_file.write(str(scores[i]) + ',' + str(ranks[i]) + ',\n')

