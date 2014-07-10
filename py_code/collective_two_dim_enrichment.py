import enrichment
import Bioisosteres

# get enrichment over all targets

# TODO: add 3D enrichment

def get_enrichment(args):
	# get a list of all the ranked 2D pairs
	all_fragment_pairs = []
	all_bioisosteres = []

	for target in args: 

		# get ranked list of 3D targets which includes 2D similarity scores
		target_ranks = enrichment.rank_targets_test(True,False, target) # True for pair objects

		# identify bioisosteres from the target
		bioisostere_pairs = [x.group for x in  Bioisosteres.get_bioisosteres(target,test=True,return_pairs=True)]

		# add to overall lists
		all_fragment_pairs.extend(target_ranks)
		all_bioisosteres.extend(bioisostere_pairs)
		
	# sort list of fragment pairs into order
	
	for i in range(len(all_fragment_pairs)):
		test_pair = all_fragment_pairs[i].frags
		for ranked_pair in all_bioisosteres:
			if (test_pair[0].are_similar(ranked_pair[0],1) and test_pair[1].are_similar(ranked_pair[1],1)) or \
					(test_pair[0].are_similar(ranked_pair[0],1) and test_pair[1].are_similar(ranked_pair[1],1)):
						test_pair.active = 1
	
	actives = len(all_bioisosteres)
	total_pairs = len(all_fragment_pairs)

	print "total actives = " + str(actives)
	print "total inactives = " + str(total_pairs-actives)

	f = open("2D_3D_activity_15.csv",'w')
	f.write("2Dscore,3Dscore,active,\n")

	for i in range(len(all_fragment_pairs)):
		f.write(str(all_fragment_pairs[i].twoDim)+','+str(all_fragment_pairs[i].threeDim)+','+str(all_fragment_pairs[i].active)+',\n')



	




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
get_enrichment(targets[:15])


for i in range(len(ranks)):
	stats_file.write(str(scores[i]) + ',' + str(ranks[i]) + ',\n')

#get_enrichment(targets[:75])
