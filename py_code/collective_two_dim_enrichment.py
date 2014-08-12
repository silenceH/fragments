import enrichment
import Bioisosteres
from rdkit import Chem

# get enrichment over all targets

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
		
	print str(len(all_fragment_pairs)) + " pairs"

	# deduplicate all_fragment_pairs list
	print "start deduplication"
	all_fragment_pairs = enrichment.deduplicate_list(all_fragment_pairs)
	print "finish deduplication"
	
	print str(len(all_fragment_pairs)) + " unique pairs"

	for i in range(len(all_fragment_pairs)):
		test_pair = all_fragment_pairs[i]
		for ranked_pair in all_bioisosteres:
			ranked_pair_obj = enrichment.Pair(ranked_pair)
			if enrichment.pairs_are_the_same(test_pair,ranked_pair_obj):
				test_pair.active = 1
	
	actives = sum([1 for x in all_fragment_pairs if x.active == 1])
	total_pairs = len(all_fragment_pairs)

	print "total actives = " + str(actives)
	print "total inactives = " + str(total_pairs-actives)

	f = open("2D_3D_activity_" + str(len(args)) + ".csv",'w')
	f.write("2Dscore,3Dscore,active,\n")

	for i in range(len(all_fragment_pairs)):
		f.write(str(all_fragment_pairs[i].twoDim)+','+str(all_fragment_pairs[i].threeDim)+','+str(all_fragment_pairs[i].active)+',\n')
	
	directory = '../test_output/weird_pairs/' 

	try: 
		os.makedirs(directory)
		print "created new directory: " + directory
	except OSError:
		print directory + " already exists" 
			
	#weird_stats = open(directory+'weird_stats.csv','w')
	#weird_stats.write('pair_num,2D,3D,\n')
	#count = 0 
	#for pair in all_fragment_pairs: 
	#	if abs(pair.twoDim - pair.threeDim) > 0.6:
	#		count += 1
	#		w = Chem.SDWriter(directory+'pair_'+str(all_fragment_pairs.index(pair))+'.sdf')
	#		w.write(pair.frags[0].frag)
	#		w.write(pair.frags[1].frag)
	#		w.flush()
	#		weird_stats.write(str(all_fragment_pairs.index(pair))+','  + \
	#			str(pair.twoDim)+','+str(pair.threeDim)+',\n')
	#		
	#print "num with difference > 0.6: " + str(count)		
	#weird_stats.close()


	




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


#get_enrichment(targets[:75])
