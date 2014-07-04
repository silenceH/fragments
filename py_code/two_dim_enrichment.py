import enrichment
import Bioisosteres

def get_enrichment(target,*args):
	# get a list of all the ranked 2D pairs
	all_fragments = enrichment.rank_targets_test(target)[0]
	bioisostere_pairs = [x.group for x in  Bioisosteres.get_bioisosteres(target,test=True,return_pairs=True)]
	ranked_list = [0 for x in range(len(all_fragments))]

	for i in range(len(all_fragments)):
		test_pair = all_fragments[i]
		for ranked_pair in bioisostere_pairs:
			if (test_pair[0].are_similar(ranked_pair[0],1) and test_pair[1].are_similar(ranked_pair[1],1)) or \
					(test_pair[0].are_similar(ranked_pair[0],1) and test_pair[1].are_similar(ranked_pair[1],1)):
				ranked_list[i] = 1
	
	# assert sum(ranked_list) == len(bioisostere_pairs)


	actives = len(bioisostere_pairs)
	total_pairs = len(all_fragments)

	print "total actives = " + str(actives)
	print "total inactives = " + str(total_pairs-actives)
	print "enrichment factor threshold" 

	efs = []

	for threshold in args:
		try:
			n = int(total_pairs * threshold)
			pick_n = ranked_list[:n]
			t_p = sum(pick_n)
			ef = (t_p/float(n))/(actives/float(total_pairs))
			print str(threshold) +'\t' + str(ef)
			efs.append(ef)
		except ZeroDivisionError:
			print "zero division"
			efs.append('nan')

	return efs

	




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
f = open('full_enrichment_stats.csv','w')
f.write('target,0.01,0.05,0.13,0.2,0.5,\n')
for target in targets:
	efs = get_enrichment(target,0.01,0.05,0.13,0.2,0.5)
	f.write(str(target)+',')
	for ef in efs:
		f.write(str(ef) + ',')
	f.write('\n')
	

