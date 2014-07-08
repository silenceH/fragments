import enrichment
import Bioisosteres

def get_enrichment(threeDim,target,*args):
	# get a list of all the ranked 2D pairs
	all_fragments = enrichment.rank_targets_test(False,True, target) # True for pair objects
	threeDim_frags = enrichment.rank_targets_test(True,False,target)
	bioisostere_pairs = [x.group for x in  Bioisosteres.get_bioisosteres(target,test=True,return_pairs=True)]
	ranked_list = [0 for x in range(len(all_fragments))]
	three_Dim_ranked_list = [0 for x in range(len(all_fragments))]

	for i in range(len(all_fragments)):
		test_pair = all_fragments[i].frags
		test_pair_3D = threeDim_frags[i].frags
		for ranked_pair in bioisostere_pairs:
			if (test_pair[0].are_similar(ranked_pair[0],1) and test_pair[1].are_similar(ranked_pair[1],1)) or \
					(test_pair[0].are_similar(ranked_pair[0],1) and test_pair[1].are_similar(ranked_pair[1],1)):
				ranked_list[i] = 1
			if threeDim:
				if (test_pair_3D[0].are_similar(ranked_pair[0],1) and test_pair_3D[1].are_similar(ranked_pair[1],1)) or \
					(test_pair_3D[0].are_similar(ranked_pair[0],1) and test_pair_3D[1].are_similar(ranked_pair[1],1)):
					three_Dim_ranked_list[i] = 1
	
	actives = len(bioisostere_pairs)
	total_pairs = len(all_fragments)

	print "total actives = " + str(actives)
	print "total inactives = " + str(total_pairs-actives)
	print "enrichment factor threshold" 

	efs = []
	efs_3D = []

	for threshold in args:
		try:
			n = int(round(total_pairs * threshold))
			pick_n = ranked_list[:n]
			t_p = sum(pick_n)
			ef = (float(t_p)/float(n))/(float(actives)/float(total_pairs))
			pick_n_3D = three_Dim_ranked_list[:n]
			t_p_3D = sum(pick_n_3D)
			ef_3D = (float(t_p_3D)/float(n))/(float(actives)/float(total_pairs))
			print str(threshold) +'\t' + str(ef) + '\t' + str(ef_3D)
			efs.append(ef)
			efs_3D.append(ef_3D)
		except ZeroDivisionError:
			print "zero division"
			print total_pairs * threshold
			efs.append('nan')
	if threeDim:
		return efs,efs_3D
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
f.write('target,0.01,0.05,0.13,0.01,0.05,0.13,0.2,0.5,\n')
for target in targets[:10]:
	efs,efs_3D = get_enrichment(True,target,0.01,0.05,0.13)
	f.write(str(target)+',')
	for ef in efs:
		f.write(str(ef) + ',')

	for ef_3D in efs_3D:
		f.write(str(ef_3D) + ',')
	f.write('\n')
	
f.close()
