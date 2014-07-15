from Bioisosteres import * 
from rdkit import Chem

## Programme to find number of bioisosteres that are present in different targets

def count_cross_target_pairs(args):
	""" A function that collects the frequencies of the pairs over a number of ligand overlays """
	coll = []
	for data_file in args:
		pairs = get_bioisosteres(data_file,return_pairs=True)
		coll.extend(pairs)

	coll = [x.group for x in coll]
	coll[0]=(coll[0][0],coll[0][1])		
	final_dict = {coll[0]:[coll[0][0].ligand]}
	for frag0,frag1 in coll[1:]:
		not_in_dict = True
		for final_pair in final_dict:
			if frag0.are_similar(final_pair[0],1.0) or frag1.are_similar(final_pair[0],1): 
				if frag0.are_similar(final_pair[1],1.0) or frag1.are_similar(final_pair[1],1): 
					not_in_dict = False 
					pair_target = frag0.ligand
					if pair_target not in final_dict[final_pair]:
						final_dict[final_pair].append(pair_target)
		if not_in_dict:
			final_dict[(frag0,frag1)] = [frag0.ligand]

	
	## make dirs by count number 
	for frag_pair in final_dict:
		count = len(final_dict[frag_pair])
		if count > 1:
			directory = '../test_output/common_pairs/'+str(count)+'/' 
			try: 
				os.makedirs(directory)
				print "created new directory: " + directory
			except OSError:
				print directory + " already exists" 
			
			w = Chem.SDWriter(directory+str(frag_pair)+'.sdf')
			for mol in frag_pair: 
				w.write(mol.frag)
			w.flush()
			
	stats_file = '../test_output/common_pairs/common_pair_stats.csv' 
	stats = open(stats_file,'w')
	for x in final_dict.values():
		for y in [len(x) for x in final_dict.values()]:
			if y>1:
				stats.write(str(x) + ',' + str(y) + ',\n')
	stats.close()



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

count_cross_target_pairs(targets[:75])


#get_enrichment(targets[:75])
