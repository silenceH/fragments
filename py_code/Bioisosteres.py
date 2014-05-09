## KENNEWELL DEV

import os
import time
import matplotlib.pyplot as plt
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, BRICS, Draw, rdShapeHelpers, Descriptors
from RDMols import * 
from Fragment import *
from FragmentGroup import *

def get_overlapping_fragments(mol):
	## returbs a list of overlapping fragments using a SMARTS pattern for a particular type of bond
	## bond_smarts = Chem.MolFromSmarts("[*]!@!#!=[*]")	## single bonds
	## bond_smarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]') ## simple rotatable bond smarts
	bond_smarts = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]') ## rotatable bonds
	if mol.GetSubstructMatches(bond_smarts):
		num_bonds = len(mol.GetSubstructMatches(bond_smarts))
		final_frags = Group()
		for i in range(num_bonds):
			test_mol_frags = get_i_fragment_from_smarts(mol,bond_smarts,i)
			if test_mol_frags is not None:
				final_frags.merge(test_mol_frags)
				while(test_mol_frags.size()>0):
					temp_frags = Group()
					for temp_mol in test_mol_frags.group:
						temp = get_i_fragment_from_smarts(temp_mol.frag,bond_smarts,0)
						if temp is not None:
							temp_frags.merge(temp)
					final_frags.merge(temp_frags)
					test_mol_frags=temp_frags
		return final_frags.group
	return []

def get_i_fragment_from_smarts(mol,smarts,i):
	## fragments a molecule at first occurence of smarts 
	## note that smarts is a RDMol Object
	if mol.GetSubstructMatches(smarts):
		bonds = mol.GetSubstructMatches(smarts)
		bonds = [((x,y),(0,0)) for x,y in bonds]
		frags = list(Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol,bonds=[bonds[i]]),asMols=True))
		frag_group = Group()
		for tmp_frag in frags:
			frag_group.add(Fragment(tmp_frag))
		return frag_group

## TODO :: update this method for Fragment object
def get_section_set(sections):
	## function to collect different sections of candidate bioisosteres
	grouped = []
	for section in sections:
		
		qs = [x for x in sections if x is not mols]
		for q in qs:
			q = set(q)
			if m.intersection(q):
				m.update(q)
		grouped.append(list(m))
	return grouped
		

def score_pairs_TD(atom1,atom2):
	## scores two fragments using Tanimoto overlap of shape
	## Note: shape calculated from 3D grid represnentation of mol
	score = rdShapeHelpers.ShapeTanimotoDist(atom1.frag,atom2.frag)
	if score < 0.3: 
		return True
	else:
		return False

def merge_frag_groups(ref_group,merge_group):
	## take two groups of Frag objects and merge them 
	## so that there are no duplicates
	add_group = []
	for mol in merge_group:
		unique = True
		for ref in ref_group:
			if ref.are_similar(mol,1):
				unique = False
				break
		if unique:
			add_group.append(mol)
	return ref_group.extend(add_group)

def get_bioisosteres(data_file,noHs=True,brics=True, kennewell=True,overlap=False,test=False,return_pairs=False,count_pairs=False,debug=False):
	## gets bioisosteric pairs from an sdf file
	## data_file : string of data file without .sdf
	## noHs : specifies whether Hs are to be removed (boolean)
	## brics : specifies whether brics bonds are broken at fragmentation (boolean)
	## kennewell : specifies whether kennewell scoring is to be used (boolean)
	## overlap : specifies whether fragments should be overlapping (boolean)
	## test : if true, does not write files 

	## import mols from data_file 
	mols = get_mols_from_sdf_file(data_file,noHs)
	# list of candidate bioisosteric pairs to view individual pairs if needed commented out
	candidate_pairs = []
	grouped = []
	pairs = []
	count = 0

	## create dictionary for counting frequency of pairs
	if count_pairs:
		pair_frequency = {}
	
	## create a list of fragment objects here
	mols = [get_fragments(mol,brics,data_file) for mol in mols]	
	## for a reference molecule mol in mols
	for i in range(len(mols)-1):
		## set query data set to be the molecules that are not the reference
		ref_frag = mols[i] 			# a list of Fragment objects
		query_set = mols[i+1:]
		## for each remaining ligand make it the query ligand
		for q in query_set:
			frags = q
			for s in ref_frag:
				# create fragment group
				section_group = Group(s)
				for f in frags:
					candidate = True
					if kennewell:  
						candidate = s.score_pairs_kennewell(f) 	## Kennewell scores
					else: 
						candidate = score_pairs_TD(s,f)		## Tanimoto distance
					if candidate and not s.are_similar(f,1.0): 
						# For dedugging purposes print details
						if debug:
							print "pair: " + str(count+1)
							print "Ref: " + str(i) + "\tfrag: " + str(ref_frag.index(s))
							print "Query: " + str(i+query_set.index(q)) + "\tfrag: " + str(frags.index(f))
						# code for returning a list of pairs
						if return_pairs:
							pair = Group(s)
							pair.add(f)
							candidate_pairs.append(pair)

						# add pair to the dictionary if count_pairs
						if count_pairs:
							# look for test pair in dict pairs 
							not_in_dict = True
							for pair in pair_frequency:
								if s.are_similar(pair[0],1.0) or f.are_similar(pair[0],1): 
									if s.are_similar(pair[1],1.0) or f.are_similar(pair[1],1): 
										not_in_dict = False
										pair_frequency[pair] += 1
							# add pair
							if not_in_dict:
								pair_frequency[(s,f)] = 1

						section_group.add(f)
						count += 1
				if section_group.size() > 1:
					in_grouped = False
					for x in grouped:
						if section_group.get_mol(0) in x.group:
							in_grouped = True
							x.merge(section_group) 
					if not in_grouped:
						grouped.append(section_group)

	if noHs: data_file += "_noHs"
	if brics: data_file += "_brics"
	if overlap:  data_file += "_overlap"
	if kennewell: data_file += "_KN"

	directory = '../test_output/'+data_file+'_candidate_bioisosteres/' 
	try: 
		os.makedirs(directory)
		print "created new directory: " + directory
	except OSError:
		print directory + " already exists."	
	if not test: 
		if return_pairs:
			write_mols_to_file(candidate_pairs,'pair',directory)
		# write groups to file
		write_mols_to_file(grouped,'final_group',directory)
		
		if return_pairs:
			draw_mols_to_png(candidate_pairs,'pair',directory)
		# draw mols to png
		draw_mols_to_png(grouped,'final_group',directory)
	
	## produce statistics and write to file ## 
	f = open(directory+"stats",'w')
	## number of fragments with overlapping scheme?? 

	total_sim_of_groups = 0
	f.write("SUMMARY STATISTICS FOR OVERLAY " + data_file + "\n")
	f.write("Tanimoto cut off for similarity is 1.\n")
	if overlap: f.write("Overlapping fragments.\n")
	if noHs: f.write("Hydrogen's removed.\n")
	if kennewell: f.write("Kennewell scoring\n")
	f.write("Number of fragments: " + str(sum([len(m) for m in mols])) + "\n")
	f.write("Number of identified pairs : " + str(count) + "\n")
	f.write("Number of overlaid sections: " + str(len(grouped)) + "\n")
	f.write("group \tfrags \tav similarity\n")
	for i in range(len(grouped)):
		av_sim = grouped[i].get_av_similarity()
		f.write(str(i+1) + "\t" + str(grouped[i].size()) + "\t" + str(av_sim)+"\n")
		total_sim_of_groups += av_sim
	f.write("\naverage similarity over the group: " + str(total_sim_of_groups/len(grouped)) + "\n\n")
	print "no pairs: " + str(count)
	if return_pairs:
		if len(candidate_pairs) == count:
			print "PAIRS WORKED"
		return candidate_pairs
	if count_pairs:
		return pair_frequency
	return grouped

def get_pair_frequency(write_stats,*args):
	""" A function that collects the frequencies of the pairs over a number of ligand overlays """
	coll = [get_bioisosteres(data_file,count_pairs=True) for data_file in args]
	## NOTE:: this gives a list of dictionaries

	final_dict = coll.pop()
	for merge_dict in coll:
		for merge_pair in merge_dict:
			not_in_dict = True
			for final_pair in final_dict:
				if merge_pair[0].are_similar(final_pair[0],1.0) or merge_pair[1].are_similar(final_pair[0],1): 
					if merge_pair[0].are_similar(final_pair[1],1.0) or merge_pair[1].are_similar(final_pair[1],1): 
						not_in_dict = False 
						final_dict[final_pair] += merge_dict[merge_pair]
			if not_in_dict:
				final_dict[merge_pair] = merge_dict[merge_pair]

	## create a dictionary of frequencies
	frequency_stats = {}
	for pair in final_dict:
		val = final_dict[pair]
		if val in frequency_stats:
			frequency_stats[val] += 1
		else:
			frequency_stats[val] = 1
	
	print frequency_stats
	
	home= os.environ['HOME']		## get data env
	directory = home+'/Dropbox/test_output/pairs_with_frequency_'

	frequencies = frequency_stats.keys()
	pairs_by_frequency = []
	for i in frequencies:
		pairs_with_freqency_i = [val for val,key in final_dict.items() if key == i]
		pairs_by_frequency.append(pairs_with_freqency_i)
		print i
		print pairs_with_freqency_i
		
		try: 
			os.makedirs(directory+str(i))
		except OSError:
			print directory + str(i) + " already exists."	
		for p in pairs_with_freqency_i:
			draw = [p[0],p[1]]
			draw_mols_to_png([draw],"/pair_"+str(pairs_with_freqency_i.index(p)),directory)

	

	bar = []
	for x in range(max(frequency_stats)):
		if (x+1) in frequency_stats:
			bar.append(frequency_stats[x+1])
		else:
			bar.append(0)
	
	print bar

	if write_stats:
		plt.bar(range(len(frequency_stats)),frequency_stats.values(),align='center')
		plt.xticks(range(len(frequency_stats)),frequency_stats.keys())

		plt.show()
	f = open("stats.csv",'w')
	f.write("SUMMARY STATISTICS,\n")
	f.write("number of occurences,number of pairs,\n")
	for i in bar:
		f.write(str(bar.index(i) + 1)+','+str(i)+',\n')
	f.close()




def collect_bioisosteres(*args):
	coll = [get_bioisosteres(data_file,True,True,True,False,True) for data_file in args]
	## NOTE:: this gives a list of Group objects
	collection = [coll[i][j] for i in range(len(coll)) for j in range(len(coll[i]))]
	final_collection = [collection[0]] 
	collection= collection[1:]
	print " number of groups to compare: " + str(len(collection)+ 1)
	print "comparing..."
	while len(collection) > 0: 
		count = 0
		## compares each fragment in the group to those in the final collection 
		not_extended = True
		q_frags = collection[0]
		for ref_frags in final_collection:
			match = False
			for ref,q in ((m1,m2) for m1 in ref_frags.group for m2 in q_frags.group):
				count += 1
				if ref.are_similar(q,1.):
					match = True
					break
			if match:
				# need to merge so that mols are not duplicated
				ref_frags.merge(q_frags)
				not_extended = False
		if not_extended:
			final_collection.append(q_frags)
		collection= collection[1:]
	lig_per_group = [set([frag.ligand for frag in frag_group.group]) for frag_group in final_collection]
	for i in range(len(final_collection)):
		print str(i) + "\t" + str(final_collection[i].size()) + "\t" + str(list(lig_per_group[i]))
	directory = '../test_output/compared_results_NO_SMILES/' 
	try: 
		os.makedirs(directory)
		print "created new directory: " + directory
	except OSError:
		print directory + " already exists."	
	write_mols_to_file(final_collection,'final_collection',directory)
	f = open(directory+"stats.csv",'w')
	f.write("files,\n")
	for i in args:
		f.write(i + ",\n")
	f.write("total," +  str(len(final_collection))+",\n")
	f.write("group,number,ligand_sets,av_sim,\n")
	for i in range(len(final_collection)):
		f.write(str(i+1) + "," + str(final_collection[i].size())+ "," + str(len(lig_per_group[i]))+ "," + \
				str(final_collection[i].get_av_similarity()) + ",\n")
	f.close()
	print "statistics written"
	draw_mols_to_png(final_collection,'final_collection',directory)
	return final_collection

def collect_bioisosteres_greedy(from_pairs,*args):
	coll = [get_bioisosteres(data_file,True,True,True,False,True) for data_file in args]
	if from_pairs:
		coll = [get_bioisosteres(data_file,return_pairs=True) for data_file in args]
	## NOTE:: this gives a list of Group objects
	## reduce to give a one dimesional list of groups
	collection = [coll[i][j] for i in range(len(coll)) for j in range(len(coll[i]))]
	final_collection = [collection[0]] 
	collection= collection[1:]
	print " number of groups to compare: " + str(len(collection)+ 1)
	finished = False
	while not finished:
		tmp_changes = 0
		original_len = len(collection) + 1
		print "comparing..."
		while len(collection) > 0: 
			## compares each fragment in the group to those in the final collection 
			not_extended = True
			q_frags = collection[0]
			for ref_frags in final_collection:
				match = False
				for ref,q in ((m1,m2) for m1 in ref_frags.group for m2 in q_frags.group):
					if ref.are_similar(q,1.):
						match = True
						break
				if match:
					# need to merge so that mols are not duplicated
					ref_frags.merge(q_frags)
					not_extended = False
			if not_extended:
				final_collection.append(q_frags)
			collection= collection[1:]
		new_len = len(final_collection)
		changes = original_len - new_len
		if changes == 0:
			finished = True
		else:
			collection = final_collection
			final_collection = [collection[0]] 
			collection= collection[1:]

	lig_per_group = [set([frag.ligand for frag in frag_group.group]) for frag_group in final_collection]
	for i in range(len(final_collection)):
		print str(i) + "\t" + str(final_collection[i].size()) + "\t" + str(list(lig_per_group[i]))
	directory = '../test_output/compared_results_NO_SMILES/greedy/' 
	if from_pairs:
		directory += 'from_pairs/'
	try: 
		os.makedirs(directory)
		print "created new directory: " + directory
	except OSError:
		print directory + " already exists."	
	write_mols_to_file(final_collection,'final_collection',directory)
	f = open(directory+"stats.csv",'w')
	f.write("files,\n")
	for i in args:
		f.write(i + ",\n")
	f.write("total," +  str(len(final_collection))+",\n")
	f.write("group,number,ligand_sets,av_sim,\n")
	for i in range(len(final_collection)):
		f.write(str(i+1) + "," + str(final_collection[i].size())+ "," + str(len(lig_per_group[i]))+ "," + \
				str(final_collection[i].get_av_similarity()) + ",\n")
	f.close()
	print "statistics written"
	draw_mols_to_png(final_collection,'final_collection',directory)
	return final_collection

def debug_bioisosteres(file_name):
	mols = get_mols_from_sdf_file(file_name,True)
	ref_mol = get_fragments(mols[0],True,False)
	query_mol = get_fragments(mols[1],True,False)
	print "num ref frags: " + str(len(ref_mol))
	print "num query frags: " + str(len(query_mol))
	for r in ref_mol: 
		print "i am reference fragment: " + str(ref_mol.index(r))
		for q in query_mol:
			print "\ti am query fragment: " + str(query_mol.index(q))
			candidate = r.score_pairs_kennewell(q,test=True) 	## Kennewell scores




def collect_bioisosteres_by_smiles(*args):
	coll = [get_bioisosteres(data_file,True,True,True,False,False) for data_file in args]
	collection = [coll[i][j] for i in range(len(coll)) for j in range(len(coll[i]))]
	# create a dictionary of the smiles with mol objects as values
	mols_by_smiles = dict()
	for mol in (collection[i].get_mol(j) for i in range(len(collection)) for j in range(collection[i].size())):
		mol.smiles = Chem.MolToSmiles(mol.frag)
		if mol.smiles in mols_by_smiles:
			mols_by_smiles[mol.smiles].append(mol)
		else: 
			mols_by_smiles[mol.smiles] = [mol]
	smiles = [[collection[i].get_mol(j).smiles for j in range(collection[i].size())] for i in range(len(collection))]
	final_collection = [set(smiles[0])]
	number_of_files = [1]
	smiles= smiles[1:]
	print " number of groups to compare: " + str(len(collection)+ 1)
	print "comparing..."
	it = 0
	while len(smiles) > 0: 
		count = 0
		q_frags = set(smiles[0])
		matched = False
		for ref_frags in final_collection:
			if len(ref_frags.intersection(q_frags)) > 0:
				ref_frags.update(q_frags)
				number_of_files[final_collection.index(ref_frags)] += 1
				count += 1
				matched = True
		if not matched:
			final_collection.append(q_frags)
			number_of_files.append(1)
		smiles= smiles[1:]
		it += 1
		## testing code 
		## print "groups extended: " + str(count)
		## print "length of final collection: " + str(len(final_collection))
		## print "number of iterations: " + str(it)
		## print "comparing... \n"
	directory = '../test_output/compared_results/' 
	try: 
		os.makedirs(directory)
		print "created new directory: " + directory
	except OSError:
		print directory + " already exists."	
	## write to file
	f = open(directory+"stats.csv",'w')
	f.write("files,\n")
	for i in args:
		f.write(i + ",\n")
	f.write("total," +  str(len(final_collection))+",\n")
	f.write("group,number,\n")
	for i in range(len(final_collection)):
		f.write(str(i+1) + "," + str(len(final_collection[i]))+",\n")
	f.close()
	print "statistics written"
	final_collection = [[Fragment(Chem.MolFromSmiles(mol),'final') for mol in list(smi)] for smi in final_collection]
	draw_mols_to_png(final_collection,'final_collection',directory)
	print "pictures drawn"
	return final_collection

def two_dim_similars(data_file,threshold,full_screen=False):
	mols = get_mols_from_sdf_file(data_file,True)
	fragments = [] 
	## obtain 1 dimensional list of fragments
	for m in mols:
		fragments.extend(get_fragments(m,True,data_file))
	
	## obtain average fragment size and statistics
	num_atoms = [f.frag.GetNumAtoms() for f in fragments]
	min_num_atoms = min(num_atoms)
	max_num_atoms = max(num_atoms)
	av_num_atoms = sum(num_atoms)/len(num_atoms)
	
	## obtain upper triangular boolean matrix based on similarity test
	sim_matrix = []
	for i in range(len(fragments)):
		## are similar at a threshold but not identical
		row = [fragments[i].are_similar(fragments[j],threshold) and not fragments[i].are_similar(fragments[j],1) for j in range(len(fragments)) if j > i]
		sim_matrix.append(row)
	pairs = [(x,y+1+x) for x in range(len(sim_matrix)) for y in range(len(sim_matrix[x])) if sim_matrix[x][y]]
	## Output data about the fragments and the pairs
	print data_file
	print "number of fragments: " + str(len(fragments))
	print "largest fragment: " + str(max_num_atoms)
	print "smallest fragment: " + str(min_num_atoms)
	print "average number size: " + str(av_num_atoms)
	print "there are " + str(len(pairs)) + " pairs at threshold " + str(threshold) + " that are not identical."
	final = []
	for pair in pairs:
		for i in range(2):
			final.append(fragments[pair[i]])
	try:	
		unique_pairs = [final[0]]
		for ref in final:
			ref_is_unique = True
			for q in unique_pairs:
				if ref.are_similar(q,1):
					ref_is_unique = False
			if ref_is_unique:
				unique_pairs.append(ref)
		print "of which there are " + str(len(unique_pairs)) + " unique."
		print "The smiles are: "
		for mol in unique_pairs:
			if mol.smiles is None:
				mol.smiles = Chem.MolToSmiles(mol.frag)
			print mol.smiles
		directory = '../test_output/compared_results/two_dim_pairs/' 
		try: 
			os.makedirs(directory)
			print "created new directory: " + directory
		except OSError:
			print directory + " already exists."	
		draw_mols_to_png([unique_pairs],data_file,directory)
	except IndexError:
		print "EXIT PROGRAMME!"
		return
	## compare fragments against full screen 
	if full_screen:
		## get pairs from files 1 to 10 
		file_1 = 'P39900' 
		file_2 = 'P56817'
		file_3 = 'P35557'
		file_4 = 'Q92731'
		file_5 = 'P25440'
		file_6 = 'P00918'
		file_7 = 'P0AE18'
		file_8 = 'P43235'
		file_9 = 'Q00511'
		file_10 = 'P16184'
		test_set = collect_bioisosteres(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
		valid_groups = []
		for i in xrange(len(test_set)):
			test_group = test_set[i]
			all_in_group = True
			count = 0	
			for q in unique_pairs:
				mol_in_group = False
				for ref in test_group.group:
					if ref.are_similar(q,1):
						mol_in_group = True
				all_in_group = all_in_group and mol_in_group
			if all_in_group:
					valid_groups.append(i+1)
		f = open(directory+data_file+"_valid groups",'w')
		f.write("valid groups: \n") 
		for g in valid_groups:
			f.write(str(g)+"\n")
		print "\n\n\n"

