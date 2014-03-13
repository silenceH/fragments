## KENNEWELL DEV

import os
import time
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
	test_mol_frags = get_first_fragment_from_smarts(mol,bond_smarts)
	final_frags = []
	final_frags.extend(test_mol_frags)
	while(len(test_mol_frags)>0):
		temp_frags = []
		for mol in test_mol_frags:
			temp = get_first_fragment_from_smarts(mol,bond_smarts)
			if temp is not None:
				temp_frags.extend(temp)
		temp_frags = [x for x in temp_frags if x is not None]
		final_frags.extend(temp_frags)
		test_mol_frags=temp_frags
	return final_frags

## TODO :: update this method for Fragment object
def get_first_fragment_from_smarts(mol,smarts):
	## fragments a molecule at first occurence of smarts 
	## note that smarts is a RDMol Object
	if mol.GetSubstructMatches(smarts):
		bonds = mol.GetSubstructMatches(smarts)
		bonds = [((x,y),(0,0)) for x,y in bonds]
		frags = list(Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol,bonds=[bonds[0]]),asMols=True))
		return frags

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

def get_bioisosteres(data_file,noHs,brics, kennewell,overlap,test):
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
	# candidate_pairs = []
	grouped = []
	count = 0
	## create a list of fragment objects here
	mols = [get_fragments(mol,brics,data_file) for mol in mols]	
	## for a reference molecule mol in mols
	for i in range(len(mols)):
		## set query data set to be the molecules that are not the reference
		ref_frag = mols[i] 			# a list of Fragment objects
		query_set = mols[i:]
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
						count += 1
					else: 
						candidate = score_pairs_TD(s,f)		## Tanimoto distance
					if candidate and not s.are_similar(f,1.0): 
						section_group.add(f)
				if section_group.size() > 1:
					#candidate_pairs.append(section_group)
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
		# write groups to file
		write_mols_to_file(grouped,'final_group',directory)
		
		# draw mols to png
		draw_mols_to_png(grouped,'final_group',directory)
	
	## produce statistics and write to file ## 
	f = open(directory+"stats",'w')
	## potential stats
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
	return grouped

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
	#for i in range(len(final_collection)):
		#print str(i) + "\t" + str(len(final_collection[i])) + "\t" + str(list(lig_per_group[i]))
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
	f.write("group,number,\n")
	for i in range(len(final_collection)):
		f.write(str(i+1) + "," + str(final_collection[i].size())+",\n")
	f.close()
	print "statistics written"
	draw_mols_to_png(final_collection,'final_collection',directory)
	return final_collection


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

def two_dim_similars(data_file,threshold):
	mols = get_mols_from_sdf_file(data_file,True)
	fragments = [] 
	## obtain 1 dimensional list of fragments
	for m in mols:
		fragments.extend(get_fragments(m,True,data_file))
	
	## obtain upper triangular boolean matrix based on similarity test
	sim_matrix = []
	for i in range(len(fragments)):
		## are similar at a threshold but not identical
		row = [fragments[i].are_similar(fragments[j],threshold) and not fragments[i].are_similar(fragments[j],1) for j in range(len(fragments)) if j > i]
		sim_matrix.append(row)
	pairs = [(x,y+1+x) for x in range(len(sim_matrix)) for y in range(len(sim_matrix[x])) if sim_matrix[x][y]]
	print data_file
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
	## get pairs from files 1 to 10 
	test_set = collect_bioisosteres(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
	valid_groups = []
	for i in xrange(len(test_set)):
		group = test_set[i]
		all_in_group = True
		count = 0	
		for q in unique_pairs:
			mol_in_group = False
			for ref in group:
				if ref.are_similar(q,1):
					mol_in_group = True
			all_in_group = all_in_group and mol_in_group
		if all_in_group:
				valid_groups.append(i+1)
	f = open(directory+data_file+"_valid groups",'w')
	f.write("valid groups: \n") 
	for group in valid_groups:
		f.write(str(group)+"\n")
	print "\n\n\n"

