## KENNEWELL DEV

import os
import time
from math import sqrt, exp
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, BRICS, Draw, rdShapeHelpers, Descriptors

class Fragment(object):
	def __init__(self, frag,ligand, coords = None, fp=None,smiles=None):
		self.frag = frag

		self.ligand = ligand

		if coords is None:
			self.coords = None
		
		if fp is None:
			self.fp = None
		self.fp = fp

		if smiles is None:
			self.smiles=None
		self.smiles = smiles


def get_mols_from_sdf_file(data_file, noHs):
	## Method to return molecules for a given sdf data file
	try:
		data = os.environ['DATA']		## get data env
		print "found data environment: " + data
	except KeyError:
		print "cannot find data environment variable"

	if noHs:
		suppl = Chem.SDMolSupplier(data + '/validation_overlays/'+data_file+'.sdf')
	else: 
		suppl = Chem.SDMolSupplier(data + '/validation_overlays/'+data_file+'.sdf',removeHs=False)
	return [x for x in suppl if x is not None]

def atom_coords(molecule,atom_no):
	## returns a tuple of atom coordinates
	pos = molecule.GetAtomPosition(atom_no)
	return (pos.x,pos.y,pos.z)

def get_distance(coords1,coords2):
	## returns Euclidean distance between two 3D coordinates
	return sqrt(pow((coords1[0]-coords2[0]),2)+pow((coords1[1]-coords2[1]),2)+pow((coords1[2]-coords2[2]),2))

def are_similar(frag1,frag2,threshold):
	## returns False if Tanimoto similarity is greater than threshold
	if frag1.fp is None:
		frag1.fp = AllChem.GetMorganFingerprint(frag1.frag,2)
	if frag2.fp is None:
		frag2.fp = AllChem.GetMorganFingerprint(frag2.frag,2)
	tanimoto = DataStructs.TanimotoSimilarity(frag1.fp,frag2.fp)
	if tanimoto >= threshold:
		return True 
	else: 
		return False

def remove_2D_equivalents(mols):
	## remove duplicate 2D mols
	u = []
	for mol in mols:
		include = True
		for m in u:
			if are_similar(mol,m,1.0):
				include = False
		if include:
			u.append(mol)
	return u

def write_mols_to_file(mols,title,directory):
	## write to file
	for i in mols:
		w = Chem.SDWriter(directory+title+str(mols.index(i))+'.sdf')
		for mol in i: 
			w.write(mol.frag)
		w.flush()

def draw_mols_to_png(mols,title,directory):
	## test for Fragment object
	if isinstance(mols[0][0],Fragment):
		mols = [[mol.frag for mol in group] for group in mols]
	## draw image
	for i in mols:
		for mol in i:
			tmp = AllChem.Compute2DCoords(mol)
		img = Draw.MolsToGridImage(i,legends=[str(x+1) for x in range(len(i))])
		img.save(directory+title+str(mols.index(i))+'.png')

def get_all_coords(mol):
	## return a list of coords for the 3D shape of a mol 
	# define a constructor of position objects
	tmp = (mol.GetConformer().GetAtomPosition(i) for i in range(mol.GetNumAtoms()))
	
	# return the list of the coordinates from those positions
	return [(atom.x,atom.y,atom.z) for atom in tmp]

def get_fragments(mol,brics,data_file):
	## returns a list of non overlapping fragments.
	## if brics == True then returns the BRICS fragments
	if brics: 
		return [Fragment(x,data_file) for x in Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol),asMols=True)]	## BRICS bonds
	else: 
		## bond_smarts = Chem.MolFromSmarts("[*]!@!#!=[*]")	## single bonds
		## bond_smarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]') ## simple rotatable bond smarts
		bond_smarts = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]') ## rotatable bonds
		bonds = mol.GetSubstructMatches(bond_smarts)
		bonds = [((x,y),(0,0)) for x,y in bonds]
		return [Fragment(x,data_file) for x in Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol,bonds=bonds),asMols=True)]

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

def score_pairs_kennewell(mol1,mol2):
	## get number of atoms for each molecule
	frag1, frag2 = mol1.frag,mol2.frag
	atoms_ref = frag1.GetNumAtoms()
	atoms_f = frag2.GetNumAtoms()
	## get coordinates for the reference molecule
	ref_atoms = get_all_coords(frag1)
	section_score = []
	for section_atom in ref_atoms:
		dist = [get_distance(section_atom,frag_atom) for frag_atom in get_all_coords(frag2)]
		section_score.append(sum([exp(-pow(d,2)) for d in dist]))
	av_score = sum(section_score)*(2./(atoms_f+atoms_ref))
	if av_score>0.7:
		return True 
	else: 
		return False

def get_av_similarity(mols):
	total_sim = 0
	num_mols = len(mols)
	for i in range(num_mols):
		row = (DataStructs.TanimotoSimilarity(mols[i].fp,mols[j].fp) for j in range(num_mols) if j > i)
		total_sim += sum(row)
	return total_sim/num_mols

def merge_frag_groups(ref_group,merge_group):
	## take two groups of Frag objects and merge them 
	## so that there are no duplicates
	add_group = []
	for mol in merge_group:
		unique = True
		for ref in ref_group:
			if are_similar(ref,mol,1):
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
	candidate_pairs = []
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
				section_pairs = [s]
				for f in frags:
					candidate = True
					if kennewell:  
						candidate = score_pairs_kennewell(s,f) 	## Kennewell scores
						count += 1
					else: 
						candidate = score_pairs_TD(s,f)		## Tanimoto distance
					if candidate and not are_similar(s,f,1.0): 
						section_pairs.append(f)
				if len(section_pairs) > 1:
					candidate_pairs.append(section_pairs)
					in_grouped = False
					for group in grouped:
						if section_pairs[0] in group:
							in_grouped = True
							group += section_pairs
					if not in_grouped:
						grouped.append(section_pairs)
	# depuplicate groups
	for group in grouped:
		u = []
		for mol in group:
			if mol not in u:
				u.append(mol)	
		group = u
	# remove 2D equivalent mols
	final_group = [remove_2D_equivalents(g) for g in grouped]

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
		write_mols_to_file(final_group,'final_group',directory)
		
		# draw mols to png
		draw_mols_to_png(final_group,'final_group',directory)
	
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
	f.write("Number of overlaid sections: " + str(len(final_group)) + "\n")
	f.write("group \tfrags \tav similarity\n")
	for i in range(len(final_group)):
		av_sim = get_av_similarity(final_group[i])
		f.write(str(i+1) + "\t" + str(len(final_group[i])) + "\t" + str(av_sim)+"\n")
		total_sim_of_groups += av_sim
	f.write("\naverage similarity over the group: " + str(total_sim_of_groups/len(final_group)) + "\n\n")
	return final_group

def collect_bioisosteres(*args):
	coll = [get_bioisosteres(data_file,True,True,True,False,True) for data_file in args]
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
			for ref,q in ((m1,m2) for m1 in ref_frags for m2 in q_frags):
				count += 1
				if are_similar(ref,q,1.):
					match = True
					break
			if match:
				# need to merge so that mols are not duplicated
				merge_frag_groups(ref_frags,q_frags)
				not_extended = False
		if not_extended:
			final_collection.append(q_frags)
		collection= collection[1:]
	lig_per_group = [set([frag.ligand for frag in group]) for group in final_collection]
	#for i in range(len(final_collection)):
		#print str(i) + "\t" + str(len(final_collection[i])) + "\t" + str(list(lig_per_group[i]))
	directory = '../test_output/compared_results_NO_SMILES/' 
	try: 
		os.makedirs(directory)
		print "created new directory: " + directory
	except OSError:
		print directory + " already exists."	
	write_mols_to_file(final_collection,'final_collection',directory)
	draw_mols_to_png(final_collection,'final_collection',directory)
	return final_collection


def collect_bioisosteres_by_smiles(*args):
	coll = [get_bioisosteres(data_file,True,True,True,False,False) for data_file in args]
	collection = [coll[i][j] for i in range(len(coll)) for j in range(len(coll[i]))]
	# create a dictionary of the smiles with mol objects as values
	mols_by_smiles = dict()
	for mol in [collection[i][j] for i in range(len(collection)) for j in range(len(collection[i]))]:
		mol.smiles = Chem.MolToSmiles(mol.frag)
		if mol.smiles in mols_by_smiles:
			mols_by_smiles[mol.smiles].append(mol)
		else: 
			mols_by_smiles[mol.smiles] = [mol]
	smiles = [[collection[i][j].smiles for j in range(len(collection[i]))] for i in range(len(collection))]
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
		row = [are_similar(fragments[i],fragments[j],threshold) and not are_similar(fragments[i],fragments[j],1) for j in range(len(fragments)) if j > i]
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
				if are_similar(ref,q,1):
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
				if are_similar(ref,q,1):
					mol_in_group = True
			all_in_group = all_in_group and mol_in_group
		if all_in_group:
				valid_groups.append(i+1)
	f = open(directory+data_file+"_valid groups",'w')
	f.write("valid groups: \n") 
	for group in valid_groups:
		f.write(str(group)+"\n")
	print "\n\n\n"

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

#get_bioisosteres(file_1, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_1, noHs=True, brics=False, kennewell = True, overlap = True, test = False)
#get_bioisosteres(file_1, noHs=False, brics=True, kennewell=False, overlap = True, test = False)
#get_bioisosteres(file_1, noHs=False, brics=True, kennewell=False, overlap = False, test = False)
#get_bioisosteres(file_2, noHs=True, brics=True, kennewell=True, overlap = False, test = False)
#get_bioisosteres(file_1, noHs=False, brics=False, kennewell=False, overlap = True, test = False)
#get_bioisosteres(file_1, noHs=False, brics=False, kennewell=False, overlap = False, test = False)
#get_bioisosteres(file_2, noHs=False, brics=False, kennewell=False, overlap = True, test = False)
#get_bioisosteres(file_3, noHs=True, brics=False, kennewell=False, overlap = True, test = False)
#get_bioisosteres(file_2, noHs=False, brics=False, kennewell=True, overlap = True, test = False)
#get_bioisosteres(file_3, noHs=True, brics=False, kennewell=True, overlap = True, test = False)
t1 = []
t2 = []
for i in range(5):
	start1 = time.time()
	collect_bioisosteres_by_smiles(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
	t1.append(time.time()-start1)
	start2 = time.time()
	collect_bioisosteres(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
	t2.append(time.time()-start2)

print "smiles : " + str(t1)
print "non-smiles : " + str(t2)

print "smiles max = " + str(max(t1))
print "non-smiles max = " + str(max(t2))

print "smiles min = " + str(min(t1))
print "non-smiles min = " + str(min(t2))
#two_dim_similars(file_1, 0.7)
#two_dim_similars(file_2, 0.7)
#two_dim_similars(file_3, 0.7)
#two_dim_similars(file_4, 0.7)
#two_dim_similars(file_5, 0.7)
#two_dim_similars(file_6, 0.7)
#two_dim_similars(file_7, 0.7)
#two_dim_similars(file_8, 0.7)
#two_dim_similars(file_9, 0.7)
#two_dim_similars(file_10, 0.7)

## TODO:: ARE THE SMILES OR TANIMOTO EFFECTED BY THE DUMMY ATOM FROM THE FRAGMENTATION???
## TODO:: SPLIT OVER FILES RELATING TO TASK AND LEAVE ONE TEST FILE TO TEST CODE
