## KENNEWELL DEV
import os
from math import sqrt, exp
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, BRICS, Draw, rdShapeHelpers, PyMol, Descriptors


def atom_coords(molecule,atom_no):
	## returns a tuple of atom coordinates
	pos = molecule.GetAtomPosition(atom_no)
	return (pos.x,pos.y,pos.z)

def get_distance(coords1,coords2):
	## returns Euclidean distance between two 3D coordinates
	return sqrt(pow((coords1[0]-coords2[0]),2)+pow((coords1[1]-coords2[1]),2)+pow((coords1[2]-coords2[2]),2))

def are_similar(mol1,mol2,threshold):
	## returns False if Tanimoto similarity is greater than 0.8
	fp1 = AllChem.GetMorganFingerprint(mol1,2)
	fp2 = AllChem.GetMorganFingerprint(mol2,2)
	tanimoto = DataStructs.TanimotoSimilarity(fp1,fp2)
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
			w.write(mol)
		w.flush()

def draw_mols_to_png(mols,title,directory):
	## draw image
	for i in mols:
		for mol in i:
			tmp = AllChem.Compute2DCoords(mol)
		img = Draw.MolsToGridImage(i)
		img.save(directory+title+str(mols.index(i))+'.png')


### FUNCTION NOT USED ###
def are_appropriate(mol1,mol2):
	## check to make sure that the fragments are not similar and obey Ro3
	flag = True
	if Descriptors.MolWt(mol1) > 300. or Descriptors.MolWt(mol2) > 300.:
		flag = False
	return flag and not are_similar(mol1,mol2,0.8)
##########################

def get_all_coords(mol):
	## return a list of coords for the 3D shape of a mol 
	# define a constructor of position objects
	tmp = (mol.GetConformer().GetAtomPosition(i) for i in range(mol.GetNumAtoms()))
	
	# return the list of the coordinates from those positions
	return [(atom.x,atom.y,atom.z) for atom in tmp]

def get_fragments(mol,brics):
	## returns a list of non overlapping fragments.
	## if brics == True then returns the BRICS fragments
	if brics: 
		return Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol),asMols=True)	## BRICS bonds
	else: 
		## bond_smarts = Chem.MolFromSmarts("[*]!@!#!=[*]")	## single bonds
		## bond_smarts = Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]') ## simple rotatable bond smarts
		bond_smarts = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]') ## rotatable bonds
		bonds = mol.GetSubstructMatches(bond_smarts)
		bonds = [((x,y),(0,0)) for x,y in bonds]
		return Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol,bonds=bonds),asMols=True)

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

def get_first_fragment_from_smarts(mol,smarts):
	## fragments a molecule at first occurence of smarts 
	## note that smarts is a RDMol Object
	if mol.GetSubstructMatches(smarts):
		bonds = mol.GetSubstructMatches(smarts)
		bonds = [((x,y),(0,0)) for x,y in bonds]
		frags = list(Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol,bonds=[bonds[0]]),asMols=True))
		return frags

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
	score = rdShapeHelpers.ShapeTanimotoDist(atom1,atom2)
	if score < 0.3: 
		return True
	else:
		return False

def score_pairs_kennewell(mol1,mol2):
	## get number of atoms for each molecule
	atoms_ref = mol1.GetNumAtoms()
	atoms_f = mol2.GetNumAtoms()
	## get coordinates for the reference molecule
	ref_atoms = get_all_coords(mol1)
	section_score = []
	######################
	## visual debugging ##
	#v = PyMol.MolViewer()
	#v.DeleteAll()
	#v.ShowMol(mol1)
	#v.ShowMol(mol2,showOnly=False,name='other')
	######################
	## create a list of all the distances from atom xyz to the fragment atoms
	for section_atom in ref_atoms:
		dist = [get_distance(section_atom,frag_atom) for frag_atom in get_all_coords(mol2)]
		section_score.append(sum([exp(-pow(d,2)) for d in dist]))
	av_score = sum(section_score)*(2./(atoms_f+atoms_ref))
	if av_score>0.7:
		return True 
	else: 
		return False

def get_bioisosteres(data_file,noHs,brics, kennewell,overlap):
	## gets bioisosteric pairs from an sdf file
	## data_file : string of data file without .sdf
	## noHs : specifies whether Hs are to be removed (boolean)
	## brics : specifies whether brics bonds are broken at fragmentation (boolean)
	## kennewell : specifies whether kennewell scoring is to be used (boolean)
	## overlap : specifies whether fragments should be overlapping (boolean)

	## import mols from data_file 
	if noHs:
		suppl = Chem.SDMolSupplier('../data/validation_overlays/'+data_file+'.sdf')
	else: 
		suppl = Chem.SDMolSupplier('../data/validation_overlays/'+data_file+'.sdf',removeHs=False)

	mols = [x for x in suppl if x is not None]
	candidate_pairs = []
	grouped = []
	count = 0
	## for a reference molecule mol in mols
	for mol in mols:
		## set query data set to be the molecules that are not the reference
		query_set = [x for x in mols if x is not mol]
		## split the reference molecule by non ring single bonds
		ref_frag = get_fragments(mol, brics)
		## for each remaining ligand make it the query ligand
		for q in query_set:
			if overlap:
				frags = get_overlapping_fragments(q)
			else:
				frags = get_fragments(q, brics) 
			for s in ref_frag:
				section_pairs = [s]
				for f in frags:
					candidate = True
					if kennewell:  
						candidate = score_pairs_kennewell(s,f) 	## Kennewell scores
					else: 
						candidate = score_pairs_TD(s,f)		## Tanimoto distance
					if candidate and not are_similar(s,f,1.0): 
						section_pairs.append(f)
						count += 1
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
	
	# write groups to file
	write_mols_to_file(grouped,'group',directory)
	write_mols_to_file(final_group,'final_group',directory)
	
	# draw mols to png
	draw_mols_to_png(grouped,'group',directory)
	draw_mols_to_png(final_group,'final_group',directory)
	
	print "Tanimoto cut off is 1."
	if overlap: print "Overlapping fragments."
	if noHs: print "Hydrogen's removed."
	if kennewell: print "Kennewell scoring"
	print "Number of groups: " + str(len(grouped))
	print "Number of sections: " + str(len(candidate_pairs))					
	print "Total pairs : " + str(count) + "\n"
	return final_group

def collect_bioisosteres(*args):
	final_collection = []
	collection = [get_bioisosteres(data,True,True,True,False) for data in args]
	print "comparing..."
	compared = 0
	count = 0
	for i in range(len(collection)-1):
		grouped = collection[i]
		compare = collection[i+1:]
		for compare_group in compare:
			## obtain groups of molecules to be compared
			for reference_group,query_group in ((x,y) for x in grouped for y in compare_group):
				## obtain the molecules from the groups to be compared
				compared += 1
				for reference,query in ((m1,m2) for m1 in reference_group for m2 in query_group):
					if are_similar(reference,query,1.0):
						new_group = reference_group + query_group
						final_collection.append(new_group)
						count += 1
						break
	print "groups compared: " + str(compared)
	print "groups united: " + str(count)

	directory = '../test_output/compared_results' 
	try: 
		os.makedirs(directory)
		print "created new directory: " + directory
	except OSError:
		print directory + " already exists."	
	write_mols_to_file(final_collection,'final_collection',directory)
	draw_mols_to_png(final_collection,'final_collection',directory)

file_1 = 'P39900'
file_2 = 'P56817'
file_3 = 'O14757'

#get_bioisosteres(file_1, noHs=True, brics=True, kennewell = False, overlap = True)
#get_bioisosteres(file_1, noHs=True, brics=False, kennewell = True, overlap = True)
#get_bioisosteres(file_1, noHs=False, brics=True, kennewell=False, overlap = True)
#get_bioisosteres(file_1, noHs=False, brics=True, kennewell=False, overlap = False)
#get_bioisosteres(file_2, noHs=True, brics=True, kennewell=True, overlap = False)
#get_bioisosteres(file_1, noHs=False, brics=False, kennewell=False, overlap = True)
#get_bioisosteres(file_1, noHs=False, brics=False, kennewell=False, overlap = False)
##get_bioisosteres(file_2, noHs=False, brics=False, kennewell=False, overlap = True)
##get_bioisosteres(file_3, noHs=True, brics=False, kennewell=False, overlap = True)
##get_bioisosteres(file_2, noHs=False, brics=False, kennewell=True, overlap = True)
##get_bioisosteres(file_3, noHs=True, brics=False, kennewell=True, overlap = True)
collect_bioisosteres(file_1,file_2,file_3,'P12758')

