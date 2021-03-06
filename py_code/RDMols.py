import os
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, BRICS, Draw, rdShapeHelpers, Descriptors
from Fragment import * 

def get_mols_from_sdf_file(data_file, noHs=True,debug=False):
	## Method to return molecules for a given sdf data file
	try:
		data = os.environ['DATA']		## get data env
		if debug:
			print "found data environment: " + data
	except KeyError:
		print "cannot find data environment variable"

	if noHs:
		suppl = Chem.SDMolSupplier(data + '/validation_overlays/'+data_file+'.sdf')
	else: 
		suppl = Chem.SDMolSupplier(data + '/validation_overlays/'+data_file+'.sdf',removeHs=False)
	return [x for x in suppl if x is not None]

def get_fragments(mol,brics,data_file=""):
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

def write_mols_to_file(group_list,title,directory):
	## write to file
	for i in group_list:
		w = Chem.SDWriter(directory+title+str(group_list.index(i))+'.sdf')
		for mol in i.group: 
			w.write(mol.frag)
		w.flush()

def draw_mols_to_png(group_list,title,directory):
	## test for Fragment object
	try: 
		mols = [[mol.frag for mol in frag_group.group] for frag_group in group_list]
	except AttributeError:
		mols = [[mol.frag for mol in frag_group] for frag_group in group_list]
	## draw image
	for i in mols:
		for mol in i:
			tmp = AllChem.Compute2DCoords(mol)
		img = Draw.MolsToGridImage(i,legends=[str(x+1) for x in range(len(i))])
		img.save(directory+title+str(mols.index(i))+'.png')
