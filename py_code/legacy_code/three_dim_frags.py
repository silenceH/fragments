import os
from rdkit import Chem
from rdkit.Chem import AllChem, BRICS, Draw, PyMol

data_file = 'P39900'

def get_mols(data_file):
	## function that takes a file name and retuns a list of mols from the file
	suppl = Chem.SDMolSupplier('../data/validation_overlays/'+data_file+'.sdf',removeHs=False)
	return [x for x in suppl if x is not None]

def get_brics_fragments(mol):
	## function that returns a list of brics fragments as Mol objects
	brics = BRICS.BreakBRICSBonds(mol)
	return Chem.GetMolFrags(brics,asMols=True)

def show_frags(data_file, mol_number):
	## displays mol and frags in PyMol
	#os.system('pymol -R') 	## start PyMol as server
	v = PyMol.MolViewer() 	## create rdkit PyMol object
	
	## get mol and fragments
	mols = get_mols(data_file)
	mol = mols[mol_number]
	fragments = get_brics_fragments(mols[mol_number])

	## show mol
	v.ShowMol(mol,name=mol.GetProp('_Name'))

	## then show all frags on the screen
	for i in range(len(fragments)):
		v.ShowMol(fragments[i],showOnly=False,name='frag_'+str(i+1))



show_frags(data_file,0)
show_frags(data_file,1)


## TODO note that for efficiency, the mols should be produced outside of the show_frags function
