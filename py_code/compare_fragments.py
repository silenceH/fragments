## this is a programme for comparing some fragment functions

#!/usr/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import BRICS, Recap
import os

smarts_RB = "insert smarts"

def get_mols(data_file):
	suppl = Chem.SDMolSupplier('../data/validation_overlays/' + data_file + '.sdf')
	try:
		os.mkdir(data_file)
	except OSError:
		print data_file + " directory already exists."
	return [x for x in suppl if x is not None]

def compare_frags(mol, n, data_file):
	brics_frag = [AllChem.MolFromSmiles(frag) for frag in brics_frags(mol)]
	recap_frag = [AllChem.MolFromSmiles(frag) for frag in recap_frags(mol)]
	rb_frag = [AllChem.MolFromSmiles(frag) for frag in rb_frags(mol)]
	single_bond_frag_brics = [AllChem.MolFromSmiles(frag) for frag in single_bond_BRICS(mol)]
	
	brics_img = Draw.MolsToGridImage(brics_frag)
	recap_img = Draw.MolsToGridImage(recap_frag)
	rb_img = Draw.MolsToGridImage(rb_frag)
	single_bond_brics_img = Draw.MolsToGridImage(single_bond_frag_brics) 

	brics_img.save(data_file+"/brics_"+data_file+"_"+str(n)+".png")
	recap_img.save(data_file+"/recap"+data_file+"_"+str(n)+".png")
	rb_img.save(data_file+"/rb_"+data_file+"_"+str(n)+".png")
	single_bond_brics_img.save(data_file+"/single_bond_brics_"+data_file+"_"+str(n)+".png")

def brics_frags(mol):
	## for a given molecule return a list of BRICS fragments as SMILES strings
	frags = list(BRICS.BRICSDecompose(mol))
	return frags

def recap_frags(mol):
	frags = list(Recap.RecapDecompose(mol).GetAllChildren())
	return frags
	
def rb_frags(mol):
	rb_smarts = Chem.MolFromSmarts('[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]')

	#find the rotatable bonds
	bonds = mol.GetSubstructMatches(rb_smarts)

	#we can use BRICS package in this way
	bonds = [((x,y),(0,0)) for x,y in bonds]
	pieces = BRICS.BreakBRICSBonds(mol,bonds=bonds)
	return [Chem.MolToSmiles(x,True) for x in Chem.GetMolFrags(pieces,asMols=True)]

def single_bond_BRICS(mol):
	single_bond_smarts = Chem.MolFromSmarts("[*]!@!#!=[*]")
	bonds = mol.GetSubstructMatches(single_bond_smarts)
	bonds = [((x,y),(0,0)) for x,y in bonds]
	pieces = BRICS.BreakBRICSBonds(mol,bonds=bonds)
	return [Chem.MolToSmiles(x,True) for x in Chem.GetMolFrags(pieces,asMols=True)]
	
def single_bond_trad(mol):
	single_bond_smarts = Chem.MolFromSmarts("[*]!@!#!=[*]")
	# create list of tuples showing atoms with single bonds
	bonds = mol.GetSubstructMatches(single_bond_smarts)
	
	# create an editable molecule and break the bonds
	tmp = Chem.EditableMol(mol)



		

data_files = ['O14757', 'Q00511', 'P39900'] ## ['rings', 'disimilarity', 'standard']
df='O14757'

mols = get_mols(df)
compare_frags(mols[23],24,df)
compare_frags(mols[12],13,df)
compare_frags(mols[24],25,df)
compare_frags(mols[2],3,df)
compare_frags(mols[30],31,df)
mols = get_mols(data_files[2])
compare_frags(mols[16],17,data_files[2])
compare_frags(mols[4],5,data_files[2])
compare_frags(mols[13],14,data_files[2])
compare_frags(mols[2],3,data_files[2])
compare_frags(mols[16],17,data_files[2])
