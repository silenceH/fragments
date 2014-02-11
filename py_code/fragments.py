import os
from rdkit import Chem


def getFragments(sdf):
	## get a set of the fragments for each file in the list
	return fragments

def fragmentMol(mol):
	## for a given mol apply a fragmentation scheme
	## returns a set of fragments using BRICS
	smiles = Chem.MolToSmiles(mol)
	frags = BRICS.BRICSDecompose(smiles)
	return frags

directory = "$DATA/"

for f in directory:
	##	newFrags = getFragments(f)
	## union of fragment sets

## TODO : implement as a class called fragments
