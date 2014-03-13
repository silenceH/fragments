#!/usr/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit import DataStructs
from rdkit.Chem import BRICS, Recap

def get_mols (data_file): 

	suppl = Chem.SDMolSupplier(data_file)

	ms = [x for x in suppl if x is not None]


	return ms


def decompose_mol(mols):
	num_recap = 0
	num_brics = 0
	num_same = 0
	for n in range(len(mols)):
		mol = mols[n]
		brics = list(BRICS.BRICSDecompose(mol))
	 	recap = Recap.RecapDecompose(mol)

		print str(n)
	 	print "BRICS: " + str(len(brics))
	 	print "Recap: " + str(len(recap.GetAllChildren())) + "\n"
	 	
		if len(brics) > len(recap.GetAllChildren()):
			num_brics += 1
		elif  len(brics) == len(recap.GetAllChildren()):
			num_same += 1
		else:
			num_recap += 1

	print "TOTAL GREATER"
	print "BRICS : " + str(num_brics)
	print "RECAP: " + str(num_recap)
	print "SAME: " + str(num_same)

def general_info (mols):
	num = len(mols)
	print "there are " + str(num) + " molecules in the file."
	fps = [AllChem.GetMorganFingerprint(m,2) for m in mols]

	sim = [[DataStructs.TanimotoSimilarity(fps[i],fps[j]) for i in range(len(fps))] for j in range(len(fps))]

	sumT = 0;
	count = 0;
	for i in range(len(fps)):
		for j in range(len(fps)):
			sumT += sim[i][j]
			count += 1

	average = (sumT-num)/(count-num)
	#print sumT
	print "largest Tanimoto is: " + str(max(max(filter(lambda x: x<1.0, row)) for row in sim))
	print "average tanimoto similarity is " + str(average)



data_files = ['O14757', 'Q00511', 'P39900'] ## ['rings', 'disimilarity', 'standard']

for data in data_files: 
	mols = get_mols('../data/validation_overlays/' + data + '.sdf')
	print data
	general_info(mols)
	decompose_mol(mols)
	print "\n"
