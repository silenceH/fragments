#!/usr/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from rdkit import DataStructs
from kennewell_dev import collect_bioisosteres_by_smiles
import os

def av_sim(group):
	
	sim = [[DataStructs.TanimotoSimilarity(group[i].fp,group[j].fp) for i in range(len(group))] for j in range(len(group))]
 	#print sim
	sumT = 0;
	count = 0;
	for i in range(len(group)):
		for j in range(len(group)):
			sumT += sim[i][j]
			count += 1

	average = (sumT-len(group))/(count-len(group))
	return average

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

groups = collect_bioisosteres_by_smiles(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
#for group in groups:
#	print av_sim(group)
