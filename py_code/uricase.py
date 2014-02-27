#!/usr/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from rdkit import DataStructs
import os

def av_sim(data_file):
	try:
		data = os.environ['DATA']		## get data env
		print "found data environment: " + data
	except KeyError:
		print "cannot find data environment variable"

	suppl = Chem.SDMolSupplier(data+'/validation_overlays/' + data_file + '.sdf')

	ms = [x for x in suppl if x is not None]


	fps = [Chem.GetMorganFingerprint(m,2) for m in ms]

	sim = [[DataStructs.TanimotoSimilarity(fps[i],fps[j]) for i in range(len(fps))] for j in range(len(fps))]

	sumT = 0;
	count = 0;
	for i in range(len(fps)):
		for j in range(len(fps)):
			sumT += sim[i][j]
			count += 1

	average = (sumT)/(count)
	print "average tanimoto similarity is " + str(average)
	return average

p00730 = av_sim('P00730')
p42574 = av_sim('P42574')
print "my attempt : " + str(p42574/p00730)

print "paper : " + str(0.440/0.588)
