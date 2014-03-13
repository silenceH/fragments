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
 	#print sim
	sumT = 0;
	count = 0;
	for i in range(len(fps)):
		for j in range(len(fps)):
			sumT += sim[i][j]
			count += 1

	average = (sumT-len(fps))/(count-len(fps))
	print "average:  " + str(average) + "\tnumber: " + str(len(fps)) 
	return average

p00730 = av_sim('P00730')
print "paper: " + str(0.588) + "\t number: " + str(8)
p27487 = av_sim('P27487')
print "paper: " + str(0.128) + "\t number: " + str(39)
p00760 = av_sim('P00760')
print "paper: " + str(0.317) + "\t number: " + str(22)
p35968 = av_sim('P35968')
print "paper: " + str(0.151) + "\t number: " + str(8)
q04771 = av_sim('Q04771')
print "paper : " + str(0.255) + "\t number: " + str(5)
