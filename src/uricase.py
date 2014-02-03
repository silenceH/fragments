#!/usr/bin/python

from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
from rdkit import DataStructs

suppl = Chem.SDMolSupplier('../data/validation_overlays/P39900.sdf')

ms = [x for x in suppl if x is not None]

## for m in ms: tmp=Chem.Compute2DCoords(m)
## 
## img = Draw.MolsToGridImage(ms[:7], molsPerRow=4,legends=[m.GetProp("_Name") for m in ms])
## 
## img.save('uricase.png')

fps = [Chem.GetMorganFingerprint(m,2) for m in ms]

sim = [[DataStructs.TanimotoSimilarity(fps[i],fps[j]) for i in range(len(fps))] for j in range(len(fps))]

sumT = 0;
count = 0;
for i in range(len(fps)):
	for j in range(len(fps)):
		sumT += sim[i][j]
		count += 1

average = (sumT)/(count)
#print sumT
#print sim
print "average tanimoto similarity is " + str(average)

