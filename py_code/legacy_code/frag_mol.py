from kennewell_dev import get_mols_from_sdf_file
from rdkit import Chem
from rdkit.Chem import AllChem, BRICS, Draw, rdShapeHelpers, Descriptors


mols = get_mols_from_sdf_file('P39900',True)

frag_list = [[x for x in Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol),asMols=True)] for mol in mols]	## BRICS bonds

for i in xrange(len(frag_list)):
	w = Chem.SDWriter(str(i)+'.sdf')
	for mol in frag_list[i]: 
		w.write(mol)
	w.flush()

