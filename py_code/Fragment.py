from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, BRICS, Draw, rdShapeHelpers, Descriptors

class Fragment(object):
	##
	## define a fragment class that holds references to
	## different representations of the fragment.
	##
	## self : the ROMol fragment (boost::shared_ptr in c++)
	## fp : the RDKit Morgan 2 fingerprint for Tanimoto comparison
	## smiles : canonical smiles representation
	##
	def __init__(self, frag,ligand, coords = None, fp=None,smiles=None):
		self.frag = frag

		self.ligand = ligand

		if coords is None:
			self.coords = None
		
		if fp is None:
			self.fp = None
		self.fp = fp

		if smiles is None:
			self.smiles=None
		self.smiles = smiles

	def are_similar(self,frag2,threshold):
		## returns False if Tanimoto similarity is greater than threshold
		if self.fp is None:
			self.fp = AllChem.GetMorganFingerprint(self.frag,2)
		if frag2.fp is None:
			frag2.fp = AllChem.GetMorganFingerprint(frag2.frag,2)
		tanimoto = DataStructs.TanimotoSimilarity(self.fp,frag2.fp)
		if tanimoto >= threshold:
			return True 
		else: 
			return False

	def score_pairs_kennewell(self,mol2):
		## get number of atoms for each molecule
		frag1, frag2 = self.frag,mol2.frag
		atoms_ref = frag1.GetNumAtoms()
		atoms_f = frag2.GetNumAtoms()
		## get coordinates for the reference molecule
		ref_atoms = get_all_coords(frag1)
		section_score = []
		for section_atom in ref_atoms:
			dist = [get_distance(section_atom,frag_atom) for frag_atom in get_all_coords(frag2)]
			section_score.append(sum([exp(-pow(d,2)) for d in dist]))
		av_score = sum(section_score)*(2./(atoms_f+atoms_ref))
		if av_score>0.7:
			return True 
		else: 
			return False
