from math import sqrt, exp
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

	def get_all_coords(self):
			## return a list of coords for the 3D shape of a mol 
		if self.coords is None:
			# define a constructor of position objects
			tmp = (self.frag.GetConformer().GetAtomPosition(i) for i in range(self.frag.GetNumAtoms()))
			
			# return the list of the coordinates from those positions
			self.coords = [(atom.x,atom.y,atom.z) for atom in tmp]

	def score_pairs_kennewell(self,mol2):
		## get number of atoms for each molecule
		atoms_ref = self.frag.GetNumAtoms()
		atoms_f = mol2.frag.GetNumAtoms()
		## get coordinates for the reference molecule
		self.get_all_coords()
		mol2.get_all_coords()
		section_score = []
		for section_atom in self.coords:
			dist = [sqrt(pow((section_atom[0]-x[0]),2)+pow((section_atom[1]-x[1]),2)+pow((section_atom[2]-x[2]),2)) for x in mol2.coords]
			section_score.append(sum([exp(-pow(d,2)) for d in dist]))
		av_score = sum(section_score)*(2./(atoms_f+atoms_ref))
		if av_score>0.7:
			return True 
		else: 
			return False
