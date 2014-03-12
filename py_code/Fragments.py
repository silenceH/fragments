import os
from rdkit import Chem

class Fragment(object):
	##
	## define a fragment class that holds references to
	## different representations of the fragment.
	##
	## self : the ROMol fragment (boost::shared_ptr in c++)
	## fp : the RDKit Morgan 2 fingerprint for Tanimoto comparison
	## smiles : canonical smiles representation
	##
	def __init__(self,frag,fp=None,smiles=None):
		self.frag = frag

		if fp is None:
			self.fp = None
		self.fp = fp

		if smiles is None:
			self.smiles=None
		self.smiles = smiles

	
	## an instance method that sets the instance variables when needed	

