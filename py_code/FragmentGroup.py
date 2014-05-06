from Fragment import *

class Group(object):
	## a group of fragments as a list in order to make the 
	## management of fragment groups easier
	def __init__(self,first_frag=None):
		if first_frag is None:
			self.group = []
		else:
			self.group = [first_frag]
	
	def add(self,frag):
		self.group.append(frag)

	def join(self,another):
		self.group.extend(another)

	def size(self):
		return len(self.group)

	def get_mol(self,index):
		return self.group[index]

	def merge(self,anotherGroup):
		q_group = anotherGroup.group
		add_group = []
		for mol in anotherGroup.group:
			unique = True
			for ref in self.group:
				if ref.are_similar(mol,1):
					unique = False
					break
			if unique:
				add_group.append(mol)
		self.group.extend(add_group)

	def show(self):
		return self.group

	def remove_2D_equivalents(self):
		## remove duplicate 2D mols
		u = []
		for mol in self.group:
			include = True
			for m in u:
				#if mol.are_similar(m,1.0):
				if mol == m:
					include = False
			if include:
				u.append(mol)
		self.group = u

	def get_av_similarity(self):
		total_sim = 0
		num_mols = self.size() 
		for i in range(num_mols):
			row = (DataStructs.TanimotoSimilarity(self.get_mol(i).fp,self.get_mol(j).fp) for j in range(num_mols) if j > i)
			total_sim += sum(row)
		return 2*total_sim/(num_mols**2 - num_mols)
