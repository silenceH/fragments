"""
	This programme takes a target from the validation overlays and writes each fragment to an sdf file
"""
#general 
import os 

# rdkit libraries
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, BRICS, Draw, rdShapeHelpers, Descriptors

# local modules
import Bioisosteres  

# define target
target = 'P39900'

# get fragment list
fragments = Bioisosteres.get_unique_fragments_from_files(target)

# write each fragment to sdf file
try:
	data = os.environ['DATA']
except KeyError:
	print "data env not found"

directory = data + 'fragment_files/' + target + '/'

try:
	os.makedirs(directory)
except OSError:
	print directory + " already exists."
	
for i in xrange(len(fragments)):
	f = Chem.SDWriter(directory+'f'+str(i)+'.sdf')
	f.write(fragments[i].frag)
	f.flush()
