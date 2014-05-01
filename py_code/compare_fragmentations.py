from Bioisosteres import *
## TEST FRAGMENTATION CODE
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

# define methods 
def write_to_file(fragged_mol,title,directory):
	## write to file
	for i in fragged_mol:
		w = Chem.SDWriter(directory+title+str(fragged_mol.index(i))+'.sdf')
		for mol in i: 
			w.write(mol.frag)
		w.flush()

## get mols
mols = get_mols_from_sdf_file(file_1)

#new_frags = get_overlapping_fragments(mols[14])
#
#print "total overlapping fragments: " + str(len(new_frags))
#
#for mol in new_frags:
#	print mol.GetNumAtoms()

## get fragments

brics = [get_fragments(x,brics=True) for x in mols]
kennewell = [get_fragments(x,brics=False) for x in mols]
kennewell_overlap = [get_overlapping_fragments(x) for x in mols]

kennewell_overlap_frags = list(map((lambda x: Fragment(x,"")),kennewell_overlap))
kennewell_overlap_frags = [[Fragment(x,"") for x in mol] for mol in kennewell_overlap]

print "brics: " + str(sum([len(x) for x in brics]))
print "kennewell: " + str(sum([len(x) for x in kennewell]))
print "kennewell_overlap: " + str(sum([len(x) for x in kennewell_overlap]))

## write fragments from overlapping fragmentation schemes
directory = '/usr/users/people/matts/Dropbox/test_output/FragmentationSchemes/' 
try: 
	os.makedirs(directory)
	print "created new directory: " + directory
except OSError:
	print directory + " already exists."	
	
write_to_file(brics,"brics",directory)
write_to_file(kennewell,"kennewell",directory)
write_to_file(kennewell_overlap_frags,"kennewell_overlap",directory)

draw_mols_to_png(brics,"brics",directory)
draw_mols_to_png(kennewell,"kennewell",directory)
draw_mols_to_png(kennewell_overlap_frags,"kennewell_overlap_frags",directory)
