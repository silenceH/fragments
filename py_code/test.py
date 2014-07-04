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

#get_bioisosteres(file_1,return_pairs=True)

#collect_bioisosteres_greedy(True,file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)

print get_bioisosteres(file_1)
#print get_bioisosteres(file_2,count_pairs=True)
#print get_bioisosteres(file_3,count_pairs=True)
#print get_bioisosteres(file_4,count_pairs=True)
#print get_bioisosteres(file_5,count_pairs=True)
#print get_bioisosteres(file_6,count_pairs=True)
#print get_bioisosteres(file_7,count_pairs=True)
#print get_bioisosteres(file_8,count_pairs=True)
#print get_bioisosteres(file_9,count_pairs=True)
#print get_bioisosteres(file_10,count_pairs=True)
#
#collect_bioisosteres_greedy(True,file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
coll = collect_bioisosteres_greedy(False,file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
#get_pair_frequency(False,"test1",file_1,file_1)

frags = coll[0].group[30:39]
mols = [mol.frag for mol in frags]
## draw image
for i in mols:
	tmp = AllChem.Compute2DCoords(i)
img = Draw.MolsToGridImage(mols,legends=[str(x+31) for x in range(len(mols))])
img.save("illustration"+str(mols.index(i))+'.png')
