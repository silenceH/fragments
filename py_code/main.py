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

#get_bioisosteres(file_1, noHs=True, brics=True, kennewell = True, overlap = False, test = False)
#collect_bioisosteres(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
#collect_bioisosteres_by_smiles(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
#get_bioisosteres(file_2, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_3, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_4, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_5, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_6, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_7, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_8, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_9, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_10, noHs=True, brics=True, kennewell = True, overlap = False, test = True)
#get_bioisosteres(file_1, noHs=True, brics=False, kennewell = True, overlap = True, test = False)
#get_bioisosteres(file_1, noHs=False, brics=True, kennewell=False, overlap = True, test = False)
#get_bioisosteres(file_1, noHs=False, brics=True, kennewell=False, overlap = False, test = False)
#get_bioisosteres(file_2, noHs=True, brics=True, kennewell=True, overlap = False, test = False)
#get_bioisosteres(file_1, noHs=False, brics=False, kennewell=False, overlap = True, test = False)
#get_bioisosteres(file_1, noHs=False, brics=False, kennewell=False, overlap = False, test = False)
#get_bioisosteres(file_2, noHs=False, brics=False, kennewell=False, overlap = True, test = False)
#get_bioisosteres(file_3, noHs=True, brics=False, kennewell=False, overlap = True, test = False)
#get_bioisosteres(file_2, noHs=False, brics=False, kennewell=True, overlap = True, test = False)
#get_bioisosteres(file_3, noHs=True, brics=False, kennewell=True, overlap = True, test = False)
t1 = []
t2 = []
for i in range(5):
	start1 = time.time()
	collect_bioisosteres_by_smiles(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
	t1.append(time.time()-start1)
	start2 = time.time()
	collect_bioisosteres(file_1,file_2,file_3,file_4,file_5,file_6,file_7,file_8,file_9,file_10)
	t2.append(time.time()-start2)

print "smiles : " + str(t1)
print "non-smiles : " + str(t2)

print "smiles max = " + str(max(t1))
print "non-smiles max = " + str(max(t2))

print "smiles min = " + str(min(t1))
print "non-smiles min = " + str(min(t2))
#two_dim_similars(file_1, 0.7)
#two_dim_similars(file_2, 0.7)
#two_dim_similars(file_3, 0.7)
#two_dim_similars(file_4, 0.7)
#two_dim_similars(file_5, 0.7)
#two_dim_similars(file_6, 0.7)
#two_dim_similars(file_7, 0.7)
#two_dim_similars(file_8, 0.7)
#two_dim_similars(file_9, 0.7)
#two_dim_similars(file_10, 0.7)

## TODO:: ARE THE SMILES OR TANIMOTO EFFECTED BY THE DUMMY ATOM FROM THE FRAGMENTATION???
## TODO:: SPLIT OVER FILES RELATING TO TASK AND LEAVE ONE TEST FILE TO TEST CODE
