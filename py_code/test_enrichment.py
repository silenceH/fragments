import enrichment 
from Bioisosteres import *

## test to see whether the list finds all pairs that have the query 


def test_1():
	pairs = get_bioisosteres("P39900",return_pairs=True)
	query = pairs[0].group[0]
	test_list = enrichment.pairs_with_query(pairs,query)

	## test to see whether the list is the same as an exhaustive search

	count=0
	for pair in pairs:
		if query.are_similar(pair.group[0],1.0):
			count += 1 
		elif query.are_similar(pair.group[1],1.0):
			count += 1
	return count == len(test_list)

def test_2():
	## test to make sure that the get all fragments function works

	frags_1 = get_fragments_from_files("P39900")

	part_1 = len(frags_1) == 78

	if not part_1:
		print "test 2 : failed one file test"
	
	frags_2 = get_fragments_from_files("P39900","Q00511")

	part_2 = len(frags_2) == 88

	if not part_2:
		print "test 2 : failed two file test"
	
	return part_1 and part_2
	
def test_3():
	## test score similarity list
	frags = get_fragments_from_files("P39900")
	query = frags[10]
	scores = enrichment.score_similarity(query,frags)
	part_1 = scores[10] == 1.0

	if not part_1:
		print "test 3 : cannot find sel similarity"
	
	part_2 = query.tanimoto_score(frags[4]) == scores[4]

	if not part_2:
		print "test 3 : cannot find arbitrary score"
	
	return part_1 and part_2

def test_4():
	frags = get_fragments_from_files("P39900")
	frags_2 = get_unique_fragments_from_files("P39900")
	query = frags[10]
	import copy
	copy_frags = copy.deepcopy(frags)
	new_frags = enrichment.rank_list_of_fragments(query,frags)
	new_frags_2 = enrichment.rank_list_of_fragments(query,frags_2)

	part_1 = new_frags[0].are_similar(query,1.0)
	if not part_1:
		print "test 4 : most similar is not self similar"

	part_2 = copy_frags != new_frags
	if not part_2:
		print "test 4 : list has not changed"

	part_3 = new_frags_2[0].are_similar(query,1.0)
	if not part_3:
		print "test 4 : unique fragments not self-similar"

	return part_1 and part_2 and part_3
	
def test_5():
	## test to see whether unique frags works
	frags_1 = get_unique_fragments_from_files("P39900")
	frags_2 = get_fragments_from_files("P39900")

	## see whether any have been removed
	part_1 = len(frags_1) < len(frags_2)

	if not part_1:
		print "test 5 : lists are still the same length"
	
	## compare to group de-depulicating method
	group_example = Group()
	for i in frags_2:
		group_example.add(i)
	
	group_example.remove_2D_equivalents()

	part_2 = len(frags_1) == group_example.size()
	if not part_2:
		print "test 5 : group method and this method return different length lists"
	
	return part_1 and part_2



####################################################################
##                        MAIN                                    ##
####################################################################

if test_1():
	print "test_1 passed" 

if test_2():
	print "test_2 passed" 

if test_3():
	print "test_3 passed" 

if test_4():
	print "test_4 passed"

if test_5():
	print "test_5 passed"

enrichment.pair_enrichment(0,"P39900","P39900","Q00511","P00918")
