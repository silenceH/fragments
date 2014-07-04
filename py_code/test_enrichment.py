import enrichment 
from Bioisosteres import *
import time

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

def test_6():
	## test to see whether get pairs from group is well behaved 

	target_pairs = get_bioisosteres("P39900",return_pairs=True)
	target_groups = get_bioisosteres("P39900")

	## test to see whether a two member group returns two pairs
	two_frag_group = target_pairs[0]

	first_pairs = two_frag_group.get_pairs_from_group()

	part_1 = len(first_pairs) == 2

	if not part_1:
		print "test 6 : group with two fragments does not return two pairs"

	part_2 = first_pairs[0] != first_pairs[1]

	if not part_2:
		print "test 6 : pairs are not different"

	second_group = target_groups[0]
	second_group_size = second_group.size()
	second_pairs = second_group.get_pairs_from_group()

	part_3 = len(second_pairs) == (second_group_size**2) - second_group_size

	if not part_3:
		print "test 6 : pairs from larger group and not double the size"

	return part_1 and part_2 and part_3


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

if test_6():
	print "test_6 passed"

#enrichment.pair_enrichment(0,"P39900","P39900","Q00511","P00918")
#enrichment.group_enrichment("P39900","P39900","Q00511","P00918")
import os
import fnmatch
try:
	data = os.environ['DATA']		## get data env
except KeyError:
	print "cannot find data environment variable"
data_dir = data + '/validation_overlays/'
files = []
for file in os.listdir(data_dir):
	if fnmatch.fnmatch(file, '*.sdf'):
		files.append(file)

targets = [target[:-4] for target in files]


#enrichment.group_enrichment("P39900",targets)
t0 = time.time()

frags,scores =  enrichment.rank_targets(targets[0])
assert len(frags) == len(scores)
print "all frags ranked from one target " + str(len(scores))

frags,scores =  enrichment.rank_targets(targets[:2])
assert len(frags) == len(scores)
print "all frags ranked from two targets " + str(len(scores))

frags,scores =  enrichment.rank_targets(targets[:3])
assert len(frags) == len(scores)
print "all frags ranked from three targets " + str(len(scores))

frags,scores =  enrichment.rank_targets(targets[:5])
assert len(frags) == len(scores)
print "all frags ranked from two targets " + str(len(scores))
t1= time.time()

print "old way took " + str(t1-t0)
#frags,scores =  enrichment.rank_targets(targets[:len(targets)/2])
#assert len(frags) == len(scores)
#print "all frags ranked from half targets" + str(len(scores))

t3 = time.time()
frags,scores =  enrichment.rank_targets_test(targets[0])
assert len(frags) == len(scores)
print "all frags ranked from one target " + str(len(scores))

frags,scores =  enrichment.rank_targets_test(targets[:2])
assert len(frags) == len(scores)
print "all frags ranked from two targets " + str(len(scores))

frags,scores =  enrichment.rank_targets_test(targets[:3])
assert len(frags) == len(scores)
print "all frags ranked from three targets "  + str(len(scores))

frags,scores =  enrichment.rank_targets_test(targets[:5])
assert len(frags) == len(scores)
print "all frags ranked from two targets " + str(len(scores))
t4=time.time()

print "new method took " + str(t4-t3)

print "improvement = " + str((t4-t3)/(t1-t0)*100)

print "now testing with half database. Last time it took 5 hours" 

frags,scores =  enrichment.rank_targets_test(targets[:len(targets)/2])
assert len(frags) == len(scores)
print "all frags ranked from half targets" + str(len(scores))
