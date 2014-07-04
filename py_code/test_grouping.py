from Bioisosteres import * 
import grouping

## tests for grouping algorithm 

target = 'P39900'

def test_1():
	## test to see how many groups with algorithm 1 
	grouping.get_groups_1(target)
	grouping.get_groups_2(target)

test_1()


groups = get_bioisosteres(target)

count = 0
for grp in groups:
	for frg in grp.group:
		test = [x for x in grp.group if x is not frg]
		for frg2 in test:
			if not frg.score_pairs_kennewell(frg2):
				count += 1

print count
