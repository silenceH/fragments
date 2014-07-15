## programme to extract all the tanimoto similarities in shape-it output

import os
import fnmatch

scores = []

class Pair:
	def __init__(self,pair,score):
		self.pair = pair
		self.score = score

	

for file in os.listdir('.'):

	## if file is an output file read results
	if fnmatch.fnmatch(file, '*.out'):
		f = open(file,'r')
		results = f.readlines()[1]
		first_tab = results.find('\t')
		second_tab = results.find('\t',first_tab+1)
		third_tab = results.find('\t',second_tab+1)

		score = float(results[second_tab:third_tab])
		name = results[:first_tab]

		new_pair = Pair(name,score)

		if score != 1:
			scores.append(new_pair)




scores.sort(key=lambda x: x.score, reverse=True)
print [x.pair for x in scores]


