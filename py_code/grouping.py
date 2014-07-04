from Bioisosteres import *

def get_groups_1(data_file,noHs=True,return_pairs=False,debug=False):
	## gets bioisosteric pairs from an sdf file
	## data_file : string of data file without .sdf
	## noHs : specifies whether Hs are to be removed (boolean)

	## import mols from data_file 
	mols = get_mols_from_sdf_file(data_file,noHs)
	# list of candidate bioisosteric pairs to view individual pairs if needed commented out
	candidate_pairs = []
	grouped = []
	pairs = []
	count = 0

	## create a list of fragment objects here
	mols = [get_fragments(mol,True,data_file) for mol in mols]	
	## for a reference molecule mol in mols
	for i in range(len(mols)-1):
		## set query data set to be the molecules that are not the reference
		ref_frag = mols[i] 			# a list of Fragment objects
		query_set = mols[i+1:]
		## for each remaining ligand make it the query ligand
		for q in query_set:
			frags = q
			for s in ref_frag:
				# create fragment group
				section_group = Group(s)
				for f in frags:
					candidate = True
					## note we are using the Kennewell score
					candidate = s.score_pairs_kennewell(f) 	## Kennewell scores
					if candidate and not s.are_similar(f,1.0): 
						# For dedugging purposes print details
						if debug:
							print "pair: " + str(count+1)
							print "Ref: " + str(i) + "\tfrag: " + str(ref_frag.index(s))
							print "Query: " + str(i+query_set.index(q)) + "\tfrag: " + str(frags.index(f))
						# code for returning a list of pairs
						if return_pairs:
							pair = Group(s)
							pair.add(f)
							candidate_pairs.append(pair)

						section_group.add(f)
						count += 1
				if section_group.size() > 1:
					if debug:
						print "size of section group: " + str(section_group.size())
					in_grouped = False
					for x in grouped:
						if section_group.get_mol(0) in x.group:
							in_grouped = True
							x.merge(section_group) 
					if not in_grouped:
						grouped.append(section_group)

	print "SUMMARY STATISTICS FOR 1"
	print "Number of fragments: " + str(sum([len(m) for m in mols])) 
	print"Number of identified pairs : " + str(count) 
	print"Number of overlaid sections: " + str(len(grouped))
	print"group \tfrags \tav similarity\n"
##	for i in range(len(grouped)):
##		av_sim = grouped[i].get_av_similarity()
##		print(str(i+1) + "\t" + str(grouped[i].size()) + "\t" + str(av_sim)+"\n")
##	if debug:
##		print "no pairs: " + str(count)
##	if return_pairs:
##		if debug:
##			if len(candidate_pairs) == count:
##				print "PAIRS WORKED"
##		return candidate_pairs
	return grouped

def get_groups_2(data_file,noHs=True,return_pairs=False,debug=False):
	## gets bioisosteric pairs from an sdf file
	## data_file : string of data file without .sdf
	## noHs : specifies whether Hs are to be removed (boolean)

	## import mols from data_file 
	mols = get_mols_from_sdf_file(data_file,noHs)
	# list of candidate bioisosteric pairs to view individual pairs if needed commented out
	candidate_pairs = []
	grouped = []
	pairs = []
	count = 0

	## create a list of fragment objects here
	mols = [get_fragments(mol,True,data_file) for mol in mols]	
	## for a reference molecule mol in mols
	for i in range(len(mols)):
		## set query data set to be the molecules that are not the reference
		ref_frag = mols[i] 			# a list of Fragment objects
		query_set = mols[i+1:]
		## query_set = [x for x in mols if x is not ref_frag]
		## for each remaining ligand make it the query ligand
		for q in query_set:
			frags = q
			for s in ref_frag:
				# create fragment group
				section_group = Group(s)
				for f in frags:
					candidate = True
					## note we are using the Kennewell score
					candidate = s.score_pairs_kennewell(f) 	## Kennewell scores
					if candidate and not s.are_similar(f,1.0): 
						# For dedugging purposes print details
						if debug:
							print "pair: " + str(count+1)
							print "Ref: " + str(i) + "\tfrag: " + str(ref_frag.index(s))
							print "Query: " + str(i+query_set.index(q)) + "\tfrag: " + str(frags.index(f))
						# code for returning a list of pairs
						if return_pairs:
							pair = Group(s)
							pair.add(f)
							candidate_pairs.append(pair)

						section_group.add(f)
						count += 1
				if section_group.size() > 1:
					if debug:
						print "size of section group: " + str(section_group.size())
					in_grouped = False
					for x in grouped:
						for bio in section_group.group:
							if bio in x.group:
								in_grouped = True
								x.merge(section_group) 
					if not in_grouped:
						grouped.append(section_group)

	print "SUMMARY STATISTICS FOR 2" 
	print "Number of fragments: " + str(sum([len(m) for m in mols])) 
	print"Number of identified pairs : " + str(count) 
	print"Number of overlaid sections: " + str(len(grouped))
	print"group \tfrags \tav similarity\n"
	for i in range(len(grouped)):
		av_sim = grouped[i].get_av_similarity()
		print(str(i+1) + "\t" + str(grouped[i].size()) + "\t" + str(av_sim)+"\n")
	if debug:
		print "no pairs: " + str(count)
	if return_pairs:
		if debug:
			if len(candidate_pairs) == count:
				print "PAIRS WORKED"
		return candidate_pairs
	return grouped
