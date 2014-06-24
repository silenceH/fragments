from Bioisosteres import *

## define return list of frags that have the query as a bioisosteric pair
def pairs_with_query(list_of_pairs, query):
	## test for group objects
	try: 
		## take pairs from Group data type 
		list_of_pairs = [x.group for x in list_of_pairs]
	except AttributeError:
		list_of_pairs = list_of_pairs
	## get the first pair if the second is equal to the query
	final_list_a = [x[0] for x in list_of_pairs if x[1].are_similar(query,1.0)]	
	final_list_b = [x[1] for x in list_of_pairs if x[0].are_similar(query,1.0)]	

	## combine two lists
	final_list = final_list_a + final_list_b
	return final_list

def score_similarity(query,list_of_fragments):
	## returna  list of all the fragments in a set of overlays
	similarity_scores = [x.tanimoto_score(query) for x in list_of_fragments]
	return similarity_scores

def rank_list_of_fragments(query,list_of_fragments):
	## rank a list of fragments by the similarity to a query
	scores = score_similarity(query,list_of_fragments)
	changed = True
	while changed:
		changed = False
		for i in range(len(list_of_fragments)-1):
			if scores[i] < scores[i+1]:
				scores[i], scores[i+1] = scores[i+1], scores[i]
				list_of_fragments[i], list_of_fragments[i+1] = list_of_fragments[i+1], list_of_fragments[i]
				changed = True
	return list_of_fragments
				
