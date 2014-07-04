import enrichment 
from Bioisosteres import *



all_fragments = enrichment.get_unique_fragments_from_files("P39900")

frags_1,scores_1 = enrichment.rank_all_frags(all_fragments)
frags_2,scores_2 = enrichment.ranked_list_of_all_fragments(all_fragments)

frags_3,scores_3 = enrichment.mergeSorted([3,2,1],[3,2,1],[4,3,2],[4,3,2])

print scores_3

print scores_1[-1]
print scores_2[-1]

print frags_1[0]
print frags_2[0]

