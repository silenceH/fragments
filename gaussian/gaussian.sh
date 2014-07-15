# write overlap scores of each fragment pair 

FILES=$DATA/fragment_files/P39900/*

count=0



for i in $FILES; do 
	for j in $FILES; do
		declare prefix_$suffix=$count.out
		score=prefix_$suffix
		shape-it -r $i -d $j -s ${!score}
		count=$((count+1))
		
	done
done
