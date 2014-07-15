## test of dynamic names


i=0
declare prefix_$t=$i.txt

varname=prefix_$t

echo ${!varname}

i=$((i+1))

declare prefix_$t=$i.txt

echo ${!varname}
