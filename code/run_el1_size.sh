rm -f el1_size.csv
for i in 102 122 142 162 182 202 222 242 262 282 302 322 342 362 382 402
do
	for j in 0 1 2 3 4 5 6 7 8 9
	do
	./el1_size $i $j >> el1_size.csv
	done
done
