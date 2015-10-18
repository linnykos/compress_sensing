rm -f el1sparse.csv
for i in 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 62 64 66 68 70 72 74 76 78 80 82 84 86 88 90 92 94 96 98 100 125 150
do
	for j in 0 1 2 3 4 5 6 7 8 9
	do
		./el1sparse $i $j >> el1sparse.csv
	done
done
