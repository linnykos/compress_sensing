for num in 2 20 50 70 100 150
do
	for j in 0 1 2 3 4 5 6 7 8 9
	do
		rm -f loqo_res_${num}_${j}_el1.out
		./printfiles_el1 $num $j
		time ampl loqo_el1.mod > loqo_res_${num}_${j}_el1.out
	done
done
