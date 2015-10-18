for num in 2 10 20 30 40 50 60 70 80 90 100 125 150
do
	for j in 0 1
	do
		rm -f loqo_res_${num}_${j}_el1sparse.out
		./printfiles_el1sparse $num $j
		time ampl loqo_el1sparse.mod > loqo_res_${num}_${j}_el1sparse.out
	done
done
