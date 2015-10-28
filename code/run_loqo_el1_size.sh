for num in 102 202 302 402
do
	for j in 0 1 2 3 4 5 6 7 8 9
	do 
		rm -f loqo_res_size_${num}_${j}_el1.out
		./printfiles_el1_size $num $j
		time ampl loqo_el1_size.mod > loqo_res_size_${num}_${j}_el1.out
	done
done
