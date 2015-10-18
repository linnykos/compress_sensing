for num in 102 202 302 402
do 
	for j in 0 1 2 3 4
	do
		rm -f loqo_res_size_${num}_${j}_el1sparse.out
		./printfiles_el1sparse_size $num $j
		time ampl loqo_el1sparse_size.mod > loqo_res_size_${num}_${j}_el1sparse.out
	done
done
