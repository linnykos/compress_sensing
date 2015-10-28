#remove residue
sh clean.sh

#after moving to directory, compile all necessary files
make -f makefile_el1.u
make -f makefile_el1sparse.u
make -f makefile_el1_size.u
make -f makefile_el1sparse_size.u
make -f makefile_print_el1.u
make -f makefile_print_el1sparse.u
make -f makefile_print_el1_size.u
make -f makefile_print_el1sparse_size.u


#start with simulation set for changing sparsity
sh run_el1.sh #simplex, generates el1.csv
sh run_el1sparse.sh #simplex KCS, generates el1sparse.csv

sh run_loqo_el1.sh #generates res_loqo_*_el1.out files
sh run_loqo_el1sparse.sh #generates res_loqo_*_el1sparse.out files

matlab -r run_fhtp #generates fhtp_res.csv
matlab -r run_fhtp_agnostic #generates fhtp_agnostic_res.csv
matlab -r run_l1ls #generats l1ls_res.csv and l2diff.csv
matlab -r run_MirrorProx #requires l2diff.csv, generates mirrorprox_res.csv

#####################################################
#####################################################

#next do the simulation set for changing dimension
sh run_el1_size.sh 
sh run_el1sparse_size.sh

sh run_loqo_el1_size.sh 
sh run_loqo_el1sparse_size.sh

matlab -r run_fhtp_size
matlab -r run_l1ls_size
matlab -r run_MirrorProx_size

#####################################################
#####################################################

R CMD BATCH read_ampl.R
R CMD BATCH process_loqo.R
R CMD BATCH process_loqo_size.R
R CMD BATCH rmUnbounded.R
R CMD BATCH rmUnbounded_size.R
R CMD BATCH plotter_xtable.R
R CMD BATCH plotter_xtable_size.R
matlab -r plotres
matlab -r plotres_size
