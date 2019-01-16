#!/bin/bash

#SBATCH -N 4
#SBATCH -n 32
#SBATCH -p n4c32

SRCDIR="/home/jaka/Extended-Hubbard/source_codeVV" #MODIFY THIS #SOURCE DIRECTORY
for l in 0.1 0.3 0.5 #eps
    do
    for k in 0.1 0.4 0.7 1. #t
        do
        for i in 2. 3. 4. #U
	        do
	        for j in 0.1 0.3 0.5 0.7 0.9 1.2 1.4 #V
		        do
		        #set and make work directory
		        WORKDIR="/home/jaka/Extended-Hubbard/DataXH/eps$l-t$k-U$i-V$j" #MODIFY THIS #HOME DIRECTORY
		        mkdir -p $WORKDIR

		        #files list
		        INFILE="input_main"
		        FILES="global.f90 allocation.f90 deallocation.f90 main.f90 fermi.f90 rtbis.f90 Optical_Response.f90 input_opt_response $INFILE"
		        
		        #copy file and change directory
		        cp $FILES $WORKDIR
		        cd $WORKDIR

                
		        #change input_file
		        str_e="s/EPSILON/$l/g"
		        sed -i -e $str_e input_main
		        str_t="s/HOPPING/$k/g"
		        sed -i -e $str_t input_main
		        str_u="s/COULOMB/$i/g"
		        sed -i -e $str_u input_main
		        str_v="s/INTERSITE/$j/g"
		        sed -i -e $str_v input_main
	        
		        #compile and run
		        mpif90 -g -o XH.exe global.f90 allocation.f90 deallocation.f90 main.f90 fermi.f90 rtbis.f90 Optical_Response.f90 -llapack

		        mpirun -np $SLURM_NTASKS ./XH.exe < $INFILE > logfile.dat
		        
		        #back to source_file
		        cd $SRCDIR
		         
		        done
	        done
        done
    done
