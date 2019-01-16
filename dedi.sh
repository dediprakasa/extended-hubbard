#!/bin/bash

#number of processors
np=4

SRCDIR="/media/deprak/743C8CAD0AC915D3/TCMP/Get_It_Done/Extended-Hubbard/source_code" #MODIFY THIS #SOURCE DIRECTORY
for i in 1. 2. 3. 4. 5. #MODIFY THIS #VARIASI U (JANGAN LUPA PAKE TITIK DI BELAKANGNYA)
	do
	for j in 0. #0.2 0.4 0.6 0.8 1. 1.2 1.4 1.6 1.8 2. 2.2 #MODIFY THIS #VARIASI TEMPERATURE
		do
		#set and make work directory
		WORKDIR="/media/deprak/743C8CAD0AC915D3/TCMP/Get_It_Done/Extended-Hubbard/DataXH/U$i-V$j" #MODIFY THIS #HOME DIRECTORY
		mkdir -p $WORKDIR

		#files list
		INFILE="input_main"
		FILES="global.f90 allocation.f90 deallocation.f90 main.f90 fermi.f90 rtbis.f90 $INFILE"
		
		#copy file and change directory
		cp $FILES $WORKDIR
		cd $WORKDIR
        
        if [ $i -eq 1. ]
        then
            wmin=-7.5
            wmax=8.5
    		str_wmin="s/WMIN/$wmin/g"
		    sed -i -e $str_wmin input_main
		    str_wmax="s/WMAX/$wmax/g"
		    sed -i -e $str_wmax input_main
		elif [ $i -eq 2. ]    
        then
            wmin=-7.
            wmax=9.
    		str_wmin="s/WMIN/$wmin/g"
		    sed -i -e $str_wmin input_main
		    str_wmax="s/WMAX/$wmax/g"
		    sed -i -e $str_wmax input_main
		elif [ $i -eq 3. ]    
        then
            wmin=-6.5
            wmax=9.5
    		str_wmin="s/WMIN/$wmin/g"
		    sed -i -e $str_wmin input_main
		    str_wmax="s/WMAX/$wmax/g"
		    sed -i -e $str_wmax input_main
		elif [ $i -eq 4. ]    
        then
            wmin=-6.
            wmax=10.
    		str_wmin="s/WMIN/$wmin/g"
		    sed -i -e $str_wmin input_main
		    str_wmax="s/WMAX/$wmax/g"
		    sed -i -e $str_wmax input_main
		else [ $i -eq 5. ]    
            wmin=-5.5
            wmax=10.5
    		str_wmin="s/WMIN/$wmin/g"
		    sed -i -e $str_wmin input_main
		    str_wmax="s/WMAX/$wmax/g"
		    sed -i -e $str_wmax input_main
	    fi
		    		    
		#change input_file
		str_u="s/COULOMB/$i/g"
		sed -i -e $str_u input_main
		str_v="s/INTERSITE/$j/g"
		sed -i -e $str_v input_main
	
		#compile and run
		mpif90 -g -o XH.exe global.f90 allocation.f90 deallocation.f90 main.f90 fermi.f90 rtbis.f90 -llapack

		mpirun -np $np ./XH.exe < $INFILE > logfile.dat
		
		#back to source_file
		cd $SRCDIR
		 
		done
	done
