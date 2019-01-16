program main
	use mpi
	use global	
	implicit none
	
! For use in MPI environment
	integer, parameter :: nprocmax=200
    integer ierr,myid
	integer numtasks,numextra
	integer sourceid,destid
    integer nstart(0:nprocmax-1),nfinish(0:nprocmax-1),numdata(0:nprocmax-1)	
	integer status(MPI_STATUS_SIZE)

! Declaration of variables
	real*8, external :: fermi, rtbis, integrand
	integer ::i, j, l, m, n, lp, mp, np, counter, iw, pw, is, Niter, iconv
	real*8 :: time_start,time_finish,time_elapsed
	character*50 :: DOSupfile,DOSdownfile,DOStotfile,mufile,nnnfile

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)     

	if (myid==0) then

        time_start=real(MPI_Wtime(),4)

!		OPEN(unit=10, file="input_main", status="old")
		read(*,*) N1D
		read(*,*) eps
		read(*,*) t1
		read(*,*) t2
		read(*,*) t3
		read(*,*) Nw
		read(*,*) wmin
		read(*,*) wmax
		read(*,*) nfilling
		read(*,*) U
		read(*,*) V
		read(*,*) T
		read(*,*) eta
		read(*,*) alpha
		read(*,*) Niter
		read(*,*) tolmu
		read(*,*) tolsigma
		read(*,*) DOSupfile
		read(*,*) DOSdownfile
		read(*,*) DOStotfile
		read(*,*) mufile
		read(*,*) nnnfile
!		close(10)
	end if
	
    call MPI_BCAST(N1D,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(eps,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	call MPI_BCAST(t1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(t2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(t3,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Nw,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(wmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(wmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nfilling,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(U,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(T,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(eta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(alpha,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Niter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tolmu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(tolsigma,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(DOSupfile,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(DOSdownfile,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(DOStotfile,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(mufile,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(nnnfile,1,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

	Nsites = N1D**3
	nfilling = nfilling*real(2*Nsites, 8)

	if (myid==0) write(*,*) "Finished reading input_main."

	call allocation

!	Simpson composite weighting
	dw = (wmax-wmin)/real(Nw-1,8)

	do i=1,Nw 
		w(i) = wmin + dw*(i-1)
	end do	

	wfreq(1) = dw/3.D0
	do i = 2, Nw-1, 2
		wfreq(i) = 4.D0*dw/3.D0
	end do

	do i=3, Nw-2, 2
		wfreq(i) = 2.D0*dw/3.D0
	end do

	wfreq(Nw) = dw/3.D0
	
!	Identity matrix
	IdMat(:,:)=0.D0
	do i = 1, 2*Nsites
		IdMat(i,i) = 1.D0
	end do
	
!	Generate Hamiltonian matrix
	counter = 0
	do l = 0, N1D-1
	do m = 0, N1D-1
	do n = 0, N1D-1
		ind(l,m,n) = counter + 1
		counter = counter + 1
	end do
	end do
	end do

	H(:,:)=0.D0
	do l=0,N1D-1
	do m=0,N1D-1
	do n=0,N1D-1
		do lp=0,N1D-1
		do mp=0,N1D-1
		do np=0,N1D-1
			i=ind(l,m,n)
			j=ind(lp,mp,np)
			
			if(i==j) then
				H(i,j)=eps
			else
				if((abs(l-lp)==1 .and. abs(m-mp)==0 .and. abs(n-np)==0) .or. &
					&(abs(l-lp)==0 .and. abs(m-mp)==1 .and. abs(n-np)==0) .OR. &
					&(abs(l-lp)==0 .and. abs(m-mp)==0 .and. abs(n-np)==1)) then
					H(i,j)=-t1
				else if((abs(l-lp)==1 .and. abs(m-mp)==1 .and. abs(n-np)==0) .or. &
					&(abs(l-lp)==1 .and. abs(m-mp)==0 .and. abs(n-np)==1) .OR. &
					&(abs(l-lp)==0 .and. abs(m-mp)==1 .and. abs(n-np)==1)) then
					H(i,j)=-t2
				else if((abs(l-lp)==1 .and. abs(m-mp)==1 .and. abs(n-np)==1)) then
					H(i,j)=-t3
				else
					H(i,j)=0.d0
				end if
			end if
		end do
		end do
		end do
	end do
	end do
	end do
	
	H(Nsites+1:2*Nsites,Nsites+1:2*Nsites)=H(1:Nsites,1:Nsites)
	
! Initial sigma matrix
	sigma(:,:) = 0.d0
	
! For U>0, setup initial occupation number with AFM configuration
	n_up(ind(0,0,0)) = 1.d0
	n_down(ind(0,0,0)) = abs(1.d0 - n_up(ind(0,0,0)))
	do n = 1, N1D-1
		n_up(ind(0,0,n)) = abs(1.d0 - n_up(ind(0,0,n-1)))
		n_down(ind(0,0,n)) = abs(1.d0 - n_up(ind(0,0,n)))
	end do
	do m = 1, N1D-1
	do n = 0, N1D-1
		n_up(ind(0,m,n)) = abs(1.d0 - n_up(ind(0,m-1,n)))
		n_down(ind(0,m,n)) = abs(1.d0 - n_up(ind(0,m,n)))
	end do
	end do
	do l = 1, N1D-1
	do m = 0, N1D-1
	do n = 0, N1D-1
		n_up(ind(l,m,n)) = abs(1.d0 - n_up(ind(l-1,m,n)))
		n_down(ind(l,m,n)) = abs(1.d0 - n_up(ind(l,m,n)))
	end do
	end do
	end do

! For U=0, setup initial occupation number with PM configuration
	if (U == 0.d0) then
		n_up = 0.5d0
		n_down = 0.5d0
	end if
	
	Vm(:,:) = 0.d0
	counter = 0
	do l=0,N1D-1
	do m=0,N1D-1
	do n=0,N1D-1
		ind(l,m,n) = counter + 1
		i = ind(l,m,n)
		
	    ! ========================= Sudut ======================== !
		if (l == 0 .and. m == 0 .and. n == 0) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                 &   + n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                 &   + n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V
        elseif (l == 0 .and. m == 0 .and. n == N1D-1) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                 &   + n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                 &   + n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V      
        elseif (l == 0 .and. m == N1D-1 .and. n == 0) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                 &   + n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                 &   + n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V        
	    elseif (l == 0 .and. m == N1D-1 .and. n == N1D-1) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                 &   + n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                 &   + n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V	
	    elseif (l == N1D-1 .and. m == 0 .and. n == 0) then
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                 &   + n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                 &   + n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V
	    elseif (l == N1D-1 .and. m == 0 .and. n == N1D-1) then
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                 &   + n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                 &   + n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V	
	    elseif (l == N1D-1 .and. m == N1D-1 .and. n == 0) then
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                 &   + n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                 &   + n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V	
	    elseif (l == N1D-1 .and. m == N1D-1 .and. n == N1D-1) then
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                 &   + n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                 &   + n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        ! ======================================================== !
        
        ! ========================= Rusuk ======================== !
        elseif (l == 0     .and. m == 0    .and. (n > 0) .and. (n < N1D-1)) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
	            &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        elseif (l == 0     .and. m == N1D-1    .and. (n > 0) .and. (n < N1D-1)) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
	            &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        elseif (l == N1D-1     .and. m == 0    .and. (n > 0) .and. (n < N1D-1)) then
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
	            &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        elseif (l == N1D-1     .and. m == N1D-1    .and. (n > 0) .and. (n < N1D-1)) then
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
	            &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V     
        elseif (l == 0     .and. (m > 0) .and. (m < N1D-1)   .and. n == 0) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
	            &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V
        elseif (l == 0     .and. (m > 0) .and. (m < N1D-1)   .and. n == N1D-1) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
	            &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V                
        elseif (l == N1D-1     .and. (m > 0) .and. (m < N1D-1)   .and. n == 0) then
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
	            &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V                
        elseif (l == N1D-1     .and. (m > 0) .and. (m < N1D-1)   .and. n == N1D-1) then
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
	            &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        elseif ((l > 0) .and. (l < N1D-1) .and. m == 0 .and. n == 0) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
	            &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V                
        elseif ((l > 0) .and. (l < N1D-1) .and. m == 0 .and. n == N1D-1) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
	            &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V                
        elseif ((l > 0) .and. (l < N1D-1) .and. m == N1D-1 .and. n == 0) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
	            &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V                 
        elseif ((l > 0) .and. (l < N1D-1) .and. m == N1D-1 .and. n == N1D-1) then
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
	            &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        ! ======================================================== !
        
        ! ========================= Sisi ========================= !
        elseif (l == 0 .and. (m > 0) .and. (m < N1D-1) .and. (n > 0) .and. (n < N1D-1)) then                        
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &   
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V                                
        elseif (l == N1D-1 .and. (m > 0) .and. (m < N1D-1) .and. (n > 0) .and. (n < N1D-1)) then                        
            Vm(i,i) = (n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &   
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        elseif ((l > 0) .and. (l < N1D-1) .and. m == 0 .and. (n > 0) .and. (n < N1D-1)) then                        
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &   
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V                
        elseif ((l > 0) .and. (l < N1D-1) .and. m == N1D-1 .and. (n > 0) .and. (n < N1D-1)) then                        
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &   
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V                
         elseif ((l > 0) .and. (l < N1D-1) .and. (m > 0) .and. (m < N1D-1) .and. n == 0) then                        
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &   
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)))*V               
         elseif ((l > 0) .and. (l < N1D-1) .and. (m > 0) .and. (m < N1D-1) .and. n == N1D-1) then                        
            Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &   
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        ! ======================================================== !
        
        ! ========================= Tengah ========================= !                        
		else
		    Vm(i,i) = (n_up(ind(l+1,m,n)) + n_down(ind(l+1,m,n)) &
                &   +  n_up(ind(l-1,m,n)) + n_down(ind(l-1,m,n)) &
                &   +  n_up(ind(l,m+1,n)) + n_down(ind(l,m+1,n)) &
                &   +  n_up(ind(l,m-1,n)) + n_down(ind(l,m-1,n)) &   
                &   +  n_up(ind(l,m,n+1)) + n_down(ind(l,m,n+1)) &
                &   +  n_up(ind(l,m,n-1)) + n_down(ind(l,m,n-1)))*V
        end if
        counter = counter + 1
    end do
    end do
    end do                
    
    Vm(Nsites+1:2*Nsites,Nsites+1:2*Nsites) = Vm(1:Nsites,1:Nsites)
                                            
! Filling sigma elements
	do j = 1, Nsites
		sigma(j,j) = U*n_down(j) + Vm(j,j)
		sigma(j+Nsites,j+Nsites) = U*n_up(j) + Vm(j+Nsites,j+Nsites)
	end do
	
	

! Arrange the job sharing by distributing the real frequency jobs to all processors 

    numtasks = nw/nproc
    numextra = nw - numtasks*nproc
             
    if (numextra .gt. 0) then
		do i = 0, nproc-1-numextra
			nstart(i) = i*numtasks + 1
			nfinish(i) = nstart(i) + numtasks - 1
		end do
		
		do i = nproc-numextra, nproc-1
			nstart(i) = nfinish(i-1) + 1
			nfinish(i) = nstart(i) + numtasks
		end do	
    else
		do i = 0, nproc-1
			nstart(i) = i*numtasks + 1
			nfinish(i) = nstart(i) + numtasks-1
		end do
    end if

    do i=0,nproc-1
        numdata(i)=nfinish(i)-nstart(i)+1
    end do 
        

    if (myid == 0) then
        open(unit=555, file='check_jobs_dist', status='unknown')
        write(555,*) 'For the self-consistency process with Nw =',Nw
        do i=0,nproc-1 
           write(555,*) 'id =',i,'nstrt=',nstart(i),'nfin =',&
          &             nfinish(i)
        end do
        close(555)
        write(*,*) 'Finished writing check_jobs_dist for the Self-Consistency process!'
    end if

   	call MPI_Barrier(MPI_COMM_WORLD,ierr)

! Start iteration
	
	if (U == 0.d0) Niter=1

	iconv=0

	do is = 1, Niter	
		if (myid == 0) then
			write(*,*)
			write(*,*) 'Iteration number', is
		end if		

		do iw = nstart(myid),nfinish(myid)		
			Gdummy(:,:) = (w(iw) + ii*eta)*IdMat(:,:) - H(:,:)- sigma(:,:)	
		
			!	Matrix Inversion
			call zgetrf(2*Nsites, 2*Nsites, Gdummy, 2*Nsites, ipiv, info)
			call zgetri(2*Nsites, Gdummy, 2*Nsites, ipiv, work, 2*2*Nsites, info)

        	G(:,:,iw) = Gdummy(:,:)
			
			do j = 1, 2*Nsites
				PDOS(j,iw) = -(1.d0/pi)*aimag(G(j,j,iw))
			end do	
		end do

		if (myid .ne. 0) then
			call MPI_SSEND(PDOS(1,nstart(myid)),2*Nsites*numdata(myid), &
		&          MPI_DOUBLE_PRECISION,0,111,MPI_COMM_WORLD,ierr)
		else
			do sourceid=1, nproc-1
				call MPI_RECV(PDOS(1,nstart(sourceid)),2*Nsites*numdata(sourceid), &
		&       	MPI_DOUBLE_PRECISION,sourceid,111,MPI_COMM_WORLD,status,ierr)
			end do
		end if

		if (myid == 0) then
			norm = 0.d0	
			do iw = 1, Nw
				DOSup(iw) = 0.d0
				DOSdown(iw) = 0.d0
				do j = 1, Nsites
					DOSup(iw) = DOSup(iw) + PDOS(j,iw)	
					DOSdown(iw) = DOSdown(iw) + PDOS(j+Nsites,iw)
				end do
				DOS(iw) = DOSup(iw) + DOSdown(iw)
				norm = norm + DOS(iw)*wfreq(iw)
			end do

			mu = rtbis(integrand,wmin,wmax,tolmu)
			write(*,*) "Chemical potential:",mu

			!	Calculate <n>
			nn(:) = 0.d0
			do j = 1, 2*Nsites	
				do iw = 1,Nw
					nn(j) = nn(j) + PDOS(j,iw)*fermi(mu,T,w(iw))*wfreq(iw)			
				end do	
			end do
		end if

    	call MPI_BCAST(nn(1),2*Nsites,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

		sigmac(:,:) = 0.d0
		parterrorsigma = 0.d0
!		do iw = nstart(myid), nfinish(myid)
		do j = 1, Nsites
			sigmac(j,j) = nn(j+Nsites)*U
			sigmac(j+Nsites,j+Nsites) = nn(j)*U
		end do

			! Error check
		rmat(:,:) = abs(sigmac(:,:) - sigma(:,:))
		parterrorsigma = max(parterrorsigma, maxval(rmat))
!		end do

        call MPI_REDUCE(parterrorsigma, errorsigma, 1, MPI_DOUBLE_PRECISION,&
     	&                   MPI_MAX, 0, MPI_COMM_WORLD, ierr )

		
		if (myid == 0) then
			write(*,*) 'Error', errorsigma
			write(*,*) 'Norm', norm
			open(unit=11, file=DOStotfile, status='unknown')
			open(unit=12, file=DOSupfile, status='unknown')
			open(unit=13, file=DOSdownfile, status='unknown')
			open(unit=14, file=mufile, status='unknown')
			open(unit=15,file=nnnfile,status='replace')	
				do iw = 1, Nw
					write(11,*)w(iw), DoS(iw)
					write(12,*)w(iw), DoSup(iw)
					write(13,*)w(iw), DOSdown(iw)
				end do
					do j = 1,Nsites
					write(15,*)j,(sigma(j+Nsites,j+Nsites)-Vm(j+Nsites,j+Nsites))/U,(sigma(j,j)-Vm(j,j))/U
		 		end do
				write(14,*)mu,0.D0
				write(14,*)mu,Nsites	
			close(11)
			close(12)
			close(13)
			close(14)
			close(15)

			if (errorsigma .lt. tolsigma) then 
				write(*,*) 'Convergent!'
				iconv = 1
			end if

		end if
		
			call MPI_BCAST(iconv,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

		if (iconv==0) then
			sigma(:,:) = alpha*sigmac(:,:) + (1-alpha)*sigma(:,:)
		else
			EXIT
		end if


	end do ! iter

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

	! 	Reconstructing the Green Function using the already converged <n(j)> by all processors
	do iw = 1, Nw		
		sigmaf(:,:) = 0.d0
		do j = 1,Nsites
			sigmaf(j,j) = nn(j+Nsites)*U + Vm(j+Nsites,j+Nsites)
			sigmaf(j+Nsites,j+Nsites) = nn(j)*U + Vm(j,j)
		end do

		Gdummy(:,:) = (w(iw)+ii*eta)*IdMat(:,:) - H(:,:) - sigmaf(:,:)	
		
		!	   Matrix Inversion
		call zgetrf(2*Nsites, 2*Nsites, Gdummy, 2*Nsites,ipiv,info)
		call zgetri(2*Nsites, Gdummy, 2*Nsites, ipiv,work, 2*2*Nsites,info)

        G(:,:,iw)=Gdummy(:,:)			
	end do

	call Optical_Response(myid)

!	call deallocation

	if (myid == 0) then
			time_finish = real(MPI_Wtime(),4)	
			time_elapsed = time_finish - time_start
			write(*,*) 'All calculations done !!!'
			write(*,*)
			write(*,*) 'Elapsed time (seconds) =', time_elapsed
			write(*,*) 'Elapsed time (minutes) =', time_elapsed/60.0
			write(*,*) 'Elapsed time (hours) =', time_elapsed/3600.0
			write(*,*)
	end if

	call MPI_FINALIZE(ierr)

	stop
    
end program
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	FUNCTION integrand(x)
	
	USE global
	IMPLICIT NONE
	INTEGER :: i
	REAL*8::integrand,x,tot
	REAL*8,EXTERNAL::fermi

	tot=0.D0
	DO i=1,Nw
		tot=tot+DOS(i)*fermi(x,T,w(i))*wfreq(i)
	END DO
	integrand=nfilling-(tot*2*Nsites/norm)
	
	RETURN

	
	END FUNCTION integrand
