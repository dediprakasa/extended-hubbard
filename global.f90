	MODULE global
	
!	Mathematical and physical constants
	COMPLEX*16,PARAMETER::ii=CMPLX(0.D0,1.D0)
	REAL*8,PARAMETER::PI=2.D0*ASIN(1.D0)

	
!	Variables
	REAL*8:: eps,t1,t2,t3,T,wmin,wmax,dw,dws,nfilling,norm,mu,V
	REAL*8 :: delta_sigma,U,alpha,traceup,tracedown,tolmu,tolsigma,eta,errorsigma,parterrorsigma
	INTEGER :: Nw, N1D, Nsites, nproc		
	
!	Matrices
	COMPLEX*16,ALLOCATABLE::G(:,:,:), Gdummy(:,:), Ginv(:,:,:)
	REAL*8,ALLOCATABLE::IdMat(:,:),H(:,:),w(:),DoS(:),PDOS(:,:),wfreq(:),sigma(:,:),sigmac(:,:)
	REAL*8,ALLOCATABLE::sigma_down(:,:,:),sigmaf(:,:),sigma_up(:,:,:),cek_sigma(:,:,:),nn(:), n_up(:), n_down(:)
	REAL*8,ALLOCATABLE::DOSup(:),DOSdown(:),rmat(:,:),Vm(:,:)
	integer, allocatable::ind(:,:,:)

!	Arrays and variables for LAPACK
	complex*16, allocatable :: work(:) !dimension(2*2*Nsites)
	integer, allocatable :: ipiv(:) !dimension(2*Nsites)
	integer::info
		
	END MODULE global

