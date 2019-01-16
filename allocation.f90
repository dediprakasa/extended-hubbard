	subroutine allocation
	
	use global
	implicit none

	allocate(work(2*2*Nsites))
	allocate(ipiv(2*Nsites))	
	allocate(Gdummy(2*Nsites,2*Nsites))
	allocate(G(2*Nsites,2*Nsites,Nw))
	allocate(IdMat(2*Nsites,2*Nsites))
	allocate(H(2*Nsites,2*Nsites))
	allocate(rmat(2*Nsites,2*Nsites))
	allocate(w(Nw))
	allocate(DoS(Nw))
	allocate(PDOS(2*Nsites,Nw))
	allocate(wfreq(Nw))
	allocate(sigma(2*Nsites,2*Nsites))
	allocate(sigmac(2*Nsites,2*Nsites))
	allocate(sigmaf(2*Nsites,2*Nsites))
	allocate(sigma_up(Nsites,Nsites,Nw))
	allocate(sigma_down(Nsites,Nsites,Nw))
	allocate(nn(2*Nsites))
	allocate(n_up(Nsites))
	allocate(n_down(Nsites))
	allocate(DOSup(Nw))
	allocate(DOSdown(Nw))
	allocate(ind(0:N1D-1,0:N1D-1,0:N1D-1))
	allocate(Vm(2*Nsites,2*Nsites))

	return
	
	end subroutine allocation

