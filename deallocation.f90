	subroutine deallocation
	
	use global
	implicit none

	deallocate(work)
	deallocate(ipiv)	
	deallocate(Gdummy)
	deallocate(G)
	deallocate(IdMat)
	deallocate(H)
	deallocate(rmat)
	deallocate(w)
	deallocate(DoS)
	deallocate(PDOS)
	deallocate(wfreq)
	deallocate(sigma)
	deallocate(sigmac)
	deallocate(sigmaf)
	deallocate(sigma_up)
	deallocate(sigma_down)
	deallocate(cek_sigma)
	deallocate(nn)
	deallocate(n_up)
	deallocate(n_down)
	deallocate(DOSup)
	deallocate(DOSdown)
	deallocate(ind)

	return
	
	end subroutine deallocation
