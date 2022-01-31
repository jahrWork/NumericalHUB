	PROGRAM xigray
!	driver for routine igray
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: jp,n,ng,nmax,nmin,nni,nxor
	do
		write(*,*) 'INPUT NMIN,NMAX:'
		read(*,*,END=999) nmin,nmax
		jp=max(1,(nmax-nmin)/11)
		write(*,*) 'n, Gray(n), Gray(Gray(n)), Gray(n) .xor. Gray(n+1)'
		do n=nmin,nmax
			ng=igray(n,1)
			nni=igray(ng,-1)
			if (nni /= n) write(*,*) 'WRONG ! AT ',n,ng,nni
			if (mod(n-nmin,jp) == 0) then
				nxor=ieor(ng,igray(n+1,1))
				write(*,*) n,ng,nni,nxor
			end if
		end do
	end do
999	write(*,*) 'NORMAL COMPLETION'
	END PROGRAM xigray
