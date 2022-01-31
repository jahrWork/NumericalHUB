	PROGRAM xpade
!	driver for routine pade
	USE nrtype; USE nrutil
	USE nr, ONLY : pade,ratval
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NMAX=100
	INTEGER(I4B) :: i,n,n2p
	INTEGER(I4B), DIMENSION(NMAX) :: itmp
	REAL(SP) :: resid
	REAL(DP) :: b,d,fac,x
	REAL(DP), DIMENSION(NMAX) :: c,cc
	itmp=arth(1,1,NMAX)
	do
		write(*,*) 'Enter n for PADE routine:'
		read(*,*,END=999) n
		n2p=2*n+1
		fac=1
		c(1:n2p)=1/real(itmp(1:n2p),dp)
		where (mod(itmp(1:n2p),2) == 0) c(1:n2p)=-c(1:n2p)
		cc(1:n2p)=c(1:n2p)
		call pade(c(1:n2p),resid)
		write(*,'(1x,a,1p,d16.8)') 'Norm of residual vector=',resid
		write(*,*) 'point, func. value, pade series, power series'
		do i=1,21
			x=(i-1)*0.25_dp
			b=poly(x,cc(1:n2p))
			d=ratval(x,c,n,n)
			write(*,'(1p,4d16.8)') x,fn(x),d,b
		end do
	end do
	CONTAINS
!BL
	FUNCTION fn(x)
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: fn
	if (x == 0.0) then
		fn=1.0
	else
		fn=log(1.0_dp+x)/x
	end if
	END FUNCTION fn
999	END PROGRAM xpade
