	PROGRAM xfleg
!	driver for routine fleg
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAL=5,NPOLY=5
	REAL(SP), PARAMETER :: DX=0.2_sp
	INTEGER(I4B) :: i,j
	REAL(SP) :: x
	REAL(SP), DIMENSION(NPOLY) :: afunc
	write(*,'(/1x,t25,a)') 'Legendre Polynomials'
	write(*,'(/1x,t8,a,t18,a,t28,a,t38,a,t48,a)')&
		'N=1','N=2','N=3','N=4','N=5'
	do i=1,NVAL
		x=i*DX
		afunc(1:NPOLY)=fleg(x,NPOLY)
		write(*,'(1x,a,f6.2)') 'X =',x
		write(*,'(1x,5f10.4,a)') (afunc(j),j=1,NPOLY),'  routine FLEG'
		write(*,'(1x,5f10.4,a/)') (plgndr(j-1,0,x),j=1,NPOLY),&
			'  routine PLGNDR'
	end do
	END PROGRAM xfleg
