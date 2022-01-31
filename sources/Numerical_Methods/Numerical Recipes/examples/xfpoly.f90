	PROGRAM xfpoly
!	driver for routine fpoly
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAL=15,NPOLY=5
	REAL(SP), PARAMETER :: DX=0.1_sp
	INTEGER(I4B) :: i,j
	REAL(SP) :: x
	REAL(SP), DIMENSION(NPOLY) :: afunc
	write(*,'(/1x,t29,a)') 'Powers of X'
	write(*,'(/1x,t9,a,t17,a,t27,a,t37,a,t47,a,t57,a)') 'X','X**0',&
		'X**1','X**2','X**3','X**4'
	do i=1,NVAL
		x=i*DX
		afunc(1:NPOLY)=fpoly(x,NPOLY)
		write(*,'(1x,6f10.4)') x,(afunc(j),j=1,NPOLY)
	end do
	END PROGRAM xfpoly
