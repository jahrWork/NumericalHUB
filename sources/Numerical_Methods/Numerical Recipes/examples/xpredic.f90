	PROGRAM xpredic
!	driver for routine predic
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPTS=500,NPOLES=10,NFUT=20
	INTEGER(I4B) :: i
	REAL(SP) :: dum
	REAL(SP), DIMENSION(NFUT) :: future
	REAL(SP), DIMENSION(NPOLES) :: d
	REAL(SP), DIMENSION(NPTS) :: data
	do i=1,NPTS
		data(i)=func(i)
	end do
	call memcof(data,dum,d)
	call fixrts(d)
	future(1:NFUT)=predic(data,d,NFUT)
	write(*,'(6x,a,t13,a,t25,a)') 'I','Actual','PREDIC'
	do i=1,NFUT
		write(*,'(1x,i6,2f12.6)') i,func(i+NPTS),future(i)
	end do
	CONTAINS
!BL
	FUNCTION func(n)
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: n
	REAL(SP) :: func
	func=exp(-1.0_sp*n/NPTS)*sin(2.0_sp*PI*n/50.0_sp)+&
		exp(-2.0_sp*n/NPTS)*sin(2.2_sp*PI*n/50.0_sp)
	END FUNCTION func
	END PROGRAM xpredic
