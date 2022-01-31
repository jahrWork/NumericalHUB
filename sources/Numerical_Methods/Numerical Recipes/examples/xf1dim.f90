	PROGRAM xf1dim
!	driver for routine f1dim
	USE nrtype
	USE nr
	USE f1dim_mod
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDIM=3
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(NDIM), TARGET :: xi,p = (/ 0.0_sp,0.0_sp,0.0_sp /)
	ncom=NDIM
	pcom=>p
	xicom=>xi
	write(*,'(/1x,a)') 'Enter vector direction along which to'
	write(*,'(1x,a)') 'plot the function. Minimum is in the'
	write(*,'(1x,a)') 'direction 1.0,1.0,1.0 - Enter X,Y,Z:'
	read(*,*) (xi(i),i=1,3)
	call scrsho(f1dim)
	END PROGRAM xf1dim

	FUNCTION func(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP) :: func
	func=sum((x(:)-1.0_sp)**2)
	END FUNCTION func
