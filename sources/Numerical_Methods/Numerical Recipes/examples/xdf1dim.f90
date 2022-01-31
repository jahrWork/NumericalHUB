	PROGRAM xdf1dim
!	driver for routine df1dim
	USE nrtype; USE nrutil
	USE nr
	USE df1dim_mod
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
	call scrsho(df1dim)
	END PROGRAM xdf1dim

	FUNCTION dfunc(x)
	USE nrtype; USE nrutil
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: dfunc
	dfunc(:)=(x(:)-1.0_sp)**2
	END FUNCTION dfunc

	FUNCTION func(x)
	USE nrtype; USE nrutil
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP) :: func
	call nrerror('xdf1dim: dummy function, should never be called')
	func=sum(x)
!	just to quiet compiler warnings
	END FUNCTION func
