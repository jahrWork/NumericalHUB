	PROGRAM xtrapzd
!	driver for routine trapzd
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NMAX=14
	INTEGER(I4B) :: i
	REAL(SP) :: a,b,s
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
!BL
		FUNCTION fint(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: fint
		END FUNCTION fint
	END INTERFACE
	a=0.0
	b=PIO2
	write(*,'(1x,a)') 'Integral of FUNC with 2^(n-1) points'
	write(*,'(1x,a,f10.6)') 'Actual value of integral is',&
		fint(b)-fint(a)
	write(*,'(1x,t7,a,t16,a)') 'n','Approx. Integral'
	do i=1,NMAX
		call trapzd(func,a,b,s,i)
		write(*,'(1x,i6,f20.6)') i,s
	end do
	END PROGRAM xtrapzd

	FUNCTION func(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: func
	func=(x**2)*(x**2-2.0_sp)*sin(x)
	END FUNCTION func

	FUNCTION fint(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: fint
	fint=4.0_sp*x*((x**2)-7.0_sp)*sin(x)-((x**4)-14.0_sp*(x**2)+28.0_sp)*cos(x)
	END FUNCTION fint
