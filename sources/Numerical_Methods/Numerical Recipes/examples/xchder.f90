	PROGRAM xchder
!	driver for routine chder
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAL=40
	INTEGER(I4B) :: i,mval
	REAL(SP) :: a,b,x
	REAL(SP), DIMENSION(NVAL) :: c,cder
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
!BL
		FUNCTION fder(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: fder
		END FUNCTION fder
	END INTERFACE
	a=-PIO2/2.0_sp
	b=PIO2
	c(1:NVAL)=chebft(a,b,NVAL,func)
!	test derivative
	do
		write(*,*) 'How many terms in Chebyshev evaluation?'
		write(*,'(1x,a,i2,a)') 'Enter n between 6 and ',NVAL,&
			'. Enter n=0 to END.'
		read(*,*) mval
		if ((mval <= 0) .or. (mval > NVAL)) exit
		cder(1:mval)=chder(a,b,c(1:mval))
		write(*,'(1x,t10,a,t19,a,t28,a)') 'X','Actual','Cheby. Deriv.'
		do i=-5,10
			x=i*PIO2/10.0_sp
			write(*,'(1x,3f12.6)') x,fder(x),chebev(a,b,cder(1:mval),x)
		end do
	end do
	END PROGRAM xchder

	FUNCTION func(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: func
	func=(x**2)*(x**2-3.0_sp*x-2.0_sp)*sin(x)
	END FUNCTION func

	FUNCTION fder(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: fder
!	derivative of func
	fder=x*(4.0_sp*x**2-9.0_sp*x-4.0_sp)*sin(x)+(x**2)*(x**2-3.0_sp*x-2.0_sp)*cos(x)
	END FUNCTION fder
