	PROGRAM xchint
!	driver for routine chint
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAL=40
	REAL(SP), PARAMETER :: PIO4=PIO2/2.0_sp
	INTEGER(I4B) :: i,mval
	REAL(SP) :: a,b,x
	REAL(SP), DIMENSION(NVAL) :: c,cint
	INTERFACE
		FUNCTION fint(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: fint
		END FUNCTION fint
!BL
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	a=-PIO4
	b=PIO2
	c(1:NVAL)=chebft(a,b,NVAL,func)
!	test integral
	do
		write(*,*) 'How many terms in Chebyshev evaluation?'
		write(*,'(1x,a,i2,a)') 'Enter n between 6 and ',NVAL,&
			'. Enter n=0 to END.'
		read(*,*) mval
		if ((mval <= 0) .or. (mval > NVAL)) exit
		cint(1:mval)=chint(a,b,c(1:mval))
		write(*,'(1x,t10,a,t19,a,t29,a)') 'X','Actual','Cheby. Integ.'
		do i=-5,10
			x=i*PIO2/10.0_sp
			write(*,'(1x,3f12.6)') x,fint(x)-fint(-PIO4),&
				chebev(a,b,cint(1:mval),x)
		end do
	end do
	END PROGRAM xchint

	FUNCTION func(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: func
	func=(x**2)*(x**2-3.0_sp*x-2.0_sp)*sin(x)
	END FUNCTION func

	FUNCTION fint(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: fint
!	integral of func
	fint=(4.0_sp*x**3-9.0_sp*x**2-28.0_sp*x+18.0_sp)*sin(x)-&
		(x**4-3.0_sp*x**3-14.0_sp*x**2+18.0_sp*x+28.0_sp)*cos(x)
	END FUNCTION fint
