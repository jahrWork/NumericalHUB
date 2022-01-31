	PROGRAM xchebev
!	driver for routine chebev
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAL=40,MAXX=16
	INTEGER(I4B) :: i,mval
	REAL(SP) :: a,b
	REAL(SP), DIMENSION(MAXX) :: x,funcx,chebxs,chebxv
	REAL(SP), DIMENSION(NVAL) :: c
	INTERFACE
		FUNCTION func(x)
		USE nrtype; USE nrutil
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	a=-PIO2/2.0_sp
	b=PIO2
	c(1:NVAL)=chebft(a,b,NVAL,func)
	x=arth(-5,1,MAXX)*PIO2/10.0_sp
	funcx=func(x)
!	test Chebyshev evaluation routine
	do
		write(*,*) 'How many terms in Chebyshev evaluation?'
		write(*,'(1x,a,i2,a)') 'Enter n between 6 and ',NVAL,&
			'. Enter n=0 to END.'
		read(*,*) mval
		if ((mval <= 0) .or. (mval > NVAL)) exit
		do i=1,MAXX
			chebxs(i)=chebev(a,b,c(1:mval),x(i))
		end do
		chebxv(1:MAXX)=chebev(a,b,c(1:mval),x(1:MAXX))
		write(*,'(1x,t10,a,t19,a,t28,a,t42,a)') &
			'X','Actual','Chebyshev fit','Chebyshev fit'
		do i=1,MAXX
			write(*,'(1x,4f12.6)') x(i),funcx(i),chebxs(i),chebxv(i)
		end do
	end do
	END PROGRAM xchebev

	FUNCTION func(x)
	USE nrtype; USE nrutil
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: func
	func=(x**2)*(x**2-3.0_sp*x-2.0_sp)*sin(x)
	END FUNCTION func
