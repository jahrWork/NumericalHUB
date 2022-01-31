	PROGRAM xchebft
!	driver for routine chebft
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAL=40,NSHOW=16
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	INTEGER(I4B) :: i,j,mval
	REAL(SP) :: a,b,dum,f,t0,t1,term
	REAL(SP), DIMENSION(NVAL) :: c
	REAL(SP), DIMENSION(NSHOW) :: x,y,fx
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
!	test result
	do
		write(*,*) 'How many terms in Chebyshev evaluation?'
		write(*,'(1x,a,i2,a)') 'Enter n between 6 and ',NVAL,&
			'. Enter n=0 to END.'
		read(*,*) mval
		if ((mval <= 0) .or. (mval > NVAL)) exit
		write(*,'(1x,t10,a,t19,a,t28,a)') 'X','Actual','Chebyshev fit'
		x(1:NSHOW)=arth(-5,1,NSHOW)*PIO2/10.0_sp
		fx(1:NSHOW)=func(x(1:NSHOW))
		y(1:NSHOW)=(x(1:NSHOW)-0.5_sp*(b+a))/(0.5_sp*(b-a))
		do i=1,NSHOW
!	evaluate Chebyshev polynomial without using routine CHEBEV
			t0=1.0
			t1=y(i)
			f=c(2)*t1+c(1)*0.5_sp
			do j=3,mval
				dum=t1
				t1=2.0_sp*y(i)*t1-t0
				t0=dum
				term=c(j)*t1
				f=f+term
			end do
			write(*,'(1x,3f12.6)') x(i),fx(i),f
		end do
	end do
	END PROGRAM xchebft

	FUNCTION func(x)
	USE nrtype; USE nrutil
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: func
	func=(x**2)*(x**2-3.0_sp*x-2.0_sp)*sin(x)
	END FUNCTION func
