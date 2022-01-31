	PROGRAM xchebpc
!	driver for routine chebpc
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAL=40,NSHOW=16
	INTEGER(I4B) :: i,mval
	REAL(SP) :: a,b
	REAL(SP), DIMENSION(NVAL) :: c,d
	REAL(SP), DIMENSION(NSHOW) :: x,y,fx,ply
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
	do
		write(*,*) 'How many terms in Chebyshev evaluation?'
		write(*,'(1x,a,i2,a)') 'Enter n between 6 and ',NVAL,&
			'. Enter n=0 to END.'
		read(*,*) mval
		if ((mval <= 0) .or. (mval > NVAL)) exit
		d(1:mval)=chebpc(c(1:mval))
!	test polynomial
		write(*,'(1x,t10,a,t19,a,t29,a)') 'X','Actual','Polynomial'
		x(1:NSHOW)=arth(-5,1,NSHOW)*PIO2/10.0_sp
		fx(1:NSHOW)=func(x(1:NSHOW))
		y(1:NSHOW)=(x(1:NSHOW)-0.5_sp*(b+a))/(0.5_sp*(b-a))
		ply(1:NSHOW)=poly(y(1:NSHOW),d(1:mval))
		do i=1,NSHOW
			write(*,'(1x,3f12.6)') x(i),fx(i),ply(i)
		end do
	end do
	END PROGRAM xchebpc

	FUNCTION func(x)
	USE nrtype; USE nrutil
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: func
	func=(x**2)*(x**2-3.0_sp*x-2.0_sp)*sin(x)
	END FUNCTION func
