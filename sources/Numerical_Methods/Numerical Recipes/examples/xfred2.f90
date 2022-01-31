	PROGRAM xfred2
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
!	driver for routine fred2
	INTEGER(I4B), PARAMETER :: N=8
	INTEGER(I4B) :: i
	REAL(SP) :: a,b
	REAL(SP), DIMENSION(N) :: t,f,w
	INTERFACE
		FUNCTION g(t)
		USE nrtype; USE nrutil
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: t
		REAL(SP), DIMENSION(size(t)) :: g
		END FUNCTION g
!BL
		FUNCTION ak(t,s)
		USE nrtype; USE nrutil
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
		REAL(SP), DIMENSION(size(t),size(s)) :: ak
		END FUNCTION ak
	END INTERFACE
	a=0.0
	b=PIO2
	call fred2(a,b,t,f,w,g,ak)
!	compare with exact solution
	write(*,*) 'Abscissa, Calc soln, True soln'
	do i=1,N
		write(*,'(1x,3f10.6)') t(i),f(i),sqrt(t(i))
	end do
	END PROGRAM xfred2

	FUNCTION g(t)
	USE nrtype; USE nrutil
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: t
	REAL(SP), DIMENSION(size(t)) :: g
	g=sqrt(t)-PIO2**2.25_sp*t**0.75_sp/2.25_sp
	END FUNCTION g

	FUNCTION ak(t,s)
	USE nrtype; USE nrutil
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: t,s
	REAL(SP), DIMENSION(size(t),size(s)) :: ak
	ak=outerprod(t,s)**0.75_sp
	END FUNCTION ak
