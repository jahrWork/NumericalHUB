	PROGRAM xfredin
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
!	driver for routine fredin
	INTEGER(I4B), PARAMETER :: N=8,MMAX=100
	INTEGER(I4B) :: m=1,mm
	REAL(SP) :: a,b
	REAL(SP), DIMENSION(N) :: t,f,w
	REAL(SP), DIMENSION(MMAX) :: ans,x
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
	do
		write(*,*) 'Enter T between 0 and PI/2'
		read(*,*,END=99) x(m)
		m=m+1
		if (m > mmax) exit
	end do
99	m=m-1
	ans(1:m)=fredin(x(1:m),a,b,t,f,w,g,ak)
	write(*,*) 'T, Calculated answer, True answer'
	do mm=1,m
		write(*,'(1x,3f10.6)') x(mm),ans(mm),sqrt(x(mm))
	end do
	END PROGRAM xfredin

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
