	PROGRAM xvoltra
	USE nrtype
	USE nr
	IMPLICIT NONE
!	driver for routine voltra
	INTERFACE
		FUNCTION g(t)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: t
		REAL(SP), DIMENSION(:), POINTER :: g
		END FUNCTION g
!BL
		FUNCTION ak(t,s)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: t,s
		REAL(SP), DIMENSION(:,:), POINTER :: ak
		END FUNCTION ak
	END INTERFACE
	REAL(SP), PARAMETER :: H=0.05_sp
	INTEGER(I4B), PARAMETER :: N=30,M=2
	INTEGER(I4B) :: i
	REAL(SP) :: t0
	REAL(SP), DIMENSION(N) :: t
	REAL(SP), DIMENSION(M,N) :: f
	t0=0.0
	call voltra(t0,H,t,f,g,ak)
!	exact soln is f(1)=exp(-t), f(2)=2sin(t)
	write(*,*)&
		'abscissa,voltra answer1,real answer1,voltra answer2,real answer2'
	do i=1,N
		write(*,'(1x,5f10.6)') t(i),f(1,i),exp(-t(i)),f(2,i),2.0_sp*sin(t(i))
	end do
	END PROGRAM xvoltra

	FUNCTION g(t)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: t
	REAL(SP), DIMENSION(2), TARGET, SAVE :: gg
	REAL(SP), DIMENSION(:), POINTER :: g
	g=>gg
	g(1)=cosh(t)+t*sin(t)
	g(2)=2.0_sp*sin(t)+t*(sin(t)**2+exp(t))
	END FUNCTION g

	FUNCTION ak(t,s)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: t,s
	REAL(SP), DIMENSION(2,2), TARGET, SAVE :: akk
	REAL(SP), DIMENSION(:,:), POINTER :: ak
	ak=>akk
	ak(1,1)=-exp(t-s)
	ak(1,2)=-cos(t-s)
	ak(2,1)=-exp(t+s)
	ak(2,2)=-t*cos(s)
	END FUNCTION ak
