	PROGRAM xmidpnt
!	driver for routine midpnt
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NMAX=10
	INTEGER(I4B) :: i
	REAL(SP) :: a,b,s
	INTERFACE
		FUNCTION func2(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func2
		END FUNCTION func2
!BL
		FUNCTION fint2(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: fint2
		END FUNCTION fint2
	END INTERFACE
	a=0.0
	b=1.0
	write(*,*) 'Integral of FUNC2 computed with MIDPNT'
	write(*,*) 'Actual value of integral is',fint2(b)-fint2(a)
	write(*,'(1x,t7,a,t20,a)') 'n','Approx. Integral'
	do i=1,NMAX
		call midpnt(func2,a,b,s,i)
		write(*,'(1x,i6,f24.6)') i,s
	end do
	END PROGRAM xmidpnt

	FUNCTION func2(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: func2
	func2=1.0_sp/sqrt(x)
	END FUNCTION func2

	FUNCTION fint2(x)
!	integral of FUNC2
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: fint2
	fint2=2.0_sp*sqrt(x)
	END FUNCTION fint2
