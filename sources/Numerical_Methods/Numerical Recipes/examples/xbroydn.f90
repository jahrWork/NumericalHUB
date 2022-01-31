	PROGRAM xbroydn
!	driver for routine broydn
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTERFACE
		FUNCTION funcv(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funcv
		END FUNCTION funcv
	END INTERFACE
	INTEGER(I4B), PARAMETER :: N=2
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(N) :: x,f
	LOGICAL(LGT) :: check
	x(1)=2.0
	x(2)=0.5_sp
	call broydn(x,check)
	f=funcv(x)
	if (check) write(*,*) 'Convergence problems.'
	write(*,'(1x,a5,t10,a1,t22,a1)') 'Index','x','f'
	do i=1,N
		write(*,'(1x,i2,2x,2f12.6)') i,x(i),f(i)
	end do
	END PROGRAM xbroydn

	FUNCTION funcv(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: funcv
	funcv(1)=x(1)**2+x(2)**2-2.0_sp
	funcv(2)=exp(x(1)-1.0_sp)+x(2)**3-2.0_sp
	END FUNCTION funcv
