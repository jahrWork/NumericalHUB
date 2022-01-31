	PROGRAM xgauher
!	driver for routine gauher
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=64
	REAL(SP), PARAMETER :: SQRTPI=1.7724539_sp
	INTEGER(I4B) :: i,n
	REAL(SP) :: check,xx
	REAL(SP), DIMENSION(NP) :: x,w
	do
		write(*,*) 'Enter N'
		read(*,*,END=99) n
		call gauher(x(1:n),w(1:n))
		write(*,'(/1x,t3,a,t10,a,t22,a/)') '#','X(I)','W(I)'
		do i=1,n
			write(*,'(1x,i2,2e14.6)') i,x(i),w(i)
		end do
		check=sum(w(1:n))
		write(*,'(/1x,a,e15.7,a,e15.7)') 'Check value:',check,&
			'  should be:',SQRTPI
!	demonstrate the use of GAUHER for an integral
		xx=0.0
		do i=1,n
			xx=xx+w(i)*func(x(i))
		end do
		write(*,'(/1x,a,f12.6)') 'Integral from GAUHER:',xx
		write(*,'(1x,a,f12.6)') 'Actual value:        ',SQRTPI*exp(-0.25_sp)
	end do
	CONTAINS
!BL
	FUNCTION func(x)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: func
	func=cos(x)
	END FUNCTION func
99	END PROGRAM xgauher
