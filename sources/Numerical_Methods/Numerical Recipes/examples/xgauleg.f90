	PROGRAM xgauleg
!	driver for routine gauleg
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPOINT=10
	REAL(SP), PARAMETER :: X1=0.0_sp,X2=1.0_sp,X3=10.0
	INTEGER(I4B) :: i
	REAL(SP) :: xx
	REAL(SP), DIMENSION(NPOINT) :: x,w
	call gauleg(X1,X2,x,w)
	write(*,'(/1x,t3,a,t10,a,t22,a/)') '#','X(I)','W(I)'
	do i=1,NPOINT
		write(*,'(1x,i2,2f12.6)') i,x(i),w(i)
	end do
!	demonstrate the use of GAULEG for an integral
	call gauleg(X1,X3,x,w)
	xx=0.0
	do i=1,NPOINT
		xx=xx+w(i)*func(x(i))
	end do
	write(*,'(/1x,a,f12.6)') 'Integral from GAULEG:',xx
	write(*,'(1x,a,f12.6)') 'Actual value:',1.0-(1.0+X3)*exp(-X3)
	CONTAINS
!BL
	FUNCTION func(x)
	IMPLICIT NONE
	REAL(SP) :: func,x
	func=x*exp(-x)
	END FUNCTION func
	END PROGRAM xgauleg
