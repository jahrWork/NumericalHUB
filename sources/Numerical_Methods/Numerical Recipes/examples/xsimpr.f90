	PROGRAM xsimpr
!	driver for routine simpr
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAR=3
	REAL(SP), PARAMETER :: HTOT=50.0_sp,X1=0.0
	INTEGER(I4B) :: i
	REAL(SP) :: a1,a2,a3
	REAL(SP), DIMENSION(NVAR) :: y,yout,dfdx,dydx
	REAL(SP), DIMENSION(NVAR,NVAR) :: dfdy
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	y(1)=1.0
	y(2)=1.0
	y(3)=0.0
	a1=0.5976_sp
	a2=1.4023_sp
	a3=0.0
	call derivs(X1,y,dydx)
	call jacobn(X1,y,dfdx,dfdy)
	write(*,'(1x,a/)') 'Test Problem'
	do i=5,50,5
		call simpr(y,dydx,dfdx,dfdy,X1,HTOT,i,yout,derivs)
		write(*,'(1x,a,f6.4,a,f7.4,a,i2,a)') 'X = ',X1,&
			' to ',X1+HTOT,' in ',i,' steps'
		write(*,'(1x,t5,a,t20,a)') 'Integration','Answer'
		write(*,'(1x,2f12.6)') yout(1),a1
		write(*,'(1x,2f12.6)') yout(2),a2
		write(*,'(1x,2f12.6)') yout(3),a3
	end do
	END PROGRAM xsimpr
