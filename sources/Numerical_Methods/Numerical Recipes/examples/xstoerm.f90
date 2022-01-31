	PROGRAM xstoerm
!	driver for routine stoerm
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NVAR=4
	REAL(SP), PARAMETER :: HTOT=1.570796_sp,X1=0.0
	INTEGER(I4B) :: i
	REAL(SP) :: a1,a2,x,xf,d2y1,d2y2
	REAL(SP), DIMENSION(NVAR) :: y,yout,d2y
	d2y1(x)=x+sin(x)
	d2y2(x)=x**2+cos(x)-2
	y(1)=0.0
	y(2)=-1.0
	y(3)=2.0
	y(4)=0.0
	call derivs(X1,y,d2y)
	xf=X1+HTOT
	a1=d2y1(xf)
	a2=d2y2(xf)
	write(*,'(1x,a/)') 'Stoermer''s Rule:'
	do i=5,45,10
		call stoerm(y,d2y,X1,HTOT,i,yout,derivs)
		write(*,'(1x,a,f6.4,a,f6.4,a,i2,a)') 'X = ',X1,&
			' to ',X1+HTOT,' in ',i,' steps'
		write(*,'(1x,t5,a,t20,a)') 'Integration','Answer'
		write(*,'(1x,2f12.6)') yout(1),a1
		write(*,'(1x,2f12.6)') yout(2),a2
	end do
	END PROGRAM xstoerm

	SUBROUTINE derivs(x,y,d2y)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: d2y
	d2y(1)=x-y(1)
	d2y(2)=x*x-y(2)
	END SUBROUTINE derivs
