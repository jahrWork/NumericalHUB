	PROGRAM xrk4
!	driver for routine rk4
	USE nrtype
	USE nr
	INTEGER(I4B), PARAMETER :: N=4
	INTEGER(I4B) :: i,j
	REAL(SP) :: h,x
	REAL(SP), DIMENSION(N) :: y,dydx,yout
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	x=1.0
	y(1)=bessj0(x)
	y(2)=bessj1(x)
	y(3)=bessj(2,x)
	y(4)=bessj(3,x)
	call derivs(x,y,dydx)
	write(*,'(/1x,a,t19,a,t31,a,t43,a,t55,a)')&
		'Bessel Function:','J0','J1','J3','J4'
	do i=1,5
		h=0.2*i
		call rk4(y,dydx,x,h,yout,derivs)
		write(*,'(/1x,a,f6.2)') 'For a step size of:',h
		write(*,'(1x,a10,4f12.6)') 'RK4:',(yout(j),j=1,4)
		write(*,'(1x,a10,4f12.6)') 'Actual:',bessj0(x+h),&
			bessj1(x+h),bessj(2,x+h),bessj(3,x+h)
	end do
	END PROGRAM xrk4

	SUBROUTINE derivs(x,y,dydx)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
	dydx(1)=-y(2)
	dydx(2)=y(1)-(1.0_sp/x)*y(2)
	dydx(3)=y(2)-(2.0_sp/x)*y(3)
	dydx(4)=y(3)-(3.0_sp/x)*y(4)
	END SUBROUTINE derivs
