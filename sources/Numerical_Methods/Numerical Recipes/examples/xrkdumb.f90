	PROGRAM xrkdumb
!	driver for routine rkdumb
	USE nrtype
	USE nr
	USE rkdumb_path
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NSTEP=1500
	INTEGER(I4B) :: i,j
	REAL(SP) :: x1,x2
	REAL(SP), DIMENSION(4) :: vstart
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	x1=1.0
	vstart(1)=bessj0(x1)
	vstart(2)=bessj1(x1)
	vstart(3)=bessj(2,x1)
	vstart(4)=bessj(3,x1)
	x2=20.0
	call rkdumb(vstart,x1,x2,NSTEP,derivs)
	write(*,'(/1x,t9,a,t17,a,t31,a/)') 'X','Integrated','BESSJ3'
	do i=1,(NSTEP/10)
		j=10*i
		write(*,'(1x,f10.4,2x,2f12.6)') xx(j),y(4,j),bessj(3,xx(j))
	end do
	END PROGRAM xrkdumb

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
