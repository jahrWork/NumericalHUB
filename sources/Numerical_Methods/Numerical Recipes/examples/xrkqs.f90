	PROGRAM xrkqs
!	driver for routine rkqs
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=4
	INTEGER(I4B) :: i
	REAL(SP) :: eps,hdid,hnext,htry,x
	REAL(SP), DIMENSION(N) :: y,dydx,dysav,ysav,yscal
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
	ysav(1)=bessj0(x)
	ysav(2)=bessj1(x)
	ysav(3)=bessj(2,x)
	ysav(4)=bessj(3,x)
	call derivs(x,ysav,dysav)
	yscal(:)=1.0
	htry=0.6_sp
	write(*,'(/1x,t8,a,t19,a,t31,a,t43,a)')&
		'eps','htry','hdid','hnext'
	do i=1,15
		eps=exp(-real(i,sp))
		x=1.0
		y(:)=ysav(:)
		dydx(:)=dysav(:)
		call rkqs(y,dydx,x,htry,eps,yscal,hdid,hnext,derivs)
		write(*,'(2x,e12.4,f8.2,2x,2f12.6)') eps,htry,hdid,hnext
	end do
	END PROGRAM xrkqs

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
