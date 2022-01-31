	PROGRAM xbsstep
!	driver for routine bsstep
	USE nrtype
	USE nr
	USE ode_path
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAR=4
	INTEGER(I4B) :: i,nrhs
	REAL(SP) :: eps,h1,hmin,x1,x2
	REAL(SP), DIMENSION(NVAR) :: ystart
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	COMMON nrhs
	nrhs=0
	save_steps=.true.
	x1=1.0
	x2=10.0
	ystart(1)=bessj0(x1)
	ystart(2)=bessj1(x1)
	ystart(3)=bessj(2,x1)
	ystart(4)=bessj(3,x1)
	eps=1.0e-4_sp
	h1=0.1_sp
	hmin=0.0
	dxsav=(x2-x1)/20.0_sp
	call odeint(ystart,x1,x2,eps,h1,hmin,derivs,bsstep)
	write(*,'(/1x,a,t30,i3)') 'Successful steps:',nok
	write(*,'(1x,a,t30,i3)') 'Bad steps:',nbad
	write(*,'(1x,a,t30,i3)') 'Function evaluations:',nrhs
	write(*,'(1x,a,t30,i3)') 'Stored intermediate values:',kount
	write(*,'(/1x,t9,a,t20,a,t33,a)') 'X','Integral','BESSJ(3,X)'
	do i=1,kount
		write(*,'(1x,f10.4,2x,2f14.6)') xp(i),yp(4,i),bessj(3,xp(i))
	end do
	END PROGRAM xbsstep

	SUBROUTINE derivs(x,y,dydx)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
	INTEGER(I4B) :: nrhs
	COMMON nrhs
	nrhs=nrhs+1
	dydx(1)=-y(2)
	dydx(2)=y(1)-(1.0_sp/x)*y(2)
	dydx(3)=y(2)-(2.0_sp/x)*y(3)
	dydx(4)=y(3)-(3.0_sp/x)*y(4)
	END SUBROUTINE derivs
