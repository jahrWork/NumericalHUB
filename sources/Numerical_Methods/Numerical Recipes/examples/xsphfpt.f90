!	auxiliary routine for sphfpt
	SUBROUTINE derivs(x,y,dydx)
	USE nrtype
	USE sphfpt_data
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(IN) :: y
	REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
	dydx(1)=y(2)
	dydx(2)=(2.0_sp*x*(m+1.0_sp)*y(2)-(y(3)-c2*x*x)*y(1))/(1.0_sp-x*x)
	dydx(3)=0.0
	END SUBROUTINE derivs
