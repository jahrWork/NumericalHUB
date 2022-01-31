	PROGRAM xquad3d
!	driver for routine quad3d
	USE nrtype
	USE quad3d_qgaus_mod
	USE quad3d_qromb_mod
	IMPLICIT NONE
	COMMON xmax
	REAL(SP) :: xmax
	INTEGER(I4B), PARAMETER :: NVAL=10
	INTEGER(I4B) :: i
	REAL(SP) :: s,xmin
	write(*,*) 'Integral of r^2 over a spherical volume (from qgaus)'
	write(*,'(/4x,a,t14,a,t24,a)') 'Radius','QUAD3D','Actual'
	do i=1,NVAL
		xmax=0.1_sp*i
		xmin=-xmax
		call quad3d_qgaus(xmin,xmax,s)
		write(*,'(1x,f8.2,2f10.4)') xmax,s,4.0*PI*(xmax**5)/5.0
	end do
	write(*,*)
	write(*,*) 'Integral of r^2 over a spherical volume (from qromb)'
	write(*,'(/4x,a,t14,a,t24,a)') 'Radius','QUAD3D','Actual'
	do i=1,NVAL
		xmax=0.1_sp*i
		xmin=-xmax
		call quad3d_qromb(xmin,xmax,s)
		write(*,'(1x,f8.2,2f10.4)') xmax,s,4.0*PI*(xmax**5)/5.0
	end do
	END PROGRAM xquad3d

	FUNCTION func(x,y,z)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y
	REAL(SP), DIMENSION(:), INTENT(IN) :: z
	REAL(SP), DIMENSION(size(z)) :: func
	func=x**2+y**2+z(:)**2
	END FUNCTION func

	FUNCTION z1(x,y)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y
	REAL(SP) :: z1
	COMMON xmax
	REAL(SP) :: xmax
	z1=-sqrt(abs(xmax**2-x**2-y**2))
	END FUNCTION z1

	FUNCTION z2(x,y)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,y
	REAL(SP) :: z2
	COMMON xmax
	REAL(SP) :: xmax
	z2=sqrt(abs(xmax**2-x**2-y**2))
	END FUNCTION z2

	FUNCTION y1(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: y1
	COMMON xmax
	REAL(SP) :: xmax
	y1=-sqrt(abs(xmax**2-x**2))
	END FUNCTION y1

	FUNCTION y2(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: y2
	COMMON xmax
	REAL(SP) :: xmax
	y2=sqrt(abs(xmax**2-x**2))
	END FUNCTION y2
