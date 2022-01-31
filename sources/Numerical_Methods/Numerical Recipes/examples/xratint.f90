	PROGRAM xratint
!	driver for routine ratint
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPT=6
	REAL(SP), PARAMETER :: EPSSQ=1.0
	INTEGER(I4B) :: i
	REAL(SP) :: dyy,xx,yexp,yy
	REAL(SP), DIMENSION(NPT) :: x,y
	x(1:NPT)=arth(1,1,NPT)*2.0_sp/NPT
	do i=1,NPT
		y(i)=f(x(i))
	end do
	write(*,'(/1x,a/)') 'Diagonal rational function interpolation'
	write(*,'(1x,t6,a,t13,a,t26,a,t40,a)')&
		'x','interp.','accuracy','actual'
	do i=1,10
		xx=0.2_sp*i
		call ratint(x,y,xx,yy,dyy)
		yexp=f(xx)
		write(*,'(1x,f6.2,f12.6,e15.4,f12.6)') xx,yy,dyy,yexp
	end do
	CONTAINS
!BL
	FUNCTION f(x)
	IMPLICIT NONE
	REAL(SP) :: f,x
	f=x*exp(-x)/((x-1.0_sp)**2+EPSSQ)
	END FUNCTION f
	END PROGRAM xratint
