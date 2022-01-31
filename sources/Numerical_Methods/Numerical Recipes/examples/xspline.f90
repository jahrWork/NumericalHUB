	PROGRAM xspline
!	driver for routine spline
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=20
	INTEGER(I4B) :: i
	REAL(SP) :: yp1,ypn
	REAL(SP), DIMENSION(N) :: x,y,y2
	write(*,*) 'Second-derivatives for sin(x) from 0 to PI'
!	generate array for interpolation
	x(1:N)=arth(1,1,N)*PI/N
	y(:)=sin(x(:))
!	calculate 2nd derivative with SPLINE
	yp1=cos(x(1))
	ypn=cos(x(N))
	call spline(x,y,yp1,ypn,y2)
!	test result
	write(*,'(t19,a,t35,a)') 'spline','actual'
	write(*,'(t6,a,t17,a,t33,a)') 'angle','2nd deriv','2nd deriv'
	do i=1,N
		write(*,'(1x,f8.2,2f16.6)') x(i),y2(i),-sin(x(i))
	end do
	END PROGRAM xspline
