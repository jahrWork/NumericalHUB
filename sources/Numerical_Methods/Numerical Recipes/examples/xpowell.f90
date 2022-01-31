	PROGRAM xpowell
!	driver for routine powell
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDIM=3
	REAL(SP), PARAMETER :: FTOL=1.0e-6_sp
	INTEGER(I4B) :: i,iter
	REAL(SP) :: fret
	REAL(SP), DIMENSION(NDIM) :: p = (/ 1.5_sp,1.5_sp,2.5_sp /)
	REAL(SP), DIMENSION(NDIM,NDIM) :: xi = reshape( (/ &
	1.0_sp,0.0_sp,0.0_sp,0.0_sp,1.0_sp,0.0_sp,0.0_sp,0.0_sp,1.0_sp /),&
	(/ NDIM,NDIM /) )
	call powell(p,xi,FTOL,iter,fret)
	write(*,'(/1x,a,i3)') 'Iterations:',iter
	write(*,'(/1x,a/1x,3f12.6)') 'Minimum found at: ',(p(i),i=1,NDIM)
	write(*,'(/1x,a,f12.6)') 'Minimum function value =',fret
	write(*,'(/1x,a)') 'True minimum of function is at:'
	write(*,'(1x,3f12.6/)') 1.0,2.0,3.0
	END PROGRAM xpowell

	FUNCTION func(x)
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP) :: func
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
!	func=1.0_sp-product(bessj0(x(:)-0.5_sp))
	func=0.5_sp-bessj0((x(1)-1.0_sp)**2+(x(2)-2.0_sp)**2+(x(3)-3.0_sp)**2)
	END FUNCTION func
