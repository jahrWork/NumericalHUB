	PROGRAM xsplie2
!	driver for routine splie2
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: M=10,N=10
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(M) :: x1
	REAL(SP), DIMENSION(M,N) :: y,y2
	REAL(SP), DIMENSION(N) :: x2
	x1(1:M)=0.2_sp*arth(1,1,M)
	x2(:)=x1(:)
	y(:,:)=outerprod(x1(:),x2(:))**2
	call splie2(x1,x2,y,y2)
	write(*,'(/1x,a)') 'Second derivatives from SPLIE2'
	write(*,'(1x,a/)') 'Natural spline assumed'
	do i=1,5
		write(*,'(1x,5f12.6)') (y2(i,j),j=1,5)
	end do
	write(*,'(/1x,a/)') 'Actual second derivatives'
	y2(1:5,1:5)=2.0_sp*(spread(x1(1:5),dim=2,ncopies=5)**2)
	do i=1,5
		write(*,'(1x,5f12.6)') (y2(i,j),j=1,5)
	end do
	END PROGRAM xsplie2
