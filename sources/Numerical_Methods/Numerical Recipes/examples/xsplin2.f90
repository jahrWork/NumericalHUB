	PROGRAM xsplin2
!	driver for routine splin2
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: M=10,N=10
	INTEGER(I4B) :: i
	REAL(SP) :: f,ff,x1x2,xx1,xx2
	REAL(SP), DIMENSION(M) :: x1
	REAL(SP), DIMENSION(M,N) :: y,y2
	REAL(SP), DIMENSION(N) :: x2
	x1(1:M)=0.2_sp*arth(1,1,M)
	x2(:)=x1(:)
	y(:,:)=outerprod(x1(:),x2(:))
	y(:,:)=y(:,:)*exp(-y(:,:))
	call splie2(x1,x2,y,y2)
	write(*,'(/1x,t9,a,t21,a,t31,a,t43,a)')	'x1','x2','splin2','actual'
	do i=1,10
		xx1=0.1_sp*i
		xx2=xx1**2
		f=splin2(x1,x2,y,y2,xx1,xx2)
		x1x2=xx1*xx2
		ff=x1x2*exp(-x1x2)
		write(*,'(1x,4f12.6)') xx1,xx2,f,ff
	end do
	END PROGRAM xsplin2
