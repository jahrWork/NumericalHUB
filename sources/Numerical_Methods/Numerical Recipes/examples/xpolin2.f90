	PROGRAM xpolin2
!	driver for routine polin2
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=5
	INTEGER(I4B) :: i,j
	REAL(SP) :: dy,f,x1,x2,y
	REAL(SP), DIMENSION(N) :: x1a,x2a
	REAL(SP), DIMENSION(N,N) :: ya
	x1a(1:N)=arth(1,1,N)*PI/N
	x2a(1:N)=arth(1.0_sp,1.0_sp,N)/N
	ya(:,:)=outerprod(sin(x1a(:)),exp(x2a(:)))
!	test 2-dimensional interpolation
	write(*,'(t9,a,t21,a,t32,a,t40,a,t58,a)')&
		'x1','x2','f(x)','interpolated','error'
	do i=1,N-1
		x1=(-0.1_sp+i/5.0_sp)*PI
		do j=1,N-1
			x2=-0.1_sp+j/5.0_sp
			f=sin(x1)*exp(x2)
			call polin2(x1a,x2a,ya,x1,x2,y,dy)
			write(*,'(1x,4f12.6,f14.6)') x1,x2,f,y,dy
		end do
		write(*,*) '***********************************'
	end do
	END PROGRAM xpolin2
