	PROGRAM xtoeplz
!	driver for routine toeplz
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=5,N2M1=2*N-1
	INTEGER(I4B) :: i
	REAL(SP) :: sm
	REAL(SP), DIMENSION(N) :: x,y
	REAL(SP), DIMENSION(N2M1) :: r
	y(:)=arth(0.1_sp,0.1_sp,N)
	r(:)=1.0_sp/arth(1,1,N2M1)
	x(1:N)=toeplz(r,y)
	write(*,*) 'Solution vector:'
	do i=1,N
		write(*,'(5x,a2,i1,a4,e13.6)') 'X(',i,') = ',x(i)
	end do
	write(*,'(/1x,a)') 'Test of solution:'
	write(*,'(1x,t6,a,t19,a)') 'mtrx*soln','original'
	do i=1,N
		sm=dot_product(r(N-1+i:i:-1),x(1:N))
		write(*,'(1x,2f12.4)') sm,y(i)
	end do
	END PROGRAM xtoeplz
