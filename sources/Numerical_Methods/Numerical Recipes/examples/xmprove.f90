	PROGRAM xmprove
!	driver for routine mprove
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=5
	INTEGER(I4B) :: i,idum
	INTEGER(I4B), DIMENSION(N) :: indx
	REAL(SP) :: d
	REAL(SP), DIMENSION(N) :: b=(/ 1.0_sp,1.0_sp,1.0_sp,1.0_sp,1.0_sp /),x
	REAL(SP), DIMENSION(N,N) :: aa,a = reshape( (/ &
		1.0_sp,2.0_sp,1.0_sp,4.0_sp,5.0_sp,&
		2.0_sp,3.0_sp,1.0_sp,5.0_sp,1.0_sp,&
		0.0_sp,4.0_sp,1.0_sp,1.0_sp,2.0_sp,&
		4.0_sp,5.0_sp,1.0_sp,2.0_sp,3.0_sp,&
		0.0_sp,1.0_sp,1.0_sp,3.0_sp,4.0_sp /), (/ N,N /) )
	x(:)=b(:)
	aa(:,:)=a(:,:)
	call ludcmp(aa,indx,d)
	call lubksb(aa,indx,x)
	write(*,'(/1x,a)') 'Solution vector for the equations:'
	write(*,'(1x,5f12.6)') (x(i),i=1,N)
!	now phoney up x and let mprove fix it
	idum=-123
	do i=1,N
		x(i)=x(i)*(1.0_sp+0.2_sp*ran(idum))
	end do
	write(*,'(/1x,a)') 'Solution vector with noise added:'
	write(*,'(1x,5f12.6)') (x(i),i=1,N)
	call mprove(a,aa,indx,b,x)
	write(*,'(/1x,a)') 'Solution vector recovered by MPROVE:'
	write(*,'(1x,5f12.6)') (x(i),i=1,N)
	END PROGRAM xmprove
