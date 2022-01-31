	PROGRAM xtred2
!	driver for routine tred2
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=10
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(NP) :: d,e
	REAL(SP), DIMENSION(NP,NP) :: a,f,c = reshape( (/ &
		5.0_sp,4.3_sp,3.0_sp,2.0_sp,1.0_sp,&
		0.0_sp,-1.0_sp,-2.0_sp,-3.0_sp,-4.0_sp,&
		4.3_sp,5.1_sp,4.0_sp,3.0_sp,2.0_sp,&
		1.0_sp,0.0_sp,-1.0_sp,-2.0_sp,-3.0_sp,&
		3.0_sp,4.0_sp,5.0_sp,4.0_sp,3.0_sp,&
		2.0_sp,1.0_sp,0.0_sp,-1.0_sp,-2.0_sp,&
		2.0_sp,3.0_sp,4.0_sp,5.0_sp,4.0_sp,&
		3.0_sp,2.0_sp,1.0_sp,0.0_sp,-1.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,&
		4.0_sp,3.0_sp,2.0_sp,1.0_sp,0.0_sp,&
		0.0_sp,1.0_sp,2.0_sp,3.0_sp,4.0_sp,&
		5.0_sp,4.0_sp,3.0_sp,2.0_sp,1.0_sp,&
		-1.0_sp,0.0_sp,1.0_sp,2.0_sp,3.0_sp,&
		4.0_sp,5.0_sp,4.0_sp,3.0_sp,2.0_sp,&
		-2.0_sp,-1.0_sp,0.0_sp,1.0_sp,2.0_sp,&
		3.0_sp,4.0_sp,5.0_sp,4.0_sp,3.0_sp,&
		-3.0_sp,-2.0_sp,-1.0_sp,0.0_sp,1.0_sp,&
		2.0_sp,3.0_sp,4.0_sp,5.0_sp,4.0_sp,&
		-4.0_sp,-3.0_sp,-2.0_sp,-1.0_sp,0.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp /),&
		(/ NP,NP /) )
	a(:,:)=c(:,:)
	call tred2(a,d,e)
	write(*,'(/1x,a)') 'Diagonal elements'
	write(*,'(1x,5f12.6)') (d(i),i=1,NP)
	write(*,'(/1x,a)') 'Off-diagonal elements'
	write(*,'(1x,5f12.6)') (e(i),i=2,NP)
!	check transformation matrix
	f=matmul(c,a)
	f=matmul(transpose(f),a)
!	how does it look?
	write(*,'(/1x,a)') 'Tridiagonal matrix'
	do i=1,NP
		write(*,'(1x,10f7.2)') (f(i,j),j=1,NP)
	end do
	END PROGRAM xtred2
