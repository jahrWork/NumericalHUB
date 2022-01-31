	PROGRAM xeigsrt
!	driver for routine eigsrt
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=10
	INTEGER(I4B) :: i,j,nrot
	REAL(SP), DIMENSION(NP) :: d
	REAL(SP), DIMENSION(NP,NP) :: v,c = reshape( (/ &
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
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp /),	(/ NP,NP /) )
	call jacobi(c,d,v,nrot)
	write(*,*) 'Unsorted Eigenvectors:'
	do i=1,NP
		write(*,'(/1x,a,i3,a,f12.6)') 'Eigenvalue',i,' =',d(i)
		write(*,*) 'Eigenvector:'
		write(*,'(10x,5f12.6)') (v(j,i),j=1,NP)
	end do
	write(*,'(//,1x,a,//)') '****** sorting ******'
	call eigsrt(d,v)
	write(*,*) 'Sorted Eigenvectors:'
	do i=1,NP
		write(*,'(/1x,a,i3,a,f12.6)') 'Eigenvalue',i,' =',d(i)
		write(*,*) 'Eigenvector:'
		write(*,'(10x,5f12.6)') (v(j,i),j=1,NP)
	end do
	END PROGRAM xeigsrt
