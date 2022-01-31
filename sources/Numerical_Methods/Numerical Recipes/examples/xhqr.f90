	PROGRAM xhqr
!	driver for routine hqr
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP1=8,NP2=4,NP3=5
	INTEGER(I4B) :: i,j,n
	LOGICAL(LGT) :: hessenberg
	INTEGER(I4B), DIMENSION(:), ALLOCATABLE :: indx
	REAL(SP), DIMENSION(:), ALLOCATABLE :: wr,wi
	REAL(SP), DIMENSION(:,:), ALLOCATABLE :: a
	REAL(SP), DIMENSION(NP1,NP1) :: b = reshape( (/ &
		3.0_sp,2.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,&
		2.0_sp,1.0_sp,3.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,&
		1.0_sp,3.0_sp,1.0_sp,1.0_sp,0.0_sp,0.0_sp,0.0_sp,0.0_sp,&
		2.0_sp,1.0_sp,2.0_sp,1.0_sp,1.0e-6_sp,0.0_sp,0.0_sp,0.0_sp,&
		1.0_sp,2.0_sp,1.0_sp,2.0_sp,3.0_sp,1.0e-6_sp,0.0_sp,0.0_sp,&
		4.0_sp,2.0_sp,2.0_sp,1.0_sp,1.0_sp,2.0_sp,1.0_sp,0.0_sp,&
		1.0_sp,1.0_sp,1.0_sp,3.0_sp,4.0_sp,1.0_sp,2.0_sp,3.0_sp,&
		2.0_sp,4.0_sp,3.0_sp,1.0_sp,2.0_sp,4.0_sp,3.0_sp,2.0_sp /),&
		(/ NP1,NP1 /) )
	REAL(SP), DIMENSION(NP2,NP2) :: c = reshape( (/ &
		3.0_sp,2.0_sp,3.0_sp,4.0_sp,1.0_sp,1.0_sp,1.0_sp,1.0_sp,&
		2.0_sp,3.0_sp,2.0_sp,3.0_sp,5.0_sp,7.0_sp,4.0_sp,2.0_sp /),&
		(/ NP2,NP2 /) )
	REAL(SP), DIMENSION(NP3,NP3) :: d = reshape( (/ &
		1.0_sp,-2.0_sp,3.0_sp,-4.0_sp,-5.0_sp,&
		2.0_sp,3.0_sp,4.0_sp,5.0_sp,6.0_sp,&
		0.0_sp,0.0_sp,50.0_sp,-60.0_sp,-70.0_sp,&
		0.0_sp,0.0_sp,0.0_sp,7.0_sp,8.0_sp,&
		0.0_sp,0.0_sp,0.0_sp,0.0_sp,-9.0_sp /), (/ NP3,NP3 /) )
! test exceptional shifts (matrix on p367 of Handbook)
	n=8
	allocate(a(n,n),wr(n),wi(n),indx(n))
	call unit_matrix(a)
	a=cshift(a,1,dim=2)
	hessenberg=.true.
	call doit
! test splitting criterion for 2 small consecutive subdiag. elements (p368,Hbk)
	a=b
	call doit
! well-conditioned, full matrix: Table 8, p394 Handbook
	deallocate(a,wr,wi,indx)
	n=4
	allocate(a(n,n),wr(n),wi(n),indx(n))
	a=c
	hessenberg=.false.
	call doit
! lower Hessenberg input
	deallocate(a,wr,wi,indx)
	n=5
	allocate(a(n,n),wr(n),wi(n),indx(n))
	a=d
	call doit
	deallocate(a,wr,wi,indx)
	CONTAINS
!BL
	SUBROUTINE doit
	write(*,'(/1x,a)') 'Matrix:'
	do i=1,n
		write(*,'(1x,8f9.2)') (a(i,j),j=1,n)
	end do
	if (.not. hessenberg) then
		call balanc(a)
		call elmhes(a)
		write(*,'(/1x,a)') 'Transformed Matrix:'
		do i=1,n
			write(*,'(1x,8f9.2)') (a(i,j),j=1,n)
		end do
	end if
	call hqr(a,wr,wi)
	call indexx(wr,indx)
	wr=wr(indx)
	wi=wi(indx)
	write(*,'(/1x,a)') 'Eigenvalues:'
	write(*,'(/1x,t9,a,t24,a/)') 'Real','Imag.'
	do i=1,n
		write(*,'(1x,2e15.6)') wr(i),wi(i)
	end do
	END SUBROUTINE doit
	END PROGRAM xhqr
