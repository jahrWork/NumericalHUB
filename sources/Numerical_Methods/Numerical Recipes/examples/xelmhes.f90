	PROGRAM xelmhes
!	driver for routine elmhes
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=5
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(NP,NP) :: a = reshape( (/ &
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,&
		2.0_sp,3.0_sp,4.0_sp,5.0_sp,6.0_sp,&
		300.0_sp,400.0_sp,5.0_sp,600.0_sp,700.0_sp,&
		4.0_sp,5.0_sp,6.0_sp,7.0_sp,8.0_sp,&
		5.0_sp,6.0_sp,7.0_sp,8.0_sp,9.0_sp /), (/ NP,NP /) )
	write(*,'(/1x,a/)') '***** Original Matrix *****'
	do i=1,NP
		write(*,'(1x,5f12.2)') (a(i,j),j=1,NP)
	end do
	write(*,'(/1x,a/)') '***** Balance Matrix *****'
	call balanc(a)
	do i=1,NP
		write(*,'(1x,5f12.2)') (a(i,j),j=1,NP)
	end do
	write(*,'(/1x,a/)') '***** Reduce to Hessenberg Form *****'
	call elmhes(a)
	do j=1,NP-2
		a(j+2:NP,j)=0.0
	end do
	do i=1,NP
		write(*,'(1x,5e12.4)') (a(i,j),j=1,NP)
	end do
	END PROGRAM xelmhes
