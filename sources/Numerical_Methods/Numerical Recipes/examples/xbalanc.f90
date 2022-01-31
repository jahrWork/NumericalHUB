	PROGRAM xbalanc
!	driver for routine balanc
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=5
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(NP) :: r,c
	REAL(SP), DIMENSION(NP,NP) :: aa,a = reshape( (/ &
		1.0_sp,1.0_sp,1.0_sp,1.0_sp,1.0_sp,&
		100.0_sp,1.0_sp,100.0_sp,1.0_sp,100.0_sp,&
		1.0_sp,1.0_sp,1.0_sp,1.0_sp,1.0_sp,&
		100.0_sp,1.0_sp,100.0_sp,1.0_sp,100.0_sp,&
		1.0_sp,1.0_sp,1.0_sp,1.0_sp,1.0_sp /), (/ NP,NP /) )
!	print norms
	aa=abs(a)
	r(:)=sum(aa,dim=2)
	c(:)=sum(aa,dim=1)
	write(*,*) 'Rows:'
	write(*,*) (r(i),i=1,NP)
	write(*,*) 'Columns:'
	write(*,*) (c(i),i=1,NP)
	write(*,'(/1x,a/)') '***** Balancing Matrix *****'
	call balanc(a)
!	print norms
	aa=abs(a)
	r(:)=sum(aa,dim=2)
	c(:)=sum(aa,dim=1)
	write(*,*) 'Rows:'
	write(*,*) (r(i),i=1,NP)
	write(*,*) 'Columns:'
	write(*,*) (c(i),i=1,NP)
	END PROGRAM xbalanc
