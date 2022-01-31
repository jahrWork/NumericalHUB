	PROGRAM xsort2
!	driver for routine sort2
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=10,N2=N*N
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(N2) :: a,b
	open(7,file='TARRAY.DAT',status='OLD')
	read(7,*) (a(i),i=1,N2)
	close(7)
!	generate B-array
	b(1:N2)=arth(1,1,N2)
!	sort A and mix B
	call sort2(a,b)
	write(*,*) 'After sorting A and mixing B, array A is:'
	do i=1,N
		write(*,'(1x,10f7.2)') (a(N*(i-1)+j), j=1,N)
	end do
	write(*,*) '...and array B is:'
	do i=1,N
		write(*,'(1x,10f7.2)') (b(N*(i-1)+j), j=1,N)
	end do
	write(*,*) 'press RETURN to continue...'
	read(*,*)
!	sort B and mix A
	call sort2(b,a)
	write(*,*) 'After sorting B and mixing A, array A is:'
	do i=1,N
		write(*,'(1x,10f7.2)') (a(N*(i-1)+j), j=1,N)
	end do
	write(*,*) '...and array B is:'
	do i=1,N
		write(*,'(1x,10f7.2)') (b(N*(i-1)+j), j=1,N)
	end do
	END PROGRAM xsort2
