	PROGRAM xrank
!	driver for routine rank
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=10,N2=N*N
	INTEGER(I4B) :: i,j,k,l
	INTEGER(I4B), DIMENSION(N2) :: indx,irankv
	REAL(SP), DIMENSION(N) :: b
	REAL(SP), DIMENSION(N2) :: a
	open(7,file='TARRAY.DAT',status='OLD')
	read(7,*) (a(i),i=1,100)
	close(7)
	call indexx(a,indx)
	irankv=rank(indx)
	write(*,*) 'Original array is:'
	do i=1,N
		write(*,'(1x,10f7.2)') (a(N*(i-1)+j), j=1,N)
	end do
	write(*,*) 'Table of ranks is:'
	do i=1,N
		write(*,'(1x,10i6)') (irankv(N*(i-1)+j), j=1,N)
	end do
	write(*,*) 'press RETURN to continue...'
	read(*,*)
	write(*,*) 'Array sorted according to rank table:'
	do i=1,N
		do j=1,N
			k=N*(i-1)+j
				do l=1,N2
					if (irankv(l) == k) b(j)=a(l)
				end do
		end do
		write(*,'(1x,10f7.2)') (b(j),j=1,N)
	end do
	END PROGRAM xrank
