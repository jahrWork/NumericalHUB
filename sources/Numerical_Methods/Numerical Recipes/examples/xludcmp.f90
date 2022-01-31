	PROGRAM xludcmp
!	driver for routine ludcmp
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=20
	INTEGER(I4B) :: i,j,k,m,n
	INTEGER(I4B), DIMENSION(NP) :: indx,jndx
	REAL(SP) :: d
	REAL(SP), DIMENSION(NP,NP) :: a,xl,xu,x
	CHARACTER(3) :: txt
	open(7,file='MATRX1.DAT',status='old')
	read(7,*)
	do
		read(7,*)
		read(7,*) n,m
		read(7,*)
		read(7,*) ((a(i,j), j=1,n), i=1,n)
		read(7,*)
		read(7,*) ((x(i,j), i=1,n), j=1,m)
!	print out a-matrix for comparison with product of lower
!	and upper decomposition matrices.
		write(*,*) 'Original matrix:'
		do i=1,n
			write(*,'(1x,6f12.6)') (a(i,j), j=1,n)
		end do
!	perform the decomposition
		call ludcmp(a(1:n,1:n),indx(1:n),d)
!	compose separately the lower and upper matrices
		do i=1,n
			xu(i,i+1:n)=a(i,i+1:n)
			xl(i,i+1:n)=0.0
			xu(i,1:i-1)=0.0
			xl(i,1:i-1)=a(i,1:i-1)
			xu(i,i)=a(i,i)
			xl(i,i)=1.0
		end do
!	compute product of lower and upper matrices for
!	comparison with original matrix.
		jndx(1:n)=arth(1,1,n)
		x(1:n,1:n)=matmul(xl(1:n,1:n),xu(1:n,1:n))
		write(*,*) 'Product of lower and upper matrices (unscrambled):'
		do i=1,n
			call swap(jndx(indx(i)),jndx(i))
		end do
		do k=1,n
			do i=1,n
				if (jndx(i) == k) then
					write(*,'(1x,6f12.6)') (x(i,j), j=1,n)
				end if
			end do
		end do
		write(*,*) 'Lower matrix of the decomposition:'
		do i=1,n
			write(*,'(1x,6f12.6)') (xl(i,j), j=1,n)
		end do
		write(*,*) 'Upper matrix of the decomposition:'
		do i=1,n
			write(*,'(1x,6f12.6)') (xu(i,j), j=1,n)
		end do
		write(*,*) '***********************************'
		write(*,*) 'Press RETURN for next problem:'
		read(*,*)
		read(7,'(a3)') txt
		if (txt == 'END') exit
	end do
	close(7)
	END PROGRAM xludcmp
