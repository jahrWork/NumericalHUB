	PROGRAM xsvdcmp
!	driver for routine svdcmp
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: MP=20,NP=20
	INTEGER(I4B) :: i,j,m,n
	REAL(SP), DIMENSION(MP,NP) :: a,u
	REAL(SP), DIMENSION(NP) :: w
	REAL(SP), DIMENSION(NP,NP) :: v
	CHARACTER(3) :: dummy
	open(7,file='MATRX3.DAT',status='old')
	do
		read(7,'(a)') dummy
		if (dummy == 'END') exit
		read(7,*)
		read(7,*) m,n
		read(7,*)
!	copy original matrix into u
		do i=1,m
			read(7,*) (a(i,j), j=1,n)
			u(i,1:n)=a(i,1:n)
		end do
!	perform decomposition
		call svdcmp(u(1:m,1:n),w(1:n),v(1:n,1:n))
!	print results
		write(*,*) 'Decomposition Matrices:'
		write(*,*) 'Matrix U'
		do i=1,m
			write(*,'(1x,6f12.6)') (u(i,j),j=1,n)
		end do
		write(*,*) 'Diagonal of Matrix W'
		write(*,'(1x,6f12.6)') (w(i),i=1,n)
		write(*,*) 'Matrix V-Transpose'
		do i=1,n
			write(*,'(1x,6f12.6)') (v(j,i),j=1,n)
		end do
		write(*,*) 'Check product against original matrix:'
		write(*,*) 'Original Matrix:'
		do i=1,m
			write(*,'(1x,6f12.6)') (a(i,j),j=1,n)
		end do
		write(*,*) 'Product U*W*(V-Transpose):'

		do i=1,m
			a(i,1:n)=matmul(v(1:n,1:n),u(i,1:n)*w(1:n))
			write(*,'(1x,6f12.6)') (a(i,j),j=1,n)
		end do
		write(*,*) '***********************************'
		write(*,*) 'Press RETURN for next problem'
		read(*,*)
	end do
	close(7)
	END PROGRAM xsvdcmp
