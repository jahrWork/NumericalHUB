	PROGRAM xlubksb
!	driver for routine lubksb
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=20
	REAL(SP) :: p
	REAL(SP), DIMENSION(NP) :: x
	REAL(SP), DIMENSION(NP,NP) :: a,b,c
	INTEGER(I4B) :: i,j,m,n
	INTEGER(I4B), DIMENSION(NP) :: indx
	CHARACTER(3) :: txt
	open(7,file='MATRX1.DAT',status='old')
	read(7,*)
	do
		read(7,*)
		read(7,*) n,m
		read(7,*)
		read(7,*) ((a(i,j), j=1,n), i=1,n)
		read(7,*)
		read(7,*) ((b(i,j), i=1,n), j=1,m)
!	save matrix a for later testing
		c(1:n,1:n)=a(1:n,1:n)
!	do LU decomposition
		call ludcmp(c(1:n,1:n),indx(1:n),p)
!	solve equations for each right-hand vector
		do j=1,m
			x(1:n)=b(1:n,j)
			call lubksb(c(1:n,1:n),indx(1:n),x(1:n))
!	test results with original matrix
			write(*,*) 'Right-hand side vector:'
			write(*,'(1x,6f12.6)') (b(i,j), i=1,n)
			write(*,*) 'Result of matrix applied to sol''n vector'
			b(1:n,j)=matmul(a(1:n,1:n),x(1:n))
			write(*,'(1x,6f12.6)') (b(i,j), i=1,n)
			write(*,*) '***********************************'
		end do
		write(*,*) 'Press RETURN for next problem:'
		read(*,*)
		read(7,'(a3)') txt
		if (txt == 'END') exit
	end do
	close(7)
	END PROGRAM xlubksb
