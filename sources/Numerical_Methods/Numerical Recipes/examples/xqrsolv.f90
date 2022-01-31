	PROGRAM xqrsolv
!	driver for routine qrsolv
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=20
	REAL(SP), DIMENSION(NP) :: c,d,x
	REAL(SP), DIMENSION(NP,NP) :: a,ai,b
	INTEGER(I4B) :: i,j,m,n
	CHARACTER(3) :: txt
	LOGICAL(LGT) :: sing
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
		ai(1:n,1:n)=a(1:n,1:n)
!	do qr decomposition
		call qrdcmp(a(1:n,1:n),c(1:n),d(1:n),sing)
		if (sing) write(*,*) 'Singularity in QR decomposition.'
!	solve equations for each right-hand vector
		do j=1,m
			x(1:n)=b(1:n,j)
			call qrsolv(a(1:n,1:n),c(1:n),d(1:n),x(1:n))
!	test results with original matrix
			write(*,*) 'Right-hand side vector:'
			write(*,'(1x,6f12.6)') (b(i,j), i=1,n)
			write(*,*) 'Result of matrix applied to sol''n vector'
			b(1:n,j)=matmul(ai(1:n,1:n),x(1:n))
			write(*,'(1x,6f12.6)') (b(i,j), i=1,n)
			write(*,*) '***********************************'
		end do
		write(*,*) 'Press RETURN for next problem:'
		read(*,*)
		read(7,'(a3)') txt
		if (txt == 'END') exit
	end do
	close(7)
	END PROGRAM xqrsolv
