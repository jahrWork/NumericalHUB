	PROGRAM xsvbksb
!	driver for routine svbksb, which calls routine svdcmp
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: MP=20,NP=20
	INTEGER(I4B) :: i,j,m,n
	REAL(SP), DIMENSION(NP) :: w
	REAL(SP), DIMENSION(NP,MP) :: b
	REAL(SP), DIMENSION(NP,NP) :: a,u
	REAL(SP), DIMENSION(NP) :: c,x
	REAL(SP), DIMENSION(NP,NP) :: v
	REAL(SP) :: wmax,wmin
	CHARACTER(3) :: dummy
	open(7,file='MATRX1.DAT',status='old')
	do
		read(7,'(a)') dummy
		if (dummy == 'END') exit
		read(7,*)
		read(7,*) n,m
		read(7,*)
		read(7,*) ((a(i,j), j=1,n), i=1,n)
		read(7,*)
		read(7,*) ((b(i,j), i=1,n), j=1,m)
!	copy a into u
		u(1:n,1:n)=a(1:n,1:n)
!	decompose matrix a
		call svdcmp(u(1:n,1:n),w(1:n),v(1:n,1:n))
!	find maximum singular value
		wmax=max(maxval(w(1:n)),0.0_sp)
!	define "small"
		wmin=wmax*(1.0e-6)
!	zero the "small" singular values
		where (w(1:n) < wmin) w(1:n)=0.0
!	backsubstitute for each right-hand side vector
		do j=1,m
			write(*,'(1x,a,i2)') 'Vector number ',j
			c(1:n)=b(1:n,j)
			call svbksb(u(1:n,1:n),w(1:n),v(1:n,1:n),c(1:n),x(1:n))
			write(*,*) '    Solution vector is:'
			write(*,'(1x,6f12.6)') (x(i), i=1,n)
			write(*,*) '    Original right-hand side vector:'
			write(*,'(1x,6f12.6)') (c(i), i=1,n)
			write(*,*) '    Result of (matrix)*(sol''n vector):'
			c(1:n)=matmul(a(1:n,1:n),x(1:n))
			write(*,'(1x,6f12.6)') (c(i), i=1,n)
		end do
		write(*,*) '***********************************'
		write(*,*) 'Press RETURN for next problem'
		read(*,*)
	end do
	close(7)
	END PROGRAM xsvbksb
