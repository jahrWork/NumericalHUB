	PROGRAM xgaussj
!	driver for routine gaussj
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: MP=20,NP=20
	INTEGER(I4B) :: i,j,m,n
	REAL(SP), DIMENSION(NP,MP) :: b,x
	REAL(SP), DIMENSION(NP,NP) :: a,ai
	REAL(SP), DIMENSION(NP,MP) :: t
	REAL(SP), DIMENSION(NP,NP) :: u
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
!	save Matrices for later testing of results
		ai(1:n,1:n)=a(1:n,1:n)
		x(1:n,1:m)=b(1:n,1:m)
!	invert Matrix
		call gaussj(ai(1:n,1:n),x(1:n,1:m))
		write(*,*) 'Inverse of Matrix A : '
		do i=1,n
			write(*,'(1x,(6f12.6))') (ai(i,j), j=1,n)
		end do
!	test Results
!	check Inverse
		write(*,*) 'A times A-inverse (compare with unit matrix)'
		u(1:n,1:n)=matmul(ai(1:n,1:n),a(1:n,1:n))
		do i=1,n
			write(*,'(1x,(6f12.6))') (u(i,j), j=1,n)
		end do
!	check Vector Solutions
		write(*,*) 'Check the following vectors for equality:'
		write(*,'(t12,a8,t23,a12)') 'Original','Matrix*Sol''n'
		t(1:n,1:m)=matmul(a(1:n,1:n),x(1:n,1:m))
		do j=1,m
			write(*,'(1x,a,i2,a)') 'Vector ',j,':'
			do i=1,n
				write(*,'(8x,2f12.6)') b(i,j),t(i,j)
			end do
		end do
		write(*,*) '***********************************'
		write(*,*) 'Press RETURN for next problem:'
		read(*,*)
	end do
	close(7)
	END PROGRAM xgaussj
