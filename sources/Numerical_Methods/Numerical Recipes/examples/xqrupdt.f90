	PROGRAM xqrupdt
!	driver for routine qrdupdt
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=20
	INTEGER(I4B) :: i,j,m,n
	REAL(SP) :: con
	REAL(SP), DIMENSION(NP) :: c,d,u,v
	REAL(SP), DIMENSION(NP,NP) :: a,au,q,qt,r,s,x
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
		read(7,*) ((s(i,j), i=1,n), j=1,m)
!	print out a-matrix for comparison with product of Q and R
!	decomposition matrices.
		write(*,*) 'Original matrix:'
		do i=1,n
			write(*,'(1x,6f12.6)') (a(i,j), j=1,n)
		end do
!	updated matrix we'll use later
		au(1:n,1:n)=a(1:n,1:n)+outerprod(s(1:n,1),s(1:n,2))
!	perform the initial decomposition
		call qrdcmp(a(1:n,1:n),c(1:n),d(1:n),sing)
		if (sing) write(*,*) 'Singularity in QR decomposition.'
!	find the Q and R matrices
		do i=1,n
			r(i,i+1:n)=a(i,i+1:n)
			q(i,i+1:n)=0.0
			r(i,1:i-1)=0.0
			q(i,1:i-1)=0.0
			r(i,i)=d(i)
			q(i,i)=1.0
		end do
		do i=n-1,1,-1
			con=dot_product(a(i:n,i),a(i:n,i))/2.0_sp
			qt(i:n,i:n)=outerprod(a(i:n,i),matmul(a(i:n,i),q(i:n,i:n)))/con
			q(i:n,i:n)=q(i:n,i:n)-qt(i:n,i:n)
		end do
!	compute product of Q and R matrices for comparison with original matrix.
		x(1:n,1:n)=matmul(q(1:n,1:n),r(1:n,1:n))
		write(*,*) 'Product of Q and R matrices:'
		do i=1,n
			write(*,'(1x,6f12.6)') (x(i,j), j=1,n)
		end do
		write(*,*) 'Q matrix of the decomposition:'
		do i=1,n
			write(*,'(1x,6f12.6)') (q(i,j), j=1,n)
		end do
		write(*,*) 'R matrix of the decomposition:'
		do i=1,n
			write(*,'(1x,6f12.6)') (r(i,j), j=1,n)
		end do
!	Q transpose
		qt(1:n,1:n)=transpose(q(1:n,1:n))
		v(1:n)=s(1:n,2)
		u(1:n)=0.0
		u(1:n)=matmul(qt(1:n,1:n),s(1:n,1))
		call qrupdt(r(1:n,1:n),qt(1:n,1:n),u(1:n),v(1:n))
		x(1:n,1:n)=matmul(transpose(qt(1:n,1:n)),r(1:n,1:n))
		write(*,*) 'Updated matrix:'
		do i=1,n
			write(*,'(1x,6f12.6)') (au(i,j), j=1,n)
		end do
		write(*,*) 'Product of new Q and R matrices:'
		do i=1,n
			write(*,'(1x,6f12.6)') (x(i,j), j=1,n)
		end do
		write(*,*) 'New Q matrix'
		do i=1,n
			write(*,'(1x,6f12.6)') (qt(j,i), j=1,n)
		end do
		write(*,*) 'New R matrix'
		do i=1,n
			write(*,'(1x,6f12.6)') (r(i,j), j=1,n)
		end do
		write(*,*) '***********************************'
		write(*,*) 'Press RETURN for next problem:'
		read(*,*)
		read(7,'(a3)') txt
		if (txt == 'END') exit
	end do
	close(7)
	END PROGRAM xqrupdt
