	PROGRAM xtridag
!	driver for routine tridag
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=20
	INTEGER(I4B) :: k,n
	REAL(SP), DIMENSION(NP) :: diag,superd,subd,rhs,u
	CHARACTER :: txt*3
	open(7,file='MATRX2.DAT',status='old')
	do
		read(7,'(a3)') txt
		if (txt == 'END') exit
		read(7,*)
		read(7,*) n
		read(7,*)
		read(7,*) (diag(k), k=1,n)
		read(7,*)
		read(7,*) (superd(k), k=1,n-1)
		read(7,*)
		read(7,*) (subd(k), k=2,n)
		read(7,*)
		read(7,*) (rhs(k), k=1,n)
!	carry out solution
!	use tridag_par since it also tests tridag_ser!
		call tridag_par(subd(2:n),diag(1:n), &
			superd(1:n-1),rhs(1:n),u(1:n))
		write(*,*) 'The solution vector is:'
		write(*,'(1x,6f12.6)') (u(k), k=1,n)
!	test solution
		write(*,*) '(matrix)*(sol''n vector) should be:'
		write(*,'(1x,6f12.6)') (rhs(k), k=1,n)
		write(*,*) 'Actual result is:'
		rhs(1)=diag(1)*u(1) + superd(1)*u(2)
		rhs(2:n-1)=subd(2:n-1)*u(1:n-2) + diag(2:n-1)*u(2:n-1) +&
					superd(2:n-1)*u(3:n)
		rhs(n)=subd(n)*u(n-1) + diag(n)*u(n)
		write(*,'(1x,6f12.6)') (rhs(k), k=1,n)
		write(*,*) '***********************************'
		write(*,*) 'Press RETURN for next problem:'
		read(*,*)
	end do
	close(7)
	END PROGRAM xtridag
