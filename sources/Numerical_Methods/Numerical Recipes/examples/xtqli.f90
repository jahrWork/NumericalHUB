	PROGRAM xtqli
!	driver for routine tqli
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=10
	REAL(SP), PARAMETER :: TINY=1.0e-6_sp
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(NP) :: d,e,f
	REAL(SP), DIMENSION(NP,NP) :: a,c = reshape( (/ &
		5.0_sp,4.3_sp,3.0_sp,2.0_sp,1.0_sp,&
		0.0_sp,-1.0_sp,-2.0_sp,-3.0_sp,-4.0_sp,&
		4.3_sp,5.1_sp,4.0_sp,3.0_sp,2.0_sp,&
		1.0_sp,0.0_sp,-1.0_sp,-2.0_sp,-3.0_sp,&
		3.0_sp,4.0_sp,5.0_sp,4.0_sp,3.0_sp,&
		2.0_sp,1.0_sp,0.0_sp,-1.0_sp,-2.0_sp,&
		2.0_sp,3.0_sp,4.0_sp,5.0_sp,4.0_sp,&
		3.0_sp,2.0_sp,1.0_sp,0.0_sp,-1.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,&
		4.0_sp,3.0_sp,2.0_sp,1.0_sp,0.0_sp,&
		0.0_sp,1.0_sp,2.0_sp,3.0_sp,4.0_sp,&
		5.0_sp,4.0_sp,3.0_sp,2.0_sp,1.0_sp,&
		-1.0_sp,0.0_sp,1.0_sp,2.0_sp,3.0_sp,&
		4.0_sp,5.0_sp,4.0_sp,3.0_sp,2.0_sp,&
		-2.0_sp,-1.0_sp,0.0_sp,1.0_sp,2.0_sp,&
		3.0_sp,4.0_sp,5.0_sp,4.0_sp,3.0_sp,&
		-3.0_sp,-2.0_sp,-1.0_sp,0.0_sp,1.0_sp,&
		2.0_sp,3.0_sp,4.0_sp,5.0_sp,4.0_sp,&
		-4.0_sp,-3.0_sp,-2.0_sp,-1.0_sp,0.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp /),&
		(/ NP,NP /) )
	a(:,:)=c(:,:)
	call tred2(a,d,e)
	call tqli(d,e,a)
	write(*,'(/1x,a)') 'Eigenvectors for a real symmetric matrix'
	do i=1,NP
		f=matmul(c(:,:),a(:,i))
		write(*,'(/1x,a,i3,a,f10.6)') 'Eigenvalue',i,' =',d(i)
		write(*,'(/1x,t7,a,t17,a,t31,a)') 'Vector','Mtrx*Vect.','Ratio'
		do j=1,NP
			if (abs(a(j,i)) < TINY) then
				write(*,'(1x,2f12.6,a12)') a(j,i),f(j),'div. by 0'
			else
				write(*,'(1x,2f12.6,e14.6)') a(j,i),f(j),&
					f(j)/a(j,i)
			end if
		end do
		write(*,'(/1x,a)') 'press ENTER to continue...'
		read(*,*)
	end do
	END PROGRAM xtqli
