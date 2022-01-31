	PROGRAM xjacobi
!	driver for routine jacobi
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=10,NMAT=3
	INTEGER(I4B) :: i,j,k,kk,l,ll,nrot
	INTEGER(I4B), DIMENSION(NMAT) :: num = (/ 3,5,10 /)
	REAL(SP) :: ratio
	REAL(SP), DIMENSION(NP) :: d,r
	REAL(SP), DIMENSION(NP,NP) :: v
	REAL(SP), DIMENSION(NP,NP) :: e
	REAL(SP), DIMENSION(3,3) :: a = reshape( (/ &
	1.0_sp,2.0_sp,3.0_sp,2.0_sp,2.0_sp,3.0_sp,3.0_sp,3.0_sp,3.0_sp/),&
	(/ 3,3 /) )
	REAL(SP), DIMENSION(5,5) :: b = reshape( (/ &
		-2.0_sp,-1.0_sp,0.0_sp,1.0_sp,2.0_sp,&
		-1.0_sp,-1.0_sp,0.0_sp,1.0_sp,2.0_sp,&
		0.0_sp,0.0_sp,0.0_sp,1.0_sp,2.0_sp,&
		1.0_sp,1.0_sp,1.0_sp,1.0_sp,2.0_sp,&
		2.0_sp,2.0_sp,2.0_sp,2.0_sp,2.0_sp /), (/ 5,5 /) )
	REAL(SP), DIMENSION(10,10) :: c = reshape( (/ &
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
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp /), (/ 10,10 /) )
	do i=1,NMAT
		if (i == 1) then
			e(1:num(i),1:num(i))=a(1:num(i),1:num(i))
		else if (i == 2) then
			e(1:num(i),1:num(i))=b(1:num(i),1:num(i))
		else if (i == 3) then
			e(1:num(i),1:num(i))=c(1:num(i),1:num(i))
		end if
		call jacobi(e(1:num(i),1:num(i)),d(1:num(i)),v(1:num(i),1:num(i)),nrot)
		write(*,'(/1x,a,i2)') 'Matrix Number ',i
		write(*,'(1x,a,i3)') 'Number of JACOBI rotations: ',nrot
		write(*,'(/1x,a)') 'Eigenvalues:'
		do j=1,num(i)
			write(*,'(1x,5f12.6)') d(j)
		end do
		write(*,'(/1x,a)') 'Eigenvectors:'
		do j=1,num(i)
			write(*,'(1x,t5,a,i3)') 'Number',j
			write(*,'(1x,5f12.6)') (v(k,j),k=1,num(i))
		end do
!	eigenvector test
		write(*,'(/1x,a)') 'Eigenvector Test'
		do j=1,num(i)
			do l=1,num(i)
				r(l)=0.0
				do k=1,num(i)
					if (k > l) then
						kk=l
						ll=k
					else
						kk=k
						ll=l
					end if
					if (i == 1) then
						r(l)=r(l)+a(ll,kk)*v(k,j)
					else if (i == 2) then
						r(l)=r(l)+b(ll,kk)*v(k,j)
					else if (i == 3) then
						r(l)=r(l)+c(ll,kk)*v(k,j)
					end if
				end do
			end do
			write(*,'(/1x,a,i3)') 'Vector Number',j
			write(*,'(/1x,t7,a,t18,a,t31,a)')&
				'Vector','Mtrx*Vec.','Ratio'
			do l=1,num(i)
				ratio=r(l)/v(l,j)
				write(*,'(1x,3f12.6)') v(l,j),r(l),ratio
			end do
		end do
		write(*,*) 'press RETURN to continue...'
		read(*,*)
	end do
	END PROGRAM xjacobi
