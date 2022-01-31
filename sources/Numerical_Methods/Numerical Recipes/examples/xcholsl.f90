	PROGRAM xcholsl
!	driver for routine cholsl
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=3
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(N) :: p,x,b = (/ 0.4_sp,0.02_sp,99.0_sp /)
	REAL(SP), DIMENSION(N,N) :: aorig,atest,chol,a = reshape( (/ &
		100.0_sp,15.0_sp,0.01_sp,15.0_sp,2.3_sp,0.01_sp,&
		0.01_sp,0.01_sp,1.0_sp /), (/ N,N /) )
	aorig(:,:)=a(:,:)
	call choldc(a,p)
	do i=1,N
		chol(i,1:i-1)=a(i,1:i-1)
		chol(i,i)=p(i)
		chol(i,i+1:N)=0.0
	end do
	atest(:,:)=matmul(chol,transpose(chol))
	write(*,*) 'Original matrix:'
	write(*,'(1p,3e16.6)') ((aorig(i,j),j=1,N),i=1,N)
	write(*,*)
	write(*,*) 'Product of Cholesky factors:'
	write(*,'(1p,3e16.6)') ((atest(i,j),j=1,N),i=1,N)
	write(*,*)
	call cholsl(a,p,b,x)
	p=matmul(aorig,x)
	write(*,*) 'Check solution vector:'
	write(*,'(1p,2e16.6)') (p(i),b(i),i=1,N)
	END PROGRAM xcholsl
