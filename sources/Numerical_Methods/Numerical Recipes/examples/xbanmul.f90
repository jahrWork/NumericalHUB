	PROGRAM xbanmul
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
!	driver for routine banmul
	INTEGER(I4B), PARAMETER :: NP=7,M1=2,M2=1,MP=M1+1+M2
	INTEGER(I4B) :: i,j,k
	INTEGER(I4B), DIMENSION(NP) :: itmp
	REAL(SP), DIMENSION(NP) :: ax,b,x
	REAL(SP), DIMENSION(NP,MP) :: a
	REAL(SP), DIMENSION(NP,NP) :: aa
	itmp(1:NP)=arth(1,1,NP)
	a(1:NP,1:M1)=outersum(10.0_sp*itmp(1:NP),arth(1.0_sp,1.0_sp,M1))
!	lower band
	a(1:NP,M1+1)=itmp(1:NP)
!	diagonal
	a(1:NP,M1+2:M1+1+M2)=outersum(0.1_sp*itmp(1:NP),arth(1.0_sp,1.0_sp,M2))
!	upper band
	aa(:,:)=0.0
	do i=1,NP
		do j=1,NP
			k=i-M1-1
			if (j >= max(1,1+k) .and. j <= min(M1+M2+1+k,NP)) &
				aa(i,j)=a(i,j-k)
		end do
	end do
	x(1:NP)=itmp(1:NP)/10.0_sp
	call banmul(a(1:NP,1:M1+M2+1),M1,M2,x,b)
	ax=matmul(aa,x)
	write(*,'(t8,a,t32,a)') 'Reference vector','banmul vector'
	do i=1,NP
		write(*,'(t8,f12.4,t32,f12.4)') ax(i),b(i)
	end do
	END PROGRAM xbanmul
