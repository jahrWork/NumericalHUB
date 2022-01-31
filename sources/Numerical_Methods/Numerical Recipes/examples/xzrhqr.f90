	PROGRAM xzrhqr
!	driver for routine zrhqr
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: M=4,MP1=M+1,NTRY=21
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(M) :: rtr,rti
	REAL(SP), DIMENSION(MP1) :: a = (/-1.0_sp,0.0_sp,0.0_sp,0.0_sp,1.0_sp/)
	write(*,'(/1x,a)') 'Roots of polynomial x^4-1'
	write(*,'(/1x,t16,a,t29,a/)') 'Real','Complex'
	call zrhqr(a,rtr,rti)
	do i=1,M
		write(*,'(1x,i5,2f15.6)') i,rtr(i),rti(i)
	end do
	END PROGRAM xzrhqr
