	PROGRAM xpzextr
!	driver for routine pzextr
	USE nrtype
	USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NV=4
	INTEGER(I4B) :: i,iest,j
	REAL(SP) :: xest
	REAL(SP), DIMENSION(NV) :: yest,yz,dy
	do i=1,10
		iest=i
		xest=1.0_sp/real(i,sp)
		yest(1:NV)=(1.0_sp-xest+xest**3)/geop(xest+1.0_sp,xest+1.0_sp,NV)
		call pzextr(iest,xest,yest,yz,dy)
		write(*,'(/1x,a,i2)') 'I = ',i
		write(*,'(1x,a,4f12.6)') 'Extrap. function: ',(yz(j),j=1,NV)
		write(*,'(1x,a,4f12.6)') 'Estimated error:  ',(dy(j),j=1,NV)
	end do
	write(*,'(/1x,a,4f12.6)') 'Actual values:    ',1.0,1.0,1.0,1.0
	END PROGRAM xpzextr
