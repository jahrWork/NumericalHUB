	PROGRAM xbandec
!	driver for routine bandec
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP) :: d
	REAL(SP), DIMENSION(7) :: x,b
	REAL(SP), DIMENSION(7,2) :: al
	REAL(SP), DIMENSION(7,4) :: a
	INTEGER(I4B), DIMENSION(7) :: indx
	INTEGER(I4B) :: i,idum,j
	idum=-1
	do i=1,7
		x(i)=ran(idum)
		do j=1,4
			a(i,j)=ran(idum)
		end do
	end do
	call banmul(a(1:7,1:4),2,1,x,b)
	do i=1,7
		write(*,*) i,b(i),x(i)
	end do
	call bandec(a(1:7,1:4),2,1,al(1:7,1:2),indx(1:7),d)
	call banbks(a(1:7,1:4),2,1,al(1:7,1:2),indx(1:7),b)
	do i=1,7
		write(*,*) i,b(i),x(i)
	end do
	END PROGRAM xbandec
