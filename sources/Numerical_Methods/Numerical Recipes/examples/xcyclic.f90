	PROGRAM xcyclic
!	driver for routine cyclic
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=20
	REAL(SP) :: alpha,betavar,d
	REAL(SP), DIMENSION(N) :: a,b,c,r,x
	REAL(SP), DIMENSION(N,N) :: aa
	INTEGER(I4B), DIMENSION(N) :: indx
	INTEGER(I4B) :: i,idum
	idum= -23
	aa(:,:)=0.0
	do i=1,N
		b(i)=ran(idum)
		r(i)=ran(idum)
		aa(i,i)=b(i)
	end do
	do i=1,N-1
		a(i+1)=ran(idum)
		aa(i+1,i)=a(i+1)
		c(i)=ran(idum)
		aa(i,i+1)=c(i)
	end do
	alpha=ran(idum)
	aa(N,1)=alpha
	betavar=ran(idum)
	aa(1,N)=betavar
	call cyclic(a,b,c,alpha,betavar,r,x)
	call ludcmp(aa,indx,d)
	call lubksb(aa,indx,r)
	do i=1,N
		write(*,*) i,(x(i)-r(i))/(x(i)+r(i))
	end do
	END PROGRAM xcyclic
