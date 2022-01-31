	PROGRAM xbnldev
!	driver for routine bnldev
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=20,NPTS=1000,NN1=20,NN2=100
	INTEGER(I4B) :: i,j
	REAL(SP) :: pp,xm
	REAL(SP), DIMENSION(N+1) :: dist
	call ran_seed(sequence=133)
	do
		write(*,*)&
		'Mean of binomial distribution (0 to 20) (Negative to END)'
		read(*,*) xm
		if (xm < 0) exit
! call short branch of code (n < 25)
		pp=xm/NN1
		dist(:)=0.0
		do i=1,NPTS
			j=int(bnldev(pp,NN1))
			if ((j >= 0) .and. (j <= N)) dist(j+1)=dist(j+1)+1
		end do
		call graphout(dist)
! call other branch of code (n >= 25)
		pp=xm/NN2
		dist(:)=0.0
		do i=1,NPTS
			j=int(bnldev(pp,NN2))
			if ((j >= 0) .and. (j <= N)) dist(j+1)=dist(j+1)+1
		end do
		call graphout(dist)
	end do
	CONTAINS
!BL
	SUBROUTINE graphout(dist)
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: dist
	INTEGER(I4B) :: i,k,n,klim
	INTEGER(I4B), PARAMETER :: ISCAL=200,LLEN=50
	CHARACTER(1), DIMENSION(LLEN) :: text
	n=size(dist)-1
	write(*,'(1x,t5,a,t10,a,t18,a)') 'x','p(x)','graph:'
	dist(1:n)=dist(1:n)/NPTS
	do i=1,n
		text(:)=' '
		text(1)='*'
		klim=min(int(ISCAL*dist(i)),LLEN)
		text(1:klim)='*'
		write(*,'(1x,f5.1,f8.4,3x,50a1)') real(i-1,sp),dist(i),&
			(text(k),k=1,LLEN)
	end do
	END SUBROUTINE graphout
	END PROGRAM xbnldev
