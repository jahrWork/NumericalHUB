	PROGRAM xgamdev
!	driver for routine gamdev
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=20,NPTS=10000,ISCAL=200,LLEN=50
	INTEGER(I4B) :: i,ia,j,k,klim
	REAL(SP), DIMENSION(N+1) :: dist
	CHARACTER(1), DIMENSION(50) :: text
	call ran_seed(sequence=1776)
	do
		do
			dist(:)=0.0
			write(*,*) 'Order of Gamma distribution (n=1..20); -1 to END.'
			read(*,*) ia
			if (ia <= 20) exit
		end do
		if (ia <= 0) exit
		do i=1,NPTS
			j=int(gamdev(ia))+1
			if ((j >= 1) .and. (j <= N+1)) dist(j)=dist(j)+1
		end do
		write(*,'(1x,a,i2,a,i6,a)') 'Gamma-distribution deviate, order ',&
			ia,' of ',NPTS,' points'
		write(*,'(1x,t6,a,t14,a,t23,a)') 'x','p(x)','graph:'
		dist(1:N)=dist(1:N)/NPTS
		do i=1,N
			text(:)=' '
			klim=min(int(ISCAL*dist(i)),LLEN)
			text(1:klim)='*'
			write(*,'(1x,f7.2,f10.4,4x,50a1)')&
				real(i,sp),dist(i),(text(k),k=1,50)
		end do
	end do
	END PROGRAM xgamdev
