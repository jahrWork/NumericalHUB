	PROGRAM xpoidev
!	driver for routine poidev
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=20,NPTS=10000,ISCAL=200,LLEN=50
	INTEGER(I4B) :: i,j,k,klim
	REAL(SP) :: xm
	REAL(SP), DIMENSION(N+1) :: dist
	CHARACTER(1), DIMENSION(50) :: text
	call ran_seed(sequence=15)
	do
		do
			dist(:)=0.0
			write(*,*) 'Mean of Poisson distrib. (x=0 to 20); neg. to END'
			read(*,*) xm
			if (xm <= 20.0) exit
		end do
		if (xm < 0.0) exit
		do i=1,NPTS
			j=int(poidev(xm))+1
			if ((j >= 1) .and. (j <= N+1)) dist(j)=dist(j)+1
		end do
		write(*,'(1x,a,f5.2,a,i6,a)') &
			'Poisson-distributed deviate, mean ',&
			xm,' of ',NPTS,' points'
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
	END PROGRAM xpoidev
