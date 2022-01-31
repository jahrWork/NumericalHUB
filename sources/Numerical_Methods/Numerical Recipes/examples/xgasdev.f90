	PROGRAM xgasdev
!	driver for routine gasdev
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: &
		N=20,NP1=N+1,NOVER2=N/2,NPTS=10000,ISCAL=400,LLEN=50
	INTEGER(I4B) :: i,j,klim
	REAL(SP) :: harvest,x
	REAL(SP), DIMENSION(NP1) :: dist
	CHARACTER text(50)*1
	call ran_seed(sequence=1976)
	dist(:)=0.0
	do i=1,NPTS
		call gasdev(harvest)
		j=nint(0.25_sp*N*harvest)+NOVER2+1
		if ((j >= 1) .and. (j <= NP1)) dist(j)=dist(j)+1
	end do
	write(*,'(1x,a,i6,a)') &
		'Normally distributed deviate of ',NPTS,' points'
	write(*,'(1x,t6,a,t14,a,t23,a)') 'x','p(x)','graph:'
	dist(:)=dist(:)/NPTS
	do j=1,NP1
		text(1:50)=' '
		klim=int(ISCAL*dist(j))
		if (klim > LLEN) klim=LLEN
		text(1:klim)='*'
		x=float(j)/(0.25_sp*N)
		write(*,'(1x,f7.2,f10.4,4x,50a1)') &
			x,dist(j),(text(i),i=1,50)
	end do
	END PROGRAM xgasdev
