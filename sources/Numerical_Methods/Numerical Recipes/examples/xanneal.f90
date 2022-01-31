	PROGRAM xanneal
!	driver for routine anneal
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NCITY=10
	INTEGER(I4B) :: i,ii
	INTEGER(I4B), DIMENSION(NCITY) :: iorder
	REAL(SP) :: harvest
	REAL(SP), DIMENSION(NCITY) :: x,y
!	create points of sale
	call ran_seed(sequence=113)
	do i=1,NCITY
		call ran1(harvest)
		x(i)=harvest
		call ran1(harvest)
		y(i)=harvest
	end do
	iorder(1:NCITY)=arth(1,1,NCITY)
	call anneal(x,y,iorder)
	write(*,*) '*** System Frozen ***'
	write(*,*) 'Final path:'
	write(*,'(1x,t3,a,t13,a,t23,a)') 'city','x','y'
	do i=1,NCITY
		ii=iorder(i)
		write(*,'(1x,i4,2f10.4)') ii,x(ii),y(ii)
	end do
	END PROGRAM xanneal
