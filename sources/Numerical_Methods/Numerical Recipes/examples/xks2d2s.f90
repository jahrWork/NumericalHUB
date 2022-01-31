	PROGRAM xks2d2s
!	driver for routine ks2d2s
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NMAX=1000
	INTEGER(I4B) :: j,jtrial,n1,n2,ntrial
	REAL(SP) :: d,harvest,prob,shrink,u,v
	REAL(SP), DIMENSION(NMAX) :: x1,y1,x2,y2
	do
		write(*,*) 'INPUT N1,N2'
		read(*,*,END=999) n1,n2
		if (n1 > NMAX) then
			write(*,*) 'n1 too large.'
			cycle
		end if
		if (n2 > NMAX) then
			write(*,*) 'n2 too large.'
			cycle
		end if
		write(*,*) 'WHAT SHRINKAGE?'
		read(*,*,END=999) shrink
		write(*,*) 'HOW MANY TRIALS?'
		read(*,*,END=999) ntrial
		if (ntrial > NMAX) then
			write(*,*) 'Too many trials.'
			cycle
		end if
		call ran_seed(sequence=287+ntrial+n1+n2)
		do jtrial=1,ntrial
			do j=1,n1
				call gasdev(harvest)
				u=harvest
				call gasdev(harvest)
				v=harvest*shrink
				x1(j)=u+v
				y1(j)=u-v
			end do
			do j=1,n2
				call gasdev(harvest)
				u=harvest*shrink
				call gasdev(harvest)
				v=harvest
				x2(j)=u+v
				y2(j)=u-v
			end do
			call ks2d2s(x1(1:n1),y1(1:n1),x2(1:n2),y2(1:n2),d,prob)
			write(*,'(1x,a7,2f12.6)') 'D,PROB= ',d,prob
		end do
	end do
999	write(*,*) 'NORMAL COMPLETION'
	END PROGRAM xks2d2s
