	PROGRAM xks2d1s
!	driver for routine ks2d1s
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NMAX=1000
	INTEGER(I4B) :: idum,j,jtrial,n1,ntrial
	REAL(SP) :: d,prob,factor,u,v
	REAL(SP), DIMENSION(NMAX) :: x1,y1
	do
		write(*,*) 'HOW MANY POINTS?'
		read(*,*,END=999) n1
		if (n1 > NMAX) then
			write(*,*) 'n1 too large.'
			cycle
		end if
		do
			write(*,*) 'WHAT FACTOR NONLINEARITY (0 to 1)?'
			read(*,*,END=999) factor
			if (factor < 0.0) then
				write(*,*) 'factor less than 0'
				cycle
			end if
			if (factor > 1.0) then
				write(*,*) 'factor greater than 1'
				cycle
			end if
			write(*,*) 'HOW MANY TRIALS?'
			read(*,*,END=999) ntrial
			idum=-289-ntrial-n1
			do jtrial=1,ntrial
				do j=1,n1
					u=ran(idum)
					u=u*((1.0_sp-factor)+u*factor)
					x1(j)=2.0_sp*u-1.0_sp
					v=ran(idum)
					v=v*((1.0_sp-factor)+v*factor)
					y1(j)=2.0_sp*v-1.0_sp
				end do
				call ks2d1s(x1(1:n1),y1(1:n1),quadvl,d,prob)
				write(*,'(1x,a7,2f12.6)') 'D,PROB= ',d,prob
			end do
		end do
	end do
999	write(*,*) 'NORMAL COMPLETION'
	END PROGRAM xks2d1s
