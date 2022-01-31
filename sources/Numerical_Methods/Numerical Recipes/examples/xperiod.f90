	PROGRAM xperiod
!	driver for routine period
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=90,NPR=11
	INTEGER(I4B) :: j,jmax,n,nout
	REAL(SP) :: harvest,prob
	REAL(SP), DIMENSION(:), POINTER :: px,py
	REAL(SP), DIMENSION(NP) :: x,y
	call ran_seed(sequence=2099)
	j=0
	do n=1,NP+10
		if (n /= 3 .and. n /= 4 .and. n /= 6 .and. n /= 21 .and. &
			n /= 38 .and. n /= 51 .and. n /= 67 .and. &
			n /= 68 .and. n /= 83 .and. n /= 93) then
			j=j+1
			x(j)=n
			call gasdev(harvest)
			y(j)=0.75_sp*cos(0.6_sp*x(j))+harvest
		end if
	end do
	call period(x(1:j),y(1:j),4.0_sp,1.0_sp,px,py,jmax,prob)
	nout=size(px)
	write(*,*) 'PERIOD results for test signal (cos(0.6x) + noise):'
	write(*,*) 'NOUT,JMAX,PROB=',nout,jmax,prob
	do n=max(1,jmax-NPR/2),min(nout,jmax+NPR/2)
		write(*,*) n,TWOPI*px(n),py(n)
	end do
	END PROGRAM xperiod
