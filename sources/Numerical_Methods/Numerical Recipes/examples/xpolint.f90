	PROGRAM xpolint
!	driver for routine polint
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=10
	INTEGER(I4B) :: i,n,nfunc
	REAL(SP) :: dy,f,x,y
	REAL(SP), DIMENSION(NP) :: xa,ya
	write(*,*) 'Generation of interpolation tables'
	write(*,*) ' ... sin(x)    0<x<PI'
	write(*,*) ' ... exp(x)    0<x<1 '
	write(*,*) 'How many entries go in these tables? (note: N<10)'
	read(*,*) n
	do nfunc=1,2
		if (nfunc == 1) then
			write(*,*) 'sine function from 0 to PI'
			xa(1:NP)=arth(1,1,NP)*PI/NP
			ya(:)=sin(xa(:))
		else if (nfunc == 2) then
			write(*,*) 'exponential function from 0 to 1'
			xa(1:NP)=arth(1.0_sp,1.0_sp,NP)/NP
			ya(:)=exp(xa(:))
		else
			stop
		end if
		write(*,'(t10,a1,t20,a4,t28,a12,t46,a5)')&
			'x','f(x)','interpolated','error'
		do i=1,10
			if (nfunc == 1) then
				x=(-0.05_sp+i/10.0_sp)*PI
				f=sin(x)
			else if (nfunc == 2) then
				x=(-0.05_sp+i/10.0_sp)
				f=exp(x)
			end if
			call polint(xa,ya,x,y,dy)
			write(*,'(1x,3f12.6,e15.4)') x,f,y,dy
		end do
		write(*,*) '***********************************'
		write(*,*) 'Press RETURN'
		read(*,*)
	end do
	END PROGRAM xpolint
