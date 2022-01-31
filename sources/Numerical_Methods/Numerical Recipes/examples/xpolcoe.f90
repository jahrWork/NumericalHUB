	PROGRAM xpolcoe
!	driver for routine polcoe
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=5
	INTEGER(I4B) :: i,nfunc
	REAL(SP) :: f,x
	REAL(SP), DIMENSION(NP) :: xa,ya,coeff
	do nfunc=1,2
		if (nfunc == 1) then
			write(*,*) 'Sine function from 0 to PI'
			xa(1:NP)=arth(1,1,NP)*PI/NP
			ya(:)=sin(xa(:))
		else if (nfunc == 2) then
			write(*,*) 'Exponential function from 0 to 1'
			xa(1:NP)=arth(1.0_sp,1.0_sp,NP)/NP
			ya(:)=exp(xa(:))
		else
			stop
		end if
		coeff(:)=polcoe(xa(:),ya(:))
		write(*,*) '    coefficients'
		write(*,'(1x,6f12.6)') (coeff(i),i=1,NP)
		write(*,'(1x,t10,a1,t20,a4,t29,a10)') 'x','f(x)','polynomial'
		do i=1,10
			if (nfunc == 1) then
				x=(-0.05_sp+i/10.0_sp)*PI
				f=sin(x)
			else if (nfunc == 2) then
				x=-0.05_sp+i/10.0_sp
				f=exp(x)
			end if
			write(*,'(1x,3f12.6)') x,f,poly(x,coeff)
		end do
		write(*,*) '***********************************'
		write(*,*) 'Press RETURN'
		read(*,*)
	end do
	END PROGRAM xpolcoe
