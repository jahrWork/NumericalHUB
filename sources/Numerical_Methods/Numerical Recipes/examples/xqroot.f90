	PROGRAM xqroot
!	driver for routine qroot
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=7,NTRY=10
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp,TINY=1.0e-5_sp
	INTEGER(I4B) :: i,nroot
	REAL(SP), DIMENSION(N) :: p = &
		(/ 10.0_sp,-18.0_sp,25.0_sp,-24.0_sp,16.0_sp,-6.0_sp,1.0_sp /)
	REAL(SP), DIMENSION(NTRY) :: b,c
	write(*,'(/1x,a)') 'P(x)=x^6-6x^5+16x^4-24x^3+25x^2-18x+10'
	write(*,'(1x,a)') 'Quadratic factors x^2+Bx+C'
	write(*,'(/1x,a,t15,a,t27,a/)') 'Factor','B','C'
	nroot=0
	c(1:NTRY)=0.5_sp*arth(1,1,NTRY)
	b(1:NTRY)=-c(1:NTRY)
	do i=1,NTRY
		call qroot(p,b(i),c(i),EPS)
		if (nroot == 0) then
			write(*,'(1x,i3,2x,2f12.6)') nroot,b(i),c(i)
			nroot=1
		else
			if (any(abs(b(i)-b(1:nroot)) < TINY .and. abs(c(i)-c(1:nroot)) < TINY)) cycle
			write(*,'(1x,i3,2x,2f12.6)') nroot,b(i),c(i)
			nroot=nroot+1
		end if
	end do
	END PROGRAM xqroot
