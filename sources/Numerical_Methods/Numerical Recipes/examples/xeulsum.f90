	PROGRAM xeulsum
!	driver for routine eulsum
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NVAL=40
	INTEGER(I4B) :: i,j,mval
	REAL(SP) :: sersum,term,x,xpower
!	evaluate ln(1+x)=x-x^2/2+x^3/3-x^4/4.0_sp..  for -1<x<1
	do
		write(*,*) 'How many terms in polynomial?'
		write(*,'(1x,a,i2,a)') 'Enter n between 1 and ',NVAL,&
			'. Enter n=0 to END.'
		read(*,*) mval
		if ((mval <= 0) .or. (mval > NVAL)) exit
		write(*,'(1x,t9,a1,t18,a6,t28,a10)') 'X','Actual','Polynomial'
		do i=-8,8,1
			x=i/10.0_sp
			sersum=0.0
			xpower=-1
			do j=1,mval
				xpower=-x*xpower
				term=xpower/j
				call eulsum(sersum,term,j)
			end do
			write(*,'(3f12.6)') x,log(1.0+x),sersum
		end do
	end do
	END PROGRAM xeulsum
