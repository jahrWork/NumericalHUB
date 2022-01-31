	PROGRAM xhypgeo
!	driver for routine hypgeo
	USE nrtype
	USE nr
	IMPLICIT NONE
	COMPLEX(SPC) :: a,b,c,z,zi,q1,q2,q3,q4
	REAL(SP) :: x,y
	a=0.5_sp
	b=1.0
	c=1.5_sp
	do
		write(*,*) 'INPUT X,Y OF COMPLEX ARGUMENT:'
		read(*,*,END=999) x,y
		z=cmplx(x,y,kind=spc)
		q1=hypgeo(a,b,c,z*z)
		q2=0.5_sp*log((1.0_sp+z)/(1.0_sp-z))/z
		q3=hypgeo(a,b,c,-z*z)
		zi=cmplx(-y,x,kind=spc)
		q4=0.5_sp*log((1.0_sp+zi)/(1.0_sp-zi))/zi
		write(*,*) '2F1(0.5,1.0,1.5;Z**2) =',q1
		write(*,*) 'check using log form   ',q2
		write(*,*) '2F1(0.5,1.0,1.5;-Z**2)=',q3
		write(*,*) 'check using log form   ',q4
	end do
999	END PROGRAM xhypgeo
