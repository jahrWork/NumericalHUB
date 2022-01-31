	PROGRAM xbessik
!	driver for routine bessik
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i,nval
	REAL(SP) :: ri,rk,rip,rkp,x,xnu,xri,xrk,xrip,xrkp
	CHARACTER(25) :: text
	open(7,file='FNCVAL.DAT',status='OLD')
	do
		read(7,'(a)') text
		if (text == 'Modified Bessel Functions') exit
	end do
	read(7,*) nval
	write(*,*) text
	write(*,'(1x,t3,a3,t8,a1)') 'XNU','X'
	write(*,'(1x,t5,a2,t21,a2,t37,a3,t53,a3)') 'RI','RK','RIP','RKP'
	write(*,'(1x,t5,a3,t21,a3,t37,a4,t53,a4)') 'XRI','XRK','XRIP','XRKP'
	do i=1,nval
		read(7,*) xnu,x,ri,rk,rip,rkp
		call bessik(x,xnu,xri,xrk,xrip,xrkp)
		write(*,'(2f6.2,/,1p,4e16.6,/,1p,4e16.6)')&
			xnu,x,ri,rk,rip,rkp,xri,xrk,xrip,xrkp
	end do
	close(7)
	END PROGRAM xbessik
