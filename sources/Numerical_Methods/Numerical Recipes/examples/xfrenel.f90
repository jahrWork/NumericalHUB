	PROGRAM xfrenel
!	driver for routine frenel
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i,nval
	REAL(SP) :: c,s,x,xc,xs
	CHARACTER(17) :: text
	open(7,file='FNCVAL.DAT',status='OLD')
	do
		read(7,'(a)') text
		if (text == 'Fresnel Integrals') exit
	end do
	read(7,*) nval
	write(*,*) text
	write(*,'(1x,t5,a1,t12,a6,t28,a4,t42,a6,t58,a4)')&
		'X','Actual','S(X)','Actual','C(X)'
	do i=1,nval
		read(7,*) x,xs,xc
		call frenel(x,s,c)
		write(*,'(f6.2,4e15.6)') x,xs,s,xc,c
	end do
	close(7)
	END PROGRAM xfrenel
