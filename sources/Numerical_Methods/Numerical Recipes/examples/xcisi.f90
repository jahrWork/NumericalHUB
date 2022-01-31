	PROGRAM xcisi
!	driver for routine cisi
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i,nval
	REAL(SP) :: ci,si,x,xci,xsi
	CHARACTER(25) :: text
	open(7,file='FNCVAL.DAT',status='OLD')
	do
		read(7,'(a)') text
		if (text == 'Cosine and Sine Integrals') exit
	end do
	read(7,*) nval
	write(*,*) text
	write(*,'(1x,t5,a1,t13,a6,t25,a5,t37,a6,t48,a5)')&
		'X','Actual','CI(X)','Actual','SI(X)'
	do i=1,nval
		read(7,*) x,xci,xsi
		call cisi(x,ci,si)
		write(*,'(f6.2,4f12.6)') x,xci,ci,xsi,si
	end do
	close(7)
	END PROGRAM xcisi
