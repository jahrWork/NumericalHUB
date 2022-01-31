	PROGRAM xairy
!	driver for routine airy
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i,nval
	REAL(SP) :: ai,bi,aip,bip,x,xai,xbi,xaip,xbip
	CHARACTER(14) :: text
	open(7,file='FNCVAL.DAT',status='OLD')
	do
		read(7,'(a)') text
		if (text == 'Airy Functions') exit
	end do
	read(7,*) nval
	write(*,*) text
	write(*,'(1x,t3,a1)') 'X'
	write(*,'(1x,t5,a2,t21,a2,t37,a3,t53,a3)') 'AI','BI','AIP','BIP'
	write(*,'(1x,t5,a3,t21,a3,t37,a4,t53,a4)')&
		'XAI','XBI','XAIP','XBIP'
	do i=1,nval
		read(7,*) x,ai,bi,aip,bip
		call airy(x,xai,xbi,xaip,xbip)
		write(*,'(f6.2,/,1p,4e16.6,/,1p,4e16.6)')&
			x,ai,bi,aip,bip,xai,xbi,xaip,xbip
	end do
	close(7)
	END PROGRAM xairy
