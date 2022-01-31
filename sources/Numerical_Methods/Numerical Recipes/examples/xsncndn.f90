	PROGRAM xsncndn
!	driver for routine sncndn
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i,nval
	REAL(SP) :: cn,dn,em,emmc,resul1,resul2,sn,uu,value
	CHARACTER(26) :: text
	open(7,file='FNCVAL.DAT',status='OLD')
	do
		read(7,'(a)') text
		if (text == 'Jacobian Elliptic Function') exit
	end do
	read(7,*) nval
	write(*,*) text
	write(*,'(1x,t4,a,t13,a,t21,a,t38,a,t49,a,t60,a)')&
		'Mc','U','Actual','SN','SN^2+CN^2','(Mc)*(SN^2)+DN^2'
	do i=1,nval
		read(7,*) em,uu,value
		emmc=1.0_sp-em
		call sncndn(uu,emmc,sn,cn,dn)
		resul1=sn*sn+cn*cn
		resul2=em*sn*sn+dn*dn
		write(*,'(1x,f5.2,f8.2,2e15.5,f12.5,f14.5)')&
			emmc,uu,value,sn,resul1,resul2
	end do
	close(7)
	END PROGRAM xsncndn
