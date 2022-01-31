	PROGRAM xcntab1
!	driver for routine cntab1
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDAT=9,NMON=12
	INTEGER(I4B) :: i,j
	INTEGER(I4B), DIMENSION(NDAT,NMON) :: nmbr
	REAL(SP) :: ccc,chisq,cramrv,df,prob
	CHARACTER(64) :: text
	CHARACTER(15), DIMENSION(NDAT) :: fate
	CHARACTER(5), DIMENSION(NMON) :: mon
	open(7,file='TABLE.DAT',status='OLD')
	read(7,*)
	read(7,'(a)') text
	read(7,'(15x,12a5/)') (mon(i),i=1,12)
	do i=1,NDAT
		read(7,'(a15,12i5)') fate(i),(nmbr(i,j),j=1,12)
	end do
	close(7)
	write(*,'(/1x,a/)') text
	write(*,'(1x,15x,12a5)') (mon(i),i=1,12)
	do i=1,NDAT
		write(*,'(1x,a,12i5)') fate(i),(nmbr(i,j),j=1,12)
	end do
	call cntab1(nmbr,chisq,df,prob,cramrv,ccc)
	write(*,'(/1x,a,t20,f20.2)') 'Chi-squared',chisq
	write(*,'(1x,a,t20,f20.2)') 'Degrees of Freedom',df
	write(*,'(1x,a,t20,f20.4)') 'Probability',prob
	write(*,'(1x,a,t20,f20.4)') 'Cramer-V',cramrv
	write(*,'(1x,a,t20,f20.4)') 'Contingency Coeff.',ccc
	END PROGRAM xcntab1
