	PROGRAM xmemcof
!	driver for routine memcof
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=1000,M=10
	INTEGER(I4B) :: i
	REAL(SP) :: pm
	REAL(SP), DIMENSION(M) :: cof
	REAL(SP), DIMENSION(N) :: data
	open(7,file='SPCTRL.DAT',status='OLD')
	read(7,*) (data(i),i=1,N)
	close(7)
	call memcof(data,pm,cof)
	write(*,'(/1x,a/)') 'Coeff. for spectral estim. of SPCTRL.DAT'
	do i=1,M
		write(*,'(1x,a,i2,a,f12.6)') 'a[',i,'] =',cof(i)
	end do
	write(*,'(/1x,a,f12.6/)') 'a0 =',pm
	END PROGRAM xmemcof
