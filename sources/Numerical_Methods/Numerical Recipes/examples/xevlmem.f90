	PROGRAM xevlmem
!	driver for routine evlmem
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=1000,M=10,NFDT=16
	INTEGER(I4B) :: i
	REAL(SP) :: fdt,pm
	REAL(SP), DIMENSION(M) :: cof
	REAL(SP), DIMENSION(N) :: data
	open(7,file='SPCTRL.DAT',status='old')
	read(7,*) (data(i),i=1,N)
	close(7)
	call memcof(data,pm,cof)
	write(*,*) 'Power spectrum estimate of DATA in SPCTRL.DAT'
	write(*,'(1x,t6,a,t20,a)') 'f*delta','power'
	do i=0,NFDT
		fdt=0.5_sp*i/NFDT
		write(*,'(1x,2f12.6)') fdt,evlmem(fdt,cof,pm)
	end do
	END PROGRAM xevlmem
