	PROGRAM xavevar
!	driver for routine avevar
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPTS=1000
	REAL(SP), PARAMETER :: EPS=0.1_sp
	INTEGER(I4B) :: i
	REAL(SP) :: ave,shift,var
	REAL(SP), DIMENSION(NPTS) :: data
!	generate Gaussian distributed data
	call ran_seed(sequence=1)
	write(*,'(1x,t4,a,t14,a,t26,a)') 'Shift','Average','Variance'
	do i=1,11
		shift=(i-1)*EPS
		call gasdev(data)
		data(:)=shift+i*data(:)
		call avevar(data,ave,var)
		write(*,'(1x,f6.2,2f12.2)') shift,ave,var
	end do
	END PROGRAM xavevar
