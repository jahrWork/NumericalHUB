	PROGRAM xscrsho
!	driver for routine scrsho
	USE nrtype
	USE nr
	IMPLICIT NONE
	write(*,*) 'Graph of the Bessel Function J0:'
	call scrsho(bessj0_s)
	END PROGRAM xscrsho
