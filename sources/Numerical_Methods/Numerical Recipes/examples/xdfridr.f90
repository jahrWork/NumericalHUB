	PROGRAM xdfridr
!	driver for routine dfridr
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP) :: dx,err,h,x
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	do
		write(*,*) 'input x,h '
		read(*,*,END=999) x,h
		dx=dfridr(func,x,h,err)
		write(*,*) 'DFRIDR=',dx,1./cos(x)**2,err
	end do
999	END PROGRAM xdfridr

	FUNCTION func(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: func
	func=tan(x)
	END FUNCTION func
