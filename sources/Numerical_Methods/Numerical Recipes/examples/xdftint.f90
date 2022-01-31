MODULE xdftint_parms
	USE nrtype
	IMPLICIT NONE
	REAL(SP) :: c,d
END MODULE xdftint_parms

	PROGRAM xdftint
!	driver for routine dftint
	USE nrtype
	USE nr
	USE xdftint_parms
	IMPLICIT NONE
	REAL(SP) :: a,b,cans,cosint,sans,sinint,w
	INTERFACE
		FUNCTION coscxd(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: coscxd
		END FUNCTION coscxd
	END INTERFACE
	write(*,'(2x,a,t10,a,t37,a,t44,a,t70,a)') 'Omega',&
		'Integral cosine*test func','Err','Integral sine*test func','Err'
	w=999.0
1	do
		if (w >= 0.0) then
			write(*,*) 'INPUT C,D: '
			read(*,*) c,d
		end if
		write(*,*) 'INPUT A,B: '
		read(*,*) a,b
		if (a /= b) exit
	end do
	do
		write(*,*) 'INPUT W: '
		read(*,*,END=999) w
		if (w < 0.0) goto 1
		call dftint(coscxd,a,b,w,cosint,sinint)
		call getans(w,a,b,cans,sans)
		write(*,'(1x,1p,5e15.6)') w,cans,cosint-cans,sans,sinint-sans
	end do
999	write(*,*) 'Normal completion'
	CONTAINS
!BL
	SUBROUTINE getans(w,a,b,cans,sans)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: w,a,b
	REAL(SP), INTENT(OUT) :: cans,sans
	cans=ci(b,w)-ci(a,w)
	sans=si(b,w)-si(a,w)
	END SUBROUTINE getans
!BL
	FUNCTION ci(x,w)
	USE nrtype
	USE xdftint_parms
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,w
	REAL(SP) :: ci
	ci=sin((w-c)*x-d)/(2.0_sp*(w-c))+sin((w+c)*x+d)/(2.0_sp*(w+c))
	END FUNCTION ci
!BL
	FUNCTION si(x,w)
	USE nrtype
	USE xdftint_parms
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x,w
	REAL(SP) :: si
	si=-cos((w-c)*x-d)/(2.0_sp*(w-c))-cos((w+c)*x+d)/(2.0_sp*(w+c))
	END FUNCTION si
	END PROGRAM xdftint

	FUNCTION coscxd(x)
	USE nrtype
	USE xdftint_parms
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: coscxd
	coscxd=cos(c*x+d)
	END FUNCTION coscxd
