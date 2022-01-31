	PROGRAM xstiff
!	driver for routine stiff
	USE nrtype
	USE nr
	USE ode_path
	IMPLICIT NONE
	INTERFACE
		SUBROUTINE derivs(x,y,dydx)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(IN) :: y
		REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
		END SUBROUTINE derivs
	END INTERFACE
	REAL(SP) :: eps,hstart,x1,x2
	REAL(SP), DIMENSION(3) :: y
	do
		write(*,*) 'ENTER EPS,HSTART:'
		read(*,*,END=999) eps,hstart
		save_steps=.true.
		x1=0.0
		x2=50.0
		y(1)=1.0
		y(2)=1.0
		y(3)=0.0
		call odeint(y,x1,x2,eps,hstart,0.0_sp,derivs,stiff)
		write(*,'(/1x,a,t30,i3)') 'Successful steps:',nok
		write(*,'(1x,a,t30,i3)') 'Bad steps:',nbad
		write(*,*) 'Y(END)=',y(1),y(2),y(3)
	end do
999	write(*,*) 'NORMAL COMPLETION'
	END PROGRAM xstiff
