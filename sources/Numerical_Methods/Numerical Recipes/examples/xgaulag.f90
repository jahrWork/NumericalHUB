	PROGRAM xgaulag
!	driver for routine gaulag
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=64
	INTEGER(I4B) :: i,n
	REAL(SP) :: alf,checkw,checkx,xx
	REAL(SP), DIMENSION(NP) :: x,w
	alf=1.0
	do
		write(*,*) 'Enter N'
		read(*,*,END=99) n
		call gaulag(x(1:n),w(1:n),alf)
		write(*,'(/1x,t3,a,t10,a,t22,a/)') '#','X(I)','W(I)'
		do i=1,n
			write(*,'(1x,i2,2e14.6)') i,x(i),w(i)
		end do
		checkx=sum(x(1:n))
		checkw=sum(w(1:n))
		write(*,'(/1x,a,e15.7,a,e15.7)') 'Check value:',checkx,&
			'  should be:',n*(n+alf)
		write(*,'(/1x,a,e15.7,a,e15.7)') 'Check value:',checkw,&
			'  should be:',exp(gammln(1.0_sp+alf))
!	demonstrate the use of GAULAG for an integral
		xx=0.0
		do i=1,n
			xx=xx+w(i)*func(x(i))
		end do
		write(*,'(/1x,a,f12.6)') 'Integral from GAULAG:',xx
		write(*,'(1x,a,f12.6)') 'Actual value:        ',&
			1.0_sp/(2.0_sp*sqrt(2.0_sp))
	end do
	CONTAINS
!BL
	FUNCTION func(x)
	IMPLICIT NONE
	REAL(SP) :: func,x
	func=bessj0(x)
	END FUNCTION func
99	END PROGRAM xgaulag
