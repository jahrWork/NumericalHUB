	PROGRAM xgaujac
!	driver for routine gaujac
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=64
	INTEGER(I4B) :: i,n
	REAL(SP) :: ak,alf,bet,checkw,checkx,xx
	REAL(SP), DIMENSION(NP) :: x,w
	alf=-0.5_sp
	bet=-0.5_sp
	do
		write(*,*) 'Enter N'
		read(*,*,END=99) n
		call gaujac(x(1:n),w(1:n),alf,bet)
		write(*,'(/1x,t3,a,t10,a,t22,a/)') '#','X(I)','W(I)'
		do i=1,n
			write(*,'(1x,i2,2e14.6)') i,x(i),w(i)
		end do
		checkx=sum(x(1:n))
		checkw=sum(w(1:n))
		write(*,'(/1x,a,e15.7,a,e15.7)') 'Check value:',checkx,&
			'  should be:',n*(bet-alf)/(alf+bet+2*n)
		write(*,'(/1x,a,e15.7,a,e15.7)') 'Check value:',checkw,&
			'  should be:',exp(gammln(1.0_sp+alf)+gammln(1.0_sp+bet)-&
			gammln(2.0_sp+alf+bet))*2**(alf+bet+1.0_sp)
!	demonstrate the use of GAUJAC for an integral
		ak=0.5_sp
		xx=0.0
		do i=1,n
			xx=xx+w(i)*func(ak,x(i))
		end do
		write(*,'(/1x,a,f12.6)') 'Integral from GAUJAC:',xx
		write(*,'(1x,a,f12.6)') 'Actual value:        ',2.0_sp*ellf(PIO2,ak)
	end do
	CONTAINS
!BL
	FUNCTION func(ak,x)
	IMPLICIT NONE
	REAL(SP) :: func,ak,x
	func=1.0_sp/sqrt(1.0_sp-ak**2*(1.0_sp+x)/2.0_sp)
	END FUNCTION func
99	END PROGRAM xgaujac
