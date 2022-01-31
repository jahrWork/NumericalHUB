	PROGRAM xorthog
!	driver for routine orthog
	USE nrtype; USE nrutil
	USE nr, ONLY : orthog,gaucof
	INTEGER(I4B), PARAMETER :: NP=64
	INTEGER(I4B) :: i,n
	INTEGER(I4B), DIMENSION(2*NP-1) :: itmp
	REAL(SP) :: amu0,check,xx
	REAL(SP), DIMENSION(2*NP) :: anu
	REAL(SP), DIMENSION(2*NP-1) :: alpha,beta
	REAL(SP), DIMENSION(NP) :: a,b,x,w
	itmp(1:2*NP-1)=arth(1,1,2*NP-1)
	do
		write(*,*) 'Enter N'
		read(*,*,END=999) n
		alpha(1:2*n-1)=0.5_sp
		beta(1)=1.0
		beta(2:2*n-1)=1.0_sp/(4.0_sp*(4.0_sp-1.0_sp/(itmp(1:2*n-2)**2)))
		anu(1)=1.0
		anu(2)=-0.25
		do i=2,2*n-1
			anu(i+1)=-anu(i)*i*(i-1)/(2.0_sp*(i+1)*(2*i-1))
		end do
		call orthog(anu(1:2*n),alpha(1:2*n-1),beta(1:2*n-1),a(1:n),b(1:n))
		amu0=1.0
		call gaucof(a(1:n),b(1:n),amu0,x(1:n),w(1:n))
		write(*,'(/1x,t3,a,t10,a,t22,a/)') '#','X(I)','W(I)'
		do i=1,n
			write(*,'(1x,i2,2e14.6)') i,x(i),w(i)
		end do
		check=sum(w(1:n))
		write(*,'(/1x,a,e15.7,a,e15.7)') 'Check value:',check,'  should be:',amu0
!	demonstrate the use of ORTHOG for an integral
		xx=0.0
		do i=1,n
			xx=xx+w(i)*func(x(i))
		end do
		write(*,'(/1x,a,f12.6)') 'Integral from ORTHOG:',xx
		write(*,'(1x,a,f12.6)') 'Actual value:        ',log(2.0_sp)
	end do
	CONTAINS
!BL
	FUNCTION func(x)
	REAL(SP) :: x,func
	func=1.0_sp/(1.0_sp+x)**2
	END FUNCTION func
999	END PROGRAM xorthog
