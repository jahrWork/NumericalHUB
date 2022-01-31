	PROGRAM xgaucof
!	driver for routine gaucof
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=64
	REAL(SP), PARAMETER :: SQRTPI=1.7724539_sp
	INTEGER(I4B) :: i,n
	REAL(SP) :: amu0,check
	REAL(SP), DIMENSION(NP) :: a,b,x,w
	do
		write(*,*) 'Enter N'
		read(*,*,END=99) n
		a(1:n)=0.0
		b(2:n)=arth(1,1,n-1)*0.5_sp
		b(1)=0.0
!	b(1) is arbitrary for call to TQLI
		amu0=SQRTPI
		call gaucof(a(1:n),b(1:n),amu0,x(1:n),w(1:n))
		write(*,'(/1x,t3,a,t10,a,t22,a/)') '#','X(I)','W(I)'
		do i=1,n
			write(*,'(1x,i2,2e14.6)') i,x(i),w(i)
		end do
		check=sum(w(1:n))
		write(*,'(/1x,a,e15.7,a,e15.7)') 'Check value:',check,&
			'  should be:',SQRTPI
	end do
99	END PROGRAM xgaucof
