	PROGRAM xexpdev
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
!	driver for routine expdev
	INTEGER(I4B), PARAMETER :: N=21,NPTS=10000
	REAL(SP), PARAMETER :: EE=2.718281828_sp
	INTEGER(I4B) :: i
	REAL(SP) :: expect,total,y
	REAL(SP), DIMENSION(N) :: trig,x
	call ran_seed(sequence=24)
	trig(1:N)=arth(0,1,N)/20.0_sp
	x(:)=0.0
	do i=1,NPTS
		call expdev(y)
		where ((y < trig(2:N)) .and. (y > trig(1:N-1))) &
			x(2:N)=x(2:N)+1.0_sp
	end do
	total=sum(x(2:N))
	write(*,'(1x,a,i6,a)') 'Exponential distribution with',NPTS,' points:'
	write(*,'(1x,t5,a,t19,a,t31,a)') 'interval','observed','expected'
	x(2:N)=x(2:N)/total
	do i=2,N
		expect=exp(-(trig(i-1)+trig(i))/2.0_sp)
		expect=expect*0.05_sp*EE/(EE-1)
		write(*,'(1x,2f6.2,2f12.4)') trig(i-1),trig(i),x(i),expect
	end do
	END PROGRAM xexpdev
