	PROGRAM xvander
!	driver for routine vander
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=5
	INTEGER(I4B) :: i
	REAL(DP) :: sm
	REAL(DP), DIMENSION(N) :: term,w,&
		q = (/ 1.0_dp,1.5_dp,2.0_dp,2.5_dp,3.0_dp /),&
		x = (/ 1.0_dp,1.5_dp,2.0_dp,2.5_dp,3.0_dp /)
	w(1:N)=vander(x,q)
	write(*,*) 'Solution vector:'
	do i=1,N
		write(*,'(5x,a2,i1,a4,e12.6)') 'W(',i,') = ',w(i)
	end do
	write(*,'(/1x,a)') 'Test of solution vector:'
	write(*,'(1x,t6,a,t19,a)') 'mtrx*sol''n','original'
	sm=0.0
	term(:)=w(:)
	sm=sum(w(:))
	write(*,'(1x,2f12.4)') sm,q(1)
	do i=2,N
		term(:)=term(:)*x(:)
		sm=sum(term(:))
		write(*,'(1x,2f12.4)') sm,q(i)
	end do
	END PROGRAM xvander
