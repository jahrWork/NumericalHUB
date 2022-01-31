	PROGRAM xkstwo
!	driver for routine kstwo
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N1=2000,N2=1000
	REAL(SP), PARAMETER :: EPS=0.1_sp
	INTEGER(I4B) :: i
	REAL(SP) :: d,factr,prob,var
	REAL(SP), DIMENSION(N1) :: data1,harvest
	REAL(SP), DIMENSION(N2) :: data2
	call ran_seed(sequence=2001)
	call gasdev(harvest(1:N1))
	data1(1:N1)=harvest(1:N1)
	write(*,'(/1x,t6,a,t26,a,t46,a/)') &
		'Variance Ratio','K-S Statistic','Probability'
	do i=1,11
		var=1.0_sp+(i-1)*EPS
		factr=sqrt(var)
		call gasdev(harvest(1:N2))
		data2(1:N2)=factr*harvest(1:N2)
		call kstwo(data1,data2,d,prob)
		write(*,'(1x,f15.6,f19.6,e20.4)') var,d,prob
	end do
	END PROGRAM xkstwo
