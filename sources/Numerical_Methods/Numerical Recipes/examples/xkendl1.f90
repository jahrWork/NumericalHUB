	PROGRAM xkendl1
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
!	driver for routine kendl1
	INTEGER(I4B), PARAMETER :: NDAT=200
	INTEGER(I4B) :: i,idum,j
	REAL(SP) :: prob,tau,z
	REAL(SP), DIMENSION(NDAT) :: data1,data2
	CHARACTER*4, DIMENSION(5) :: text = (/ &
		'RAN0','RAN1','RAN2','RAN3','RAN '/)
	write(*,'(/1x,a/)') 'Pair correlations of random number routines'
	write(*,'(2x,a,t16,a,t34,a,t50,a,/)') &
		'Program','Kendall Tau','Std. Dev.','Probability'
	do i=1,5
		idum=-1357
		call ran_seed(sequence=-idum)
		select case(i)
		case (1)
			call ran0(data1)
			call ran0(data2)
		case (2)
			call ran1(data1)
			call ran1(data2)
		case (3)
			call ran2(data1)
			call ran2(data2)
		case (4)
			call ran3(data1)
			call ran3(data2)
		case (5)
			do j=1,NDAT
				data1(j)=ran(idum)
				data2(j)=ran(idum)
			end do
		end select
		call kendl1(data1,data2,tau,z,prob)
		write(*,'(1x,t4,a,3f17.6)') text(i),tau,z,prob
	end do
	END PROGRAM xkendl1
