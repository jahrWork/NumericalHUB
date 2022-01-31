	PROGRAM xksone
!	driver for routine ksone
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPTS=1000
	REAL(SP), PARAMETER :: EPS=0.1_sp
	INTEGER(I4B) :: i
	REAL(SP) :: d,factr,prob,var
	REAL(SP), DIMENSION(NPTS) :: data,harvest
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: func
		END FUNCTION func
	END INTERFACE
	call ran_seed(sequence=1066)
	write(*,'(/1x,t5,a,t24,a,t44,a/)') &
		'Variance Ratio','K-S Statistic','Probability'
	do i=1,11
		var=1.0_sp+(i-1)*EPS
		factr=sqrt(var)
		call gasdev(harvest)
		data(:)=factr*abs(harvest(:))
		call ksone(data,func,d,prob)
		write(*,'(1x,f14.6,f18.6,e20.4)') var,d,prob
	end do
	END PROGRAM xksone

	FUNCTION func(x)
	USE nrtype
	USE nr
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: func
	func=erf(x/sqrt(2.0_sp))
	END FUNCTION func
