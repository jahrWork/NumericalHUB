	PROGRAM xftest
!	driver for routine ftest
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPTS=1000,MPTS=500,NVAL=10
	REAL(SP), PARAMETER :: EPS=0.01_sp
	INTEGER(I4B) :: i
	REAL(SP) :: f,factor,prob,var
	REAL(SP), DIMENSION(MPTS) :: data2,data3
	REAL(SP), DIMENSION(NPTS) :: data1,harvest
!	generate two Gaussian distributions with
!	different variances
	call ran_seed(sequence=1492)
	call gasdev(harvest(1:NPTS))
	data1(1:NPTS)=harvest(1:NPTS)
	call gasdev(harvest(1:MPTS))
	data2(1:MPTS)=harvest(1:MPTS)
	write(*,'(1x,t5,a,f5.2)') 'Variance 1 = ',1.0_sp
	write(*,'(1x,t5,a,t21,a,t30,a)') 'Variance 2','Ratio','Probability'
	do i=1,NVAL+1
		var=1.0_sp+(i-1)*EPS
		factor=sqrt(var)
		data3(1:MPTS)=factor*data2(1:MPTS)
		call ftest(data1,data3,f,prob)
		write(*,'(1x,f11.4,2x,2f12.4)') var,f,prob
	end do
	END PROGRAM xftest
