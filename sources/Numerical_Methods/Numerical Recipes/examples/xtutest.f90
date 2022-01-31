	PROGRAM xtutest
!	driver for routine tutest
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPTS=5000,MPTS=1000,NSHFT=10
	REAL(SP), PARAMETER :: EPS=0.02_sp,VAR1=1.0_sp,VAR2=4.0
	INTEGER(I4B) :: i
	REAL(SP) :: fctr1,fctr2,prob,shift,t
	REAL(SP), DIMENSION(NPTS) :: data1
	REAL(SP), DIMENSION(MPTS) :: data2
!	generate two Gaussian distributions of different variance
	call ran_seed(sequence=51773)
	fctr1=sqrt(VAR1)
	call gasdev(data1)
	data1(:)=fctr1*data1(:)
	fctr2=sqrt(VAR2)
	call gasdev(data2)
	data2(:)=(NSHFT/2.0_sp)*EPS+fctr2*data2(:)
	write(*,'(1x,a,f6.2)') 'Distribution #1 : variance = ',VAR1
	write(*,'(1x,a,f6.2/)') 'Distribution #2 : variance = ',VAR2
	write(*,'(1x,t4,a,t18,a,t25,a)') 'Shift','T','Probability'
	do i=1,NSHFT+1
		call tutest(data1,data2,t,prob)
		shift=(i-1)*EPS
		write(*,'(1x,f6.2,2f12.2)') shift,t,prob
		data1(:)=data1(:)+EPS
	end do
	END PROGRAM xtutest
