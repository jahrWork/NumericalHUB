	PROGRAM xttest
!	driver for routine ttest
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPTS=1024,MPTS=512,NSHFT=10
	REAL(SP), PARAMETER :: EPS=0.02_sp
	INTEGER(I4B) :: i
	REAL(SP) :: prob,shift,t
	REAL(SP), DIMENSION(NPTS) :: data1
	REAL(SP), DIMENSION(MPTS) :: data2
!	generate Gaussian distributed data
	call ran_seed(sequence=5)
	call gasdev(data1)
	call gasdev(data2)
	data2(:)=(NSHFT/2.0_sp)*EPS+data2(:)
	write(*,'(/1x,t4,a,t18,a,t25,a)') 'Shift','T','Probability'
	do i=1,NSHFT+1
		call ttest(data1,data2,t,prob)
		shift=(i-1)*EPS
		write(*,'(1x,f6.2,2f12.2)') shift,t,prob
		data1(:)=data1(:)+EPS
	end do
	END PROGRAM xttest
