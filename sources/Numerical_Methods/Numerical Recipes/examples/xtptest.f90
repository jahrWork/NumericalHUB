	PROGRAM xtptest
!	driver for routine tptest
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPTS=500,NSHFT=10
	REAL(SP), PARAMETER :: EPS=0.01_sp,ANOISE=0.3_sp
	INTEGER(I4B) :: i
	REAL(SP) :: ave1,ave2,ave3,offset,prob1,prob2,shift,t1,t2,var1,var2,var3
	REAL(SP), DIMENSION(NPTS) :: data1,data2,data3,datatmp
	write(*,'(1x,t18,a,t46,a)') 'Correlated:','Uncorrelated:'
	write(*,'(1x,t4,a,t18,a,t25,a,t46,a,t53,a)') &
		'Shift','T','Probability','T','Probability'
	offset=(NSHFT/2)*EPS
	call ran_seed(sequence=8)
	call gasdev(data1)
	call gasdev(data2)
	call gasdev(data3)
	call gasdev(datatmp)
	data2(:)=data1(:)+ANOISE*data2(:)
	data3(:)=data3(:)+ANOISE*datatmp(:)
	call avevar(data1,ave1,var1)
	call avevar(data2,ave2,var2)
	call avevar(data3,ave3,var3)
	data1(:)=data1(:)-ave1+offset
	data2(:)=data2(:)-ave2
	data3(:)=data3(:)-ave3
	do i=1,NSHFT+1
		shift=i*EPS
		data2(:)=data2(:)+EPS
		data3(:)=data3(:)+EPS
		call tptest(data1,data2,t1,prob1)
		call tptest(data1,data3,t2,prob2)
		write(*,'(1x,f6.2,2x,2f12.4,4x,2f12.4)') &
			shift,t1,prob1,t2,prob2
	end do
	END PROGRAM xtptest
