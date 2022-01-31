	PROGRAM xprobks
!	driver for routine probks
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NT=50
	INTEGER(I4B) :: i,j,npts
	REAL(SP) :: alam,eps,scale,value
	CHARACTER(1), DIMENSION(NT) :: text
	write(*,*) 'Probability func. for Kolmogorov-Smirnov statistic'
	write(*,'(/1x,t3,a,t15,a,t27,a)') 'Lambda:','Value:','Graph:'
	npts=20
	eps=0.1_sp
	scale=40.0
	do i=1,npts
		alam=i*eps
		value=probks(alam)
		text(1)='*'
		text(1:NT)=merge('*',' ', arth(1,1,NT) <= nint(scale*value) )
		write(*,'(1x,f9.6,f12.6,4x,50a1)') alam,value,&
			(text(j),j=1,NT)
	end do
	END PROGRAM xprobks
