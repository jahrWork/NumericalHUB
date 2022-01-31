	PROGRAM xcosft1
!	driver for routine cosft1
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=17
	REAL(SP), PARAMETER :: EPS=1.0e-3_sp,WIDTH=30.0
	INTEGER(I4B) :: i,j,nlim
	REAL(SP) :: big,per,scal,small
	REAL(SP), DIMENSION(NP) :: data
	do
		write(*,'(1x,a,i2,a)') 'Period of cosine in channels (2-',NP,')'
		read(*,*) per
		if (per <= 0.0) exit
		data(1:NP)=cos(2.0_sp*PI*arth(0,1,NP)/per)
		call cosft1(data)
		big=maxval(data)
		small=minval(data)
		scal=WIDTH/(big-small)
		do i=1,NP
			nlim=scal*(data(i)-small)+EPS
			write(*,'(1x,i2,f6.2,1x,60a1)') i,data(i),('*',j=1,nlim+1)
		end do
		write(*,*) 'press continue ...'
		read(*,*)
		call cosft1(data)
		big=maxval(data)
		small=minval(data)
		scal=WIDTH/(big-small)
		do i=1,NP
			nlim=scal*(data(i)-small)+EPS
			write(*,'(1x,i4,1x,60a1)') i,('*',j=1,nlim+1)
		end do
	end do
	END PROGRAM xcosft1
