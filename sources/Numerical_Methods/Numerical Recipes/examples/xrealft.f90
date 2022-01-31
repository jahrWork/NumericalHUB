	PROGRAM xrealft
!	driver for routine realft
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=32
	REAL(SP), PARAMETER :: EPS=1.0e-3_sp,WIDTH=50.0
	INTEGER(I4B) :: i,j,n,nlim
	REAL(SP) :: big,per,scal,small
	REAL(SP), DIMENSION(NP) :: data,size
	n=NP/2
	do
		write(*,'(1x,a,i2,a)') 'Period of sinusoid in channels (2-',&
			NP,', OR 0 TO STOP)'
		read(*,*) per
		if (per <= 0.0) exit
		data(:)=cos(2.0_sp*PI*(arth(0,1,NP))/per)
		call realft(data,+1)
		size(1)=data(1)
		size(2:n)=sqrt(data(3:2*n-1:2)**2+data(4:2*n:2)**2)
		big=maxval(size(1:n))
		scal=WIDTH/big
		do i=1,n
			nlim=scal*size(i)+EPS
			write(*,'(1x,i4,1x,60a1)') i,('*',j=1,nlim+1)
		end do
		write(*,*) 'press continue ...'
		read(*,*)
		call realft(data,-1)
		small=minval(data(:))
		big=maxval(data(:))
		scal=WIDTH/(big-small)
		do i=1,NP
			nlim=scal*(data(i)-small)+EPS
			write(*,'(1x,i4,1x,60a1)') i,('*',j=1,nlim+1)
		end do
	end do
	END PROGRAM xrealft
