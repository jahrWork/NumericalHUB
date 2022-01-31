	PROGRAM xmoment
!	driver for routine moment
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPTS=10000,NBIN=100,NDAT=NPTS+NBIN
	INTEGER(I4B) :: i,j,nlim
	REAL(SP) :: adev,ave,curt,sdev,skew,var,x
	REAL(SP), DIMENSION(NDAT) :: data
	i=1
	do j=1,NBIN
		x=PI*j/NBIN
		nlim=nint(sin(x)*PIO2*NPTS/NBIN)
		data(i:i+nlim-1)=x
		i=i+nlim
	end do
	write(*,'(1x,a/)') 'Moments of a sinusoidal distribution'
	call moment(data(1:i-1),ave,adev,sdev,var,skew,curt)
	write(*,'(1x,t29,a,t42,a/)') 'Calculated','Expected'
	write(*,'(1x,a,t25,2f12.4)') 'Mean :',ave,PIO2
	write(*,'(1x,a,t25,2f12.4)') 'Average Deviation :',adev,0.570796_sp
	write(*,'(1x,a,t25,2f12.4)') 'Standard Deviation :',sdev,0.683667_sp
	write(*,'(1x,a,t25,2f12.4)') 'Variance :',var,0.467401_sp
	write(*,'(1x,a,t25,2f12.4)') 'Skewness :',skew,0.0_sp
	write(*,'(1x,a,t25,2f12.4)') 'Kurtosis :',curt,-0.806249_sp
	END PROGRAM xmoment
