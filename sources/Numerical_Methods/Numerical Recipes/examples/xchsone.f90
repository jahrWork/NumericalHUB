	PROGRAM xchsone
!	driver for routine chsone
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NBINS=10,NPTS=2000
	INTEGER(I4B) :: i,ibin
	REAL(SP) :: chsq,df,prob,x
	REAL(SP), DIMENSION(NBINS) :: bins,ebins
	call ran_seed(sequence=56)
	bins(:)=0.0
	do i=1,NPTS
		call expdev(x)
		ibin=x*NBINS/3.0_sp+1
		if (ibin <= NBINS) bins(ibin)=bins(ibin)+1.0_sp
	end do
	ebins(1:NBINS)=3.0_sp*NPTS/NBINS*exp(-3.0_sp*(arth(1,1,NBINS)-0.5_sp)/NBINS)
	call chsone(bins,ebins,0,df,chsq,prob)
	write(*,'(1x,t10,a,t25,a)') 'Expected','Observed'
	do i=1,NBINS
		write(*,'(1x,2f15.2)') ebins(i),bins(i)
	end do
	write(*,'(/1x,t9,a,e12.4)') 'Chi-squared:',chsq
	write(*,'(1x,t9,a,e12.4)') 'Probability:',prob
	END PROGRAM xchsone
