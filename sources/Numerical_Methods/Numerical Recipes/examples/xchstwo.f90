	PROGRAM xchstwo
!	driver for routine chstwo
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NBINS=10,NPTS=2000
	INTEGER(I4B) :: i,ibin
	REAL(SP) :: chsq,df,prob,x
	REAL(SP), DIMENSION(NBINS) :: bins1,bins2
	call ran_seed(sequence=41)
	bins1(:)=0.0
	bins2(:)=0.0
	do i=1,NPTS
		call expdev(x)
		ibin=x*NBINS/3.0_sp+1
		if (ibin <= NBINS) bins1(ibin)=bins1(ibin)+1.0_sp
		call expdev(x)
		ibin=x*NBINS/3.0_sp+1
		if (ibin <= NBINS) bins2(ibin)=bins2(ibin)+1.0_sp
	end do
	call chstwo(bins1,bins2,0,df,chsq,prob)
	write(*,'(1x,t10,a,t25,a)') 'Dataset 1','Dataset 2'
	do i=1,NBINS
		write(*,'(1x,2f15.2)') bins1(i),bins2(i)
	end do
	write(*,'(/1x,t10,a,e12.4)') 'Chi-squared:',chsq
	write(*,'(1x,t10,a,e12.4)') 'Probability:',prob
	END PROGRAM xchstwo
