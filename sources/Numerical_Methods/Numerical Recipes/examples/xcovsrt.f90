	PROGRAM xcovsrt
!	driver for routine covsrt
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: MA=10,MFIT=5
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(MA,MA) :: covar
	LOGICAL(LGT), DIMENSION(MA) :: maska
	covar(:,:)=0.0
	covar(1:MFIT,1:MFIT)=outersum(real(arth(1,1,MFIT),sp),real(arth(0,1,MFIT),sp))
	write(*,'(//2x,a)') 'Original matrix'
	do i=1,MA
		write(*,'(1x,10f4.1)') (covar(i,j),j=1,MA)
	end do
	write(*,*) ' press RETURN to continue...'
	read(*,*)
	write(*,'(/2x,a)') 'Test #1 - Full Fitting'
	maska(:)=.true.
	call covsrt(covar,maska)
	do i=1,MA
		write(*,'(1x,10f4.1)') (covar(i,j),j=1,MA)
	end do
	write(*,*) ' press RETURN to continue...'
	read(*,*)
	write(*,'(/2x,a)') 'Test #2 - Spread'
	covar(:,:)=0.0
	covar(1:MFIT,1:MFIT)=outersum(real(arth(1,1,MFIT),sp),real(arth(0,1,MFIT),sp))
	maska(1:MA:2)=.false.
	call covsrt(covar,maska)
	do i=1,MA
		write(*,'(1x,10f4.1)') (covar(i,j),j=1,MA)
	end do
	END PROGRAM xcovsrt
