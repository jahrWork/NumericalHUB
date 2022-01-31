	PROGRAM xsvdvar
!	driver for routine svdvar
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: MP=6,MA=3,NCVM=MA
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(MP) :: w = (/ 0.0_sp,1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp /)
	REAL(SP), DIMENSION(NCVM,NCVM) :: cvm
	REAL(SP), DIMENSION(MP,MP) :: v = reshape( (/ &
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,6.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,6.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,6.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,6.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,6.0_sp,&
		1.0_sp,2.0_sp,3.0_sp,4.0_sp,5.0_sp,6.0_sp /),&
		(/ MP,MP /) )
	REAL(SP), DIMENSION(MA,MA) :: tru = reshape( (/ &
		1.25_sp,2.5_sp,3.75_sp,2.5_sp,5.0_sp,7.5_sp,&
		3.75_sp,7.5_sp,11.25_sp /), (/ MA,MA /) )
	write(*,'(/1x,a)') 'Matrix V'
	do i=1,MP
		write(*,'(1x,6f12.6)') (v(i,j),j=1,MP)
	end do
	write(*,'(/1x,a)') 'Vector W'
	write(*,'(1x,6f12.6)') (w(i),i=1,MP)
	call svdvar(v(1:MA,1:MA),w(1:MA),cvm(1:NCVM,1:NCVM))
	write(*,'(/1x,a)') 'Covariance matrix from SVDVAR'
	do i=1,MA
		write(*,'(1x,3f12.6)') (cvm(i,j),j=1,MA)
	end do
	write(*,'(/1x,a)') 'Expected covariance matrix'
	do i=1,MA
		write(*,'(1x,3f12.6)') (tru(i,j),j=1,MA)
	end do
	END PROGRAM xsvdvar
