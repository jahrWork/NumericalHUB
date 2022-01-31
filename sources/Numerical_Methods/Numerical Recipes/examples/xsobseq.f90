	PROGRAM xsobseq
	USE nrtype
	USE nr
	IMPLICIT NONE
!	driver for routine sobseq
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(3) :: x
	call sobseq(x,1)
	do i=1,10
		call sobseq(x)
		write(*,'(3(1x,f10.5),1x,i5)') (x(j),j=1,3),i
	end do
	write (*,*)
!	test re-initializing
	call sobseq(x,1)
	do i=1,10
		call sobseq(x)
		write(*,'(3(1x,f10.5),1x,i5)') (x(j),j=1,3),i
	end do
	END PROGRAM xsobseq
