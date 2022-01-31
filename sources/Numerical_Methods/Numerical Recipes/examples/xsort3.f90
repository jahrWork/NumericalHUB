	PROGRAM xsort3
!	driver for routine sort3
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NLEN=64
	INTEGER(I4B) :: i
	INTEGER(I4B), DIMENSION(NLEN) :: itmp
	REAL(SP), DIMENSION(NLEN) :: a,b,c
	CHARACTER(1), DIMENSION(64) :: amsg,bmsg,cmsg
	CHARACTER(33) :: msg1 = 'I''d rather have a bottle in front'
	CHARACTER(31) :: msg2 = ' of me than a frontal lobotomy.'
	do i=1,33
		amsg(i)=msg1(i:i)
	end do
	do i=1,31
		amsg(33+i)=msg2(i:i)
	end do
	write(*,*) 'Original message:'
	write(*,'(1x,64a1,/)') (amsg(i),i=1,64)
!	read array of random numbers
	open(7,file='TARRAY.DAT',status='OLD')
	read(7,*) (a(i),i=1,NLEN)
	close(7)
!	create array B and array C
	itmp(1:NLEN)=arth(1,1,NLEN)
	b(1:NLEN)=itmp(1:NLEN)
	c(1:NLEN)=NLEN+1-itmp(1:NLEN)
!	sort array A while mixing B and C
	call sort3(a,b,c)
!	scramble message according to array B
	bmsg(1:NLEN)=amsg(int(b(1:NLEN)))
	write(*,*) 'Scrambled message:'
	write(*,'(1x,64a1,/)') (bmsg(i),i=1,64)
!	unscramble according to array C
	cmsg(int(c(1:NLEN)))=bmsg(1:NLEN)
	write(*,*) 'Mirrored message:'
	write(*,'(1x,64a1,/)') (cmsg(i),i=1,64)
	END PROGRAM xsort3
