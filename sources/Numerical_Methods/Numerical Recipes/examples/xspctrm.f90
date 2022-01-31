	PROGRAM xspctrm
!	driver for routine spctrm
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: M=16,M4=4*M
	INTEGER(I4B) :: j,k
	REAL(SP), DIMENSION(M) :: p,q
	LOGICAL(LGT) :: ovrlap
	open(9,file='SPCTRL.DAT',status='OLD')
	k=8
	ovrlap=.true.
	call spctrm(p,k,ovrlap)
	rewind(9)
	k=16
	ovrlap=.false.
	call spctrm(q,k,ovrlap)
	close(9)
	write(*,*) 'Spectrum of DATA in file SPCTRL.DAT'
	write(*,'(1x,t14,a,t29,a)') 'Overlapped','Non-Overlapped'
	do j=1,M
		write(*,'(1x,i4,2f17.6)') j,p(j),q(j)
	end do
	END PROGRAM xspctrm
