	PROGRAM xzroots
!	driver for routine zroots
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: M=4,M1=M+1
	INTEGER(I4B) :: i
	COMPLEX(SP), DIMENSION(M) :: roots
	COMPLEX(SP), DIMENSION(M1) :: a = (/ (0.0_sp,2.0_sp),(0.0_sp,0.0_sp),&
		(-1.0_sp,-2.0_sp),(0.0_sp,0.0_sp),(1.0_sp,0.0_sp) /)
	LOGICAL(LGT) :: polish
	write(*,'(/1x,a)') 'Roots of the polynomial x^4-(1+2i)*x^2+2i'
	polish=.false.
	call zroots(a,roots,polish)
	write(*,'(/1x,a)') 'Unpolished roots:'
	write(*,'(1x,t10,a,t25,a,t37,a)') 'Root #','Real','Imag.'
	do i=1,M
		write(*,'(1x,i11,5x,2f12.6)') i,roots(i)
	end do
	write(*,'(/1x,a)') 'Corrupted roots:'
	do i=1,M
		roots(i)=roots(i)*(1.0_sp+0.01_sp*i)
	end do
	write(*,'(1x,t10,a,t25,a,t37,a)') 'Root #','Real','Imag.'
	do i=1,M
		write(*,'(1x,i11,5x,2f12.6)') i,roots(i)
	end do
	polish=.true.
	call zroots(a,roots,polish)
	write(*,'(/1x,a)') 'Polished roots:'
	write(*,'(1x,t10,a,t25,a,t37,a)') 'Root #','Real','Imag.'
	do i=1,M
		write(*,'(1x,i11,5x,2f12.6)') i,roots(i)
	end do
	END PROGRAM xzroots
