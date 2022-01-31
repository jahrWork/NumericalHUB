	PROGRAM xlocate
!	driver for routine locate
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=100
	INTEGER(I4B) :: i,j
	REAL(SP) :: x
	REAL(SP), DIMENSION(N) :: xx
!	create array to be searched
	xx(1:N)=exp(arth(1,1,N)/20.0_sp)-74.0_sp
	write(*,*) 'Result of:   j=0 indicates x too small'
	write(*,*) '           j=100 indicates x too large'
	write(*,'(t5,a7,t17,a1,t24,a5,t34,a7)') 'locate ','j','xx(j)','xx(j+1)'
!	perform test
	do i=1,19
		x=-100.0_sp+200.0_sp*i/20.0_sp
		j=locate(xx,x)
		if (j == 0) then
			write(*,'(1x,f10.4,i6,a12,f12.6)') x,j,'lower lim',xx(j+1)
		else if (j == N) then
			write(*,'(1x,f10.4,i6,f12.6,a12)') x,j,xx(j),'upper lim'
		else
			write(*,'(1x,f10.4,i6,2f12.6)') x,j,xx(j),xx(j+1)
		end if
	end do
	END PROGRAM xlocate
