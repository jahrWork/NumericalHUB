	PROGRAM xhunt
!	driver for routine hunt
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=100
	INTEGER(I4B) :: i,j,ji
	REAL(SP) :: x
	REAL(SP), DIMENSION(N) :: xx
!	create array to be searched
	xx(1:N)=exp(arth(1,1,N)/20.0_sp)-74.0_sp
	write(*,*) 'Result of:    j=0 indicates x too small'
	write(*,*) '            j=100 indicates x too large'
	write(*,'(t7,a7,t17,a5,t25,a1,t32,a5,t42,a7)') 'locate:',&
		'guess','j','xx(j)','xx(j+1)'
!	perform test
	do i=1,19
		x=-100.0_sp+200.0_sp*i/20.0_sp
!	trial parameter
		ji=5*i
		j=ji
!	begin search
		call hunt(xx,x,j)
		if (j == 0) then
			write(*,'(1x,f12.6,2i6,a12,f12.6)') x,ji,j,&
				'lower lim',xx(j+1)
		else if (j == N) then
			write(*,'(1x,f12.6,2i6,f12.6,a12)') x,ji,j,&
				xx(j),'upper lim'
		else
			write(*,'(1x,f12.6,2i6,2f12.6)') x,ji,j,xx(j),xx(j+1)
		end if
	end do
	END PROGRAM xhunt
