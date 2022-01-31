	PROGRAM xmnbrak
!	driver for routine mnbrak
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i
	REAL(SP) :: ax,bx,cx,fa,fb,fc
	do i=1,10
		ax=i*0.5_sp
		bx=(i+1.0_sp)*0.5_sp
		call mnbrak(ax,bx,cx,fa,fb,fc,bessj0_s)
		write(*,'(1x,t13,a,t25,a,t37,a)') 'A','B','C'
		write(*,'(1x,a3,t5,3f12.6)') 'X',ax,bx,cx
		write(*,'(1x,a3,t5,3f12.6)') 'F',fa,fb,fc
	end do
	END PROGRAM xmnbrak
