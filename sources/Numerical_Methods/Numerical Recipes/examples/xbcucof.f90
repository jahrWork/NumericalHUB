	PROGRAM xbcucof
!	driver for routine bcucof
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i,j
	REAL(SP) :: d1,d2,ee,x1x2
	REAL(SP), DIMENSION(4) :: y,y1,y2
	REAL(SP), DIMENSION(4,4) :: c
	REAL(SP), DIMENSION(4) :: y12,&
		x1 = (/ 0.0_sp,2.0_sp,2.0_sp,0.0_sp /),&
		x2 = (/ 0.0_sp,0.0_sp,2.0_sp,2.0_sp /)
	d1=x1(2)-x1(1)
	d2=x2(4)-x2(1)
	do i=1,4
		x1x2=x1(i)*x2(i)
		ee=exp(-x1x2)
		y(i)=x1x2*ee
		y1(i)=x2(i)*(1.0_sp-x1x2)*ee
		y2(i)=x1(i)*(1.0_sp-x1x2)*ee
		y12(i)=(1.0_sp-3.0_sp*x1x2+x1x2**2)*ee
	end do
	call bcucof(y,y1,y2,y12,d1,d2,c)
	write(*,*) 'Coefficients for bicubic interpolation'
	do i=1,4
		write(*,'(1x,4e15.6)') (c(i,j),j=1,4)
	end do
	END PROGRAM xbcucof
