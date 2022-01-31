	PROGRAM xbcuint
!	driver for routine bcuint
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i
	REAL(SP) :: ansy,ansy1,ansy2,ey,ey1,ey2
	REAL(SP) :: x1,x1l,x1u,x1x2,x2,x2l,x2u,xxyy
	REAL(SP), DIMENSION(4) :: y,y1,y2,y12,&
		xx = (/ 0.0_sp,2.0_sp,2.0_sp,0.0_sp /),&
		yy = (/ 0.0_sp,0.0_sp,2.0_sp,2.0_sp /)
	x1l=xx(1)
	x1u=xx(2)
	x2l=yy(1)
	x2u=yy(4)
	do i=1,4
		xxyy=xx(i)*yy(i)
		y(i)=xxyy**2
		y1(i)=2.0_sp*yy(i)*xxyy
		y2(i)=2.0_sp*xx(i)*xxyy
		y12(i)=4.0_sp*xxyy
	end do
	write(*,'(/1x,t6,a,t14,a,t22,a,t28,a,t38,a,t44,a,t54,a,t60,a/)')&
		'X1','X2','Y','EXPECT','Y1','EXPECT','Y2','EXPECT'
	do i=1,10
		x1=0.2_sp*i
		x2=x1
		call bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
		x1x2=x1*x2
		ey=x1x2**2
		ey1=2.0_sp*x2*x1x2
		ey2=2.0_sp*x1*x1x2
		write(*,'(1x,8f8.4)') x1,x2,ansy,ey,ansy1,ey1,ansy2,ey2
	end do
	END PROGRAM xbcuint
