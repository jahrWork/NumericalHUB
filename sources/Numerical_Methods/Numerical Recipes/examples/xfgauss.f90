	PROGRAM xfgauss
!	driver for routine fgauss
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPT=3,NLIN=2,NA=3*NLIN
	INTEGER(I4B) :: i,j
	REAL(SP), DIMENSION(NPT) :: e1,e2,f,x,y
	REAL(SP), DIMENSION(NA) :: a = (/ 3.0_sp,0.2_sp,0.5_sp,1.0_sp,0.7_sp,0.3_sp /)
	REAL(SP), DIMENSION(NPT,NA) :: df,dyda
	write(*,'(/1x,t7,a,t15,a,t20,a,t28,a,t36,a,t44,a,t52,a,t60,a)')&
		'X','Y','DYDA1','DYDA2','DYDA3','DYDA4','DYDA5','DYDA6'
	x(1:NPT)=0.3_sp*arth(1,1,NPT)
	call fgauss(x,a,y,dyda)
	e1(:)=exp(-((x(:)-a(2))/a(3))**2)
	e2(:)=exp(-((x(:)-a(5))/a(6))**2)
	f(:)=a(1)*e1(:)+a(4)*e2(:)
	df(:,1)=e1(:)
	df(:,4)=e2(:)
	df(:,2)=a(1)*e1(:)*2.0_sp*(x(:)-a(2))/(a(3)**2)
	df(:,5)=a(4)*e2(:)*2.0_sp*(x(:)-a(5))/(a(6)**2)
	df(:,3)=a(1)*e1(:)*2.0_sp*((x(:)-a(2))**2)/(a(3)**3)
	df(:,6)=a(4)*e2(:)*2.0_sp*((x(:)-a(5))**2)/(a(6)**3)
	do i=1,NPT
		write(*,'(1x,a)') 'from FGAUSS'
		write(*,'(1x,8f8.4)') x(i),y(i),(dyda(i,j),j=1,6)
		write(*,'(1x,a)') 'independent calc.'
		write(*,'(1x,8f8.4)') x(i),f(i),(df(i,j),j=1,6)
		write(*,*)
	end do
	END PROGRAM xfgauss
