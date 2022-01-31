	PROGRAM xcorrel
!	driver for routine correl
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=64,NA=17
	INTEGER(I4B) :: i
	INTEGER(I4B), DIMENSION(N) :: ind
	REAL(SP), DIMENSION(NA) :: cmp
	REAL(SP), DIMENSION(N) :: data1,data2,ans
	data1(:)=0.0_sp
	data1(N/2-N/8+1:N/2+N/8-1)=1.0_sp
	data2(:)=data1(:)
!	calculate directly
	write(*,'(/1x,t4,a,t13,a,t25,a/)') 'n','CORREL','Direct Calc.'
	ind(1:N)=arth(0,1,N)
	do i=1,NA
		cmp(i)=dot_product(data1(mod(i-1+ind(1:N),N)+1),data2(1:N))
	end do
	ans(:)=correl(data1,data2)
	do i=1,NA
		write(*,'(1x,i3,3x,f12.6,f15.6)') i-1,ans(i),cmp(i)
	end do
	END PROGRAM xcorrel
