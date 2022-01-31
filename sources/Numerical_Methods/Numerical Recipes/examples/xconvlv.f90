	PROGRAM xconvlv
!	driver for routine convlv
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=16,N2=32,M=9
	INTEGER(I4B) :: i,isign
	INTEGER(I4B), DIMENSION(M/2) :: ind
	REAL(SP), DIMENSION(M) :: respns
	REAL(SP), DIMENSION(N) :: ans,cmp,data,resp
	data(:)=0.0
	data(N/2-N/8:N/2+N/8)=1.0
	respns(1:M)=0.0
	respns(3:6)=1.0
	resp(1:M)=respns(1:M)
!	compare with a direct convolution
	ind(1:M/2)=arth(1,1,M/2)
	write(*,'(/1x,t4,a,t13,a,t24,a)') 'I','CONVLV','Expected'
	do i=1,N
		cmp(i)=dot_product(data(mod(i-ind(1:M/2)-1+N,N)+1),respns(ind(1:M/2)+1))+&
			dot_product(data(mod(i+ind(1:M/2)-1,N)+1),respns(M-ind(1:M/2)+1))+&
			data(i)*respns(1)
	end do
	isign=1
	ans(1:N)=convlv(data,resp(1:M),isign)
	do i=1,N
		write(*,'(1x,i3,3x,2f12.6)') i,ans(i),cmp(i)
	end do
	END PROGRAM xconvlv
