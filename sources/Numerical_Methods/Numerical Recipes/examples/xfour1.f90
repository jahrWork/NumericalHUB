	PROGRAM xfour1
!	driver for routine four1
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NN=32,NN2=2*NN
	REAL(SP), DIMENSION(NN2) :: data,dcmp
	REAL(SP), DIMENSION(NN) :: tmp
	INTEGER(I4B) :: i,isign,j
	write(*,*) 'h(t)=real-valued even-function'
	write(*,*) 'H(n)=H(N-n) and real?'
	tmp=arth(-NN,2,NN)/real(NN,sp)
	data(1:2*NN-1:2)=1.0_sp/(tmp(:)**2+1.0_sp)
	data(2:2*NN:2)=0.0
	isign=1
	call rfour1(data,isign)
	call prntft(data)
	write(*,*) 'h(t)=imaginary-valued even-function'
	write(*,*) 'H(n)=H(N-n) and imaginary?'
	data(2:2*NN:2)=1.0_sp/(tmp(:)**2+1.0_sp)
	data(1:2*NN-1:2)=0.0
	isign=1
	call rfour1(data,isign)
	call prntft(data)
	write(*,*) 'h(t)=real-valued odd-function'
	write(*,*) 'H(n)=-H(N-n) and imaginary?'
	data(1:2*NN-1:2)=tmp(:)/(tmp(:)**2+1.0_sp)
	data(2:2*NN:2)=0.0
	data(1)=0.0
	isign=1
	call rfour1(data,isign)
	call prntft(data)
	write(*,*) 'h(t)=imaginary-valued odd-function'
	write(*,*) 'H(n)=-H(N-n) and real?'
	data(2:2*NN:2)=tmp(:)/(tmp(:)**2+1.0_sp)
	data(1:2*NN-1:2)=0.0
	data(2)=0.0
	isign=1
	call rfour1(data,isign)
	call prntft(data)
!	transform, inverse-transform test
	data(1:2*NN-1:2)=1.0_sp/((0.5_sp*tmp(:))**2+1.0_sp)
	dcmp(1:2*NN-1:2)=data(1:2*NN-1:2)
	data(2:2*NN:2)=(0.25_sp*tmp(:))*exp(-(0.5_sp*tmp(:))**2)
	dcmp(2:2*NN:2)=data(2:2*NN:2)
	isign=1
	call rfour1(data,isign)
	isign=-1
	call rfour1(data,isign)
	write(*,'(/1x,t10,a,t44,a)') 'Original Data:',&
		'Double Fourier Transform:'
	write(*,'(/1x,t5,a,t11,a,t24,a,t41,a,t53,a/)')&
		'k','Real h(k)','Imag h(k)','Real h(k)','Imag h(k)'
	do i=1,NN,2
		j=(i+1)/2
		write(*,'(1x,i4,2x,2f12.6,5x,2f12.6)') j,dcmp(i),&
			dcmp(i+1),data(i)/NN,data(i+1)/NN
	end do
	CONTAINS
!BL
	SUBROUTINE prntft(data)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data
	INTEGER(I4B) :: n,nn2,m,mm
	nn2=size(data)
	write(*,'(/1x,t5,a,t11,a,t23,a,t39,a,t52,a)')&
		'n','Real H(n)','Imag H(n)','Real H(N-n)','Imag H(N-n)'
	write(*,'(1x,i4,2x,2f12.6,5x,2f12.6)') 0,data(1),data(2),&
		data(1),data(2)
	do n=3,(nn2/2)+1,2
		m=(n-1)/2
		mm=nn2+2-n
		write(*,'(1x,i4,2x,2f12.6,5x,2f12.6)') m,data(n),&
			data(n+1),data(mm),data(mm+1)
	end do
	write(*,'(/1x,a)') ' press RETURN to continue ...'
	read(*,*)
	END SUBROUTINE prntft
!BL
	SUBROUTINE rfour1(rdata,isign)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: rdata
	INTEGER(I4B), INTENT(IN) :: isign
	COMPLEX(SPC), DIMENSION(size(rdata)/2) :: cdata
	INTEGER(I4B) :: n,m
	n=size(rdata)
	m=n/2
	cdata(1:m)=cmplx(rdata(1:n:2),rdata(2:n:2),kind=spc)
	call four1(cdata(1:m),isign)
	rdata(1:n:2)=real(cdata(1:m))
	rdata(2:n:2)=aimag(cdata(1:m))
	END SUBROUTINE rfour1
	END PROGRAM xfour1
