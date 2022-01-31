	PROGRAM xtwofft
!	driver for routine twofft
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	REAL(SP), PARAMETER :: PER=8.0
	INTEGER(I4B), PARAMETER :: N=32,N2=2*N
	INTEGER(I4B) :: isign
	REAL(SP), DIMENSION(N) :: data1,data2,x
	REAL(SP), DIMENSION(N2) :: fft1,fft2
	x(:)=2.0_sp*PI*arth(1,1,N)/PER
	data1(:)=nint(cos(x(:)))
	data2(:)=nint(sin(x(:)))
	call rtwofft(data1,data2,fft1,fft2)
	write(*,*) 'Fourier transform of first function:'
	call prntft(fft1)
	write(*,*) 'Fourier transform of second function:'
	call prntft(fft2)
!	invert transform
	isign=-1
	call rfour1(fft1,isign)
	write(*,*) 'Inverted transform = first function:'
	call prntft(fft1)
	call rfour1(fft2,isign)
	write(*,*) 'Inverted transform = second function:'
	call prntft(fft2)
	CONTAINS
!BL
	SUBROUTINE prntft(data)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data
	INTEGER(I4B) :: i,nn2,m
	INTEGER(I4B) :: n2
	n2=size(data)
	write(*,'(1x,t7,a,t13,a,t24,a,t35,a,t47,a)')&
		'n','Real(n)','Imag.(n)','Real(N-n)','Imag.(N-n)'
	write(*,'(1x,i6,4f12.6)') 0,data(1),data(2),data(1),data(2)
	do i=3,(n2/2)+1,2
		m=(i-1)/2
		nn2=n2+2-i
		write(*,'(1x,i6,4f12.6)') m,data(i),data(i+1),&
			data(nn2),data(nn2+1)
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
!BL
	SUBROUTINE rtwofft(data1,data2,fft1,fft2)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: data1,data2
	REAL(SP), DIMENSION(:), INTENT(OUT) :: fft1,fft2
	COMPLEX(SPC), DIMENSION(size(fft1)/2) :: cfft1,cfft2
	INTEGER(I4B) :: n,m
	n=size(fft1)
	m=n/2
	call twofft(data1,data2,cfft1(1:m),cfft2(1:m))
	fft1(1:n:2)=real(cfft1(1:m))
	fft1(2:n:2)=aimag(cfft1(1:m))
	fft2(1:n:2)=real(cfft2(1:m))
	fft2(2:n:2)=aimag(cfft2(1:m))
	END SUBROUTINE rtwofft
	END PROGRAM xtwofft
