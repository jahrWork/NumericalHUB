	PROGRAM xlaguer
!	driver for routine laguer
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: M=4,MP1=M+1,NTRY=21
	REAL(SP), PARAMETER :: EPS=1.0e-6_sp
	INTEGER(I4B) :: i,its,n
	COMPLEX(SPC) :: x
	COMPLEX(SPC), DIMENSION(MP1) :: a = (/ (0.0_sp,2.0_sp),(0.0_sp,0.0_sp),&
		(-1.0_sp,-2.0_sp),(0.0_sp,0.0_sp),(1.0_sp,0.0_sp) /)
	COMPLEX(SP), DIMENSION(NTRY) :: y
	write(*,'(/1x,a)') 'Roots of polynomial x^4-(1+2i)*x^2+2i'
	write(*,'(/1x,t16,a4,t29,a7,t39,a5/)') 'Real','Complex','#iter'
	n=0
	do i=1,NTRY
		x=cmplx((i-11.0_sp)/10.0_sp,(i-11.0_sp)/10.0_sp,kind=spc)
		call laguer(a(1:MP1),x,its)
		if (n == 0) then
			n=1
			y(1)=x
			write(*,'(1x,i5,2f15.6,i5)') n,x,its
		else
			if (any(abs(x-y(1:n)) <= EPS*abs(x))) cycle
			n=n+1
			y(n)=x
			write(*,'(1x,i5,2f15.6,i5)') n,x,its
		end if
	end do
	END PROGRAM xlaguer
