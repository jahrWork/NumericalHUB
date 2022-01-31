	PROGRAM xpccheb
!	driver for routine pccheb
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NCHECK=15,NFEW=13,NMANY=17
	INTEGER(I4B) :: i,j
	INTEGER(I4B), DIMENSION(NMANY) :: nser
	REAL(SP) :: a,b,f,sum,sume,py
	REAL(SP), DIMENSION(NFEW) :: fac
	REAL(SP), DIMENSION(NMANY) :: c,d,e,ee
!	put power series of cos(PI*y) + sin(PI*y) into e
	nser(1:NMANY)=arth(1,1,NMANY)
	e(1:NMANY)=1.0_sp/cumprod(real(nser,sp))
	nser(1:NMANY)=mod(nser(1:NMANY)-1,4)
	where (nser(1:NMANY) == 2 .or. nser(1:NMANY) == 3) &
		e(1:NMANY)=-e(1:NMANY)
	ee(:)=e(:)
	a=-PI
	b=PI
	call pcshft((-2.0_sp-b-a)/(b-a),(2.0_sp-b-a)/(b-a),e(1:NMANY))
!	i.e., inverse of PCSHFT(A,B,...) which we do below
	c(1:NMANY)=pccheb(e(1:NMANY))
	write(*,*) 'Index, series, Chebyshev coefficients'
	do j=1,NMANY
		write(*,'(i3,2e15.6)') j,e(j),c(j)
	end do
	d(1:NFEW)=chebpc(c(1:NFEW))
	call pcshft(a,b,d(1:NFEW))
	write(*,*) 'Index, new series, coefficient ratios'
	do j=1,NFEW
		write(*,'(i3,2e15.6)') j,d(j),d(j)/(ee(j)+1.0e-30_sp)
	end do
	write(*,'(7x,a)')&
		'Point tested, function value, error power series, error Cheb.'
	do i=0,15
		py=(a+i*(b-a)/15.0_sp)
		fac(1:NFEW)=geop(1.0_sp,py,NFEW)
		sum=dot_product(fac(1:NFEW),d(1:NFEW))
		sume=dot_product(fac(1:NFEW),ee(1:NFEW))
		f=cos(py)+sin(py)
		write(*,'(1x,a,4e15.6)') 'check:',py,f,sume-f,sum-f
	end do
	END PROGRAM xpccheb
