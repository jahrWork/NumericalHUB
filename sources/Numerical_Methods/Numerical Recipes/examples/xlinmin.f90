	PROGRAM xlinmin
!	driver for routine linmin
	USE nrtype
	USE nr
	USE f1dim_mod
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDIM=3
	INTEGER(I4B) :: i,j
	REAL(SP) :: fret,sr2,x
	REAL(SP), DIMENSION(NDIM), TARGET :: p,xi
	ncom=NDIM
	pcom=>p
	xicom=>xi
	write(*,'(/1x,a)') 'Minimum of a 3-D quadratic centered'
	write(*,'(1x,a)') 'at (1.0,1.0,1.0). Minimum is found'
	write(*,'(1x,a)') 'along a series of radials.'
	write(*,'(/1x,t10,a,t22,a,t34,a,t42,a/)') 'x','y','z','minimum'
	sr2=sqrt(2.0_sp)
	do i=0,10
		x=PIO2*i/10.0_sp
		xi(1)=sr2*cos(x)
		xi(2)=sr2*sin(x)
		xi(3)=1.0
		p(:)=0.0
		call linmin(p,xi,fret)
		write(*,'(1x,4f12.6)') (p(j),j=1,3),fret
	end do
	END PROGRAM xlinmin

	FUNCTION func(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP) :: func
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: tmp
	tmp(:)=x(:)-1.0_sp
	func=dot_product(tmp,tmp)
	END FUNCTION func
