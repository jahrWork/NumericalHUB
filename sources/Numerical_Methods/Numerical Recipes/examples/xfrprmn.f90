	PROGRAM xfrprmn
!	driver for routine frprmn
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDIM=3
	REAL(SP), PARAMETER :: FTOL=1.0e-6_sp
	INTEGER(I4B) :: iter,k
	REAL(SP) :: angl,fret
	REAL(SP), DIMENSION(NDIM) :: p
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP) :: func
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
!BL
		FUNCTION dfunc(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	write(*,'(/1x,a)') 'Program finds the minimum of a function'
	write(*,'(1x,a)') 'with different trial starting vectors.'
	write(*,'(1x,a)') 'True minimum is (0.5,0.5,0.5)'
	do k=0,4
		angl=PIO2*k/4.0_sp
		p(1)=2.0_sp*cos(angl)
		p(2)=2.0_sp*sin(angl)
		p(3)=0.0
		write(*,'(/1x,a,3(f6.4,a))') 'Starting vector: (',&
			p(1),',',p(2),',',p(3),')'
		call frprmn(p,FTOL,iter,fret)
		write(*,'(1x,a,i3)') 'Iterations:',iter
		write(*,'(1x,a,3(f6.4,a))') 'Solution vector: (',&
			p(1),',',p(2),',',p(3),')'
		write(*,'(1x,a,e14.6)') 'Func. value at solution',fret
	end do
	END PROGRAM xfrprmn

	FUNCTION func(x)
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP) :: func
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	func=1.0_sp-product(bessj0(x(:)-0.5_sp))
	END FUNCTION func

	FUNCTION dfunc(x)
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: dfunc
	dfunc(1)=bessj1(x(1)-0.5_sp)*bessj0(x(2)-0.5_sp)*bessj0(x(3)-0.5_sp)
	dfunc(2)=bessj0(x(1)-0.5_sp)*bessj1(x(2)-0.5_sp)*bessj0(x(3)-0.5_sp)
	dfunc(3)=bessj0(x(1)-0.5_sp)*bessj0(x(2)-0.5_sp)*bessj1(x(3)-0.5_sp)
	END FUNCTION dfunc
