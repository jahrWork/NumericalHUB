MODULE stats
	USE nrtype
	INTEGER(I4B) :: stats_ndfunc,stats_nfunc
END MODULE stats

	PROGRAM xdfpmin
!	driver for routine dfpmin
	USE nrtype
	USE nr
	USE stats
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDIM=2
	REAL(SP), PARAMETER :: GTOL=1.0e-4_sp
	INTEGER(I4B) :: iter
	REAL(SP) :: fret
	REAL(SP), DIMENSION(NDIM) :: p
	INTERFACE
		FUNCTION func(p)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP) :: func
		END FUNCTION func
!BL
		FUNCTION dfunc(p)
		USE nrtype
		REAL(SP), DIMENSION(:), INTENT(IN) :: p
		REAL(SP), DIMENSION(size(p)) :: dfunc
		END FUNCTION dfunc
	END INTERFACE
	write(*,'(/1x,a)') 'True minimum is at (-2.0,+-0.89442719)'
	stats_nfunc=0
	stats_ndfunc=0
	p(1:2)=(/ 0.1_sp,4.2_sp /)
	write(*,'(/1x,a,2(f7.4,a))') 'Starting vector: (',p(1),',',p(2),')'
	call dfpmin(p,GTOL,iter,fret,func,dfunc)
	write(*,'(1x,a,i3)') 'Iterations:',iter
	write(*,'(1x,a,i3)') 'Func. evals:',stats_nfunc
	write(*,'(1x,a,i3)') 'Deriv. evals:',stats_ndfunc
	write(*,'(1x,a,2(f9.6,a))') 'Solution vector: (',p(1),',',p(2),')'
	write(*,'(1x,a,e14.6)') 'Func. value at solution',fret
	END PROGRAM xdfpmin

	FUNCTION func(x)
	USE nrtype
	USE stats
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP) :: func
	stats_nfunc=stats_nfunc+1
	func=10.0_sp*(x(2)**2*(3.0_sp-x(1))-x(1)**2*(3.0_sp+x(1)))**2+&
		(2.0_sp+x(1))**2/(1.0_sp+(2.0_sp+x(1))**2)
	END FUNCTION func

	FUNCTION dfunc(x)
	USE nrtype
	USE nr
	USE stats
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: dfunc
	REAL(SP) :: x3m,x3p,x2p
	x3m=3.0_sp-x(1)
	x3p=3.0_sp+x(1)
	x2p=2.0_sp+x(1)
	stats_ndfunc=stats_ndfunc+1
	dfunc(1)=20.0_sp*(x(2)**2*x3m-x(1)**2*x3p)*(-x(2)**2-6.0_sp*&
		x(1)-3.0_sp*x(1)**2)+2.0_sp*x2p/(1.0_sp+x2p**2)-&
		2.0_sp*x2p**3/(1.0_sp+x2p**2)**2
	dfunc(2)=40.0_sp*(x(2)**2*x3m-x(1)**2*x3p)*x(2)*x3m
	END FUNCTION dfunc
