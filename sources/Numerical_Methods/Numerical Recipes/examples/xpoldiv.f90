	PROGRAM xpoldiv
!	driver for routine poldiv
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=6,NV=4
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(N) :: q,r,u = &
		(/ -1.0_sp,5.0_sp,-10.0_sp,10.0_sp,-5.0_sp,1.0_sp /)
	REAL(SP), DIMENSION(NV) :: v = (/ 1.0_sp,3.0_sp,3.0_sp,1.0_sp /)
	call poldiv(u,v,q,r)
	write(*,'(//1x,6(7x,a)/)') 'X^0','X^1','X^2','X^3','X^4','X^5'
	write(*,*) 'Quotient polynomial coefficients:'
	write(*,'(1x,6f10.2/)') (q(i),i=1,6)
	write(*,*) 'Expected quotient coefficients:'
	write(*,'(1x,6f10.2///)') 31.0_sp,-8.0_sp,1.0_sp,0.0_sp,0.0_sp,0.0_sp
	write(*,*) 'Remainder polynomial coefficients:'
	write(*,'(1x,4f10.2/)') (r(i),i=1,4)
	write(*,*) 'Expected remainder coefficients:'
	write(*,'(1x,4f10.2//)') -32.0_sp,-80.0_sp,-80.0_sp,0.0_sp
	END PROGRAM xpoldiv
