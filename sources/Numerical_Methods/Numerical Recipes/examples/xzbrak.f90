	PROGRAM xzbrak
!	driver for routine zbrak
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=100
	REAL(SP), PARAMETER :: X1=1.0_sp,X2=50.0
	INTEGER(I4B) :: i,nb
	REAL(SP), DIMENSION(:), POINTER :: xb1,xb2
	call zbrak(bessj0_s,X1,X2,N,xb1,xb2,nb)
	write(*,'(/1x,a/)') 'Brackets for roots of BESSJ0:'
	write(*,'(/1x,t17,a,t27,a,t40,a,t50,a/)') 'lower','upper',&
		'F(lower)','F(upper)'
	do i=1,nb
		write(*,'(1x,a,i2,2(4x,2f10.4))') 'Root ',i,xb1(i),xb2(i),&
			bessj0(xb1(i)),bessj0(xb2(i))
	end do
	END PROGRAM xzbrak
