	PROGRAM xddpoly
!	driver for routine ddpoly
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NC=6,NCM1=5,NP=20
	INTEGER(I4B) :: i,j
	REAL(SP) :: x
	REAL(SP), DIMENSION(NC) :: c = (/ &
		-1.0_sp,5.0_sp,-10.0_sp,10.0_sp,-5.0_sp,1.0_sp /)
	REAL(SP), DIMENSION(NCM1) :: pd
	REAL(SP), DIMENSION(NCM1,NP) :: d
	CHARACTER(15), DIMENSION(NCM1) :: a = (/ &
		'polynomial:  ','first deriv: ','second deriv:',&
		'third deriv: ','fourth deriv:' /)
	do i=1,NP
		x=0.1_sp*i
		call ddpoly(c(1:NC),x,pd(1:NC-1))
		d(1:NC-1,i)=pd(1:NC-1)
	end do
	do i=1,NC-1
		write(*,'(1x,t7,a)') a(i)
		write(*,'(1x,t13,a,t25,a,t40,a)') 'X','DDPOLY','actual'
		do j=1,NP
			x=0.1_sp*j
			write(*,'(1x,3f15.6)') x,d(i,j),&
				factrl(NC-1)/factrl(NC-i)*((x-1.0_sp)**(NC-i))
		end do
		write(*,*) 'press ENTER to continue...'
		read(*,*)
	end do
	END PROGRAM xddpoly
