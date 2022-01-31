	PROGRAM xran
!	driver for routines ran0, ran1, ran2, ran3
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTERFACE
		SUBROUTINE ran_port(harvest)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran_port
!BL
		SUBROUTINE integ_s(ran_func)
		USE nrtype
		IMPLICIT NONE
		INTERFACE
			SUBROUTINE ran_func(harvest)
			USE nrtype
			IMPLICIT NONE
			REAL(SP), INTENT(OUT) :: harvest
			END SUBROUTINE ran_func
		END INTERFACE
		END SUBROUTINE integ_s
!BL
		SUBROUTINE integ_v(ran_func)
		USE nrtype
		IMPLICIT NONE
		INTERFACE
			SUBROUTINE ran_func(harvest)
			USE nrtype
			IMPLICIT NONE
			REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
			END SUBROUTINE ran_func
		END INTERFACE
		END SUBROUTINE integ_v
	END INTERFACE
	write(*,*) 'Scalar Tests'
	write(*,*)
	write(*,*) 'Testing ran:'
	call integ_s(ran_port)
	write(*,*) 'Testing ran0:'
	call integ_s(ran0_s)
	write(*,*)
	write(*,*) 'Testing ran1:'
	call integ_s(ran1_s)
	write(*,*)
	write(*,*) 'Testing ran2:'
	call integ_s(ran2_s)
	write(*,*)
	write(*,*) 'Testing ran3:'
	call integ_s(ran3_s)
	write(*,*)
	write(*,*) 'Vector Tests'
	write(*,*)
	write(*,*) 'Testing ran0:'
	call integ_v(ran0_v)
	write(*,*)
	write(*,*) 'Testing ran1:'
	call integ_v(ran1_v)
	write(*,*)
	write(*,*) 'Testing ran2:'
	call integ_v(ran2_v)
	write(*,*)
	write(*,*) 'Testing ran3:'
	call integ_v(ran3_v)
	END PROGRAM xran

	SUBROUTINE ran_port(harvest)
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: idum=-1
	REAL(SP), INTENT(OUT) :: harvest
	harvest=ran(idum)
	END SUBROUTINE ran_port

	SUBROUTINE integ_s(ran_func)
	USE nrtype
	USE nr
	IMPLICIT NONE
!	calculates pi statistically using volume of unit n-sphere
	INTERFACE
		SUBROUTINE ran_func(harvest)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(OUT) :: harvest
		END SUBROUTINE ran_func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NMAX=15
	INTEGER(I4B) :: i,j,k
	INTEGER(I4B), DIMENSION(3) :: iy
	REAL(SP) :: x1,x2,x3,x4
	REAL(SP), DIMENSION(3) :: yprob
	iy(:)=0
	write(*,'(1x,t15,a)') 'Volume of unit n-sphere, n = 2, 3, 4'
	write(*,'(1x,/,t3,a,t17,a,t26,a,t37,a)')&
		'# points','PI','(4/3)*PI','(1/2)*PI^2'
	do j=1,NMAX
		do k=2**(j-1),0,-1
			call ran_func(x1)
			call ran_func(x2)
			call ran_func(x3)
			call ran_func(x4)
			if (fnc(x1,x2,0.0_sp,0.0_sp) < 1.0) iy(1)=iy(1)+1
			if (fnc(x1,x2,x3,0.0_sp) < 1.0) iy(2)=iy(2)+1
			if (fnc(x1,x2,x3,x4) < 1.0) iy(3)=iy(3)+1
		end do
		yprob(:)=1.0_sp*(2**(/ 2,3,4 /))*iy(:)/(2**j)
		write(*,'(1x,i8,3f12.6)') 2**j,(yprob(i),i=1,3)
	end do
	write(*,'(1x,/,t4,a,3f12.6,/)') 'actual',PI,4.0*PI/3.0,0.5*(PI**2)
	CONTAINS
!BL
	FUNCTION fnc(x1,x2,x3,x4)
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x1,x2,x3,x4
	REAL(SP) :: fnc
	fnc=sqrt(x1**2+x2**2+x3**2+x4**2)
	END FUNCTION fnc
	END SUBROUTINE integ_s

	SUBROUTINE integ_v(ran_func)
	USE nrtype
	USE nr
	IMPLICIT NONE
!	calculates pi statistically using volume of unit n-sphere
	INTERFACE
		SUBROUTINE ran_func(harvest)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(OUT) :: harvest
		END SUBROUTINE ran_func
	END INTERFACE
	INTEGER(I4B), PARAMETER :: NMAX=13
	INTEGER(I4B) :: i,j,k
	INTEGER(I4B), DIMENSION(3) :: iy
	REAL(SP), DIMENSION(2**NMAX) :: x1,x2,x3,x4
	REAL(SP), DIMENSION(3) :: yprob
	iy(:)=0
	write(*,'(1x,t15,a)') 'Volume of unit n-sphere, n = 2, 3, 4'
	write(*,'(1x,/,t3,a,t17,a,t26,a,t37,a)')&
		'# points','PI','(4/3)*PI','(1/2)*PI^2'
	do j=1,NMAX
		k=2**(j-1)+1
		call ran_func(x1(1:k))
		call ran_func(x2(1:k))
		call ran_func(x3(1:k))
		call ran_func(x4(1:k))
		iy(1)=iy(1)+count(fnc(x1(1:k),x2(1:k)) < 1.0)
		iy(2)=iy(2)+count(fnc(x1(1:k),x2(1:k),x3(1:k)) < 1.0)
		iy(3)=iy(3)+count(fnc(x1(1:k),x2(1:k),x3(1:k),x4(1:k)) < 1.0)
		yprob(:)=1.0_sp*(2**(/ 2,3,4 /))*iy(:)/(2**j)
		write(*,'(1x,i8,3f12.6)') 2**j,(yprob(i),i=1,3)
	end do
	write(*,'(1x,/,t4,a,3f12.6,/)') 'actual',PI,4.0*PI/3.0,0.5*(PI**2)
	CONTAINS
!BL
	FUNCTION fnc(x1,x2,x3,x4)
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x1,x2,x3,x4
	REAL(SP), DIMENSION(size(x1)) :: fnc
	OPTIONAL :: x3,x4
	if (present(x4)) then
		fnc=sqrt(x1**2+x2**2+x3**2+x4**2)
	else if (present(x3)) then
		fnc=sqrt(x1**2+x2**2+x3**2)
	else
		fnc=sqrt(x1**2+x2**2)
	end if
	END FUNCTION fnc
	END SUBROUTINE integ_v
