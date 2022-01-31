	PROGRAM xmnewt
!	driver for routine mnewt
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	REAL(SP), PARAMETER :: TOLX=1.0e-6_sp,TOLF=1.0e-6_sp
	INTEGER(I4B), PARAMETER :: NTRIAL=5,N=4
	INTEGER(I4B) :: i,j,k,kk
	REAL(SP) :: xx
	REAL(SP), DIMENSION(N) :: fvec,x
	REAL(SP), DIMENSION(N,N) :: fjac
	INTERFACE
		SUBROUTINE usrfun(x,fvec,fjac)
		USE nrtype; USE nrutil
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
		REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
		END SUBROUTINE usrfun
	END INTERFACE
	do kk=-1,1,2
		do k=1,3
			xx=0.2001_sp*k*kk
			write(*,'(/1x,a,i2)') 'Starting vector number',k
			x(1:N)=xx+0.2_sp*arth(1,1,N)
			do i=1,N
				write(*,'(1x,t5,a,i1,a,f5.2)') 'X(',i,') = ',x(i)
			end do
			do j=1,NTRIAL
				call mnewt(1,x(1:N),TOLX,TOLF,usrfun)
				call usrfun(x(1:n),fvec,fjac)
				write(*,'(/1x,t5,a,t14,a,t29,a/)') 'I','X(I)','F'
				do i=1,N
					write(*,'(1x,i4,2e15.6)') i,x(i),fvec(i)
				end do
				write(*,'(/1x,a)') 'press RETURN to continue...'
				read(*,*)
			end do
		end do
	end do
	END PROGRAM xmnewt

	SUBROUTINE usrfun(x,fvec,fjac)
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: fvec
	REAL(SP), DIMENSION(:,:), INTENT(OUT) :: fjac
	INTEGER(I4B) :: n
	n=size(x)
	fjac(1,1:4) = (/ -2.0_sp*x(1), -2.0_sp*x(2), -2.0_sp*x(3), 1.0_sp /)
	fjac(2,1:n)=2.0_sp*x(1:n)
	fjac(3,1:4) = (/ 1.0_sp, -1.0_sp, 0.0_sp, 0.0_sp /)
	fjac(4,1:4) = (/ 0.0_sp, 1.0_sp, -1.0_sp, 0.0_sp /)
	fvec(1)=-x(1)**2-x(2)**2-x(3)**2+x(4)
	fvec(2)=x(1)**2+x(2)**2+x(3)**2+x(4)**2-1.0_sp
	fvec(3)=x(1)-x(2)
	fvec(4)=x(2)-x(3)
	END SUBROUTINE usrfun
