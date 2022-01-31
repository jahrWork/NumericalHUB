	PROGRAM xmppi
!	driver for mp routines
	USE nr
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B) :: n
	do
		write(*,*) 'INPUT N'
		read(*,*,END=999) n
		call mpsqr2(n)
		call mppi(n)
	end do
999	write(*,*) 'NORMAL COMPLETION'
	CONTAINS
!BL
	SUBROUTINE mpsqr2(n)
	USE nrtype
	USE nr
	USE mpops, ONLY : mpsad,mpmov
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: IAOFF=48,NMAX=8192
	INTEGER(I4B), INTENT(IN) :: n
	INTEGER(I4B) :: j,m
	CHARACTER(1), DIMENSION(n) :: x,q,t
	CHARACTER(1), DIMENSION(n+1) :: y,r
	CHARACTER(1), DIMENSION(3*n) :: s
	t(1)=char(2)
	t(2:n)=char(0)
	call mpsqrt(x,x,t,n,n)
	call mpmov(y,x,n)
	write(*,*) 'SQRT(2)='
	s(1)=char(ichar(y(1))+IAOFF)
	s(2)='.'
!	caution: next step is N**2! omit it for large N
	call mp2dfr(y(2:),s(3:),n-1,m)
	write(*,'(1x,64a1)') (s(j),j=1,m+2)
	write(*,*) 'Result rounded to 1 less base-256 place:'
!	use s as scratch space
	call mpsad(s,x,n,128)
	call mpmov(y,s(2:n),n-1)
	s(1)=char(ichar(y(1))+IAOFF)
	s(2)='.'
!	caution: next step is N**2! omit it for large N
	call mp2dfr(y(2:),s(3:),n-2,m)
	write(*,'(1x,64a1)') (s(j),j=1,m+2)
	write(*,*) '2-SQRT(2)='
!	calculate this the hard way to exercise the mpdiv function
	call mpdiv(q,r,t,x,n,n)
	s(1)=char(ichar(r(1))+IAOFF)
	s(2)='.'
!	caution: next step is N**2! omit it for large N
	call mp2dfr(r(2:),s(3:),n-1,m)
	write(*,'(1x,64a1)') (s(j),j=1,m+2)
	END SUBROUTINE mpsqr2
	END PROGRAM xmppi
