	PROGRAM xfixrts
!	driver for routine fixrts
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPOLES=6,NPOL=NPOLES+1
	INTEGER(I4B) :: i
	REAL(SP), DIMENSION(NPOLES) :: d = &
		(/ 6.0_sp,-15.0_sp,20.0_sp,-15.0_sp,6.0_sp,0.0_sp /)
	COMPLEX(SPC) :: z
	COMPLEX(SPC), DIMENSION(NPOL) :: zcoef
	COMPLEX(SPC), DIMENSION(NPOLES) :: zeros
	LOGICAL(LGT) :: polish
!	finding roots of (z-1.0_sp)**6=1.0
!	first print roots
	zcoef(NPOLES+1)=cmplx(1.0_sp,0.0_sp,kind=spc)
	zcoef(NPOLES:1:-1)=cmplx(-d(1:NPOLES),0.0_sp,kind=spc)
	polish=.true.
	call zroots(zcoef,zeros,polish)
	write(*,'(/1x,a)') 'Roots of (z-1.0)^6 = 1.0'
	write(*,'(1x,t20,a,t42,a)') 'Root','(z-1.0)^6'
	call printorder(zeros)
	do i=1,NPOLES
		z=(zeros(i)-1.0_sp)**6
		write(*,'(1x,i6,4f12.6)') i,zeros(i),z
	end do
!	now fix them to lie within unit circle
	call fixrts(d)
!	check results
	zcoef(NPOLES+1)=cmplx(1.0_sp,0.0_sp,kind=spc)
	zcoef(NPOLES:1:-1)=cmplx(-d(1:NPOLES),0.0_sp,kind=spc)
	call zroots(zcoef,zeros,polish)
	call printorder(zeros)
	write(*,'(/1x,a)') 'Roots reflected in unit circle'
	write(*,'(1x,t20,a,t42,a)') 'Root','(z-1.0)^6'
	do i=1,NPOLES
		z=(zeros(i)-1.0_sp)**6
		write(*,'(1x,i6,4f12.6)') i,zeros(i),z
	end do
	CONTAINS
!BL
	SUBROUTINE printorder(zarr)
	IMPLICIT NONE
	COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: zarr
	INTEGER(I4B) :: i,n
	INTEGER(I4B), DIMENSION(size(zarr)) :: ind
	REAL(SP), PARAMETER :: EPS=1.0e-5_sp
	REAL(SP), DIMENSION(size(zarr)) :: tmp
	n=size(zarr)
	ind(1:n)=arth(1,1,n)
	tmp(:)=aimag(zarr)
	call indexx(tmp,ind)
	zarr(:)=zarr(ind)
	do i=1,n-1
		if ( abs(aimag(zarr(i)) - aimag(zarr(i+1))) < EPS .and. &
			real(zarr(i)) > real(zarr(i+1)) ) &
			call swap(zarr(i),zarr(i+1))
	end do
	END SUBROUTINE printorder
	END PROGRAM xfixrts
