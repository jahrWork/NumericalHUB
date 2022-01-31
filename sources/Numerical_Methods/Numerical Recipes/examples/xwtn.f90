	PROGRAM xwtn
!	driver for routine wtn
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NCMAX=50,NX=128,NY=256
	REAL(SP), PARAMETER :: EPS=1.0e-06_sp
	INTEGER(I4B) :: i,j,nerror,ntot
	REAL(SP), DIMENSION(NX,NY) :: a,aorg
	REAL(SP), DIMENSION(NX*NY) :: adata
	INTEGER(I4B), DIMENSION(2) :: ndim = (/ NX, NY /)
	nerror=0
	ntot=NX*NY
	do i=1,NX
		do j=1,NY
			if (i == j) then
				a(i,j)=-1.0
			else
				a(i,j)=1.0_sp/sqrt(abs(real(i-j,sp)))
			end if
		end do
	end do
	aorg(:,:)=a(:,:)
	call pwtset(12)
	adata(1:ntot)=reshape(a, (/ ntot /));
	call wtn(adata,ndim,1,pwt)
! here, one might set the smallest components to zero, encode and transmit
! the remaining components as a compressed form of the "image"
	call wtn(adata,ndim,-1,pwt)
	a(:,:)=reshape(adata, (/ NX, NY /));
	do i=1,NX
		do j=1,NY
			if (abs(aorg(i,j)-a(i,j)) >= EPS) then
			write(*,*) 'Compare Error at element ',i,j
			nerror=nerror+1
			end if
		end do
	end do
	if (nerror /= 0) then
		write(*,*) 'Number of comparison errors: ',nerror
	else
		write(*,*) 'Transform-inverse transform check OK'
	end if
	END PROGRAM xwtn
