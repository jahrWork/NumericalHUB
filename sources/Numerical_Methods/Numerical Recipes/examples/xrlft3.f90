	PROGRAM xrlft3
!	driver for routine rlft3
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	REAL(SP), PARAMETER :: EPS=3.0e-3_sp
	INTEGER(I4B), PARAMETER :: NX=16,NY=8,NZ=32
	INTEGER(I4B) :: i,j,ierr
	REAL(SP), DIMENSION(NX,NY,NZ) :: data1,data2
	COMPLEX(SP), DIMENSION(NX/2,NY,NZ) :: spec
	COMPLEX(SP), DIMENSION(NY,NZ) :: speq
	REAL(SP) :: fnorm
	INTEGER(I4B) :: nn1,nn2,nn3
	REAL(SP), DIMENSION(NZ) :: harvest
	call ran_seed(sequence=3)
	nn1=NX
	nn2=NY
	nn3=NZ
	fnorm=real(nn1,sp)*real(nn2,sp)*real(nn3,sp)/2.0_sp
	do i=1,nn1
		do j=1,nn2
			call ran1(harvest)
			data1(i,j,:)=2.0_sp*harvest(:)-1.0_sp
		end do
	end do
	data2(:,:,:)=data1(:,:,:)*fnorm
	call rlft3(data1,spec,speq,1)
!	here would be any processing in Fourier space
	call rlft3(data1,spec,speq,-1)
	ierr=icompare('data',reshape(data1, (/ nn1*nn2*nn3 /)), &
		reshape(data2, (/ nn1*nn2*nn3 /)) ,EPS)
	if (ierr == 0) then
		write(*,*) 'Data compares OK to tolerance',EPS
	else
		write(*,*) 'Comparison errors occured at tolerance',EPS
		write(*,*) 'Total number of errors is',ierr
	end if
	CONTAINS
!BL
	FUNCTION icompare(string,arr1,arr2,eps)
	INTEGER(I4B) :: icompare
	CHARACTER*(*), INTENT(IN) :: string
	REAL(SP), INTENT(IN) :: eps
	REAL(SP), DIMENSION(:), INTENT(IN) :: arr1,arr2
	INTEGER(I4B), PARAMETER :: IPRNT=20
	INTEGER(I4B) :: j,len
	len=size(arr1)
	if (len /= size(arr2)) call nrerror('xrlft3: comparision size error')
	write(*,*) string
	icompare=0
	do j=1,len
		if ((arr2(j) == 0 .and. abs(arr1(j)-arr2(j)) > eps) .or. &
			(abs((arr1(j)-arr2(j))/arr2(j)) > eps)) then
			icompare=icompare+1
			if (icompare <= IPRNT) write(*,*) j,arr1(j),arr2(j)
		end if
	end do
	END FUNCTION icompare
	END PROGRAM xrlft3
