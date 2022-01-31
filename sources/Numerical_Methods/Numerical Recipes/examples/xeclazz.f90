	PROGRAM xeclazz
!	driver for routine eclazz
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=15
	INTEGER(I4B) :: i,j,k,lclas,nclass
	INTEGER(I4B), DIMENSION(N) :: ind,nf,nflag,nsav
	LOGICAL(LGT), DIMENSION(N) :: eq
	INTERFACE
		FUNCTION equiv(i,j)
		USE nrtype; USE nrutil
		IMPLICIT NONE
		LOGICAL(LGT) :: equiv
		INTEGER(I4B), INTENT(IN) :: i,j
		END FUNCTION equiv
	END INTERFACE
	nf(1:N)=eclazz(equiv,N)
	nflag(:)=1
	ind(1:N)=arth(1,1,N)
	write(*,'(/1x,a)') 'Numbers from 1-15 divided according to'
	write(*,'(1x,a/)') 'their value modulo 4:'
	lclas=0
	do i=1,N
		nclass=nf(i)
		if (nflag(nclass) /= 0) then
			nflag(nclass)=0
			lclas=lclas+1
			eq(1:N)=(nf(1:N) == nf(i))
			k=count(eq)
			nsav(1:k)=pack(ind(1:N),eq(1:N))
			write(*,'(1x,a,i2,a,3x,5i3)') 'Class',&
				lclas,':',(nsav(j),j=1,k)
		end if
	end do
	END PROGRAM xeclazz

	FUNCTION equiv(i,j)
	USE nrtype; USE nrutil
	IMPLICIT NONE
	LOGICAL(LGT) :: equiv
	INTEGER(I4B), INTENT(IN) :: i,j
	equiv=merge(.true.,.false., mod(i,4) == mod(j,4) )
	END FUNCTION equiv
