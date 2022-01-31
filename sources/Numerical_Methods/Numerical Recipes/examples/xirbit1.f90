	PROGRAM xirbit1
!	driver for routine irbit1
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NBIN=15,NTRIES=10000
	INTEGER(I4B) :: i,idum,iflg,ipts,iseed,j,n
	REAL(SP), DIMENSION(NBIN) :: delay
	iseed=12345
	delay(:)=0.0
	ipts=0
	do i=1,NTRIES
		if (irbit1(iseed) == 1) then
			ipts=ipts+1
			iflg=0
			do j=1,NBIN
				idum=irbit1(iseed)
				if ((idum == 1) .and. (iflg == 0)) then
					iflg=1
					delay(j)=delay(j)+1.0_sp
				end if
			end do
		end if
	end do
	write(*,*) 'Distribution of runs of N zeros'
	write(*,'(1x,t7,a,t16,a,t38,a)') 'N','Probability','Expected'
	do n=1,NBIN
		write(*,'(1x,i6,f18.6,f20.6)')&
			n-1,delay(n)/ipts,1/(2.0_sp**n)
	end do
	END PROGRAM xirbit1
