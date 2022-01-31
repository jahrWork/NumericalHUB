	PROGRAM xeclass
!	driver for routine eclass
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: N=15,M=11
	INTEGER(I4B) :: i,j,k,lclas,nclass
	INTEGER(I4B), DIMENSION(N) :: nf,nflag,nsav
	INTEGER(I4B), DIMENSION(M) :: &
		lista = (/ 1,1,5,2,6,2,7,11,3,4,12 /),&
		listb = (/ 5,9,13,6,10,14,3,7,15,8,4 /)
	nf(1:N)=eclass(lista,listb,N)
	nflag(:)=1
	write(*,'(/1x,a)') 'Numbers from 1-15 divided according to'
	write(*,'(1x,a/)') 'their value modulo 4:'
	lclas=0
	do i=1,N
		nclass=nf(i)
		if (nflag(nclass) /= 0) then
			nflag(nclass)=0
			lclas=lclas+1
			k=0
			do j=i,N
				if (nf(j) == nf(i)) then
					k=k+1
					nsav(k)=j
				end if
			end do
			write(*,'(1x,a,i2,a,3x,5i3)') 'Class',&
				lclas,':',(nsav(j),j=1,k)
		end if
	end do
	END PROGRAM xeclass
