	PROGRAM xpsdes
!	driver for routine psdes
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i
	INTEGER(I4B), DIMENSION(4) :: lcalc,lcalc2,ircalc,ircalc2,&
		lword = (/ 1,1,99,99 /),irword = (/ 1,99,1,99 /),&
		lans = (/ 1615666638,-645954191,2015506589,-671910160 /),&
		irans = (/ 1352404003,-1502825446,1680869764,1505397227 /)
	write (*,*) 'scalar test'
	lcalc(:)=lword(:)
	ircalc(:)=irword(:)
	do i=1,4
		call psdes(lcalc(i),ircalc(i))
		write(*,*) 'PSDES now calculates:       ',&
			hexout(lcalc(i)),'   ',hexout(ircalc(i))
		write(*,*) 'Known correct answers are:  ',&
			hexout(lans(i)),'   ',hexout(irans(i))
	end do
	write (*,*) 'vector test'
	lcalc2(:)=lword(:)
	ircalc2(:)=irword(:)
	call psdes(lcalc2(:),ircalc2(:))
	do i=1,4
		write(*,*) 'PSDES now calculates:       ',&
			hexout(lcalc2(i)),'   ',hexout(ircalc2(i))
		write(*,*) 'Known correct answers are:  ',&
			hexout(lans(i)),'   ',hexout(irans(i))
	end do
	ircalc2(:)=abs(ircalc2(:)-ircalc(:))
	lcalc2(:)=abs(lcalc2(:)-lcalc(:))
	if (maxval(ircalc2(:)) /= 0 .or. maxval(lcalc2(:)) /= 0) then
		write (*,*) 'ERROR: psdes scalar - vector difference'
		do i=1,4
			if (ircalc2(i) /= 0 .or. lcalc2(i) /= 0) &
				write (*,*) 'Difference at ',i
		end do
	end if
	CONTAINS
!BL
	FUNCTION hexout(num)
!	utility routine for printing out hexadecimal values
!	(many compilers can do this using the nonstandard "Z" format)
	IMPLICIT NONE
	CHARACTER(10) :: hexout
	INTEGER(I4B), PARAMETER :: NCOMP=268435455
	INTEGER(I4B), INTENT(IN) :: num
	INTEGER(I4B) :: i,j,n
	CHARACTER(1), DIMENSION(16) :: hexit = (/ &
		'0','1','2','3','4','5','6','7','8','9',&
		'a','b','c','d','e','f' /)
	n=num
	if (n < 0) then
		i=mod(n,16)
		if (i < 0) i=16+i
		n=n/16
		n=-n
		n=merge(NCOMP-n,NCOMP+1-n, i /= 0 )
	else
		i=mod(n,16)
		n=n/16
	end if
	j=10
	hexout(j:j)=hexit(i+1)
	do
		if (n <= 0) exit
		i=mod(n,16)
		n=n/16
		j=j-1
		hexout(j:j)=hexit(i+1)
	end do
	j=j-1
	hexout(j:j)='x'
	j=j-1
	hexout(j:j)='0'
	do i=j-1,1,-1
		hexout(i:i)=' '
	end do
	END FUNCTION hexout
	END PROGRAM xpsdes
