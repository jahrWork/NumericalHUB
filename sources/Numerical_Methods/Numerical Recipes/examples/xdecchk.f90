	PROGRAM xdecchk
!	driver for routine decchk
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: j,k,l,n,nbad,ntot
	LOGICAL(LGT) :: iok,jok
	CHARACTER(1) :: ch,chh,tf
	CHARACTER(1), DIMENSION(128) :: lin
	CHARACTER(128) :: line
!	test all jump transpositions of the form 86jlk41
	ntot=0
	nbad=0
	do j=48,57
	do k=48,57
	do l=48,57
		if (j /= k) then
			ntot=ntot+1
			lin(1:7)=(/ '8','6',char(j),char(l),char(k),'4','1' /)
			iok=decchk(lin(1:7),ch)
			lin(8)=ch
			iok=decchk(lin(1:8),chh)
			lin(1:7)=(/ '8','6',char(k),char(l),char(j),'4','1' /)
			lin(8)=ch
			jok=decchk(lin(1:8),chh)
			if ((.not.iok) .or. jok) then
				nbad=nbad+1
			end if
		end if
	end do
	end do
	end do
	write(*,'(1x,a,t30,i3)') 'Total tries:',ntot
	write(*,'(1x,a,t30,i3)') 'Bad tries:',nbad
	write(*,'(1x,a,t29,f4.2)') &
		'Fraction good:',real(ntot-nbad,sp)/real(ntot,sp)
!	construct check digits for some user-supplied strings
	do
		write(*,*) 'enter string terminated by x:'
		read(*,'(a20)',END=999) line
		do j=1,128
			if (line(j:j) == 'x') exit
			lin(j)=line(j:j)
		end do
		n=j-1
		if (n == 0) goto 999
		iok=decchk(lin(1:n),ch)
		lin(n+1)=ch
		jok=decchk(lin(1:n+1),chh)
		tf=merge('T','F',jok)
		write(*,*) lin(1:n+1),' checks as ',tf
	end do
999 	write(*,*) 'NORMAL COMPLETION'
	END PROGRAM xdecchk
