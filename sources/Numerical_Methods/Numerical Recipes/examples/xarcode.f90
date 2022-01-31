	PROGRAM xarcode
!	driver for routines arcmak, arcode
	USE nrtype; USE nrutil
	USE nr
	USE arcode_info
	IMPLICIT NONE
	INTERFACE
		SUBROUTINE arcmak(nfreq,nradd,acode)
		USE nrtype; USE nrutil
		USE arcode_info
		IMPLICIT NONE
		INTEGER(I4B), INTENT(IN) :: nradd
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nfreq
		TYPE(ARITHCODE) :: acode
		END SUBROUTINE arcmak
!BL
		SUBROUTINE arcode(ich,codep,lcd,isign,acode)
		USE nrtype; USE nrutil
		USE arcode_info
		IMPLICIT NONE
		INTEGER(I4B), INTENT(INOUT) :: ich,lcd
		INTEGER(I4B), INTENT(IN) :: isign
		CHARACTER(1), DIMENSION(:), POINTER :: codep
		TYPE(ARITHCODE) :: acode
		END SUBROUTINE arcode
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXBUF=200,MAXLINE=80
	INTEGER(I4B) :: i,j,k,lc,n,nt,nch,nrad,ichvar
	INTEGER(I4B), DIMENSION(256) ::nfreq
	CHARACTER(1), DIMENSION(:), POINTER :: codep
	CHARACTER(80) :: lin
	CHARACTER(200) :: mess,ness
	TYPE(ARITHCODE) :: acode
	nullify(codep)
	codep=>reallocate(codep,256)
	open(unit=7,file='TEXT.DAT',status='old')
	nfreq(:)=0
	do
		do j=1,MAXLINE
			lin(j:j)=char(32)
		end do
		read(7,'(a)',END=1) lin
		do n=MAXLINE,1,-1
			if (lin(n:n) /= char(32)) exit
		end do
		do k=1,min(MAXLINE,n)
			j=ichar(lin(k:k))-31
			if (j >= 1) nfreq(j)=nfreq(j)+1
		end do
	end do
1	close(unit=7)
	nch=96
	nrad=256
!	here is the initialization that constructs the code
	call arcmak(nfreq(1:nch),nrad,acode)
!	now ready to prompt for lines to encode
	do
	write(*,*) 'ENTER A LINE:'
	do j=1,MAXLINE
		mess(j:j)=char(32)
	end do
	read(*,'(a)',END=999) mess
	do n=MAXLINE,1,-1
		if (mess(n:n) /= char(32)) exit
	end do
!	shift from 256 character alphabet to 96 printing characters
	do j=1,n
		mess(j:j)=char(ichar(mess(j:j))-32)
	end do
!	message initialization
	lc=1
	ichvar=0
!	need this to be a variable to satisfy compiler checking
	call arcode(ichvar,codep,lc,0,acode)
!	here we arithmetically encode mess(1:n)
	do j=1,n
		ichvar=ichar(mess(j:j))
		call arcode(ichvar,codep,lc,1,acode)
	end do
		call arcode(nch,codep,lc,1,acode)
!	message termination
	write(*,*) 'LENGTH OF LINE INPUT, CODED=',n,lc-1
!	here we decode the message, hopefully to get the original back
	lc=1
	ichvar=0
	call arcode(ichvar,codep,lc,0,acode)
	do j=1,MAXBUF
		call arcode(i,codep,lc,-1,acode)
		if (i == nch) exit
		ness(j:j)=char(i)
	end do
	if (j > MAXBUF) call nrerror('xarcode: length sanity check failed')
	nt=j-1
	write(*,*) 'DECODED OUTPUT:'
	write(*,'(1x,80a1)') (char(ichar(ness(j:j))+32),j=1,nt)
	if (nt /= n) write(*,*) 'ERROR ! J DECODED  /=  N INPUT',j,n
	end do
999	write(*,*) 'NORMAL COMPLETION'
	deallocate(codep)
	call arcode_deallocate(acode)
	END PROGRAM xarcode
