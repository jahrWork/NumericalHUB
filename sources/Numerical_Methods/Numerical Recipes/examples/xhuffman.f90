	PROGRAM xhuffman
!	driver for routines hufmak, hufenc, hufdec
	USE nrtype; USE nrutil
	USE nr
	USE huf_info
	IMPLICIT NONE
	INTERFACE
		SUBROUTINE hufmak(nfreq,ilong,nlong,hcode)
		USE nrtype; USE nrutil
		USE huf_info
		IMPLICIT NONE
		INTEGER(I4B), INTENT(OUT) :: ilong,nlong
		INTEGER(I4B), DIMENSION(:), INTENT(IN) :: nfreq
		TYPE(huffcode) :: hcode
		END SUBROUTINE hufmak
!BL
		SUBROUTINE hufdec(ich,code,nb,hcode)
		USE nrtype; USE nrutil
		USE huf_info
		IMPLICIT NONE
		INTEGER(I4B), INTENT(OUT) :: ich
		INTEGER(I4B), INTENT(INOUT) :: nb
		CHARACTER(1), DIMENSION(:), INTENT(IN) :: code
		TYPE(huffcode) :: hcode
		END SUBROUTINE hufdec
!BL
		SUBROUTINE hufenc(ich,codep,nb,hcode)
		USE nrtype; USE nrutil
		USE huf_info
		IMPLICIT NONE
		INTEGER(I4B), INTENT(IN) :: ich
		INTEGER(I4B), INTENT(INOUT) :: nb
		CHARACTER(1), DIMENSION(:), POINTER :: codep
		TYPE(huffcode) :: hcode
		END SUBROUTINE hufenc
	END INTERFACE
	INTEGER(I4B), PARAMETER :: MAXBUF=200,MAXLINE=80
	INTEGER(I4B) :: i,ilong,j,k,n,nb,nh,nlong,nt,nch
	INTEGER(I4B), DIMENSION(256) :: nfreq
	CHARACTER(1), DIMENSION(:), POINTER :: codep
	CHARACTER(200) :: mess,ness
	CHARACTER(80) :: lin
	TYPE(HUFFCODE) :: hcode
	nullify(codep)
	codep=>reallocate(codep,256)
!	construct a letter frequency table from the file TEXT.DAT
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
!	here is the initialization that constructs the code
	call hufmak(nfreq(1:nch),ilong,nlong,hcode)
	write(*,*) &
	'index, ','character, ','nfreq, ','bits in code, ','code int'
	do j=1,nch
		if (nfreq(j) /= 0.0) write(*,*) j,' ',char(j+31),&
			' ',nfreq(j),hcode%ncode(j),hcode%icode(j)
	end do
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
!	here we Huffman encode mess(1:n)
	nb=0
	do j=1,n
		call hufenc(ichar(mess(j:j)),codep,nb,hcode)
	end do
	nh=nb/8+1
!	message termination (encode a single long character)
	call hufenc(ilong,codep,nb,hcode)
!	here we decode the message, hopefully to get the original back
	nb=0
	do j=1,MAXLINE
		call hufdec(i,codep(1:nh),nb,hcode)
		if (i == nch) exit
		ness(j:j)=char(i)
	end do
	if (j > MAXLINE) call nrerror('HUFFMAN - NEVER GET HERE')
	nt=j-1
	write(*,*) 'LENGTH OF LINE INPUT,CODED=',n,nh
	write(*,*) 'DECODED OUTPUT:'
	write(*,'(1x,80a1)') (char(ichar(ness(k:k))+32),k=1,nt)
	if (nt /= n) write(*,*) 'ERROR ! N DECODED  /=  N INPUT'
	if (nt-n == 1) write(*,*) 'MAY BE HARMLESS SPURIOUS CHARACTER.'
	end do
999	write(*,*) 'NORMAL COMPLETION'
	deallocate(codep)
	call huff_deallocate(hcode)
	END PROGRAM xhuffman
