	PROGRAM xwt1
!	driver for routine wt1
	USE nrtype
	USE nr
	INTEGER(I4B), PARAMETER :: NCMAX=50,NMAX=512,NCEN=333,NWID=33
	INTEGER(I4B) :: i,itest,k,nused
	REAL(SP), DIMENSION(NMAX) :: u,v,w
	REAL(SP) :: frac,thresh,tmp
1	do
		write(*,*) 'Enter k (4, -4, 12, or 20) and frac (0. to 1.):'
		read(*,*,END=999) k,frac
		frac=min(1.0_sp,max(0.0_sp,frac))
		itest=merge(1,0,k == -4)
		k=abs(k)
		if (k == 4 .or. k == 12 .or. k == 20) exit
	end do
	v(:)=0.0
	do i=NCEN-NWID+1,NCEN+NWID-1
		v(i)=real(i-NCEN+NWID,sp)*real(NCEN+NWID-i,sp)/NWID**2
	end do
	w(:)=v(:)
	if (itest == 0) then
		call pwtset(k)
		call wt1(v,1,pwt)
	else
		call wt1(v,1,daub4)
	end if
	u(:)=abs(v(:))
	thresh=select(int((1.0-frac)*NMAX),u)
	where (abs(v(:)) <= thresh) v(:)=0.0
	nused=count(abs(v(:)) > thresh)
	if (itest == 0) then
		call wt1(v,-1,pwt)
	else
		call wt1(v,-1,daub4)
	end if
	thresh=0.0
	tmp=maxval(abs(v(:)-w(:)))
	if (tmp > thresh) thresh=tmp
	write(*,*) 'k,NMAX,nused=',k,NMAX,nused
	write(*,*) 'discrepancy=',thresh
	goto 1
999	END PROGRAM xwt1
