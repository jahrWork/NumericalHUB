	PROGRAM xflmoon
!	driver for routine flmoon
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP), PARAMETER :: TZONE=-5.0
	REAL(SP) :: frac,timzon
	INTEGER(I4B) :: i,im,id,iy,ifrac,istr,j1,j2,n,nph
	CHARACTER(3), DIMENSION(2) :: timstr = (/ ' AM',' PM' /)
	CHARACTER(15), DIMENSION(4) :: 	phase = (/ &
		'new moon     ','first quarter','full moon    ','last quarter ' /)
	write(*,*) 'Date of the next few phases of the moon'
	write(*,*) 'Enter today''s date (e.g. 12,15,1992)'
	timzon=TZONE/24.0_sp
	read(*,*) im,id,iy
!	approximate number of full moons since January 1900
	n=12.37_sp*(iy-1900+(im-0.5_sp)/12.0_sp)
	nph=2
	j1=julday(im,id,iy)
	call flmoon(n,nph,j2,frac)
	n=n+nint((j1-j2)/29.53_sp)
	write(*,'(/1x,t6,a,t19,a,t32,a)') 'Date','Time(EST)','Phase'
	do i=1,20
		call flmoon(n,nph,j2,frac)
		ifrac=nint(24.0_sp*(frac+timzon))
		if (ifrac < 0) then
			j2=j2-1
			ifrac=ifrac+24
		end if
		if (ifrac >= 12) then
			j2=j2+1
			ifrac=ifrac-12
		else
			ifrac=ifrac+12
		end if
		if (ifrac > 12) then
			ifrac=ifrac-12
			istr=2
		else
			istr=1
		end if
		call caldat(j2,im,id,iy)
		write(*,'(1x,2i3,i5,t20,i2,a,5x,a)') im,id,iy,&
			ifrac,timstr(istr),phase(nph+1)
		if (nph == 3) then
			nph=0
			n=n+1
		else
			nph=nph+1
		end if
	end do
	END PROGRAM xflmoon
