	PROGRAM xcaldat
!	driver for routine caldat
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i,im,imm,id,idd,iy,iyy,iycopy,j,n
	CHARACTER(10), DIMENSION(12) :: month = (/ &
		'January  ','February ','March    ','April    ',&
		'May      ','June     ','July     ','August   ',&
		'September','October  ','November ','December ' /)
!	check whether CALDAT properly undoes the operation of JULDAY
	open(7,file='DATES.DAT',status='OLD')
	read(7,*)
	read(7,*) n
	write(*,'(/1x,a,t40,a)') 'Original Date:','Reconstructed Date:'
	write(*,'(1x,a,t12,a,t17,a,t25,a,t40,a,t50,a,t55,a/)')&
		'Month','Day','Year','Julian Day','Month','Day','Year'
	do i=1,n
		read(7,'(i2,i3,i5)') im,id,iy
		iycopy=iy
		j=julday(im,id,iycopy)
		call caldat(j,imm,idd,iyy)
		write(*,'(1x,a,i3,i6,4x,i9,6x,a,i3,i6)') month(im),id,&
			iy,j,month(imm),idd,iyy
	end do
	END PROGRAM xcaldat
