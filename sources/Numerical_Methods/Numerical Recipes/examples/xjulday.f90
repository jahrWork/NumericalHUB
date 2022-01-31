	PROGRAM xjulday
!	driver for julday
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i,im,id,iy,n
	CHARACTER(40) :: txt
	CHARACTER(10), DIMENSION(12) :: month = (/ &
		'January  ','February ','March    ','April    ',&
		'May      ','June     ','July     ','August   ',&
		'September','October  ','November ','December ' /)
	open(7,file='DATES.DAT',status='OLD')
	read(7,'(a)') txt
	read(7,*) n
	write(*,'(/1x,a,t12,a,t17,a,t23,a,t37,a/)') 'Month','Day','Year',&
		'Julian Day','Event'
	do i=1,n
		read(7,'(i2,i3,i5,a)') im,id,iy,txt
		write(*,'(1x,a10,i3,i6,3x,i7,5x,a)') month(im),id,iy,&
			julday(im,id,iy),txt
	end do
	close(7)
	write(*,'(/1x,a/)') 'Month,Day,Year (e.g. 1,13,1905)'
	do i=1,20
		write(*,*) 'MM,DD,YYYY'
		read(*,*) im,id,iy
		if (im < 0) exit
		write(*,'(1x,a12,i8/)') 'Julian Day: ',julday(im,id,iy)
	end do
	END PROGRAM xjulday
