	PROGRAM xspear
!	driver for routine spear
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDAT=20,NMON=12
	INTEGER(I4B) :: i,j
	REAL(SP) :: d,probd,probrs,rs,zd
	REAL(SP), DIMENSION(NDAT) :: data1,data2
	REAL(SP), DIMENSION(NDAT,NMON) :: rays
	REAL(SP), DIMENSION(NDAT) :: ave,zlat
	CHARACTER(64) :: text
	CHARACTER(15), DIMENSION(NDAT) :: city
	CHARACTER(4), DIMENSION(NMON) :: mon
	open(7,file='TABLE2.DAT',status='OLD')
	read(7,*)
	read(7,'(a)') text
	read(7,'(15x,12a4/)') (mon(i),i=1,12)
	do i=1,NDAT
		read(7,'(a15,12f4.0,f6.0,f6.1)') &
			city(i),(rays(i,j),j=1,12),ave(i),zlat(i)
	end do
	close(7)
	write(*,*) text
	write(*,'(1x,15x,12a4)') (mon(i),i=1,12)
	do i=1,NDAT
		write(*,'(1x,a,12i4,i6,f6.1)') city(i),(nint(rays(i,j)),j=1,12)
	end do
!	check temperature correlations between different months
	write(*,'(/1x,a)')&
		'Are sunny summer places also sunny winter places?'
	write(*,'(1x,2a)') 'Check correlation of sampled U.S. solar ',&
		'radiation (july with other months)'
	write(*,'(/1x,a,t16,a,t23,a,t37,a,t49,a,t63,a/)')&
		'Month','D','St. Dev.','PROBD',&
		'Spearman R','PROBRS'
	data1(:)=rays(:,1)
	do j=1,12
		data2(:)=rays(:,j)
		call spear(data1,data2,d,zd,probd,rs,probrs)
		write(*,'(1x,a,f13.2,2f12.6,3x,2f12.6)')&
			mon(j),d,zd,probd,rs,probrs
	end do
	END PROGRAM xspear
