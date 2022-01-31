	SUBROUTINE caldat(julian,mm,id,iyyy)
	USE nrtype
	IMPLICIT NONE
	INTEGER(I4B), INTENT(IN) :: julian
	INTEGER(I4B), INTENT(OUT) :: mm,id,iyyy
	INTEGER(I4B) :: ja,jalpha,jb,jc,jd,je
	INTEGER(I4B), PARAMETER :: IGREG=2299161
	if (julian >= IGREG) then
		jalpha=int(((julian-1867216)-0.25_dp)/36524.25_dp)
		ja=julian+1+jalpha-int(0.25_dp*jalpha)
	else if (julian < 0) then
		ja=julian+36525*(1-julian/36525)
	else
		ja=julian
	end if
	jb=ja+1524
	jc=int(6680.0_dp+((jb-2439870)-122.1_dp)/365.25_dp)
	jd=365*jc+int(0.25_dp*jc)
	je=int((jb-jd)/30.6001_dp)
	id=jb-jd-int(30.6001_dp*je)
	mm=je-1
	if (mm > 12) mm=mm-12
	iyyy=jc-4715
	if (mm > 2) iyyy=iyyy-1
	if (iyyy <= 0) iyyy=iyyy-1
	if (julian < 0) iyyy=iyyy-100*(1-julian/36525)
	END SUBROUTINE caldat
