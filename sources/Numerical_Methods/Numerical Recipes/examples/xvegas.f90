MODULE vegas_fdata
	USE nrtype
	REAL(SP) :: xoff
END MODULE vegas_fdata

	PROGRAM xvegas
	USE nrtype
	USE nr
	USE vegas_fdata
	IMPLICIT NONE
!	driver for routine vegas
	INTEGER(I4B) :: init,itmax,ncall,ndim,nprn
	REAL(SP) :: avgi,chi2a,sd
	REAL(SP), DIMENSION(20) :: region
	INTERFACE
		FUNCTION fxn(pt,wgt)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: pt
		REAL(SP), INTENT(IN) :: wgt
		REAL(SP) :: fxn
		END FUNCTION fxn
	END INTERFACE
	write(*,*) 'ENTER NDIM,XOFF,NCALL,ITMAX,NPRN'
	read(*,*,END=999) ndim,xoff,ncall,itmax,nprn
	avgi=0.0
	sd=0.0
	chi2a=0.0
	region(1:ndim)=0.0
	region(1+ndim:2*ndim)=1.0
	init = -1
	call vegas(region(1:2*ndim),fxn,init,ncall,itmax,nprn,avgi,sd,chi2a)
	write(*,*) 'Number of iterations performed:',itmax
	write(*,*) 'Integral, Standard Dev., Chi-sq.',avgi,sd,chi2a
	init = 1
	call vegas(region(1:2*ndim),fxn,init,ncall,itmax,nprn,avgi,sd,chi2a)
	write(*,*) 'Additional iterations performed:',itmax
	write(*,*) 'Integral, Standard Dev., Chi-sq.',avgi,sd,chi2a
999	write(*,*) 'NORMAL COMPLETION'
	END PROGRAM xvegas

	FUNCTION fxn(pt,wgt)
	USE nrtype
	USE vegas_fdata
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: pt
	REAL(SP), INTENT(IN) :: wgt
	REAL(SP) :: fxn
	REAL(SP) :: sm
	sm=100.0_sp*sum((pt(:)-xoff)*(pt(:)-xoff))
	if (sm < 80.0) then
		fxn=exp(-sm)
	else
		fxn=0.0
	end if
	fxn=fxn*(5.64189_sp**size(pt))
	END FUNCTION fxn
