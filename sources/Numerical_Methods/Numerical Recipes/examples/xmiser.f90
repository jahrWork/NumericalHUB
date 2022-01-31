	MODULE xmiser_func_data
	USE nrtype
	REAL(SP) :: xoff
	END MODULE xmiser_func_data

	PROGRAM xmiser
!	driver for routine miser
	USE nrtype
	USE nr ! for miser interface
	USE xmiser_func_data ! to set xoff
	IMPLICIT NONE
	INTEGER(I4B) :: idum,n,ndim,nt,ntries
	REAL(SP) :: ave,dith,sumav,sumsd,var
	REAL(SP), DIMENSION(:), ALLOCATABLE :: region
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		REAL(SP) :: func
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		END FUNCTION func
	END INTERFACE
	idum=0
	write(*,*) 'IDUM='
	read(*,*) idum
	idum=-abs(idum)
	do
		write(*,*) 'ENTER N,NDIM,XOFF,DITH,NTRIES'
		read(*,*,END=999) n,ndim,xoff,dith,ntries
		allocate(region(2*ndim))
		sumav=0.0
		sumsd=0.0
		do nt=1,ntries
			region(1:ndim)=0.0
			region((ndim+1):(2*ndim))=1.0
			call miser(func,region,ndim,n,dith,ave,var)
			sumav=sumav+(ave-1.0_sp)**2
			sumsd=sumsd+sqrt(abs(var))
		end do
		sumav=sqrt(sumav/ntries)
		sumsd=sumsd/ntries
		deallocate(region)
		write(*,*) 'fractional error: actual,indicated=',sumav,sumsd
	end do
999	write(*,*) 'NORMAL COMPLETION'
	END PROGRAM xmiser

	FUNCTION func(pt)
	USE nrtype
	USE xmiser_func_data ! for xoff
	IMPLICIT NONE
	INTEGER(I4B) :: ndim
	REAL(SP) :: func,summ
	REAL(SP), DIMENSION(:), INTENT(IN) :: pt
	ndim=size(pt)
	summ=sum(100.0_sp*(pt-xoff)**2)
	if (summ < 80.0) then
		func=exp(-summ)
	else
		func=0.0
	end if
	func=func*(5.64189_sp**ndim)
	END FUNCTION func
