	PROGRAM xfit
!	driver for routine fit
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPT=100
	REAL(SP), PARAMETER :: SPREAD=0.5_sp
	INTEGER(I4B) :: mwt
	REAL(SP) :: a,b,chi2,q,siga,sigb
	REAL(SP), DIMENSION(NPT) :: harvest,sig,x,y
	call ran_seed(sequence=731)
	x(:)=arth(0.1_sp,0.1_sp,NPT)
	call gasdev(harvest)
	y(:)=-2.0_sp*x(:)+1.0_sp+SPREAD*harvest
	sig(:)=SPREAD
	do mwt=0,1
		if (mwt == 0) then
			write(*,'(//1x,a)') 'Ignoring standard deviation'
			call fit(x,y,a,b,siga,sigb,chi2,q)
		else
			write(*,'(//1x,a)') 'Including standard deviation'
			call fit(x,y,a,b,siga,sigb,chi2,q,sig)
		end if
		write(*,'(1x,t5,a,f9.6,t24,a,f9.6)') 'A = ',a,'Uncertainty: ',&
			siga
		write(*,'(1x,t5,a,f9.6,t24,a,f9.6)') 'B = ',b,'Uncertainty: ',&
			sigb
		write(*,'(1x,t5,a,4x,f10.6)') 'Chi-squared: ',chi2
		write(*,'(1x,t5,a,f10.6)') 'Goodness-of-fit: ',q
	end do
	END PROGRAM xfit
