	PROGRAM xmedfit
!	driver for routine medfit
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPT=100
	REAL(SP), PARAMETER :: SPREAD=0.1_sp
	REAL(SP) :: a,abdev,b,chi2,q,siga,sigb
	REAL(SP), DIMENSION(NPT) :: x,y,sig
	call ran_seed(sequence=1984)
	x(:)=arth(0.1_sp,0.1_sp,NPT)
	call gasdev(y(:))
	y(:)=SPREAD*y(:)
	y(:)=y(:)-2.0_sp*x(:)+1.0_sp
	sig(:)=SPREAD
	call fit(x,y,a,b,siga,sigb,chi2,q,sig)
	write(*,'(/1x,a)') 'According to routine FIT the result is:'
	write(*,'(1x,t5,a,f8.4,t20,a,f8.4)') 'A = ',a,'Uncertainty: ',siga
	write(*,'(1x,t5,a,f8.4,t20,a,f8.4)') 'B = ',b,'Uncertainty: ',sigb
	write(*,'(1x,t5,a,f8.4,a,i4,a)') 'Chi-squared: ',chi2,&
		' for ',NPT,' points'
	write(*,'(1x,t5,a,f8.4)') 'Goodness-of-fit: ',q
	write(*,'(/1x,a)') 'According to routine MEDFIT the result is:'
	call medfit(x,y,a,b,abdev)
	write(*,'(1x,t5,a,f8.4)') 'A = ',a
	write(*,'(1x,t5,a,f8.4)') 'B = ',b
	write(*,'(1x,t5,a,f8.4)') 'Absolute deviation (per DATA point): ',abdev
	write(*,'(1x,t5,a,f8.4,a)') '(note: Gaussian SPREAD is',SPREAD,')'
	END PROGRAM xmedfit
