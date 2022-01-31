	PROGRAM xlfit
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
!	driver for routine lfit
	INTEGER(I4B), PARAMETER :: NPT=100,NTERM=5
	REAL(SP), PARAMETER :: SPREAD=0.1_sp
	INTEGER(I4B) :: i,j
	INTEGER(I4B), DIMENSION(NTERM) :: nums
	REAL(SP) :: chisq
	REAL(SP), DIMENSION(NPT) :: x,y,sig
	REAL(SP), DIMENSION(NTERM) :: a
	REAL(SP), DIMENSION(NTERM,NTERM) :: covar
	LOGICAL(LGT), DIMENSION(NTERM) :: mask
	INTERFACE
		SUBROUTINE funcs(x,arr)
		USE nrtype; USE nrutil
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP), DIMENSION(:), INTENT(OUT) :: arr
		END SUBROUTINE funcs
	END INTERFACE
	call ran_seed(sequence=1999)
	x(:)=arth(0.1_sp,0.1_sp,NPT)
	call gasdev(y(:))
	y(:)=SPREAD*y(:)
	nums=arth(1,1,NTERM)
	do i=1,NPT
		call funcs(x(i),a)
		y(i)=y(i)+dot_product(nums,a)
	end do
	sig(1:NPT)=SPREAD
	mask(:)=.true.
	call lfit(x,y,sig,a,mask,covar,chisq,funcs)
	write(*,'(/1x,t4,a,t22,a)') 'Parameter','Uncertainty'
	do i=1,NTERM
		write(*,'(1x,t5,a,i1,a,f8.6,f11.6)') 'A(',i,') = ',&
			a(i),sqrt(covar(i,i))
	end do
	write(*,'(/3x,a,e12.6)') 'Chi-squared = ',chisq
	write(*,'(/3x,a)') 'Full covariance matrix'
	do i=1,NTERM
		write(*,'(1x,6e12.2)') (covar(i,j),j=1,NTERM)
	end do
	write(*,'(/1x,a)') 'press RETURN to continue...'
	read(*,*)
!	now check results of restricting fit parameters
	mask(2:NTERM:2)=.true.
	call lfit(x,y,sig,a,mask,covar,chisq,funcs)
	write(*,'(/1x,t4,a,t22,a)') 'Parameter','Uncertainty'
	do i=1,NTERM
		write(*,'(1x,t5,a,i1,a,f8.6,f11.6)') 'A(',i,') = ',&
			a(i),sqrt(covar(i,i))
	end do
	write(*,'(/3x,a,e12.6)') 'Chi-squared = ',chisq
	write(*,'(/3x,a)') 'Full covariance matrix'
	do i=1,NTERM
		write(*,'(1x,6e12.2)') (covar(i,j),j=1,NTERM)
	end do
	END PROGRAM xlfit

	SUBROUTINE funcs(x,afunc)
	USE nrtype; USE nrutil
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i
	REAL(SP), INTENT(IN) :: x
	REAL(SP), DIMENSION(:), INTENT(OUT) :: afunc
	afunc(1)=1.0
	afunc(2)=x
	do i=3,size(afunc)
		afunc(i)=sin(i*x)
	end do
	END SUBROUTINE funcs
