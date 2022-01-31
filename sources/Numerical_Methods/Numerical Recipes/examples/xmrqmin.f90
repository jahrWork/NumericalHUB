	PROGRAM xmrqmin
!	driver for routine mrqmin
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPT=100,MA=6
	REAL(SP), PARAMETER :: SPREAD=0.001_sp
	INTEGER(I4B) :: i,iter,itst,k,mfit
	REAL(SP) :: alamda,chisq,ochisq
	REAL(SP), DIMENSION(MA) :: &
		a = (/ 5.0_sp,2.0_sp,3.0_sp,2.0_sp,5.0_sp,3.0_sp /),&
		gues = (/ 4.5_sp,2.2_sp,2.8_sp,2.5_sp,4.9_sp,2.8_sp /)
	REAL(SP), DIMENSION(MA,MA) :: covar,alpha
	REAL(SP), DIMENSION(NPT) :: x,y,sig
	LOGICAL(LGT), DIMENSION(MA) :: mask
!	first try a sum of two Gaussians
	call ran_seed(sequence=500)
	call gasdev(y(:))
	y(:)=1.0_sp+SPREAD*y(:)
	x(:)=arth(0.1_sp,0.1_sp,NPT)
	do i=1,NPT
		y(i)=y(i)*dot_product(a(1:6:3),exp(-((x(i)-a(2:6:3))/a(3:6:3))**2))
	end do
	sig(:)=SPREAD*y(:)
	mfit=MA
	mask(1:mfit)=.true.
	a(1:mfit)=gues(1:mfit)
	do iter=1,2
		alamda=-1
		call mrqmin(x,y,sig,a,mask,covar,alpha,chisq,fgauss,alamda)
		k=1
		itst=0
		do
			write(*,'(/1x,a,i2,t18,a,f12.4,t45,a,e9.2)') &
				'Iteration #',k,'Chi-squared:',chisq,'ALAMDA:',alamda
			write(*,'(1x,t5,a,t13,a,t21,a,t29,a,t37,a,t45,a)') &
				'A(1)','A(2)','A(3)','A(4)','A(5)','A(6)'
			write(*,'(1x,6f8.4)') (a(i),i=1,6)
			k=k+1
			ochisq=chisq
			call mrqmin(x,y,sig,a,mask,covar,alpha,chisq,fgauss,alamda)
			if (chisq > ochisq) then
				itst=0
			else if (abs(ochisq-chisq) < 0.1_sp) then
				itst=itst+1
			end if
			if (itst >= 4) exit
		end do
		alamda=0.0
		call mrqmin(x,y,sig,a,mask,covar,alpha,chisq,fgauss,alamda)
		write(*,*) 'Uncertainties:'
		write(*,'(1x,6f8.4/)') (sqrt(covar(i,i)),i=1,6)
		write(*,'(1x,a)') 'Expected results:'
		write(*,'(1x,f7.2,5f8.2/)') 5.0,2.0,3.0,2.0,5.0,3.0
		if (iter == 1) then
			write(*,*) 'press return to continue with constraint'
			read(*,*)
			write(*,*) 'Holding a(2) and a(5) constant'
			a(1:MA)=a(1:MA)+0.1_sp
			a(2)=2.0
			mask(2)=.false.
			a(5)=5.0
			mask(5)=.false.
		end if
	end do
	END PROGRAM xmrqmin
