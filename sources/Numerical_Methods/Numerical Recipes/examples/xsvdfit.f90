	PROGRAM xsvdfit
!	driver for routine svdfit
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPT=100,NPOL=5
	REAL(SP), PARAMETER :: SPREAD=0.02_sp
	INTEGER(I4B) :: i
	REAL(SP) :: chisq
	REAL(SP), DIMENSION(NPOL) :: a,w
	REAL(SP), DIMENSION(NPT) :: x,y,sig
	REAL(SP), DIMENSION(NPOL,NPOL) :: cvm,v
!	polynomial fit
	call ran_seed(sequence=900)
	call gasdev(y(:))
	y(:)=1.0_sp+SPREAD*y(:)
	x(:)=arth(0.02_sp,0.02_sp,NPT)
	y(:)=y(:)*(1.0_sp+x(:)*(2.0_sp+x(:)*(3.0_sp+x(:)*(4.0_sp+x(:)*5.0_sp))))
	sig(:)=SPREAD*y(:)
	call svdfit(x,y,sig,a,v,w,chisq,fpoly)
	call svdvar(v,w,cvm)
	write(*,*) 'Polynomial fit:'
	do i=1,NPOL
		write(*,'(1x,f12.6,a,f10.6)') a(i),'  +-',sqrt(cvm(i,i))
	end do
	write(*,'(1x,a,f12.6/)') 'Chi-squared',chisq
	call svdfit(x,y,sig,a,v,w,chisq,fleg)
	call svdvar(v,w,cvm)
	write(*,*) 'Legendre polynomial fit'
	do i=1,NPOL
		write(*,'(1x,f12.6,a,f10.6)') a(i),'  +-',sqrt(cvm(i,i))
	end do
	write(*,'(1x,a,f12.6/)') 'Chi-squared',chisq
	END PROGRAM xsvdfit
