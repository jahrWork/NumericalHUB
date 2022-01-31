	PROGRAM xsor
!	driver for routine sor
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: JMAX=33,NSTEP=4
	INTEGER(I4B) :: i,j,midl
	REAL(DP) :: rjac
	REAL(DP), DIMENSION(JMAX,JMAX) :: a,b,c,d,e,f,u
	a(:,:)=1.0
	b(:,:)=1.0
	c(:,:)=1.0
	d(:,:)=1.0
	e(:,:)=-4.0
	f(:,:)=0.0
	u(:,:)=0.0
	midl=JMAX/2+1
	f(midl,midl)=2.0_dp/(JMAX-1)**2
	rjac=cos(PI/JMAX)
	call sor(a,b,c,d,e,f,u,rjac)
	write(*,'(1x,a)') 'SOR Solution:'
	do i=1,JMAX,NSTEP
		write(*,'(1x,9f8.4)') (u(i,j),j=1,JMAX,NSTEP)
	end do
	write(*,'(/1x,a)') 'Test that solution satisfies Difference Eqns:'
	do i=NSTEP+1,JMAX-1,NSTEP
		f(i,NSTEP+1:JMAX-1:NSTEP)=&
			u(i+1,NSTEP+1:JMAX-1:NSTEP)+u(i-1,NSTEP+1:JMAX-1:NSTEP)+&
			u(i,NSTEP+2:JMAX:NSTEP)+u(i,NSTEP:JMAX-2:NSTEP)-&
			4.0_dp*u(i,NSTEP+1:JMAX-1:NSTEP)
		write(*,'(7x,9f8.4)') (f(i,j),j=NSTEP+1,JMAX-1,NSTEP)
	end do
	END PROGRAM xsor
