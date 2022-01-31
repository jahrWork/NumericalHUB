	PROGRAM xmgfas
!	driver for routine mgfas
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: JMAX=33,NSTEP=4
	INTEGER(I4B) :: i,j,midl
	REAL(DP), DIMENSION(JMAX,JMAX) :: f,u
	u(:,:)=0.0
	midl=JMAX/2+1
	u(midl,midl)=2.0
	call mgfas(u,2)
	write(*,'(1x,a)') 'MGFAS Solution:'
	do i=1,JMAX,NSTEP
		write(*,'(1x,9f8.4)') (u(i,j),j=1,JMAX,NSTEP)
	end do
	write(*,'(/1x,a)') 'Test that solution satisfies Difference Eqns:'
	do i=NSTEP+1,JMAX-1,NSTEP
		f(i,NSTEP+1:JMAX-1:NSTEP)=&
			u(i+1,NSTEP+1:JMAX-1:NSTEP)+u(i-1,NSTEP+1:JMAX-1:NSTEP)+&
			u(i,NSTEP+2:JMAX:NSTEP)+u(i,NSTEP:JMAX-2:NSTEP)-&
			4.0_dp*u(i,NSTEP+1:JMAX-1:NSTEP)+&
			u(i,NSTEP+1:JMAX-1:NSTEP)**2/(JMAX-1)**2
		write(*,'(7x,7f8.4)')&
			(f(i,j)*(JMAX-1)*(JMAX-1),j=NSTEP+1,JMAX-1,NSTEP)
	end do
	END PROGRAM xmgfas
