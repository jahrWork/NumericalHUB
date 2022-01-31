	PROGRAM xpearsn
!	driver for routine pearsn
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: i
	REAL(SP) :: prob,r,z
	REAL(SP), DIMENSION(10) :: dose = (/ &
		56.1_sp,64.1_sp,70.0_sp,66.6_sp,82.0_sp,&
		91.3_sp,90.0_sp,99.7_sp,115.3_sp,110.0_sp /),&
		spore = (/ &
		0.11_sp,0.4_sp,0.37_sp,0.48_sp,0.75_sp,&
		0.66_sp,0.71_sp,1.2_sp,1.01_sp,0.95_sp /)
	write(*,'(1x,a)') 'Effect of Gamma Rays on Man-in-the-Moon Marigolds'
	write(*,'(1x,a,t29,a)') 'Count Rate (cpm)','Pollen Index'
	do i=1,10
		write(*,'(1x,f10.2,f25.2)') dose(i),spore(i)
	end do
	call pearsn(dose,spore,r,prob,z)
	write(*,'(/1x,t24,a,t38,a)') 'PEARSN','Expected'
	write(*,'(1x,a,t18,2e15.6)') 'Corr. Coeff.',r,0.906959_sp
	write(*,'(1x,a,t18,2e15.6)') 'Probability',prob,0.292650e-3_sp
	write(*,'(1x,a,t18,2e15.6/)') 'Fisher''s Z',z,1.51011_sp
	END PROGRAM xpearsn
