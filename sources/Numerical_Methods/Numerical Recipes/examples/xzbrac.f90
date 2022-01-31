	PROGRAM xzbrac
!	driver for routine zbrac
	USE nrtype
	USE nr
	IMPLICIT NONE
	LOGICAL(LGT) :: succes
	INTEGER(I4B) :: i
	REAL(SP) :: x1,x2
	write(*,'(/1x,t4,a,t29,a/)') 'Bracketing values:','Function values:'
	write(*,'(1x,t6,a,t16,a,t29,a,t41,a/)') 'X1','X2',&
		'BESSJ0(X1)','BESSJ0(X2)'
	do i=1,10
		x1=i
		x2=x1+1.0_sp
		call zbrac(bessj0_s,x1,x2,succes)
		if (succes) then
			write(*,'(1x,f7.2,f10.2,7x,2f12.6)') x1,x2,&
				bessj0(x1),bessj0(x2)
		end if
	end do
	END PROGRAM xzbrac
