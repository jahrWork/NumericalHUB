	PROGRAM xdbrent
!	driver for routine dbrent
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTERFACE
		FUNCTION deriv(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: deriv
		END FUNCTION deriv
	END INTERFACE
	REAL(SP), PARAMETER :: TOL=1.0e-6_sp,EQL=1.0e-4_sp
	INTEGER(I4B) :: i,iflag,nmin
	REAL(SP) :: ax,bx,cx,dbr,fa,fb,fc,xmin
	REAL(SP), DIMENSION(20) :: amin
	nmin=0
	write(*,'(/1x,a)') 'Minima of the function BESSJ0'
	write(*,'(/1x,t6,a,t19,a,t27,a,t39,a,t53,a/)') 'Min. #','X',&
		'BESSJ0(X)','BESSJ1(X)','DBRENT'
	do i=1,100
		ax=i
		bx=i+1.0_sp
		call mnbrak(ax,bx,cx,fa,fb,fc,bessj0_s)
		dbr=dbrent(ax,bx,cx,bessj0_s,deriv,TOL,xmin)
		if (nmin == 0) then
			amin(1)=xmin
			nmin=1
			write(*,'(1x,5x,i2,3x,4f12.6)') nmin,xmin,&
				bessj0(xmin),deriv(xmin),dbr
		else
			iflag=0
			if (any(abs(xmin-amin(1:nmin)) <= EQL*xmin)) iflag=1
			if (iflag == 0) then
				nmin=nmin+1
				amin(nmin)=xmin
				write(*,'(1x,5x,i2,3x,4f12.6)') nmin,xmin,&
					bessj0(xmin),deriv(xmin),dbr
			end if
		end if
	end do
	END PROGRAM xdbrent

	FUNCTION deriv(x)
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: deriv
	deriv=-bessj1(x)
	END FUNCTION deriv
