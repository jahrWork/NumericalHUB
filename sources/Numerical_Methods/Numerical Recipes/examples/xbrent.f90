	PROGRAM xbrent
!	driver for routine brent
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP), PARAMETER :: TOL=1.0e-6_sp,EQL=1.0e-4_sp
	INTEGER(I4B) :: i,nmin
	REAL(SP) :: ax,b,bx,cx,fa,fb,fc,xmin
	REAL(SP), DIMENSION(20) :: amin
	INTERFACE
		FUNCTION func(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	nmin=0
	write(*,'(/1x,a)') 'Minima of the function BESSJ0'
	write(*,'(/1x,t6,a,t19,a,t28,a,t40,a/)') 'Min. #','X',&
		'BESSJ0(X)','BESSJ1(X)'
	do i=1,100
		ax=i
		bx=i+1.0_sp
		call mnbrak(ax,bx,cx,fa,fb,fc,func)
		b=brent(ax,bx,cx,func,TOL,xmin)
		if (nmin == 0) then
			amin(1)=xmin
			nmin=1
			write(*,'(1x,5x,i2,3x,3f12.6)') nmin,xmin,&
				bessj0(xmin),bessj1(xmin)
		else
			if (any(abs(xmin-amin(1:nmin)) <= EQL*xmin)) cycle
			nmin=nmin+1
			amin(nmin)=xmin
			write(*,'(1x,5x,i2,3x,3f12.6)') nmin,&
				xmin,bessj0(xmin),bessj1(xmin)
		end if
	end do
	END PROGRAM xbrent

	FUNCTION func(x)
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP), INTENT(IN) :: x
	REAL(SP) :: func
	func=bessj0(x)
	END FUNCTION func
