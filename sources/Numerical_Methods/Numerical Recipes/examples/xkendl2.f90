	PROGRAM xkendl2
!	driver for routine kendl2
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NDAT=1000,IP=8,JP=8
	INTEGER(I4B) :: i,ifunc,iseed,j,k,l,m,n
	REAL(SP) :: prob,tau,z
	REAL(SP), DIMENSION(IP,JP) :: tab
	CHARACTER(3), DIMENSION(8) :: &
		text = (/ '000','001','010','011','100','101','110','111' /)
	write(*,*) 'Are ones followed by zeros and vice-versa?'
	i=IP
	j=JP
	do ifunc=1,2
		iseed=2468
		write(*,'(/1x,a,i1/)') 'Test of IRBIT',ifunc
		tab(1:i,1:j)=0.0
		do m=1,NDAT
			k=1
			do n=0,2
				if (ifunc == 1) then
					k=k+irbit1(iseed)*(2**n)
				else
					k=k+irbit2(iseed)*(2**n)
				end if
			end do
			l=1
			do n=0,2
				if (ifunc == 1) then
					l=l+irbit1(iseed)*(2**n)
				else
					l=l+irbit2(iseed)*(2**n)
				end if
			end do
			tab(k,l)=tab(k,l)+1.0_sp
		end do
		call kendl2(tab,tau,z,prob)
		write(*,'(4x,8a6/)') (text(n),n=1,8)
		do n=1,8
			write(*,'(1x,a,8i6)') text(n),(nint(tab(n,m)),m=1,8)
		end do
		write(*,'(/7x,a,t24,a,t38,a)') 'Kendall Tau','Std. Dev.',&
			'Probability'
		write(*,'(1x,3f15.6/)') tau,z,prob
		write(*,*) 'Press RETURN to continue ...'
		read(*,*)
	end do
	END PROGRAM xkendl2
