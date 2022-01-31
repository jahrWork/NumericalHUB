	PROGRAM xsavgol
!	driver for routine savgol
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NMAX=1000,NTEST=6
	INTEGER(I4B) :: i,m,nl,np,nrr
	INTEGER(I4B), DIMENSION(NTEST) :: nltest = (/ 2,3,4,5,4,5 /),&
		nrtest = (/ 2,1,0,5,4,5 /),mtest = (/ 2,2,2,2,4,4 /)
	REAL(SP), DIMENSION(NMAX) :: c,cout
	CHARACTER(32) :: spaces = '                                '
	CHARACTER(78), DIMENSION(NTEST) :: ans = (/ &
	'                          -0.086  0.343  0.486  0.343 -0.086                  ',&
	'                  -0.143  0.171  0.343  0.371  0.257                          ',&
	'           0.086 -0.143 -0.086  0.257  0.886                                  ',&
	'  -0.084  0.021  0.103  0.161  0.196  0.207  0.196  0.161  0.103  0.021 -0.084',&
	'           0.035 -0.128  0.070  0.315  0.417  0.315  0.070 -0.128  0.035      ',&
	'   0.042 -0.105 -0.023  0.140  0.280  0.333  0.280  0.140 -0.023 -0.105  0.042' /)
	write(*,*) 'M nl nr'
	write(*,'(t24,a)') 'Sample Savitzky-Golay Coefficients'
	do i=1,NTEST
		m=mtest(i)
		nl=nltest(i)
		nrr=nrtest(i)
		np=nl+nrr+1
		c(1:np)=savgol(nl,nrr,0,m)
		cout(1:nl+1)=c(nl+1:1:-1)
		cout(nl+2:nl+1+nrr)=c(np:np-nrr+1:-1)
		write(*,'(1x,3i2)') m,nl,nrr
		write(*,'(1x,a,11f7.3)') spaces(1:max(9*(5-nl)+nl-4,1)),cout(1:np)
		write(*,'(1x,a6,f7.3)') 'Sum = ',sum(c(1:np))
		write(*,'(1x,a12,/,1x,a)') 'Compare ans:',ans(i)
	end do
	END PROGRAM xsavgol
