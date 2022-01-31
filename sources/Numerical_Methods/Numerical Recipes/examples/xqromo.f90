	PROGRAM xqromo
!	driver for routine qromo
	USE nrtype
	USE nr
	REAL(SP), PARAMETER :: AINF=1.0e20_sp,ZERO=0.0
	REAL(SP) :: res1,res2,result
	INTERFACE
		FUNCTION funcl(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funcl
		END FUNCTION funcl
	END INTERFACE
	INTERFACE
		FUNCTION funcu(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: funcu
		END FUNCTION funcu
	END INTERFACE
	INTERFACE
		FUNCTION fncinf(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: fncinf
		END FUNCTION fncinf
	END INTERFACE
	INTERFACE
		FUNCTION fncend(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP), DIMENSION(size(x)) :: fncend
		END FUNCTION fncend
	END INTERFACE
	write(*,'(/1x,a)') 'Improper integrals:'
	result=qromo(funcl,ZERO,PIO2,midsql)
	write(*,'(/1x,a)') &
		'Function: SQRT(x)/SIN(x)       Interval: (0,pi/2)'
	write(*,'(1x,a,f8.4)') &
		'Using: MIDSQL                  Result:',result
	result=qromo(funcu,PIO2,PI,midsqu)
	write(*,'(/1x,a)') &
		'Function: SQRT(pi-x)/SIN(x)    Interval: (pi/2,pi)'
	write(*,'(1x,a,f8.4)') &
		'Using: MIDSQU                  Result:',result
	result=qromo(fncinf,PIO2,AINF,midinf)
	write(*,'(/1x,a)') &
		'Function: SIN(x)/x**2          Interval: (pi/2,infty)'
	write(*,'(1x,a,f8.4)') &
		'Using: MIDINF                  Result:',result
	result=qromo(fncinf,-AINF,-PIO2,midinf)
	write(*,'(/1x,a)') &
		'Function: SIN(x)/x**2          Interval: (-infty,-pi/2)'
	write(*,'(1x,a,f8.4)') &
		'Using: MIDINF                  Result:',result
	res1=qromo(fncend,ZERO,PIO2,midsql)
	res2=qromo(fncend,PIO2,AINF,midinf)
	write(*,'(/1x,a)') &
		'Function: EXP(-x)/SQRT(x)      Interval: (0.0,infty)'
	write(*,'(1x,a,f8.4)') &
		'Using: MIDSQL,MIDINF           Result:',res1+res2
	res2=qromo(fncend,PIO2,AINF,midexp)
	write(*,'(/1x,a)') &
		'Function: EXP(-x)/SQRT(x)      Interval: (0.0,infty)'
	write(*,'(1x,a,f8.4/)') &
		'Using: MIDSQL,MIDEXP           Result:',res1+res2
	END PROGRAM xqromo

	FUNCTION funcl(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: funcl
	funcl=sqrt(x)/sin(x)
	END FUNCTION funcl

	FUNCTION funcu(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: funcu
	funcu=sqrt(PI-x)/sin(x)
	END FUNCTION funcu

	FUNCTION fncinf(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: fncinf
	fncinf=sin(dble(x))/(x**2)
	END FUNCTION fncinf

	FUNCTION fncend(x)
	USE nrtype
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP), DIMENSION(size(x)) :: fncend
	fncend=exp(-x)/sqrt(x)
	END FUNCTION fncend
