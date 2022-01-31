	PROGRAM xmachar
!	driver for routine machar
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B) :: ibeta,iexp,irnd,it,machep,maxexp,minexp,negep,ngrd
	REAL(SP) :: eps,epsneg,xmax,xmin
	call machar(ibeta,it,irnd,ngrd,machep,negep,iexp,minexp,&
		maxexp,eps,epsneg,xmin,xmax)
	write(*,*) 'ibeta=',ibeta
	write(*,*) 'it=',it
	write(*,*) 'irnd=',irnd
	write(*,*) 'ngrd=',ngrd
	write(*,*) 'machep=',machep
	write(*,*) 'negep=',negep
	write(*,*) 'iexp=',iexp
	write(*,*) 'minexp=',minexp
	write(*,*) 'maxexp=',maxexp
	write(*,*) 'eps=',eps
	write(*,*) 'epsneg=',epsneg
	write(*,*) 'xmin=',xmin
	write(*,*) 'xmax=',xmax
	END PROGRAM xmachar
