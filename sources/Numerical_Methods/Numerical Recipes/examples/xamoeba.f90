	PROGRAM xamoeba
!	driver for routine amoeba
	USE nrtype
	USE nr
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NP=3,MP=4
	REAL(SP), PARAMETER :: FTOL=1.0e-6_sp
	INTEGER(I4B) :: i,iter,j,ndim
	REAL(SP), DIMENSION(MP,NP) :: p
	REAL(SP), DIMENSION(MP) :: y
	INTERFACE
		FUNCTION famoeb(x)
		USE nrtype
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: famoeb
		END FUNCTION famoeb
	END INTERFACE
	p(:,:)=0.0
	do j=2,MP
		p(j,j-1)=1.0
	end do
	ndim=NP
	do i=1,MP
		y(i)=famoeb(p(i,:))
	end do
	call amoeba(p(1:MP,1:ndim),y(1:MP),FTOL,famoeb,iter)
	write(*,'(/1x,a,i3)') 'Number of iterations: ',iter
	write(*,'(/1x,a)') 'Vertices of final 3-D simplex and'
	write(*,'(1x,a)') 'function values at the vertices:'
	write(*,'(/3x,a,t11,a,t23,a,t35,a,t45,a/)') 'I',&
		'X(I)','Y(I)','Z(I)','FUNCTION'
	do i=1,MP
		write(*,'(1x,i3,4f12.6)') i,(p(i,j),j=1,NP),y(i)
	end do
	write(*,'(/1x,a)') 'True minimum is at (0.5,0.6,0.7)'
	END PROGRAM xamoeba

	FUNCTION famoeb(x)
	USE nrtype
	USE nr
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(IN) :: x
	REAL(SP) :: famoeb
	famoeb=0.6_sp-bessj0((x(1)-0.5_sp)**2+(x(2)-0.6_sp)**2+(x(3)-0.7_sp)**2)
	END FUNCTION famoeb
