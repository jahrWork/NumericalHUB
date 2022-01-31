	PROGRAM xfitexy
!	driver for routine fitexy
	USE nrtype
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTEGER(I4B), PARAMETER :: NPT=30
	REAL(SP) :: a,b,chi2,harvest,q,sa,sb,siga,sigb
	REAL(SP), DIMENSION(NPT) :: x,y,dx,dy,dz
	INTEGER(I4B) :: i
	call ran_seed(sequence=1411)
	dz(:)=0.0
	do i=1,NPT
		call ran1(harvest)
		dx(i)=0.1_sp+harvest
		call ran1(harvest)
		dy(i)=0.1_sp+harvest
		call gasdev(harvest)
		x(i)=10.0_sp+10.0_sp*harvest
		call gasdev(harvest)
		y(i)=2.0_sp*x(i)-5.0_sp+dy(i)*harvest
		call gasdev(harvest)
		x(i)=x(i)+dx(i)*harvest
	end do
	write(*,*) 'Values of a,b,siga,sigb,chi2,q:'
	write(*,*) 'Fit with x and y errors gives:'
	call fitexy(x,y,dx,dy,a,b,siga,sigb,chi2,q)
	write(*,'(1x,6f12.6)') a,b,siga,sigb,chi2,q
	write(*,*)
	write(*,*) 'Setting x errors to zero gives:'
	call fitexy(x,y,dz,dy,a,b,siga,sigb,chi2,q)
	write(*,'(1x,6f12.6)') a,b,siga,sigb,chi2,q
	write(*,*) '...to be compared with fit result:'
	call fit(x,y,a,b,siga,sigb,chi2,q,dy)
	write(*,'(1x,6f12.6)') a,b,siga,sigb,chi2,q
	write(*,*)
	write(*,*) 'Setting y errors to zero gives:'
	call fitexy(x,y,dx,dz,a,b,siga,sigb,chi2,q)
	write(*,'(1x,6f12.6)') a,b,siga,sigb,chi2,q
	write(*,*) '...to be compared with fit result:'
	call fit(y,x,a,b,siga,sigb,chi2,q,dx)
	sa=sqrt(siga**2+sigb**2*(a/b)**2)/b
	sb=sigb/b**2
	write(*,'(1x,6f12.6)') -a/b,1./b,sa,sb,chi2,q
	END PROGRAM xfitexy
