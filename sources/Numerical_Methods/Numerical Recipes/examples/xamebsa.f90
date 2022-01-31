	PROGRAM xamebsa
!	driver for routine amebsa
	USE nrtype; USE nrutil
	USE nr
	USE ran_state, ONLY : ran_seed
	IMPLICIT NONE
	INTERFACE
		FUNCTION tfunk(x)
		USE nrtype; USE nrutil
		IMPLICIT NONE
		REAL(SP), DIMENSION(:), INTENT(IN) :: x
		REAL(SP) :: tfunk
		END FUNCTION tfunk
	END INTERFACE
	REAL(SP), PARAMETER  :: FTOL=1.0e-6_sp
	INTEGER(I4B), PARAMETER :: NP=4,MP=5,MPNP=20
	INTEGER(I4B) :: i,iiter,iter,j,jiter,ndim,nit
	REAL(SP) :: temptr,yb,ybb
	REAL(SP), DIMENSION(NP) :: xoff,pb
	REAL(SP), DIMENSION(MP) :: y
	REAL(SP), DIMENSION(MP,NP) :: p
	call ran_seed(sequence=64)
	p(:,:)=0.0
	xoff(:)=10.0
	ndim=NP
	do
		do j=2,MP
			p(j,j-1)=1.0
		end do
		p(:,:)=p(:,:)+spread(xoff(:),dim=1,ncopies=MP)
		do i=1,MP
			y(i)=tfunk(p(i,:))
		end do
		yb=1.0e30_sp
		write(*,*) 'Input T, IITER:'
		read(*,*,END=999) temptr,iiter
		write(*,*) 'T=',temptr,'     IITER=',iiter
		ybb=1.0e30_sp
		nit=0
		do jiter=1,100
			iter=iiter
			temptr=temptr*0.8_sp
			call amebsa(p(1:MP,1:ndim),y(1:MP),pb(1:ndim),&
				yb,FTOL,tfunk,iter,temptr)
			nit=nit+iiter-iter
			if (yb < ybb) then
				ybb=yb
				write(*,'(1x,i6,e10.3,4f11.5,e15.7)')&
					nit,temptr,(pb(j),j=1,NP),yb
			end if
			if (iter > 0) exit
		end do
		write(*,'(/1x,a)') 'Vertices of final 3-D simplex and'
		write(*,'(1x,a)') 'function values at the vertices:'
		write(*,'(/3x,a,t11,a,t23,a,t35,a,t45,a,t55,a/)') 'I',&
			'X(I)','Y(I)','Z(I)','W(I)','FUNCTION'
		do i=1,MP
			write(*,'(1x,i3,4f12.6,e15.7)') i,(p(i,j),j=1,NP),y(i)
		end do
		write(*,'(1x,i3,4f12.6,e15.7)') 99,(pb(j),j=1,NP),yb
	end do
999	write(*,*) 'NORMAL COMPLETION'
	END PROGRAM xamebsa

	FUNCTION tfunk(p)
	USE nrtype; USE nrutil
	IMPLICIT NONE
	REAL(SP), PARAMETER :: RAD=0.3_sp,AUG=2.0
	REAL(SP) :: sumd,sumr,tfunk
	REAL(SP), DIMENSION(:), INTENT(IN) :: p
	REAL(SP), DIMENSION(size(p)) :: wid,q,r
	if (size(p) /= 4) call nrerror('xamebsa: error in driver, in function')
	wid(1:4) = (/ 1.0_sp,3.0_sp,10.0_sp,30.0_sp /)
	q(:)=p(:)*wid(:)
	r(:)=nint(q(:))
	sumr=dot_product(q(:),q(:))
	q(:)=q(:)-r(:)
	sumd=dot_product(q(:),q(:))
	tfunk=sumr*(1.0_sp+merge(AUG,AUG*sumd/RAD**2, sumd > RAD**2 ))+1.0_sp
	END FUNCTION tfunk
