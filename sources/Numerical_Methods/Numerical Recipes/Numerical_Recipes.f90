!**********************************************************************************
! Envoltorios de las subrutinas de Numerical Recipes 
!**********************************************************************************
module Numerical_Recipes 

 character(len=20), private  :: Busca_ceros='Broyden'


    contains 
!**************************************************************************
!    Eigenvalues of A 
!**************************************************************************
   subroutine Eigenvalues_QR( A, lambda)
           real, intent(in):: A(:,:)        
           complex, intent(out):: lambda(:)   
				    
          integer :: N 
          real, allocatable :: lar(:), lai(:) 
           
          N = size( lambda ) 
          allocate( lar(N), lai(N) )
          
          call elmhes(A, N, N)
	      call hqr(A, N, N, lar, lai) 
          
          lambda = cmplx( lar, lai ) 

	end subroutine 


!**********************************************************************************
! Wrapper:  Broydn y Funcv 
!**********************************************************************************
subroutine Broyden (G,  x ) 
          interface 
           function G(x) 
               real, intent(in) :: x(:)
               real :: G(size(x)) 
           end function 
        end interface 

        real, intent(inout) :: x(:) 

      logical check 


       call Broydn(x, size(x), check) 

  contains  

!------------------------------------------------------------------------------------
! user suplied subroutine para Broydn
!------------------------------------------------------------------------------------
     subroutine Funcv(n, x, fvec) 
            integer n 
            real :: x(n), fvec(n) 
 
                 fvec =  G(x) 

     end subroutine 

!-------------------------------------------------------------------------------------
! subroutine de Numerical recipes sin modificar 
!-------------------------------------------------------------------------------------
 
      SUBROUTINE broydn(x,n,check)
      INTEGER n,nn,NP,MAXITS
      real x(n),fvec,EPS,TOLF,TOLMIN,TOLX,STPMX
      LOGICAL check
	  
      PARAMETER (NP=400,MAXITS=200,EPS=1.d-7,TOLF=1.d-4,TOLMIN=1.d-6, &
      TOLX=EPS,STPMX=100.d0)
      COMMON /newtv/ fvec(NP),nn
!CU    USES fdjac,fmin,lnsrch,qrdcmp,qrupdt,rsolv
      INTEGER i,its,j,k
      real den,f,fold,stpmax,sum,temp,test,c(NP),d(NP),fvcold(NP),g(NP), &
      p(NP),qt(NP,NP),r(NP,NP),s(NP),t(NP),w(NP),xold(NP) !,fmin
      LOGICAL restrt,sing,skip
  !    EXTERNAL fmin
      nn=n
      f=fmin(x)
      test=0.d0
      do 11 i=1,n
        if(abs(fvec(i)).gt.test)test=abs(fvec(i))
11    continue
      if(test.lt..01d0*TOLF)return
      sum=0.d0
      do 12 i=1,n
        sum=sum+x(i)**2
12    continue
      stpmax=STPMX*max(sqrt(sum), dble(n))
      restrt=.true.
      do 44 its=1,MAXITS
        if(restrt)then
          call fdjac(n,x,fvec,NP,r)
          call qrdcmp(r,n,NP,c,d,sing)
          if(sing) then 
		                 write(*,*) 'singular Jacobian in broydn'
						 stop 
          endif 

          do 14 i=1,n
            do 13 j=1,n
              qt(i,j)=0.d0
13          continue
            qt(i,i)=1.d0
14        continue
          do 18 k=1,n-1
            if(c(k).ne.0.d0)then
              do 17 j=1,n
                sum=0.d0
                do 15 i=k,n
                  sum=sum+r(i,k)*qt(i,j)
15              continue
                sum=sum/c(k)
                do 16 i=k,n
                  qt(i,j)=qt(i,j)-sum*r(i,k)
16              continue
17            continue
            endif
18        continue
          do 21 i=1,n
            r(i,i)=d(i)
            do 19 j=1,i-1
              r(i,j)=0.d0
19          continue
21        continue
        else
          do 22 i=1,n
            s(i)=x(i)-xold(i)
22        continue
          do 24 i=1,n
            sum=0.d0
            do 23 j=i,n
              sum=sum+r(i,j)*s(j)
23          continue
            t(i)=sum
24        continue
          skip=.true.
          do 26 i=1,n
            sum=0.d0
            do 25 j=1,n
              sum=sum+qt(j,i)*t(j)
25          continue
            w(i)=fvec(i)-fvcold(i)-sum
            if(abs(w(i)).ge.EPS*(abs(fvec(i))+abs(fvcold(i))))then
              skip=.false.
            else
              w(i)=0.d0
            endif
26        continue
          if(.not.skip)then
            do 28 i=1,n
              sum=0.d0
              do 27 j=1,n
                sum=sum+qt(i,j)*w(j)
27            continue
              t(i)=sum
28          continue
            den=0.d0
            do 29 i=1,n
              den=den+s(i)**2
29          continue
            do 31 i=1,n
              s(i)=s(i)/den
31          continue
            call qrupdt(r,qt,n,NP,t,s)
            do 32 i=1,n
              if(r(i,i).eq.0.d0) then 
			         write(*,*) 'r singular in broydn'
					 stop 
              endif  
              d(i)=r(i,i)
32          continue
          endif
        endif
        do 34 i=1,n
          sum=0.d0
          do 33 j=1,n
            sum=sum+qt(i,j)*fvec(j)
33        continue
          g(i)=sum
34      continue
        do 36 i=n,1,-1
          sum=0.d0
          do 35 j=1,i
            sum=sum+r(j,i)*g(j)
35        continue
          g(i)=sum
36      continue
        do 37 i=1,n
          xold(i)=x(i)
          fvcold(i)=fvec(i)
37      continue
        fold=f
        do 39 i=1,n
          sum=0.d0
          do 38 j=1,n
            sum=sum+qt(i,j)*fvec(j)
38        continue
          p(i)=-sum
39      continue
        call rsolv(r,n,NP,d,p)
        call lnsrch(n,xold,fold,g,p,x,f,stpmax,check,fmin)
        test=0.d0
        do 41 i=1,n
          if(abs(fvec(i)).gt.test)test=abs(fvec(i))
41      continue
        if(test.lt.TOLF)then
          check=.false.
          return
        endif
        if(check)then
          if(restrt)then
            return
          else
            test=0.d0
            den=max(f,.5d0*n)
            do 42 i=1,n
              temp=abs(g(i))*max(abs(x(i)),1.d0)/den
              if(temp.gt.test)test=temp
42          continue
            if(test.lt.TOLMIN)then
              return
            else
              restrt=.true.
            endif
          endif
        else
          restrt=.false.
          test=0.d0
          do 43 i=1,n
            temp=(abs(x(i)-xold(i)))/max(abs(x(i)),1.d0)
            if(temp.gt.test)test=temp
43        continue
          if(test.lt.TOLX)return
        endif
44    continue
      
      write(*,*) 'MAXITS exceeded in broydn'
	  stop 

      END subroutine 
!  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.


!--------------------------------------------------------------------------------
     FUNCTION fmin(x)
      INTEGER n,NP
      real fmin,x(*),fvec
      PARAMETER (NP=400)
      COMMON /newtv/ fvec(NP),n
      SAVE /newtv/
!CU    USES funcv
      INTEGER i
      real sum
      call funcv(n,x,fvec)
      sum=0.d0
      do 11 i=1,n
        sum=sum+fvec(i)**2
11    continue
      fmin=0.5d0*sum
      return
      END function 
!  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.


!-----------------------------------------------------------------------------
   SUBROUTINE fdjac(n,x,fvec,np,df)
      INTEGER n,np,NMAX
      real df(np,np),fvec(n),x(n),EPS
      PARAMETER (NMAX=400,EPS=1.d-4)
!U    USES funcv
      INTEGER i,j
      real h,temp,f(NMAX)
      do 12 j=1,n
        temp=x(j)
        h=EPS*abs(temp)
        if(h.eq.0.d0)h=EPS
        x(j)=temp+h
        h=x(j)-temp
        call funcv(n,x,f)
        x(j)=temp
        do 11 i=1,n
          df(i,j)=(f(i)-fvec(i))/h
11      continue
12    continue
      return
      END subroutine 
!C  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.

!-----------------------------------------------------------------------------------
      SUBROUTINE lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      INTEGER n
      LOGICAL check
      real f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
      PARAMETER (ALF=1.d-4,TOLX=1.d-7)
      EXTERNAL func
!U    USES func
      INTEGER i
      real a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2, & 
         slope,sum,temp,test,tmplam
      check=.false.
      sum=0.d0
      do 11 i=1,n
        sum=sum+p(i)*p(i)
11    continue
      sum=sqrt(sum)
      if(sum.gt.stpmax)then
        do 12 i=1,n
          p(i)=p(i)*stpmax/sum
12      continue
      endif
      slope=0.d0
      do 13 i=1,n
        slope=slope+g(i)*p(i)
13    continue
      test=0.d0
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.d0)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.d0
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.d0)then
            tmplam=-slope/(2.d0*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold2-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.d0)then
              tmplam=-slope/(2.d0*b)
            else
              disc=b*b-3.d0*a*slope
              tmplam=(-b+sqrt(disc))/(3.d0*a)
            endif
            if(tmplam.gt..5d0*alam)tmplam=.5d0*alam
          endif
        endif
        alam2=alam
        f2=f
        fold2=fold
        alam=max(tmplam,.1d0*alam)
      goto 1
      END subroutine 
!  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.


end subroutine
 
 
     SUBROUTINE qrdcmp(a,n,np,c,d,sing)
      INTEGER n,np
      real a(np,np),c(n),d(n)
      LOGICAL sing
      INTEGER i,j,k
      real scale,sigma,sum,tau
      sing=.false.
      scale=0.d0
      do 17 k=1,n-1
        do 11 i=k,n
          scale=max(scale,abs(a(i,k)))
11      continue
        if(scale.eq.0.d0)then
          sing=.true.
          c(k)=0.d0
          d(k)=0.d0
        else
          do 12 i=k,n
            a(i,k)=a(i,k)/scale
12        continue
          sum=0.d0
          do 13 i=k,n
            sum=sum+a(i,k)**2
13        continue
          sigma=sign(sqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do 16 j=k+1,n
            sum=0.d0
            do 14 i=k,n
              sum=sum+a(i,k)*a(i,j)
14          continue
            tau=sum/c(k)
            do 15 i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
15          continue
16        continue
        endif
17    continue
      d(n)=a(n,n)
      if(d(n).eq.0.d0)sing=.true.
      return
      END subroutine 
!C  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.


 

    SUBROUTINE qrupdt(r,qt,n,np,u,v)
      INTEGER n,np
       real r(np,np),qt(np,np),u(np),v(np)
!CU    USES rotate
      INTEGER i,j,k
      do 11 k=n,1,-1
        if(u(k).ne.0.d0)goto 1
11    continue
      k=1
1     do 12 i=k-1,1,-1
        call rotate(r,qt,n,np,i,u(i),-u(i+1))
        if(u(i).eq.0.d0)then
          u(i)=abs(u(i+1))
        else if(abs(u(i)).gt.abs(u(i+1)))then
          u(i)=abs(u(i))*sqrt(1.d0+(u(i+1)/u(i))**2)
        else
          u(i)=abs(u(i+1))*sqrt(1.d0+(u(i)/u(i+1))**2)
        endif
12    continue
      do 13 j=1,n
        r(1,j)=r(1,j)+u(1)*v(j)
13    continue
      do 14 i=1,k-1
        call rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
14    continue
      return
      end subroutine
!C  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.




      SUBROUTINE rsolv(a,n,np,d,b)
      INTEGER n,np
      real a(np,np),b(n),d(n)
      INTEGER i,j
      real sum
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
        sum=0.d0
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue
        b(i)=(b(i)-sum)/d(i)
12    continue
      return
      end subroutine
!C  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.



      SUBROUTINE rotate(r,qt,n,np,i,a,b)
      INTEGER n,np,i
      real a,b,r(np,np),qt(np,np)
      INTEGER j
      real c,fact,s,w,y
      if(a.eq.0.d0)then
        c=0.d0
        s=sign(1.d0,b)
      else if(abs(a).gt.abs(b))then
        fact=b/a
        c=sign(1.d0/sqrt(1.d0+fact**2),a)
        s=fact*c
      else
        fact=a/b
        s=sign(1.d0/sqrt(1.d0+fact**2),b)
        c=fact*s
      endif
      do 11 j=i,n
        y=r(i,j)
        w=r(i+1,j)
        r(i,j)=c*y-s*w
        r(i+1,j)=s*y+c*w
11    continue
      do 12 j=1,n
        y=qt(i,j)
        w=qt(i+1,j)
        qt(i,j)=c*y-s*w
        qt(i+1,j)=s*y+c*w
12    continue
      return
      end subroutine
!C  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.

      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      end subroutine
!C  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.

      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      end subroutine
!C  (C) Copr. 1986-92 Numerical Recipes Software ]u]w#!!0,)#!$!.


      SUBROUTINE elmhes(a,n,np)
      INTEGER n,np
      DOUBLE PRECISION a(np,np)
      INTEGER i,j,m
      DOUBLE PRECISION x,y
      do 17 m=2,n-1
        x=0.d0
        i=m
        do 11 j=m,n
          if(abs(a(j,m-1)).gt.abs(x))then
            x=a(j,m-1)
            i=j
          endif
11      continue
        if(i.ne.m)then
          do 12 j=m-1,n
            y=a(i,j)
            a(i,j)=a(m,j)
            a(m,j)=y
12        continue
          do 13 j=1,n
            y=a(j,i)
            a(j,i)=a(j,m)
            a(j,m)=y
13        continue
        endif
        if(x.ne.0.d0)then
          do 16 i=m+1,n
            y=a(i,m-1)
            if(y.ne.0.d0)then
              y=y/x
              a(i,m-1)=y
              do 14 j=m,n
                a(i,j)=a(i,j)-y*a(m,j)
14            continue
              do 15 j=1,n
                a(j,m)=a(j,m)+y*a(j,i)
15            continue
            endif
16        continue
        endif
17    continue
      return
      end subroutine
!  (C) Copr. 1986-92 Numerical Recipes Software *5=!R2.=D:i-13.d0


      SUBROUTINE hqr(a,n,np,wr,wi)
      INTEGER n,np
      DOUBLE PRECISION a(np,np),wi(np),wr(np)
      INTEGER i,its,j,k,l,m,nn
      DOUBLE PRECISION anorm,p,q,r,s,t,u,v,w,x,y,z
      anorm=abs(a(1,1))
      do 12 i=2,n
        do 11 j=i-1,n
          anorm=anorm+abs(a(i,j))
11      continue
12    continue
      nn=n
      t=0.d0
1     if(nn.ge.1)then
        its=0
2       do 13 l=nn,2,-1
          s=abs(a(l-1,l-1))+abs(a(l,l))
          if(s.eq.0.d0)s=anorm
          if(abs(a(l,l-1))+s.eq.s)goto 3
13      continue
        l=1
3       x=a(nn,nn)
        if(l.eq.nn)then
          wr(nn)=x+t
          wi(nn)=0.d0
          nn=nn-1
        else
          y=a(nn-1,nn-1)
          w=a(nn,nn-1)*a(nn-1,nn)
          if(l.eq.nn-1)then
            p=0.5d0*(y-x)
            q=p**2+w
            z=sqrt(abs(q))
            x=x+t
            if(q.ge.0.d0)then
              z=p+sign(z,p)
              wr(nn)=x+z
              wr(nn-1)=wr(nn)
              if(z.ne.0.d0)wr(nn)=x-w/z
              wi(nn)=0.d0
              wi(nn-1)=0.d0
            else
              wr(nn)=x+p
              wr(nn-1)=wr(nn)
              wi(nn)=z
              wi(nn-1)=-z
            endif
            nn=nn-2
          else
            if(its.eq.30)pause 'too many iterations in hqr'
            if(its.eq.10.or.its.eq.20)then
              t=t+x
              do 14 i=1,nn
                a(i,i)=a(i,i)-x
14            continue
              s=abs(a(nn,nn-1))+abs(a(nn-1,nn-2))
              x=0.75d0*s
              y=x
              w=-0.4375d0*s**2
            endif
            its=its+1
            do 15 m=nn-2,l,-1
              z=a(m,m)
              r=x-z
              s=y-z
              p=(r*s-w)/a(m+1,m)+a(m,m+1)
              q=a(m+1,m+1)-z-r-s
              r=a(m+2,m+1)
              s=abs(p)+abs(q)+abs(r)
              p=p/s
              q=q/s
              r=r/s
              if(m.eq.l)goto 4
              u=abs(a(m,m-1))*(abs(q)+abs(r))
              v=abs(p)*(abs(a(m-1,m-1))+abs(z)+abs(a(m+1,m+1)))
              if(u+v.eq.v)goto 4
15          continue
4           do 16 i=m+2,nn
              a(i,i-2)=0.d0
              if (i.ne.m+2) a(i,i-3)=0.d0
16          continue
            do 19 k=m,nn-1
              if(k.ne.m)then
                p=a(k,k-1)
                q=a(k+1,k-1)
                r=0.d0
                if(k.ne.nn-1)r=a(k+2,k-1)
                x=abs(p)+abs(q)+abs(r)
                if(x.ne.0.d0)then
                  p=p/x
                  q=q/x
                  r=r/x
                endif
              endif
              s=sign(sqrt(p**2+q**2+r**2),p)
              if(s.ne.0.d0)then
                if(k.eq.m)then
                  if(l.ne.m)a(k,k-1)=-a(k,k-1)
                else
                  a(k,k-1)=-s*x
                endif
                p=p+s
                x=p/s
                y=q/s
                z=r/s
                q=q/p
                r=r/p
                do 17 j=k,nn
                  p=a(k,j)+q*a(k+1,j)
                  if(k.ne.nn-1)then
                    p=p+r*a(k+2,j)
                    a(k+2,j)=a(k+2,j)-p*z
                  endif
                  a(k+1,j)=a(k+1,j)-p*y
                  a(k,j)=a(k,j)-p*x
17              continue
                do 18 i=l,min(nn,k+3)
                  p=x*a(i,k)+y*a(i,k+1)
                  if(k.ne.nn-1)then
                    p=p+z*a(i,k+2)
                    a(i,k+2)=a(i,k+2)-p*r
                  endif
                  a(i,k+1)=a(i,k+1)-p*q
                  a(i,k)=a(i,k)-p
18              continue
              endif
19          continue
            goto 2
          endif
        endif
      goto 1
      endif
      return
      end subroutine
!     (C) Copr. 1986-92 Numerical Recipes Software *5=!R2.=D:i-13.d0



      subroutine zsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)
      integer ldx,n,p,ldu,ldv,job,info
      complex*16 x(ldx,1),s(1),e(1),u(ldu,1),v(ldv,1),work(1)
!c
!c
!c     zsvdc is a subroutine to reduce a complex*16 nxp matrix x by
!c     unitary transformations u and v to diagonal form.  the
!c     diagonal elements s(i) are the singular values of x.  the
!c     columns of u are the corresponding left singular vectors,
!c     and the columns of v the right singular vectors.
!c
!c     on entry
!c
!c         x         complex*16(ldx,p), where ldx.ge.n.
!c                   x contains the matrix whose singular value
!c                   decomposition is to be computed.  x is
!c                   destroyed by zsvdc.
!c
!c         ldx       integer.
!c                   ldx is the leading dimension of the array x.
!c
!c         n         integer.
!c                   n is the number of rows of the matrix x.
!c
!c         p         integer.
!c                   p is the number of columns of the matrix x.
!c
!c         ldu       integer.
!c                   ldu is the leading dimension of the array u
!c                   (see below).
!c
!c         ldv       integer.
!c                   ldv is the leading dimension of the array v
!c                   (see below).
!c
!c         work      complex*16(n).
!c                   work is a scratch array.
!c
!c         job       integer.
!c                   job controls the computation of the singular
!c                   vectors.  it has the decimal expansion ab
!c                   with the following meaning
!c
!c                        a.eq.0    do not compute the left singular
!c                                  vectors.
!c                        a.eq.1    return the n left singular vectors
!c                                  in u.
!c                        a.ge.2    returns the first min(n,p)
!c                                  left singular vectors in u.
!c                        b.eq.0    do not compute the right singular
!c                                  vectors.
!c                        b.eq.1    return the right singular vectors
!c                                  in v.
!c
!c     on return
!c
!c         s         complex*16(mm), where mm=min(n+1,p).
!c                   the first min(n,p) entries of s contain the
!c                   singular values of x arranged in descending
!c                   order of magnitude.
!c
!c         e         complex*16(p).
!c                   e ordinarily contains zeros.  however see the
!c                   discussion of info for exceptions.
!c
!c         u         complex*16(ldu,k), where ldu.ge.n.  if joba.eq.1
!c                                   then k.eq.n, if joba.ge.2 then
!c
!c                                   k.eq.min(n,p).
!c                   u contains the matrix of left singular vectors.
!c                   u is not referenced if joba.eq.0.  if n.le.p
!c                   or if joba.gt.2, then u may be identified with x
!c                   in the subroutine call.
!c
!c         v         complex*16(ldv,p), where ldv.ge.p.
!c                   v contains the matrix of right singular vectors.
!c                   v is not referenced if jobb.eq.0.  if p.le.n,
!c                   then v may be identified whth x in the
!c                   subroutine call.
!c
!c         info      integer.
!c                   the singular values (and their corresponding
!c                   singular vectors) s(info+1),s(info+2),...,s(m)
!c                   are correct (here m=min(n,p)).  thus if
!c                   info.eq.0, all the singular values and their
!c                   vectors are correct.  in any event, the matrix
!c                   b = ctrans(u)*x*v is the bidiagonal matrix
!c                   with the elements of s on its diagonal and the
!c                   elements of e on its super-diagonal (ctrans(u)
!c                   is the conjugate-transpose of u).  thus the
!c                   singular values of x and b are the same.
!c
!c     linpack. this version dated 03/19/79 .
!c              correction to shift calculation made 2/85.
!c     g.w. stewart, university of maryland, argonne national lab.
!c
!c     zsvdc uses the following functions and subprograms.
!c
!c     external zdrot
!c     blas zaxpy,zdotc,zscal,zswap,dznrm2,drotg
!c     fortran dabs,dmax1,cdabs,dcmplx
!c     fortran dconjg,max0,min0,mod,dsqrt
!c
!c     internal variables
!c
      integer i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit,    mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
      complex*16 t,r !zdotc,
      double precision b,c,cs,el,emm1,f,g,scale,shift,sl,sm,sn,   smm1,t1,test,ztest !,dznrm2
      logical wantu,wantv
!c
      complex*16 csign,zdum,zdum1,zdum2
      double precision cabs1
	
      double precision abs1
!c      complex*16 zdumr,zdumi
!c      dreal(zdumr) = zdumr   
!c      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
!c      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
      csign(zdum1,zdum2) = cdabs(zdum1)*(zdum2/cdabs(zdum2))
	abs1(zdum) = abs(real(zdum)) + abs(imag(zdum)) 
!c
!c     set the maximum number of iterations.
!c
      maxit = 30

      write(*,*) ' maxit =', maxit 
!c
!c     determine what is to be computed.
!c
      wantu = .false.
      wantv = .false.
      jobu = mod(job,100)/10
      ncu = n
      if (jobu .gt. 1) ncu = min0(n,p)
      if (jobu .ne. 0) wantu = .true.
      if (mod(job,10) .ne. 0) wantv = .true.
!c
!c     reduce x to bidiagonal form, storing the diagonal elements
!c     in s and the super-diagonal elements in e.
!c
      info = 0
      nct = min0(n-1,p)
      nrt = max0(0,min0(p-2,n))
      lu = max0(nct,nrt)
      if (lu .lt. 1) go to 170
      do 160 l = 1, lu
         lp1 = l + 1
         if (l .gt. nct) go to 20
!c
!c           compute the transformation for the l-th column and
!c           place the l-th diagonal in s(l).
!c
            s(l) = dcmplx(dznrm2(n-l+1,x(l,l),1),0.0d0)
            if (abs1(s(l)) .eq. 0.0d0) go to 10
               if (abs1(x(l,l)) .ne. 0.0d0) s(l) = csign(s(l),x(l,l))
               call zscal(n-l+1,1.0d0/s(l),x(l,l),1)
               x(l,l) = (1.0d0,0.0d0) + x(l,l)
   10       continue
            s(l) = -s(l)
   20    continue
         if (p .lt. lp1) go to 50
         do 40 j = lp1, p
            if (l .gt. nct) go to 30
            if (abs1(s(l)) .eq. 0.0d0) go to 30
!c
!c              apply the transformation.
!c
               t = -zdotc(n-l+1, x(l,l),1, x(l,j),1) / x(l,l)
               call zaxpy(n-l+1,t,x(l,l),1,x(l,j),1)
   30       continue
!c
!c           place the l-th row of x into  e for the
!c           subsequent calculation of the row transformation.
!c
            e(j) = dconjg(x(l,j))
   40    continue
   50    continue
         if (.not.wantu .or. l .gt. nct) go to 70
!c
!c           place the transformation in u for subsequent back
!c           multiplication.
!c
            do 60 i = l, n
               u(i,l) = x(i,l)
   60       continue
   70    continue
         if (l .gt. nrt) go to 150
!c
!c           compute the l-th row transformation and place the
!c           l-th super-diagonal in e(l).
!c
            e(l) = dcmplx(dznrm2(p-l,e(lp1),1),0.0d0)
            if (abs1(e(l)) .eq. 0.0d0) go to 80
               if (abs1(e(lp1)) .ne. 0.0d0) e(l) = csign(e(l),e(lp1))
               call zscal(p-l,1.0d0/e(l),e(lp1),1)
               e(lp1) = (1.0d0,0.0d0) + e(lp1)
   80       continue
            e(l) = -dconjg(e(l))
            if (lp1 .gt. n .or. abs1(e(l)) .eq. 0.0d0) go to 120
!c
!c              apply the transformation.
!c
               do 90 i = lp1, n
                  work(i) = (0.0d0,0.0d0)
   90          continue
               do 100 j = lp1, p
                  call zaxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
  100          continue
               do 110 j = lp1, p
                  call zaxpy(n-l,dconjg(-e(j)/e(lp1)),work(lp1),1,  x(lp1,j),1)
  110          continue
  120       continue
            if (.not.wantv) go to 140
!c
!c              place the transformation in v for subsequent
!c              back multiplication.
!c
               do 130 i = lp1, p
                  v(i,l) = e(i)
  130          continue
  140       continue
  150    continue
  160 continue
  170 continue
!c
!c     set up the final bidiagonal matrix or order m.
!c
      m = min0(p,n+1)
      nctp1 = nct + 1
      nrtp1 = nrt + 1
      if (nct .lt. p) s(nctp1) = x(nctp1,nctp1)
      if (n .lt. m) s(m) = (0.0d0,0.0d0)
      if (nrtp1 .lt. m) e(nrtp1) = x(nrtp1,m)
      e(m) = (0.0d0,0.0d0)
!c
!c     if required, generate u.
!c
      if (.not.wantu) go to 300
         if (ncu .lt. nctp1) go to 200
         do 190 j = nctp1, ncu
            do 180 i = 1, n
               u(i,j) = (0.0d0,0.0d0)
  180       continue
            u(j,j) = (1.0d0,0.0d0)
  190    continue
  200    continue
         if (nct .lt. 1) go to 290
         do 280 ll = 1, nct
            l = nct - ll + 1
            if (abs1(s(l)) .eq. 0.0d0) go to 250
               lp1 = l + 1
               if (ncu .lt. lp1) go to 220
               do 210 j = lp1, ncu
                  t = -zdotc(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
                  call zaxpy(n-l+1,t,u(l,l),1,u(l,j),1)
  210          continue
  220          continue
               call zscal(n-l+1,(-1.0d0,0.0d0),u(l,l),1)
               u(l,l) = (1.0d0,0.0d0) + u(l,l)
               lm1 = l - 1
               if (lm1 .lt. 1) go to 240
               do 230 i = 1, lm1
                  u(i,l) = (0.0d0,0.0d0)
  230          continue
  240          continue
            go to 270
  250       continue
               do 260 i = 1, n
                  u(i,l) = (0.0d0,0.0d0)
  260          continue
               u(l,l) = (1.0d0,0.0d0)
  270       continue
  280    continue
  290    continue
  300 continue
!c
!c     if it is required, generate v.
!c
      if (.not.wantv) go to 350
         do 340 ll = 1, p
            l = p - ll + 1
            lp1 = l + 1
            if (l .gt. nrt) go to 320
            if (abs1(e(l)) .eq. 0.0d0) go to 320
               do 310 j = lp1, p
                  t = -zdotc(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
                  call zaxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
  310          continue
  320       continue
            do 330 i = 1, p
               v(i,l) = (0.0d0,0.0d0)
  330       continue
            v(l,l) = (1.0d0,0.0d0)
  340    continue
  350 continue
!c
!c     transform s and e so that they are double precision.
!c
      do 380 i = 1, m
         if (abs1(s(i)) .eq. 0.0d0) go to 360
            t = dcmplx(cdabs(s(i)),0.0d0)
            r = s(i)/t
            s(i) = t
            if (i .lt. m) e(i) = e(i)/r
            if (wantu) call zscal(n,r,u(1,i),1)
  360    continue
!c     ...exit
         if (i .eq. m) go to 390
         if (abs1(e(i)) .eq. 0.0d0) go to 370
            t = dcmplx(cdabs(e(i)),0.0d0)
            r = t/e(i)
            e(i) = t
            s(i+1) = s(i+1)*r
            if (wantv) call zscal(p,r,v(1,i+1),1)
  370    continue
  380 continue
  390 continue
!c
!c     main iteration loop for the singular values.
!c
      mm = m
      iter = 0
      write(*,*) ' iteraciones ' 
  400 continue
!c
!c        quit if all the singular values have been found.
!c
!c     ...exit
         if (m .eq. 0) go to 660
!c
!c        if too many iterations have been performed, set
!c        flag and return.
!c
         if (iter .lt. maxit) go to 410
            info = m
!c     ......exit
            go to 660
  410    continue
!c
!c        this section of the program inspects for
!c        negligible elements in the s and e arrays.  on
!c        completion the variables kase and l are set as follows.
!c
!c           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
!c           kase = 2     if s(l) is negligible and l.lt.m
!c           kase = 3     if e(l-1) is negligible, l.lt.m, and
!c                        s(l), ..., s(m) are not negligible (qr step).
!c           kase = 4     if e(m-1) is negligible (convergence).
!c
         do 430 ll = 1, m
            l = m - ll
!c        ...exit
            if (l .eq. 0) go to 440
            test = cdabs(s(l)) + cdabs(s(l+1))
            ztest = test + cdabs(e(l))
            if (ztest .ne. test) go to 420
               e(l) = (0.0d0,0.0d0)
!c        ......exit
               go to 440
  420       continue
  430    continue
  440    continue
         if (l .ne. m - 1) go to 450
            kase = 4
         go to 520
  450    continue
            lp1 = l + 1
            mp1 = m + 1
            do 470 lls = lp1, mp1
               ls = m - lls + lp1
!c           ...exit
               if (ls .eq. l) go to 480
               test = 0.0d0
               if (ls .ne. m) test = test + cdabs(e(ls))
               if (ls .ne. l + 1) test = test + cdabs(e(ls-1))
               ztest = test + cdabs(s(ls))
               if (ztest .ne. test) go to 460
                  s(ls) = (0.0d0,0.0d0)
!c           ......exit
                  go to 480
  460          continue
  470       continue
  480       continue
            if (ls .ne. l) go to 490
               kase = 3
            go to 510
  490       continue
            if (ls .ne. m) go to 500
               kase = 1
            go to 510
  500       continue
               kase = 2
               l = ls
  510       continue
  520    continue
         l = l + 1
!c
!c        perform the task indicated by kase.
!c
!         go to (530, 560, 580, 610), kase
!c
!c        deflate negligible s(m).
!c
  530    continue
            mm1 = m - 1
            f = dreal(e(m-1))
            e(m-1) = (0.0d0,0.0d0)
            do 550 kk = l, mm1
               k = mm1 - kk + l
               t1 = dreal(s(k))
               call drotg(t1,f,cs,sn)
               s(k) = dcmplx(t1,0.0d0)
               if (k .eq. l) go to 540
                  f = -sn*dreal(e(k-1))
                  e(k-1) = cs*e(k-1)
  540          continue
               if (wantv) call zdrot(p,v(1,k),1,v(1,m),1,cs,sn)
  550       continue
         go to 650
!c
!c        split at negligible s(l).
!c
  560    continue
            f = dreal(e(l-1))
            e(l-1) = (0.0d0,0.0d0)
            do 570 k = l, m
               t1 = dreal(s(k))
               call drotg(t1,f,cs,sn)
               s(k) = dcmplx(t1,0.0d0)
               f = -sn*dreal(e(k))
               e(k) = cs*e(k)
               if (wantu) call zdrot(n,u(1,k),1,u(1,l-1),1,cs,sn)
  570       continue
         go to 650
!c
!c        perform one qr step.
!c
  580    continue
!c
!c           calculate the shift.
!c
            scale = dmax1(cdabs(s(m)),cdabs(s(m-1)),cdabs(e(m-1)), cdabs(s(l)),cdabs(e(l)))
            sm = dreal(s(m))/scale
            smm1 = dreal(s(m-1))/scale
            emm1 = dreal(e(m-1))/scale
            sl = dreal(s(l))/scale
            el = dreal(e(l))/scale
            b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0d0
            c = (sm*emm1)**2
            shift = 0.0d0
            if (b .eq. 0.0d0 .and. c .eq. 0.0d0) go to 590
               shift = dsqrt(b**2+c)
               if (b .lt. 0.0d0) shift = -shift
               shift = c/(b + shift)
  590       continue
            f = (sl + sm)*(sl - sm) + shift
            g = sl*el
!c
!c           chase zeros.
!c
            mm1 = m - 1
            do 600 k = l, mm1
               call drotg(f,g,cs,sn)
               if (k .ne. l) e(k-1) = dcmplx(f,0.0d0)
               f = cs*dreal(s(k)) + sn*dreal(e(k))
               e(k) = cs*e(k) - sn*s(k)
               g = sn*dreal(s(k+1))
               s(k+1) = cs*s(k+1)
               if (wantv) call zdrot(p,v(1,k),1,v(1,k+1),1,cs,sn)
               call drotg(f,g,cs,sn)
               s(k) = dcmplx(f,0.0d0)
               f = cs*dreal(e(k)) + sn*dreal(s(k+1))
               s(k+1) = -sn*e(k) + cs*s(k+1)
               g = sn*dreal(e(k+1))
               e(k+1) = cs*e(k+1)
               if (wantu .and. k .lt. n)  call zdrot(n,u(1,k),1,u(1,k+1),1,cs,sn)
  600       continue
            e(m-1) = dcmplx(f,0.0d0)
            iter = iter + 1
         go to 650
!c
!c        convergence.
!c
  610    continue
!c
!c           make the singular value  positive
!c
            if (dreal(s(l)) .ge. 0.0d0) go to 620
               s(l) = -s(l)
               if (wantv) call zscal(p,(-1.0d0,0.0d0),v(1,l),1)
  620       continue
!c
!c           order the singular value.
!c
  630       if (l .eq. mm) go to 640
!c           ...exit
               if (dreal(s(l)) .ge. dreal(s(l+1))) go to 640
               t = s(l)
               s(l) = s(l+1)
               s(l+1) = t
               if (wantv .and. l .lt. p)   call zswap(p,v(1,l),1,v(1,l+1),1)
               if (wantu .and. l .lt. n)   call zswap(n,u(1,l),1,u(1,l+1),1)
               l = l + 1
            go to 630
  640       continue
            iter = 0
            m = m - 1
  650    continue
      go to 400
  660 continue
      return
      end subroutine
      





	subroutine drotg(da,db,c,s)
!c
!c     construct givens plane rotation.
!c     jack dongarra, linpack, 3/11/78.
!c
      double precision da,db,c,s,roe,scale,r,z
!c
      roe = db  
      if( dabs(da) .gt. dabs(db) ) roe = da
      scale = dabs(da) + dabs(db)
      if( scale .ne. 0.0d0 ) go to 10
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         z = 0.0d0
         go to 20
   10 r = scale*dsqrt((da/scale)**2 + (db/scale)**2)
      r = dsign(1.0d0,roe)*r
      c = da/r
      s = db/r
      z = 1.0d0
      if( dabs(da) .gt. dabs(db) ) z = s
      if( dabs(db) .ge. dabs(da) .and. c .ne. 0.0d0 ) z = 1.0d0/c
   20 da = r
      db = z
      return
      end subroutine
      




      DOUBLE PRECISION FUNCTION DZNRM2( N, X, INCX )
     
!*     .. Scalar Arguments ..
      INTEGER                           INCX, N
!*     .. Array Arguments ..
      COMPLEX*16                        X( * )
!*     ..
!*
!*  DZNRM2 returns the euclidean norm of a vector via the function
!*  name, so that
!*
!*     DZNRM2 := sqrt( conjg( x' )*x )
!*
!*
!*
!*  -- This version written on 25-October-1982.
!*     Modified on 14-October-1993 to inline the call to ZLASSQ.
!*     Sven Hammarling, Nag Ltd.
!*
!*
!*     .. Parameters ..
      DOUBLE PRECISION      ONE         , ZERO
      PARAMETER           ( ONE = 1.0D+0, ZERO = 0.0D+0 )
!*     .. Local Scalars ..
      INTEGER               IX
      DOUBLE PRECISION      NORM, SCALE, SSQ, TEMP
!*     .. Intrinsic Functions ..
      INTRINSIC             ABS, DIMAG, DBLE, SQRT
!*     ..
!*     .. Executable Statements ..
      IF( N.LT.1 .OR. INCX.LT.1 )THEN
         NORM  = ZERO
      ELSE
         SCALE = ZERO
         SSQ   = ONE
!*        The following loop is equivalent to this call to the LAPACK
!*        auxiliary routine:
!*        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
!*
         DO 10, IX = 1, 1 + ( N - 1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO )THEN
               TEMP = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP )THEN
                  SSQ   = ONE   + SSQ*( SCALE/TEMP )**2
                  SCALE = TEMP
               ELSE
                  SSQ   = SSQ   +     ( TEMP/SCALE )**2
               END IF
            END IF
   10    CONTINUE
         NORM  = SCALE * SQRT( SSQ )
      END IF

      DZNRM2 = NORM
      RETURN
!*
!*     End of DZNRM2.
!*
      end function
      
      
      subroutine zaxpy(n,za,zx,incx,zy,incy)
!c
!c     constant times a vector plus a vector.
!c     jack dongarra, 3/11/78.
!c     modified 12/3/93, array(1) declarations changed to array(*)
!c
      double complex zx(*),zy(*),za
      integer i,incx,incy,ix,iy,n
      !double precision dcabs1
      if(n.le.0)return
      if (dcabs1(za) .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1)go to 20
!c
!c        code for unequal increments or equal increments
!c          not equal to 1
!c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        zy(iy) = zy(iy) + za*zx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!c
!c        code for both increments equal to 1
!c
   20 do 30 i = 1,n
        zy(i) = zy(i) + za*zx(i)
   30 continue
      return
      end subroutine
      
      
      
      double complex function zdotc(n,zx,incx,zy,incy)
!c
!c     forms the dot product of a vector.
!c     jack dongarra, 3/11/78.
!c     modified 12/3/93, array(1) declarations changed to array(*)
!c
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n
      ztemp = (0.0d0,0.0d0)
      zdotc = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!c
!c        code for unequal increments or equal increments
!c          not equal to 1
!c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + dconjg(zx(ix))*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotc = ztemp
      return
!c
!c        code for both increments equal to 1
!c
   20 do 30 i = 1,n
        ztemp = ztemp + dconjg(zx(i))*zy(i)
   30 continue
      zdotc = ztemp
      return
      end function
      
      
      subroutine  zdrot (n,zx,incx,zy,incy,c,s)
!c
!c     applies a plane rotation, where the cos and sin (c and s) are
!c     double precision and the vectors zx and zy are double complex.
!c     jack dongarra, linpack, 3/11/78.
!c
      double complex zx(1),zy(1),ztemp
      double precision c,s
      integer i,incx,incy,ix,iy,n
!c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!c
!c       code for unequal increments or equal increments not equal
!c         to 1
!c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = c*zx(ix) + s*zy(iy)
        zy(iy) = c*zy(iy) - s*zx(ix)
        zx(ix) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!c
!c       code for both increments equal to 1
!c
   20 do 30 i = 1,n
        ztemp = c*zx(i) + s*zy(i)
        zy(i) = c*zy(i) - s*zx(i)
        zx(i) = ztemp
   30 continue
      return
      end subroutine
      subroutine  zscal(n,za,zx,incx)
!c
!c     scales a vector by a constant.
!c     jack dongarra, 3/11/78.
!c     modified 3/93 to return if incx .le. 0.
!c     modified 12/3/93, array(1) declarations changed to array(*)
!c
      double complex za,zx(*)
      integer i,incx,ix,n

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!c
!c        code for increment not equal to 1
!c
      ix = 1
      do 10 i = 1,n
        zx(ix) = za*zx(ix)
        ix = ix + incx
   10 continue
      return
!c
!c        code for increment equal to 1
!c
   20 do 30 i = 1,n
        zx(i) = za*zx(i)
   30 continue
      return
      end subroutine
      subroutine  zswap (n,zx,incx,zy,incy)
!c
!c     interchanges two vectors.
!c     jack dongarra, 3/11/78.
!c     modified 12/3/93, array(1) declarations changed to array(*)
!c
      double complex zx(*),zy(*),ztemp
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!c
!c       code for unequal increments or equal increments not equal
!c         to 1
!c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = zx(ix)
        zx(ix) = zy(iy)
        zy(iy) = ztemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!c
!c       code for both increments equal to 1
   20 do 30 i = 1,n
        ztemp = zx(i)
        zx(i) = zy(i)
        zy(i) = ztemp
   30 continue
      return
      end subroutine

      double precision function dcabs1(z)
      double complex z,zz
      double precision t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1 = dabs(t(1)) + dabs(t(2))
      return
      end function

end module  Numerical_Recipes 






