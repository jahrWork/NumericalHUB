module Burgers 
  

use Cauchy_Problem
use plots
use Legendre_points
use Lagrange_interpolation
implicit none 


interface 

 function f_R_R(U) result(F) 
     real, intent(in) :: U 
     real :: F 
 end function 
 
end interface 



contains 
    
subroutine  Burgers_upwinding

  integer, parameter :: N = 50, Nt = 100 
  real :: U(0:N+1), F(0:N+1), RF(0:N), dUdt(0:N+1), x(0:N+1)  
  integer :: i, j 
  real :: dx, dt, Uc  
  real, parameter :: PI = 4 * atan(1.) 
  
  dx = 2./(N+1); dt = 1./ Nt  
  x    = [ (dx*i, i=0, N+1 ) ]
  U = sin( PI * x ) 
  
  call plot(x, U, "Burgers") 
  
! Time evolution  
  dUdt = 0 
  do i= 0, Nt 
    F = U**2 / 2 
  
  ! Upwinding  by means of Riemann solver   
    do j=0, N 
      Uc = ( U(j+1) + U(j) )/2  
      RF(j) = ( F(j+1) + F(j) )/2 - abs(Uc) * ( U(j+1) - U(j) )/2   
    end do 
  
  ! Finite volume   
    do j=1, N 
      dUdt(j) = -( RF(j) - RF(j-1) )/dx   
    end do 
  
  ! Euler 
    U = U + dt * dUdt 
    
  end do 
  call plot(x, U, "Burgers") 

end subroutine 
    

subroutine  Burgers_Cauchy 
    integer :: i 
    integer, parameter :: Nt = 100, Nx= 50 
    real :: t0 = 0, tf = 1, Time(0:Nt), U(0:Nt, 0:Nx+1), x(0:Nx+1)
    real, parameter :: PI = 4 * atan(1.)
    integer :: j 
    real :: dx = 2./(Nx+1)
    
    Time = [ (t0 + (tf -t0 )*i/real(Nt), i=0, Nt ) ]
    x    = [ (dx*i, i=0, Nx+1 ) ]
       
    U(0,:) = sin( PI * x ) 
    
    call Cauchy_ProblemS( Time,  Burgers_eq, U )
    
    call plot(x, U(0,:), "Burgers t=0") 
    call plot(x, U(Nt,:), "Burgers t=1") 
  
contains

function Burgers_eq( U, t )  result(F) 
    real :: U(:), t 
    real :: F(size(U)) 
    
    integer :: N
    
    N = size(U)-1 
    F = Burgers_discretization(N, U )
    
     
end function

function Burgers_discretization( N, U )  result(F) 
    integer, intent(in) :: N 
    real, intent(in) :: U(0:N+1)
    real :: F(0:size(U)-1) 
    
    integer :: j
    real, allocatable :: Riemann_Flux(:) 
   
    allocate( Riemann_Flux(0:N) ) 

    do j=0, N 
      Riemann_Flux(j) = Riemann( F_Burgers, U(j), U(j+1) )  
    end do 
  
    F = 0 
    do j=1, N 
      F(j) = -( Riemann_Flux(j) - Riemann_Flux(j-1) )/dx   
    end do 
      
     
end function


end subroutine



function Riemann(Flux, UL, UR) result(F) 
    procedure (F_R_R) :: Flux 
    real, intent(in) :: UL, UR 
    real :: F 
    
    real :: Uc 
    
      Uc = ( UL + UR )/2  
      
      F = ( Flux(UL) + Flux(UR) )/2 - abs(Uc) * ( UR - UL )/2  
      
end function 

elemental function F_Burgers(U) result(F) 
    real, intent(in) :: U
    real :: F 
    
      F = U**2 / 2   

      
end function

elemental function F_Burgers2(U) result(F) 
    real, intent(in) :: U
    real :: F 
    
      F = U   

      
end function  


!***********************************
! Burgers: Discontinous Galerkin 
!************************************
subroutine  Burgers_DG 

    integer, parameter :: Nt = 10000, Ne= 150, p = 3
    real, parameter :: t0=0, tf=1 
    real :: Time(0:Nt), U(0:Nt, 0:p, 1:Ne), x(0:p, 1:Ne), Flux(0:p, 1:Ne), dUdt(0:p,1:Ne)
    real :: F_Riemann(0:Ne)
    real :: alpha(0:p)
    real, parameter :: PI = 4 * atan(1.)
    integer :: i, j, k, m 
    real :: dx = 2./Ne, dt = (tf-t0)/Nt 
    real :: x0, xf, x_LGL(0:p), L0(-1:p,0:p), L1(-1:p,0:p) 
    real :: U0, U1 
    
    
  ! Lagrange polynomials  
  !-1:p (integral, lagrange, derivatives) 
  ! 0:p ( L_0(x), L_1(x),... L_N(x)     ) 
  ! 0:p ( points where L_j(x) is evaluated  ) 
  ! 1:Ne different elements  
    real :: Lg( -1:p, 0:p, 0:p, 1:Ne)
    
   
    Time = [ (t0 + dt*i, i=0, Nt ) ]
    call Legendre_Gauss_Lobatto_points(p, x_LGL, alpha)
    
    
    do j=1, Ne
       
       x0 = dx*(j-1)
       xf = dx*j
       x(:,j) = (x0+xf)/2 + (xf-x0)/2 * x_LGL
       do m=0, p
         Lg( :, :, m, j) =  Lagrange_polynomials( x(:,j), x(m, j) )
       end do 
    end do 
    alpha = alpha *(xf-x0)/2 
    
    U(0, :, :) = sin( PI*x )
      
    dUdt = 0 
    do i= 0, Nt-1 
        
      
      U0 = sin (  6* PI * Time(i) ) 
      U1 = 0 
      
      F_Riemann(0) = Riemann( F_Burgers, U0, U(i,0,1)   ) 
      do j=1, Ne-1 
          F_Riemann(j) = Riemann( F_Burgers, U(i,p,j), U(i,0,j+1)   )
      end do
      F_Riemann(Ne) = Riemann( F_Burgers, U(i, p, Ne), U1 )  
      
      Flux(:,:) =  F_Burgers( U(i,:,:) )  
      do j=1, Ne 
          Flux(0,j) = F_Riemann(j-1)
          Flux(p,j) = F_Riemann(j)
      end do 
      
        
      do j=1, Ne 
          L1 =  Lagrange_polynomials( x(:,j), x(p, j) )
          L0 =  Lagrange_polynomials( x(:,j), x(0, j) )
          do m=0, p
             
              dUdt(m,j) = + L1(0,m) * Flux(p, j) - L0(0,m) * Flux(0,j)      &
                          - dot_product(  alpha * Flux(:,j), Lg(1, m, :, j)  ) 
              
              dUdt(m,j) = - dUdt(m,j) / alpha(m) 
          end do 
      end do
    
      U(i+1,:,:) = U(i,:,:) + dt * dUdt(:,:) 
     
          
    end do 

    call plot( x(0,:), U(0, 0, :), "Legendre t=0" )
    call plot( x(0,:), U(Nt, 0, :), "Legendre t=1" )
  
end subroutine



!***********************************
! Burgers: Discontinous Galerkin 
!************************************
subroutine  Burgers_DG_old 

    integer, parameter :: Nt = 10000, Ne= 50, p = 6 
    real, parameter :: t0=0, tf=1 
    real :: Time(0:Nt), U(0:Nt, 0:p, 1:Ne), x(0:p, 1:Ne), Flux(0:p, 1:Ne), dUdt(0:p,1:Ne)
    real :: F_Riemann(0:Ne)
    real :: alpha(0:p)
    real, parameter :: PI = 4 * atan(1.)
    integer :: i, j, k, m 
    real :: dx = 2./Ne, dt = (tf-t0)/Nt 
    real :: x0, xf, x_LGL(0:p) 
    
    
  ! Lagrange polynomials  
  !-1:p (integral, lagrange, derivatives) 
  ! 0:p ( L_0(x), L_1(x),... L_N(x)     ) 
  ! 0:p ( points where L_j(x) is evaluated  ) 
  ! 1:Ne different elements  
    real :: Lg( -1:p, 0:p, 0:p, 1:Ne)
    
   
    
    Time = [ (t0 + dt*i, i=0, Nt ) ]
    call Legendre_Gauss_Lobatto_points(p, x_LGL, alpha)
    
    
    do j=1, Ne
       
       x0 = dx*(j-1)
       xf = dx*j
       x(:,j) = (x0+xf)/2 + (xf-x0)/2 * x_LGL
       do m=0, p
         Lg( :, :, m, j) =  Lagrange_polynomials( x(:,j), x(m, j) )
       end do 
    end do 
    alpha = alpha *(xf-x0)/2 
    
    U(0, :, :) = sin( PI*x )
      
    dUdt = 0 
    do i= 0, Nt-1 
        
      Flux(:,:) =  F_Burgers( U(i,:,:) )   
      
      !do j=1, Ne
      !    if (j==1) then 
      !        Flux(0,j) = Riemann( F_Burgers, 0., U(i,0,j)   )
      !    else 
      !        Flux(0,j) = Riemann( F_Burgers, U(i,p,j-1), U(i,0,j)   )
      !    end if 
      !              
      !   if (j==Ne) then 
      !       Flux(p,j) = Riemann( F_Burgers, U(i,p,j), 0. ) 
      !   else 
      !       Flux(p,j) = Riemann( F_Burgers, U(i,p,j), U(i,0,j+1) )
      !   end if 
      !end do 
      !
      !!do j=1, Ne
      !!    do m=0, p
      !!        dUdt(m,j) =   Lg(j, 0, m, p) * Flux(p,j) & 
      !!                    - Lg(j, 0, m, 0) * Flux(0,j) & 
      !!                    - dot_product( alpha * Lg(j, 1, m, :), Flux(:,j) ) 
      !!        dUdt(m,j) =  dUdt(m,j) / alpha(m) 
      !!    end do 
      !!end do
      !
      !
      !do j=1, Ne
      !    dUdt(0,j) = - ( Flux(p,j) - Flux(0,j) )/dx
      !    do m=1, p-1 
      !        dUdt(m,j) = - dot_product( alpha * Lg(j, 1, m,:), Flux(:,j) ) / alpha(m) 
      !    end do 
      !    dUdt(p,j) = - ( Flux(p,j) - Flux(0,j) )/dx
      !end do 
      
    
      
      
      
      
      
      
      F_Riemann(0) = Riemann( F_Burgers, 0., U(i,0,1)   ) 
      do j=1, Ne-1 
          F_Riemann(j) = Riemann( F_Burgers, U(i,p,j), U(i,0,j+1)   )
      end do
      F_Riemann(Ne) = Riemann( F_Burgers, U(i, p, Ne), 0. )  
      
         
        
      do j=1, Ne
          do m=0, p
              Flux(0,j) = F_Riemann(j-1)
              Flux(p,j) = F_Riemann(j)
              dUdt(m,j) =                                                 &
                          ! DG Form II
                          !+ Lg(0, p, m, j)*( F_Riemann(j)   -Flux(p,j) )  & 
                          !- Lg(0, 0, m, j)*( F_Riemann(j-1) -Flux(0,j) )  & 
                          !+ dot_product(  alpha * Flux(:,j), Lg(1, :, m, j)  )
                  
                          ! DG Form I 
                          + Lg(0, p, m, j)*( F_Riemann(j)  )  & 
                          - Lg(0, 0, m, j)*( F_Riemann(j-1) )  & 
                          - dot_product(  alpha * Flux(:,j), Lg(1, m, :, j)  ) 
              
              dUdt(m,j) = - dUdt(m,j) / alpha(m) 
          end do 
      end do
      
  !    write(*,*) " alpha_0 =", alpha(0), Lg(0,0,0,:)  
      
      
      !call plot( x(0,:), U(i, 0, :), "U" )
      !call plot( x(0,:), Flux(0, :), "Flux" )
      !
      !call plot( x(0,:), dUdt(0, :), "dUdt" )
      
  
      U(i+1,:,:) = U(i,:,:) + dt * dUdt(:,:) 
     
          
    end do 

    call plot( x(0,:), U(0, 0, :), "Legendre t=0" )
    call plot( x(0,:), U(Nt, 0, :), "Legendre t=1" )
    call plot( x(p,:), U(Nt, p, :), "Legendre t=1" )


end subroutine


end module
    