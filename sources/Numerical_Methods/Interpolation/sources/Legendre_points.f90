module Legendre_points

    
    implicit none 
    
    
    
    
contains 
    
!**********************************************************************
! It determines N+1 Gauss Legendre Lobatto points
!    input: 
!            N (order or the polynomial)  
!    output: 
!            x(0:N) interpolation points 
!            alpha(0:N) weights of Gauss quadrature 
!                      
!  integral from -1 to +1 f(x) dx = sum from 0 to N ( alpha_j f(x_j) ) 
!
! Author: Juan A Hernandez (juanantonio.hernandez@upm.es)  June 2022
!**********************************************************************
subroutine Legendre_Gauss_Lobatto_points(N, x, alpha) 
  integer, intent(in) :: N 
  real, intent(out) :: x(0:N), alpha(0:N) 
 
 
 real, allocatable :: dx(:), y(:)
 real, allocatable :: L(:, :) ! first index: order of polynomial, second index: nodal point 
 real :: PI = 4 * atan(1.), eps = 1d-15 
 integer ::  k, M    
 
       
   M = N
   allocate(  dx(0:M), L(0:N+1,0:M) ) 
   x  = [ (-cos(PI*k/N), k=0, M) ]
   
   dx = 1 
   do while (norm2(dx)> eps)  
       
      L(0, :) = 1
      L(1, :) = x
      
   !  Three terms recurrence relation for Legendre polynomials    
      do k=1, N
          L(k+1,:) = (2*k+1)/(k+1.) * x * L(k, :) - k/(k+1.) * L(k-1, :) 
      end do 
    
   !  Newton method to determine interpolation points: 
   !   f(x) = L_{N+1}(x) - L_{N-1}(x) = 0  
   !   f'(x) = L'_{N+1}(x) - L'_{N-1}(x)  = (2N+1)  L_N(x)  
   !   x = x - f(x) / f'(x)    
      
      dx =  ( L(N+1,:) -  L(N-1,:) ) / (  (2*k+1) * L(N,:) ) 
      x  = x - dx 
      
   end do  
   
   alpha = 2./( N*(N+1) * L(N, :)**2 )
   
 
end subroutine 

!***********************************************************************************
!*  It computes dimensional factors for Chebyshev derivatives 
!*
! Authors : Juan A Hernandez (juanantonio.hernandez@upm.es)  June 2022
!***********************************************************************************
subroutine Legendre_Grid_Initialization( direction, nodes, alpha ) 
      character(len=*),  intent(in) ::  direction 
      real, intent(inout) :: nodes(0:)
      real, intent(out) :: alpha(0:) 
        
  integer :: N, i
  real :: x0, xf 
  
  
    N = size(nodes)-1
    x0 = nodes(0) 
    xf = nodes(N) 
    
    call Legendre_Gauss_Lobatto_points(N, nodes, alpha)   
    
    nodes = (x0+xf)/2 + (xf-x0)/2 * nodes 
    
end subroutine 




end module 