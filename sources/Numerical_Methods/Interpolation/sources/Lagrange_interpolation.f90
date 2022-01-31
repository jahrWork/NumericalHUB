module Lagrange_interpolation

implicit none 
public ::               & 
  Lagrange_polynomials, & ! Lagrange polynomial at xp from (x_i, y_i)
  Lebesgue_functions      ! Lebesgue function at xp from x_i    
         

contains 


!*******************************************************************************************************************************************
!                                           Lagrange polynomials
!
!  Recurrence relation  to obtain derivatives and integral values of the Lagrange polynomials at xp from a set of nodes x(0:N)
!  
!
!      lagrange_j{N+1} (x) = (x-x0)(x-x1)(x-2)........(x-xn) / (xj- x0)(xj-x1)(xj-x2).......(xj-xn), j=0...N
!      lagrange_j{N+1} (x) = (x -xn)/(xj-xn) lagrange_j{N} (x) 
!
!      d^k lagrange_j{r+1} (x) /dx^k  = ( d^k lagrange_j{r} (x) /dx^k (x-xr) + k d^(k-1) lagrange_j{r} (x) /dx^(k-1) ) / (xj -xr ), r=0...N  
!
!      integral( lagrange_j{N+1}(x), from x0 to xp ) = integral ( l(xp) + dl/dx(xp) ( x -xp) +....... d^N l/dx^N(xp) * (x- xp)**N / N! ) 
!                                                    = sum_k (  -  dl^k/dx^k (xp) ( xj - xp )**(k+1) /(k+1)! )  
! 
!      k = -1 means integral value of lagrange_j{N+1} (x) from xj to xp (xj  nearest node from xp)
!      k =  0 means value of lagrange_j{N+1} (x) at xp 
!      k > 1    means derivative value of lagrange_j{N+1} (x) at xp 
!
! It returns Lagrange_polynomials(k,j) 
!
!     k: meaning given above 
!     j: corresponding wth L_j(x) (this polynomial equals to one in x_j and vanishes at the rest of nodes) 
!
! Author : Juan A Hernandez (juanantonio.hernandez@upm.es) 
!*****************************************************************************************************************************************
pure function Lagrange_polynomials( x, xp ) 
   real, intent(in) :: x(0:), xp
   real Lagrange_polynomials(-1:size(x)-1,0:size(x)-1) 


   integer :: j   ! node 
   integer :: r   ! recursive index 
   integer :: k   ! derivative 
   integer :: Nk  ! maximum order of the derivative 
   integer :: N, j1(1), jp 

   real :: d(-1:size(x)-1) 
   real :: f 

   Nk = size(x) - 1 
   N  = size(x) - 1 

   do j = 0, N 
      d(-1:Nk) = 0 
      d(0) = 1

 ! ** k derivative of lagrange(x) at xp 
      do  r = 0, N 

         if (r/=j) then 
          do k = Nk, 0, -1
            d(k) = ( d(k) *( xp - x(r) ) + k * d(k-1) ) /( x(j) - x(r) )  
          end do 
        endif 

      enddo 

!  ** integral of lagrange(x) form x(jp) to xp 
      f = 1 
      j1 = minloc( abs(x - xp) ) - 2
      jp = max(0, j1(1)) 
     
      do k=0, Nk 
         f = f * ( k + 1 ) 
         d(-1) = d(-1) -  d(k) * ( x(jp) - xp )**(k+1) / f 
      enddo 

      Lagrange_polynomials(-1:Nk, j ) = d(-1:Nk) 

  end do 

  end function 

!******************************************************************************
!* Vector of Stencils of each interpolant 
!
!    0,0,...,0,0,1,2,3,4,5,6,..... N-Order, N-Order, ..N - Order 
!
!******************************************************************************
function Stencilv(Order, N) result(S)

 integer, intent(in) :: Order, N
 integer :: S(0:N) 

   integer :: i, N1, N2; 

   if (mod(Order,2)==0) then 
                     N1 = Order/2;       
                     N2 = Order/2; 
   else          
                  N1 = (Order-1)/2;   
                  N2 = (Order+1)/2; 
   endif 

   S(0:N1-1) = 0;  
   S(N1:N-N2) = [ ( i, i=0, N-N1-N2 ) ]; 
   S(N-N2+1:N) = N - Order;


end function 



!***************************************************************************************************************************************
!                                           Lagrange error polynomials
!
!  Recurrence relation  to obtain derivatives and integral values of the Lagrange error polynomials at xp from a set of nodes x(0:N)
!  
!
!      pi{N+1} (x) = (x-x0)(x-x1)(x-2)........(x-xn) 
!      pi{N+1} (x) = (x -xn) pi{N} (x) 
!
!      d^k pi{r+1} (x) /dx^k  = ( d^k pi{r} (x) /dx^k (x-xr) + k d^(k-1) pi{r} (x) /dx^(k-1) ), r=0...N  
!
!      k =  0 means value of pi{N+1} (x) at xp 
!      k > 1    means derivative value of pi{N+1} (x) at xp 
!
!*****************************************************************************************************************************************
function Lagrange_error_polynomial( x, xp ) 
   real, intent(in) :: x(0:), xp
   real Lagrange_error_polynomial(0:size(x)) 

  
   integer :: i   ! recursive index 
   integer :: k   ! derivative 
   integer :: Nk  ! maximum order of the derivative 
   integer :: N 
 
   real, allocatable :: PI(:,:) ! first index derivative, order or number of points 
   
   
      N  = size(x) - 1 
      Nk = N+1   
   
      allocate ( PI(-1:N+1, 0:N+1) ) 
 
   
      PI  = 0 
      PI(0, 0) = 1

 ! ** k derivative of pi(x) at xp 
      do  i = 0, N 
              do k = N+1, 0, -1
                               PI(k, i+1) =  PI(k, i) *( xp - x(i) ) + k * PI(k-1, i) 
              end do 
      enddo 
      
      Lagrange_error_polynomial(0:N+1) = PI(0:N+1, N+1) 
      
      deallocate( PI ) 

end function 

!*******************************************************************************************************************************************
!                                           Lebesgue function 
!
!      Definition: Lebesgue_function =  |Lagrange_0(x)| + |Lagrange_1(x)| + ....  + |Lagrange_N(x)|
!  
!      Inputs: 
!               x(0:N) nodes 
!               xp(0:M) points where the Lebesgue function is evaluated 
!
!       Outputs: 
!               Lebesgue_functions(k,i)
!                    k = -1 means integral value of lagrange_j{N+1} (x) from xj to xp (xj  nearest node from xp)
!                    k =  0 means value of lagrange_j{N+1} (x) at xp 
!                    k > 1    means derivative value of lagrange_j{N+1} (x) at xp 
!
!                    i: index of the xp point where lebesgue function is evaluated 
!
! Author : Juan A Hernandez (juanantonio.hernandez@upm.es) 
!*****************************************************************************************************************************************
pure function Lebesgue_functions( x, xp ) 
   real, intent(in) :: x(0:), xp(0:)
   real Lebesgue_functions(-1:size(x)-1, 0:size(xp)-1) 

   
    integer :: i, j, k, N, M    
    real, allocatable :: Lg(:, :, :) 
    
    N = size(x) - 1 
    M = size(xp) - 1  
    
    allocate ( Lg(-1:N, 0:N, 0:M) ) 
   
    
    do i=0, M   
      Lg(:, :, i) = Lagrange_polynomials( x, xp(i) )
    end do 
 
    Lebesgue_functions = 0 
    
    do k=-1, N 
      do j=0, N 
       Lebesgue_functions(k,:) =  Lebesgue_functions(k,:) + abs( Lg(k,j,:) )  
      end do 
    end do  
   
    deallocate( Lg ) 

end function 

!*******************************************************************************************************************************************
!                                           PI error polynomial 
!
!      Definition: pi(x) =  |(x-x0)| |(x-x1)| ..... |(x-xN)|
!  
!      Inputs: 
!               x(0:N) nodes 
!               xp(0:M) points where the pi error polynomial function is evaluated 
!
!       Outputs: 
!               PI_error_polynomial(k,i)
!                    k = -1 means integral value of lagrange_j{N+1} (x) from xj to xp (xj  nearest node from xp)
!                    k =  0 means value of lagrange_j{N+1} (x) at xp 
!                    k > 1    means derivative value of lagrange_j{N+1} (x) at xp 
!
!                    i: index of the xp point where error polynomial is evaluated 
!
! Author : Juan A Hernandez (juanantonio.hernandez@upm.es) 
!*****************************************************************************************************************************************
function PI_error_polynomial( x, xp ) 
   real, intent(in) :: x(0:), xp(0:)
   real PI_error_polynomial(0:size(x)-1, 0:size(xp)-1) 

   
    integer :: i, N, M    
    
    N = size(x) - 1 
    M = size(xp) - 1  
    
    
    do i=0, M  
         PI_error_polynomial(:, i) =  abs( Lagrange_error_polynomial( x, xp(i) ) ) 
         !write(*,*) " xp = ", xp(i) 
         ! write(*,*) " x = ", x
         !write(*,*) " PI = ", PI_error_polynomial(0,i)
         !stop 
         !read(*,*) 
    end do 
 
  
end function 



end module 
