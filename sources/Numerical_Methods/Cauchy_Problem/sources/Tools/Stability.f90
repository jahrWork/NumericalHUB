module Stability

   use Numerical_Recipes  
   use ODE_Interface
   implicit none

   
abstract interface 

 function F_R_RN(N, t) 
         integer, intent(in) :: N 
         real, intent(in) ::  t 
         real :: F_R_RN(N) 
      
 end function 

end interface 
   
   
   
contains

!***********************************************************************
!*  If gives the linearized operator F(U,t) in some point U0 
!***********************************************************************
function System_matrix( F, U0, t ) result(A)
             procedure (ODES) :: F 
             real, intent(in) :: U0(:), t 
             real :: A( size(U0), size(U0) ) 
          
          
          real ::  delta( size(U0) ), eps = 1d-6 
          integer :: j, N 
        
          N = size(U0) 
          
          do j=1, N 
              
              delta = 0 
              delta(j) = eps 
              A(:, j) = ( F( U0 + delta, t ) - F( U0 - delta, t ) )/(2*eps)
              
          end do 
          
             
end function  
        


function Eigenvalues_Jacobian( Differential_operator, U0, t0) result(lambda)
            procedure (ODES) :: Differential_operator 
            real, intent(in)  :: U0(:), t0  
            complex :: lambda( size(U0) ) 
          
       
          real, allocatable :: A(:, :)
          integer ::  N  
          
           
          N = size( U0 )
          allocate( A(N,N) )
          
          A = System_matrix( Differential_operator, U0, t0) 
          call Eigenvalues_QR( A, lambda )
          
             
end function 
                            
end module 

