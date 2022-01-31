module Jacobian_module

  implicit none    
    
    
interface  
    function FunctionRN_RN(x) result(F) 
      real, intent(in) :: x(:) 
      real :: F( size(x) ) 
     end function 
end interface  

interface  
    function FunctionRN_RNE(x) result(F) 
      real, target :: x(:) 
      real :: F( size(x) ) 
     end function 
end interface 

interface  
    function FunctionRN_R(x) result(F) 
      real, intent(in) :: x(:) 
      real :: F 
     end function 
end interface

interface  

       subroutine  NonLinearSolver(F, x)
          import :: FunctionRN_RN
          procedure(FunctionRN_RN) :: F
          real, intent(inout) :: x(:)
       end subroutine  

end interface  


contains 

!**********************************************************************
!*  Jacobian of a vector of functions F ( x in ) 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!**********************************************************************
function Jacobian( F, xp ) 
  procedure (FunctionRN_RN) :: F 
  real, intent(in) :: xp(:)
  real :: Jacobian( size(xp), size(xp) ) 

   integer ::  j, N  
   real :: xj( size(xp) ), Dx = 1d-3  
  
    N = size(xp) 
  
    do j = 1, N 
       xj = 0
       xj(j) = Dx
       Jacobian(:,j) =  ( F(xp + xj) - F(xp - xj) )/(2*Dx)       
    enddo 

end function 


!**********************************************************************
!*  Jacobian of a vector of functions F ( x inout )  
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!**********************************************************************
function Jacobianc( F, xp )
  procedure (FunctionRN_RNE) :: F 
  real, intent(in) :: xp(:)
  real :: Jacobianc( size(xp), size(xp) ) 

   integer ::  j, N  
   real :: xj( size(xp) ),  h = 1d-3 

    N = size(xp) 
    
    do j = 1, N 
       xj = 0
       xj(j) = h
       Jacobianc(:,j) =  ( F(xp + xj) - F(xp - xj) )/(2*h)
    enddo 
  
end function


end module 
