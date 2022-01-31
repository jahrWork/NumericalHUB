module Jacobian_module

  implicit none    
    
    
interface  
    function FunctionRN_RN(x) result(F) 
      real, intent(in) :: x(:) 
      real :: F( size(x) ) 
     end function 
end interface  



contains 

!**********************************************************************
!*  Jacobian of a vector of functions 
!*
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!**********************************************************************
function Jacobian( F, xp ) 
  procedure (FunctionRN_RN) :: F 
  real, intent(in) :: xp(:)
  real :: Jacobian( size(xp), size(xp) ) 

   integer ::  j, N  
   real :: xj( size(xp) ) 


    N = size(Xp) 

  
    do j = 1, N 
       xj = 0
       xj(j) = 1d-3
    !   Jacobian(:,j) = Directional_Gradient(F, xp, xj )  
       Jacobian(:,j) =  ( F(xp + xj) - F(xp - xj) )/norm2(2*xj);
    enddo 



end function 


!***********************************************************************
!* It computes the gradient of a function along the direction r(:) 
!***********************************************************************
function Directional_Gradient( F, xp, r ) 

     procedure (FunctionRN_RN) :: F 
     real, intent(in) :: xp(:), r(:) 
     real Directional_Gradient(size(xp)) 
  

  Directional_Gradient = ( F(xp + r) - F(xp-r) )/norm2(2*r); 

end function 




end module 
