module ODE_Interface
    
 implicit none 
  
 real :: ZERO = 1d-14
 
abstract interface 

  function ODES( U, t) 
    
         real :: U(:), t 
         real :: ODES( size(U) ) 
      
  end function 
    
end interface  

end module 
    