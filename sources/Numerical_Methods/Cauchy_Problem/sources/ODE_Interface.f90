module ODE_Interface
    
implicit none 
  
abstract interface 

  function ODES( U, t) 
    
         real :: U(:), t 
         real :: ODES( size(U) ) 
      
  end function 
    
end interface  

end module 
    