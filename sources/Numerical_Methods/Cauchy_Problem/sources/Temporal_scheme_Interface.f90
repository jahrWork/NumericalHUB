module Temporal_scheme_interface
  
 implicit none 
    
 abstract interface 
   subroutine Temporal_Scheme(F, t1, t2,  U1, U2, ierr )
                use ODE_Interface
                procedure (ODES) :: F 
                real, intent(in)    :: t1, t2
                real, intent(in) ::  U1(:)
                real, intent(out) ::  U2(:)
                integer, intent(out) :: ierr
  end subroutine 
end interface  

end module 
    