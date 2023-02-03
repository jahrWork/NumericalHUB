module Dependencies_BC  
 
    use Dependencies
    implicit none
	
  public :: FREE_BOUNDARY_CONDITION      ! to impose free BCs
  public :: PERIODIC_BOUNDARY_CONDITION  ! to impose periodic BCs
  public :: INTERFACE_CONDITION          ! to impose multiblocks 
  
  real, parameter :: FREE_BOUNDARY_CONDITION = 123456789.  
  real, parameter :: PERIODIC_BOUNDARY_CONDITION = 987654321. 
  real, parameter :: INTERFACE_CONDITION = 1122334455.
	
	
	 interface  
  
         function BC_IBVP2D_system(x, y, t, u, ux, uy) result(BC) 
            real, intent(in) :: x, y, t, u(:), ux(:), uy(:)
            real :: BC(size(u)) 
         end function 
	end interface 
	
    interface     
        function BC_IBVP2D(x, y, t, u, ux, uy) result(BC) 
            real, intent(in) :: x, y, t, u, ux, uy
            real :: BC 
        end function 
        
    end interface
	
	
    
 
contains


function BC_IBVP2D_Dependencies( Boundary_conditions, x0, xf, y0, yf ) result(D)
    procedure (BC_IBVP2D) :: Boundary_conditions ! x, y, t, u, ux, uy
    real, intent(in) :: x0, xf, y0, yf 
    logical :: D(2)
    
    integer, parameter :: Nd = 2, Nv = 1
    logical :: Dj(Nv, Nd)
    real :: BC1, BC2, z = 0., r  
    
     Dj =  is_Jacobian_zero( Nv, Nd, BC )  
     D = Dj(1, :)
     
     call random_number(r) 
     BC1 = Boundary_conditions( r, yf, z, z, z, z )
     BC2 = Boundary_conditions( xf, r, z, z, z, z )
     if ( BC1 == FREE_BOUNDARY_CONDITION .or. BC2 == FREE_BOUNDARY_CONDITION ) then 
         D = .true. 
     end if 
     
    
contains

function BC(X)
real, intent(in) :: X(:, :) 
real :: BC( size(X, dim=2) ) 
   
   real :: t0=0, U(Nv)
   U = 0 

   BC = Boundary_conditions( x0, y0, t0, U(1), X(1,1), X(1,2) ) 

end function
    
end function  

function BC_IBVP2D_Dependencies_system( Nv, Boundary_conditions, x0, xf, y0, yf ) result(D)
    integer, intent(in) :: Nv
    procedure (BC_IBVP2D_system) :: Boundary_conditions ! x, y, t, u, ux, uy
    real, intent(in) :: x0, xf, y0, yf 
    logical :: D(Nv, 2)
    
    integer :: Nd = 2  
    
     D =  is_Jacobian_zero( Nv, Nd, BC )  
    
contains

function BC(X)
real, intent(in) :: X(:, :) 
real :: BC( size(X, dim=2) ) 
   
   real :: t0=0, U(Nv)
   U = 0 

   BC = Boundary_conditions( x0, y0, t0, U(:), X(:,1), X(:,2) ) 

end function
    
end function  


 
end module 