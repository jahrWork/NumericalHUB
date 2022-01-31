module Dependencies2D  
 
    use Dependencies
    implicit none
    
    interface 
	       
      real function BVP2D(x, y, u, ux, uy, uxx, uyy, uxy)
            real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
      end function
      
        function BVP2D_system(x, y, u, ux, uy, uxx, uyy, uxy)
            real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
            real :: BVP2D_system(size(u)) 
        end function
   
    
    end interface
    
	
    
 
contains



function BVP2D_Dependencies( Differential_operator ) result(D)
    procedure (BVP2D) :: Differential_operator
    logical :: D(5)
    
    integer, parameter :: Nd = 5, Nv = 1
    logical :: Dj(Nv, Nd)
    real :: x0, y0, U(Nv)
   
    call random_number(x0) 
    call random_number(y0)
    call random_number(U) 
    
     Dj =  is_Jacobian_zero( Nv, Nd, F )  
     D = Dj(1, :)
     write(*,*) " BVP2D dependencies =",  D 
  
contains

function F(X)
real, intent(in) :: X(:, :) 
real :: F( size(X, dim=1) ) 
      
!                              x, y,    u,    ux,       uy,    uxx,   uyy,    uxy 
   F = Differential_operator( x0,  y0, U(1), X(1,1), X(1,2), X(1,3), X(1,4), X(1,5) ) 

end function
end function
  



function BVP2D_system_Dependencies( Differential_operator, Nv ) result(D)
    procedure (BVP2D_system) :: Differential_operator
    integer, intent(in) :: Nv 
    logical :: D(Nv, 5)
    
    integer :: Nd = 5, i 
    real :: x0, y0, U(Nv)
   
    call random_number(x0) 
    call random_number(y0)
    call random_number(U) 
    
       
     D =  is_Jacobian_zero( Nv, Nd, F )  
     do i=1, Nv 
       write(*,*) "variable =", i," BVP2D system dependencies =",  D(i, :)  
     end do  
  
contains

function F(X)
real, intent(in) :: X(:, :) 
real :: F( size(X, dim=1) ) 
   
  
!                              x, y,   u,  ux,       uy,    uxx,   uyy,    uxy 
   F = Differential_operator( x0,  y0, U, X(:,1), X(:,2), X(:,3), X(:,4), X(:,5) )  
   

end function
end function



 
end module 