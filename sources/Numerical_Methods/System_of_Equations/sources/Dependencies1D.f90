module Dependencies1D  
 
    use Dependencies
    implicit none
    
    interface 
	  real function BVP1D(x, u, ux, uxx) 
            real, intent(in) :: x, u, ux, uxx                 
      end function 
       
      function BVP1D_system(x, u, ux,  uxx)
            real, intent(in) :: x, u(:), ux(:), uxx(:)
            real :: BVP1D_system(size(u)) 
      end function
      
      real function BVP2D(x, y, u, ux, uy, uxx, uyy, uxy)
            real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
      end function
      
   
    
    end interface
    
	
    
 
contains

function BVP1D_Dependencies( Differential_operator ) result(D)
    procedure (BVP1D) :: Differential_operator
    logical :: D(2)
    
    integer, parameter :: Nd = 2, Nv = 1
    logical :: Dj(Nv, Nd)
       
     Dj =  is_Jacobian_zero( Nv, Nd, F )  
     D = Dj(1, :)
     write(*,*) " BVP1D dependencies =", D 
    
  
contains

function F(X)
real, intent(in) :: X(:, :) 
real :: F( size(X, dim=1) ) 
   
   real :: x0, U(Nv)
   
   call random_number(x0) 
   call random_number(U) 

   F = Differential_operator( x0,  U(1), X(1,1), X(1,2) ) 

end function
end function  

function BVP1D_system_Dependencies( Differential_operator, Nv ) result(D)
    procedure (BVP1D_system) :: Differential_operator
    integer, intent(in) :: Nv 
    logical :: D(Nv, 2)
    
      integer :: i, Nd = 2 
      real :: x0, U(Nv)
   
     call random_number(x0) 
     call random_number(U) 
           
     D =  is_Jacobian_zero( Nv, Nd, F )
     do i=1, Nv 
       write(*,*) "variable =", i, " BVP1D system dependencies =",  D(i, :)  
     end do   
  
contains

function F(X)
real, intent(in) :: X(:, :) 
real :: F( size(X, dim=1) ) 
   

   F = Differential_operator( x0,  U, X(:,1), X(:,2) ) 

end function
end function      
    
 


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
  

 
end module 