    
module Dependencies_IBVP2D 
 
    use Dependencies
    implicit none
	

	
	interface  
    
        real function F_IBVP2D(x, y, t, u, ux, uy, uxx, uyy, uxy) 
            real, intent(in) :: x, y, t, u, ux, uy, uxx, uyy, uxy
        end function  
    
        function F_IBVP2D_system(x, y, t, u, ux, uy, uxx, uyy, uxy) result(L) 
            real, intent(in) :: x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
            real :: L(size(u)) 
        end function 
        
    end interface
    
 
    contains



function IBVP2D_Dependencies( Spatial_operator ) result(D)
    procedure (F_IBVP2D) :: Spatial_operator ! x, y, t, u, ux, uy, uxx, uyy, uxy 
    logical :: D(5) 
    
    integer, parameter :: Nv = 1 
    logical :: Dj(Nv,5) 
    integer :: Nd = 5 
    real :: x0, y0, t0, U(Nv) 
    
     call random_number(x0) 
     call random_number(y0)
     call random_number(t0)
     call random_number(U) 
    
    
     Dj =  is_Jacobian_zero( Nv, Nd, G )  
     D = Dj(1,:) 
     write(*,*) " IBVP2D_Dependencies : ", D
     
    
contains 

function G(X) 
real, intent(in) :: X(:, :) 
real :: G( size(X, dim=1) ) 

   G = Spatial_operator( x0, y0, t0, U(1), X(1,1), X(1,2), X(1,3), X(1,4), X(1,5) ) 

end function 

end function

    


function IBVP2D_Dependencies_system( Nv, Spatial_operator ) result(D)
    integer, intent(in) :: Nv  
    procedure (F_IBVP2D_system) :: Spatial_operator ! x, y, t, u, ux, uy, uxx, uyy, uxy 
    logical :: D(Nv,5) 
    
    integer :: Nd = 5, i  
    real :: x0, y0, t0, U(Nv) 
    
     call random_number(x0) 
     call random_number(y0)
     call random_number(t0)
     call random_number(U) 
    
    
     D =  is_Jacobian_zero( Nv, Nd, G )  
     
     do i=1, Nv 
       write(*,*) " IBVP2D_Dependencies_system : ", D(i,:)
     end do 
     
    
contains 

function G(X) 
real, intent(in) :: X(:, :) 
real :: G( size(X, dim=1) ) 

   G = Spatial_operator( x0, y0, t0, U(:), X(:,1), X(:,2), X(:,3), X(:,4), X(:,5) ) 

end function 

end function


end module   
    
