module Boundary_layer
    
    use Cauchy_Problem
    use Temporal_Schemes
    use Linear_Systems
    use Collocation_methods 
    use plots
    use Interpolation
    use Initial_Boundary_Value_Problems
    use Chebyshev_interpolation
    use Boundary_value_problems
    implicit none

     real, parameter ::  x0 = 0, xf = 1, y0 = 0, yf = 0.2, xc = 0! xf/4 
     real ::   t0 = 0, tf =  5.
     real :: rho0 = 1.0
    
    
    contains
  
 !******************************************************************
 ! Boundary layer. Flat plate  
 ! Since this problem has a singularity such as sqrt(x), 
 ! it is not suitable for high order methods. 
 ! In this example, only second order is considered.  
 ! Higher orders become unstable giving exploding results.
 !******************************************************************   
 subroutine Compressible_boundary_layer_examples
 

 call Compressible_boundary_layer( grid="nonuniform", Nx=50, Ny=10, Nv=3, Nt=4000, qx=2, qy=2) 
 call Compressible_boundary_layer( grid="nonuniform",Nx=100, Ny=20, Nv=3, Nt=4000, qx=2, qy=2) 
 call Compressible_boundary_layer( grid="nonuniform",Nx=200, Ny=40, Nv=3, Nt=8000, qx=2, qy=2) 
                                                                                   
 call Compressible_boundary_layer( grid="chebyshev", Nx=50, Ny=10, Nv=3, Nt=64000,  qx=2, qy=2) 

 call  Compressible_boundary_layer_multigrid

end subroutine 
 
 
 
subroutine Compressible_boundary_layer( grid, Nx, Ny, Nv, Nt, qx, qy) 
   character(len=*), intent(in) :: grid 
   integer, intent(in) :: Nx, Ny, Nv, Nt, qx, qy   
   
     real, target ::  x(0:Nx), y(0:Ny), Time(0:Nt), W(0:Nt, 0:Nx, 0:Ny, 1:Nv )
     integer :: i, j
              
      Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     
      if (grid == "chebyshev" ) then 
          do i=0, Nx 
             x(i) = x0 + (xf-x0)*i**2/Nx**2
          end do 
           do j=0, Ny 
             y(j) = y0 + (yf-y0)*j**2/Ny**2
          end do 
          
          call Grid_Initialization( "unmodified", "x", x, qx )
          call Grid_Initialization( "unmodified", "y", y, qy )
           !write(*,*) x(0), x(1) 
           !stop 
      else 
          x = [ ( x0 + (xf-x0)*i/Nx, i=0, Nx) ]  
          y = [ ( y0 + (yf-y0)*j/Ny, j=0, Ny) ]
          call Grid_Initialization( "nonuniform", "x", x, qx )
          call Grid_Initialization( "nonuniform", "y", y, qy )
      end if 
      
     
     W(0, :, :, 1) = rho0
     W(0, :, :, 2:3) = 0      
     call Initial_Boundary_Value_Problem(Time, x, y, Layer, Layer_BC, W, Euler ) 
     call plotBL( x, y, W(Nt,:,:,:) ) 
    
   
end subroutine 


 
subroutine Compressible_boundary_layer_multigrid

   integer :: Nv, qx, qy, i, j  
   integer :: Nx1, Ny1, Nt1
   integer :: Nx2, Ny2, Nt2
   real, allocatable ::  x1(:), y1(:), x2(:), y2(:), Time1(:), Time2(:), W1(:,:,:,:), W2(:,:,:,:) 
        
      Nx1=50; Ny1=10;  Nv=3;  Nt1=4000;  qx=2;  qy=2; 
      allocate(  x1(0:Nx1), y1(0:Ny1), Time1(0:Nt1), W1(0:Nt1, 0:Nx1, 0:Ny1, 1:Nv ) ) 
      Time1 = [ (t0 + (tf-t0)*i/Nt1, i=0, Nt1 ) ]  
     
      x1 = [ ( x0 + (xf-x0)*i/Nx1, i=0, Nx1) ]  
      y1 = [ ( y0 + (yf-y0)*j/Ny1, j=0, Ny1) ]
      call Grid_Initialization( "nonuniform", "x", x1, qx )
      call Grid_Initialization( "nonuniform", "y", y1, qy )    
     
      W1(0, :, :, 1) = rho0
      W1(0, :, :, 2:3) = 0      
      call Initial_Boundary_Value_Problem(Time1, x1, y1, Layer, Layer_BC, W1, Euler ) 
      call plotBL( x1, y1, W1(Nt1,:,:,:) ) 
    
      Nx2=200; Ny2=40;  Nv=3;  Nt2=2000;  qx=2;  qy=2; 
      allocate(  x2(0:Nx2), y2(0:Ny2), Time2(0:Nt2), W2(0:Nt2, 0:Nx2, 0:Ny2, 1:Nv ) ) 
      tf = 5 
      Time2 = [ (t0 + (tf-t0)*i/Nt2, i=0, Nt2 ) ]  
     
      x2 = [ ( x0 + (xf-x0)*i/Nx2, i=0, Nx2) ]  
      y2 = [ ( y0 + (yf-y0)*j/Ny2, j=0, Ny2) ]
      do i=1, Nv  
        W2(0,:,:,i) = Interpolant( x1, y1, W1(Nt1,:,:,i), 4, x2, y2 )
        call plot_contour( x2, y2, W2(0,:,:,i), "$x$", "$y$" )
      end do 
          
      call Grid_Initialization( "nonuniform", "x", x2, qx )
      call Grid_Initialization( "nonuniform", "y", y2, qy )    
     
      call Initial_Boundary_Value_Problem(Time2, x2, y2, Layer, Layer_BC, W2, Euler ) 
      call plotBL( x2, y2, W2(Nt2,:,:,:) ) 
   
      
end subroutine 


 
 
function Layer(x, y, t, W, Wx, Wy, Wxx, Wyy, Wxy) result(F) 
           real,intent(in) :: x, y, t, W(:), Wx(:), Wy(:), Wxx(:), Wyy(:), Wxy(:) 
           real F( size(W) ) 
             
       real :: u, v, P, rho, rho_x, rho_y, ux, uy, vx, vy, uxx, uyy, vxx, vyy
       real :: Px, Py 
       real, parameter :: gam = 1.4, Re = 500, P0 =  0.1
           
       rho = W(1); rho_x = Wx(1); rho_y = Wy(1) 
       u = W(2); ux = Wx(2); uy = Wy(2); uxx = Wxx(2); uyy = Wyy(2)
       v = W(3); vx = Wx(3); vy = Wy(3); vxx = Wxx(3); vyy = Wyy(3) 
       
       Px = P0 * gam * rho**(gam-1) * rho_x 
       Py = P0 * gam * rho**(gam-1) * rho_y
            
       F(1) =  -( rho_x * u + rho_y * v + rho * (ux + vy) )
       F(2) = - u * ux - v * uy - Px/rho + 1/Re *( uxx + uyy ) / rho 
       F(3) = - u * vx - v * vy - Py/rho + 1/Re *( vxx + vyy ) / rho       
      
       F(1) = 0.01 * F(1) 
      
       
end function

!-------------------------------------------------------
function Layer_BC( x, y, t, W, Wx, Wy ) result (BC) 
     real, intent(in) :: x, y, t, W(:), Wx(:), Wy(:) 
     real :: BC( size(W) )

     real :: u, v, rho, uy, rho_y, rho_x, vy  
     
     rho = W(1); u = W(2); v = W(3); 
     rho_y = Wy(1); uy = Wy(2); vy = Wy(3); 
     rho_x = Wx(1);
     
  if (x==x0) then
                      BC = [rho-rho0, u-u0(y), v ]
       
  else if (y==y0 .and. x < xc) then
                      BC =[rho_y, uy, vy] 
                      
  else if (y==y0) then
                      BC =[rho_y, u, v]
                      
  else if (x==xf) then
                      BC = FREE_BOUNDARY_CONDITION
 
                      
  else if (y==yf ) then
                     BC = FREE_BOUNDARY_CONDITION
                
  else
            write(*,*)  "Error in Layer_BC (x,y) =", x, y; 
            write(*,*)  "x0, xf =", x0, xf
            write(*,*)  "y0, yf =", y0, yf
            stop
        end if
end function


real function u0(y) 
   real, intent(in) :: y 
   
   u0 = 1 ! - exp(-20*y) 
   
   
end function 

 

subroutine plotBL(x, y, W ) 
 real, intent(in) :: x(0:), y(0:), W(:,:, :) 

     
     integer :: i, j, Nx, Ny  
     real :: ig(1:size(x)-1), dx(1:size(x)-1)
     real, allocatable :: Wx(:, :) , Wy(:,:), U(:,:) 
    
     integer, parameter :: Mx=200, My=200 
     real :: xp(0:Mx), yp(0:My), Wp(0:Mx, 0:My) 
     
     
      Nx = size(x) -1; Ny = size(y)-1
     allocate( Wx(0:Nx, 0:Ny),  Wy(0:Nx, 0:Ny),  U(0:Nx, 0:Ny) ) 
     
     
     xp = [ (x0 + (xf-x0)*i/Mx, i=0, Mx ) ]
     yp = [ (y0 + (yf-y0)*i/My, i=0, My ) ]
      
     do i=1, size(x)-1
         ig(i) = i 
         dx(i) = x(i) - x(i-1)
 !        write(*,*) dx(i) 
     end do 
    
     
     do i=1, 3 
     call plot_contour( x, y, W(:,:,i), "$x$", "$y$" ) 
     !Wp = Interpolant( x, y, W(:,:,i), 4, xp, yp )
     !call plot_contour( xp, yp, Wp, "$x$", "$y$" )
     end do 
     
     !U =  W(:,:,1) * W(:,:,2) 
     !call Derivative( ["x", "y" ], 1, 1, U, Wx) 
     !U =  W(:,:,1) * W(:,:,3) 
     !call Derivative( ["x", "y" ], 2, 1, U, Wy) 
     !U = Wx + Wy 
     !call plot_contour( x, y, U, "$x$", "$y$" )
    
     
    !call plot(ig, dx, "grid" ) 
     
end subroutine 





 
  
end module










