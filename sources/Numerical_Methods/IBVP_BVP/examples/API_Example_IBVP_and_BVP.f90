module API_Example_IBVP_and_BVP 


use IBVPs_and_BVPs 
!use Finite_differences
use Collocation_methods
use Linear_systems
use plots
implicit none 


contains 


subroutine Nonlinear_Plate_Vibrations

     integer, parameter :: Nx = 10, Ny = 10, Nt = 300, Nu = 3 , Nv = 2 
     real ::  x(0:Nx), y(0:Ny), Time(0:Nt)
     real :: U(0:Nt, 0:Nx, 0:Ny, Nu), V(0:Nt, 0:Nx, 0:Ny, Nv)
     real ::  x0 = -1, xf = 1, y0 = -1, yf = 1, t0 = 0, tf = 1
     integer :: i, j, Order = 6 
    
     real :: mu=10,  levels(20) = 0
     character(len=5):: legends(4)  = [ "(a)", "(b)", "(c)", "(d)" ] 
     character(len=100) :: path(4) = [   & 
           "./doc/chapters/IBVP_BVP/figures/NLvibrationsa", &
           "./doc/chapters/IBVP_BVP/figures/NLvibrationsb", &
           "./doc/chapters/IBVP_BVP/figures/NLvibrationsc", &
           "./doc/chapters/IBVP_BVP/figures/NLvibrationsd"    ] 
     
    Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]      
    x(0) = x0; x(Nx) = xf; y(0) = y0; y(Ny) = yf  
    call Grid_Initialization( "nonuniform", "x", x, Order ) 
    call Grid_Initialization( "nonuniform", "y", y, Order ) 
    
   U(0,:,:,1) =  Tensor_product( exp(-10 * x**2), exp(-10 * y**2) )  
   U(0,:,:,2) =  0
   do i=0, Nx; do j=0, Ny 
     U(0,i,j,3) =  U(0,i,j,1) * (-40 + (20*x(i))**2 + (20*y(j))**2 )
   end do; end do 
   
!  Nonlinear Plate Vibration    
   call IBVP_and_BVP( Time, x, y, Lu, Lv, BCu, BCv, U, V ) 
     
   do i=1, 4
    call plot_contour(x, y, U(50*(i-1),:,:,1), "$x$", "$y$",           & 
                      levels, legends(i), path(i), "isolines" ) 
   end do 
      
contains 

function Lu(x, y, t, u, ux, uy, uxx, uyy, uxy , v, vx, vy, vxx, vyy, vxy) 
 real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
 real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
 real :: Lu(size(u))
            
     real :: wxx, wyy, wxy, pxx, pyy, pxy 
     real :: w2, w2xx, w2yy,  w3xx, w3yy 
     
     wxx  = uxx(1); wyy  = uyy(1); wxy = uxy(1);
     w2xx = uxx(2); w2yy = uyy(2); w2  = u(2); 
     w3xx = uxx(3); w3yy = uyy(3);  
     pxx  = vxx(1); pyy  = vyy(1); pxy = vxy(1);
     
     
     Lu(1)  =    w2 
     Lu(2)  =  - w3xx - w3yy + load(x, y, t)              & 
               + mu * B( wxx, wyy, wxy, pxx, pyy, pxy)
     Lu(3)  =    w2xx + w2yy 
                 
end function 



real function B( wxx, wyy, wxy, pxx, pyy, pxy) 
  real, intent(in) ::  wxx, wyy, wxy, pxx, pyy, pxy
 
  B = wxx * pyy + wyy * pxx - 2 * wxy * pxy                                 
  
end function

real  function load(x, y, t) 
   real, intent(in) :: x, y, t 
   
   load =  0 
      
end function


function BCu(x, y, t, u, ux, uy) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BCu( size(u) ) 
       
        if (x==x0 .or. x==xf .or. y==y0 .or. y==yf) then
                           BCu = u
        else
             write(*,*)  "Error in BC1 "; stop 
             stop 
        endif

end function

function Lv( x, y, t, v, vx, vy, vxx, vyy, vxy , u, ux, uy, uxx, uyy, uxy) 
 real, intent(in) ::  x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
 real, intent(in) ::           v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
 real :: Lv(size(v))
       
     real :: wxx, wyy, wxy, pxx, pyy, pxy, Fxx, Fyy, Fxy, F  
     
     pxx = vxx(1); pyy = vyy(1); pxy = vxy(1);
     Fxx = vxx(2); Fyy = vyy(2); Fxy = vxy(2); F = v(1);
     wxx = uxx(1); wyy = uyy(1); wxy = uxy(1);
 
     Lv(1) = pxx + pyy - F  
     Lv(2) = Fxx + Fyy + B(wxx, wyy, wxy, wxx, wyy, wxy)
                 
end function 


function BCv(x, y, t, v, vx, vy)  
  real, intent(in) :: x, y,  t, v(:), vx(:), vy(:) 
  real :: BCv( size(v) ) 
    
        if (x==x0 .or. x==xf .or. y==y0 .or. y==yf) then
                           BCv = v
        else
             write(*,*)  "Error in BC2 "; stop 
             stop 
        endif

end function
end subroutine


end module 
