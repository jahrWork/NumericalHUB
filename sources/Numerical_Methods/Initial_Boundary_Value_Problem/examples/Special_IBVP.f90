module Special_IBVP
    
use Linear_systems
use Initial_Boundary_Value_Problems
use Collocation_methods
use Temporal_Schemes
use Stability_regions
use Temporal_error 
use Stability
use plots 
use Non_Linear_Systems
use Chebyshev_interpolation
use Cauchy_Problem
use Dependencies_BC
implicit none 


contains 
 
subroutine Special_examples 


   call Two_solids_different_conductivity_2th_order
   call N_solids_with_different_conductivity
            
end subroutine 
 

subroutine Two_solids_different_conductivity_2th_order

       integer, parameter :: Nx = 20, Nt = 1000
       real ::  x1(0:Nx), x2(0:Nx)
       real :: Time(0:Nt), U1(0:Nt,0:Nx), U2(0:Nt,0:Nx)
       
       
       real ::  x0 = -1, xf = 1, xm = 0, t0 = 0, tf = 1 
       integer :: i, j
       real :: PI = 4 * atan(1.), dt, dx 
       real :: u_I, alpha = 0.1 
     
     dt =  (tf-t0)/Nt 
     Time = [ (t0 + dt*i, i=0, Nt ) ] 
     dx =  (xm-x0)/Nx 
     x1 = [ (x0 + dx*j, j=0, Nx ) ] 
     x2 = [ (xm + dx*j, j=0, Nx ) ]  
     
   ! Initial condition   
     U1(0,:) = sin(PI*x1) 
     U2(0,:) = sin(PI*x2) 
     
     do i=0, Nt-1 
        
     ! heat flux at interface (dT/dx)_1 = alpha (dT/dx)_2  
       u_I = ( 4*U1(i,Nx-1) - U1(i,Nx-2) + 4*alpha*U2(i,1) -alpha*U2(i,2)  & 
             ) / ( 3 + 3*alpha )
       
     ! BCs   
       U1(i,0:Nx:Nx) = [0., u_I]; 
       U2(i,0:Nx:Nx) = [u_I, 0.]; 
            
         
     ! Euler + 2 order FD Heat equation solid 1 and solid 2
       do j=1, Nx-1
          U1(i+1,j) = U1(i,j) +         dt/dx**2 *( U1(i,j+1) - 2 * U1(i,j) + U1(i,j-1) )
          U2(i+1,j) = U2(i,j) + alpha * dt/dx**2 *( U2(i,j+1) - 2 * U2(i,j) + U2(i,j-1) ) 
       end do 
     end do 
     
     call scrmod("reverse")
     call qplot([x1, x2], [ U1(Nt-1, :), U2(Nt-1, :)] , 2*(Nx+1) ) 
  
end subroutine 

subroutine N_solids_with_different_conductivity

   integer, parameter :: Nx = 20, Nt = 64000, Ns = 3 
   real :: Time(0:Nt), U(0:Nt, Ns*(Nx+1))
   character(len=2) :: grids(3) = ["x1", "x2", "x3"]
   real, target ::  x(Ns*(Nx+1)), Ux(Ns*(Nx+1)), Uxx(Ns*(Nx+1))  
   real, pointer:: xv(:, :), Uv(:, :), Uvx(:, :), Uvxx(:, :)
   
   real ::  x0 = -1, xf = 1, x1 = 0, x2 = 0.5
   real ::  x_interfaces(Ns+1)  
   integer :: i, i_interface, N_interfaces = Ns-1
   
   real :: t0 =0, tf = 1, t_spatial_discretization
   real :: PI = 4 * atan(1.) 
       
     
     x_interfaces = [ x0, x1, x2, xf ] 
     xv(0:Nx, 1:Ns) => x     
     Uvx(0:Nx, 1:Ns) => Ux  
     Uvxx(0:Nx, 1:Ns) => Uxx
     Time = [ ( t0 + (tf-t0)/Nt*i, i=0, Nt ) ] 
    
     do i=1, Ns 
         xv(:,i) = Chebyshev_extrema( Nx, x_interfaces(i),  x_interfaces(i+1) ) 
         call Grid_Initialization( "unmodified", grids(i), xv(:,i), Nx )
     end do  
   
     U(0,:) = 0 
     call Cauchy_ProblemS( Time, F, U, Euler )
   
     call plots 
  
contains

function F( U, t ) result(Fs) 
    real :: U(:), t 
    real :: Fs(size(U)) 
    
    t_spatial_discretization = t 
    Fs = Space_discretization( U, t) 
   
  
end function 
 

function Space_discretization( U, t ) result(F) 
                      real, target ::  U(:), t, F(size(U)) 
                      
     integer :: j, k  
     real :: u_I(1), U_b(2) 
   
       Uv(0:Nx, 1:Ns) => U 
     
     ! heat flux at interface (dT/dx)_1 = alpha (dT/dx)_2  
       do i_interface = 1, N_interfaces 
          u_I =  Uv(Nx, i_interface) 
          call Newton( Fluxes_at_interfaces, u_I) 
          Uv(Nx, i_interface)   = u_I(1)
          Uv(0,  i_interface+1) = u_I(1) 
       end do   
      
    ! Boundary conditions at x=x0 and x=xf    
      U_b = [  Uv(0, 1), Uv(Nx, Ns) ] 
      call Newton( Values_at_boundaries, U_b ) 
      Uv(0, 1)   = U_b(1)  
      Uv(Nx, Ns) = U_b(2) 
    
     ! calculate second derivative in every solid   
       do k=1, Ns      
         call Derivative( grids(k), 2, Uv(:,k), Uvxx(:,k) ) 
       end do 
     
     ! builds the heat equation ( Uvxx is pointed to Uxx )   
       do j  = 1, size(U)
           F(j) = Heat_equation( x(j), t, U(j), Ux(j), Uxx(j) ) 
       enddo  
       
       !write(*,*) "maxval F =", maxval(F) 
       
end function 

function Values_at_boundaries(x) result(F) 
        real, intent(in) ::  x(:)
        real :: F( size(x) ) 
            
        integer :: i
        
     Uv( 0, 1)    = x(1)   
     Uv( Nx,  Ns) = x(2)  
     call Derivative( grids(1),   1, Uv(:,1),   Uvx(:,1),   0 ) 
     call Derivative( grids(Ns),  1, Uv(:,Ns),  Uvx(:,Ns), Nx  )  
                              
     F(1) =  BC_heat_equation( x0, t_spatial_discretization, Uv(0,1),   Uvx(0,1) ) 
     F(2) =  BC_heat_equation( xf, t_spatial_discretization, Uv(Nx,Ns), Uvx(Nx,Ns) ) 
              
end function 


function Fluxes_at_interfaces(x) result(F) 
        real, intent(in) ::  x(:)
        real :: F( size(x) ) 
            
        integer :: i
        
     i = i_interface 
     
      Uv( Nx, i)   = x(1)   
      Uv( 0,  i+1) = x(1)  
      call Derivative( grids(i),   1, Uv(:,i),   Uvx(:,i),   Nx ) 
      call Derivative( grids(i+1), 1, Uv(:,i+1), Uvx(:,i+1), 0  ) 
       
      F(1) =   Heat_flux( x0, Uvx(Nx, i) )  & 
             - Heat_flux( xf, Uvx(0, i+1) )   
              
end function 

real function Heat_equation( x, t, u, ux, uxx) 
        real, intent(in) ::  x, t, u, ux, uxx
    
    Heat_equation =   alpha(x) * uxx
                 
end function 


real function BC_heat_equation( x, t, u, ux) 
        real, intent(in) ::  x, t, u, ux
  
   if (x==x0) then 
                   BC_heat_equation = u - 10 * ( 1 + sin( 8* PI * t ) ) 
                   
   else if (x==xf) then 
                   BC_heat_equation = u - 20
   else 
       write(*,*) " Error in BCs " 
       stop 
   end if 
                           
end function 

real function Heat_flux( x, ux) 
        real, intent(in) ::  x, ux
  
   Heat_flux =   -alpha(x) * ux  
              
                   
end function 

real function alpha( x ) 
        real, intent(in) ::  x
        
  
    if (x<0) then 
               alpha =   2
               
    else if (x<0.5) then 
               alpha =   1         
    else 
               alpha  = 0.1 
    end if 
            
                 
end function 

     
subroutine plots
     
     real :: xmax, xmin, ymax, ymin 

     xmax = xf; xmin = x0; 
     ymax = 20. ; ymin = 5. 
     call metafl("xwin")
     call scrmod("reverse")
    
     call disini 
     call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
     call incmrk(1) 
     call marker(21);
     call curve( x, U(Nt-1, :), Ns*(Nx+1) ) 
     call disfin

end subroutine 

end subroutine 


end module 
    
