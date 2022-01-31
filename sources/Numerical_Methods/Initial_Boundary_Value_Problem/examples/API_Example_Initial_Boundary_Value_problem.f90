module API_Example_Initial_Boundary_Value_Problem 

use API_Example_Error_IVBP
use Linear_systems
use Initial_Boundary_Value_Problems
!use Finite_differences
use Collocation_methods
use Temporal_Schemes
use Stability_regions
use Temporal_error 
use Stability
use plots 
use Non_Linear_Systems
implicit none 


contains 
 
subroutine IBVP_examples 

              call Heat_equation_1D
              call Advanced_Heat_equation_1D
              call Advection_Diffusion_1D
              call Wave_equation_1D
              
              call Heat_equation_2D
              call Advection_Diffusion_2D
              call Wave_equation_2D
end subroutine 
 

subroutine Heat_equation_1D

       integer, parameter :: Nx = 20, Nt = 1000
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = -1, xf = 1, t0 = 0, tf = 1 
       integer :: i, j, k, q = 6  
       integer, parameter :: Nl = 5 
       character(len=10) :: legends(0:Nl) 
       character(len=100) :: path(2) = [   & 
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Heat1Da", &
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Heat1Db"  ]
   
       
     write (*, '(A50)') 'Time solution of the 1D heat equation'
     
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ] 
     do i=0, Nl; legends(i) = "t="//float_to_str(time(Nt/Nl*i)) ; enddo 
         
     x(0) = x0; x(Nx) = xf
     
!    Heat equation 1D  
     call Grid_Initialization( "nonuniform", "x", x, q )
       
     U(0, :)  =  exp(-25*x**2 )
     call Initial_Boundary_Value_Problem(                              & 
                       Time_Domain = Time, x_nodes = x,                & 
                       Differential_operator =  Heat_equation1D,       & 
                       Boundary_conditions   =  Heat_BC1D,             & 
                       Solution = U ) 
     
     call plot_parametrics(x, transpose(U(0:Nt/2-1:200,:)),           & 
                   legends(0:Nl/2), "$x$", "$u(x,t)$", "(a)", path(1)) 
     call plot_parametrics(x, transpose(U(Nt/2+100:Nt:100,:)),        &
                   legends(Nl/2+1:Nl), "$x$", "$u(x,t)$", "(b)", path(2)) 
           
contains 


real function Heat_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u, ux, uxx
            
            F =   uxx
end function 
!-------------------------------------------------------
real function Heat_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u, ux 
  
        if (x==x0 .or. x==xf) then
                            BC = u 
        else
             write(*,*)  "Error in Heat_BC1D"; stop 
        endif
end function


end subroutine 



subroutine Advanced_Heat_equation_1D

    call Stability_Heat_equation_1D
    call Error_Heat_equation_1D
   ! call 
    

end subroutine 

subroutine Stability_Heat_equation_1D

   integer, parameter :: Nx = 20, Nt = 200
   real ::  x(0:Nx)
   real :: Time(0:Nt), U(0:Nt,0:Nx), dt   
   integer, parameter :: Nl = 5 
   character(len=10) :: legends(0:Nl) 
   
   integer, parameter ::N = 50
   real :: x_R(0:N), x_I(0:N), dx_R = 5./N, dx_I = 8./N , Region(0:N,0:N)
   integer, parameter :: N_levels = 9 
   real :: levels(0:N_levels)  
   real ::  A(0:Nx, 0:Nx)
   complex :: lambda(0:Nx) 
          
       real ::  x0 = -1, xf = 1, t0 = 0, tf = 1 
       integer :: i, j, k, q = 4  
    
      x_R = [ ( -4 + dx_R * i  , i=0,N )]
      x_I = [ ( -4 + dx_I * j  , j=0,N )] 
      levels = [ ( j/real(N_levels)  , j=0, N_levels )]
    
     dt = (tf-t0)/Nt 
     Time = [ (t0 + dt*i, i=0, Nt ) ] 
     x(0) = x0; x(Nx) = xf
     do i=0, Nl; legends(i) = "t="//float_to_str(Time(Nt-(Nl-i)*10)) ; enddo 
     
!    Heat equation 1D  
     call Grid_Initialization( "nonuniform", "x", x, q )
       
     U(0, :)  =  exp(-25*x**2 )
     call Initial_Boundary_Value_Problem(                              & 
                       Time_Domain = Time, x_nodes = x,                & 
                       Differential_operator =  Heat_equation1D,       & 
                       Boundary_conditions   =  Heat_BC1D,             & 
                       Solution = U ) 
     
     call plot_parametrics(x, transpose(U(Nt-Nl*10:Nt:10,:)), legends, "$x$", "$u(x,t)$") 
   
     call Absolute_Stability_Region(Runge_Kutta4, x_R, x_I, Region) 
     A =  System_matrix( F_Cauchy, U(0,:), 0. )
     call Eigenvalues_QR(A, lambda) 
     call plot_stability_region 
  
           
contains 
     
real function Heat_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u, ux, uxx
            
            F =   uxx
end function 

real function Heat_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u, ux 
    
    real :: h = 1 
    
        if (x==x0) then
                            BC = ux + h*( u -sin(8*t) )  
        else if (x==xf) then
                            BC = ux                   
        else
             write(*,*)  "Error in Heat_BC1D"; stop 
        endif
end function

function F_Cauchy( U, t ) result(F) 
   real ::  U(0:), t, F(0:size(U)-1) 
           
  F = Spatial_discretization( Heat_equation1D, Heat_BC1D, x, U, t )
         
end function 

subroutine plot_stability_region

   real :: xmin=-4, xmax=1, ymin=-4, ymax=4 
   real :: xi(0:Nx), yi(0:Nx) 
        
    call scrmod("reverse")
    call page(3000, 3000) 
    call metafl("xwin")
    call disini  
    call color("white")
   
    call graf(xmin, xmax, xmin, (xmax-xmin)/5,  & 
               ymin, ymax, ymin, (ymax-ymin)/5 )  
      do i=0, N_levels 
       call contur( x_R, size(x_R), x_I, size(x_I), Region, levels(i) ) 
      end do 
      
      xi = real( dt * lambda ) 
      yi = imag( dt * lambda ) 
      do i=0, Nx 
       call color("red");  call marker(21);  call incmrk(-1) 
       call curve( [xi(i), xi(i)], [yi(i), yi(i)], 2 )
      end do 
         
    call disfin 

       

end subroutine 

end subroutine 




!********************************************************************************************
subroutine Two_solids_different_conductivity

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

       integer, parameter :: Nx = 20, Nt = 80000, Ns = 3 
       real :: Time(0:Nt)
       real ::  U(0:Nt, Ns*(Nx+1))
       
       character(len=2) :: grids(3) = ["x1", "x2", "x3"]
       real, target ::  x(Ns*(Nx+1)), Ux(Ns*(Nx+1)), Uxx(Ns*(Nx+1))  
       real, pointer:: xv(:, :), Uv(:, :), Uvx(:, :), Uvxx(:, :)
            
       real ::  x0 = -1, xf = 1,  t0 = 0, tf = 1 
       integer :: i, it, j, q 
       real :: PI = 4 * atan(1.), dt, dx 
       real :: xmax, xmin, ymax, ymin 
           
       integer :: i_interface, N_interfaces = Ns-1 
       real :: x1 = 0, x2 = 0.5, x_interfaces(4) 
     
     x_interfaces = [ x0, x1, x2, xf ] 
     xv(0:Nx, 1:Ns) => x     
     Uvx(0:Nx, 1:Ns) => Ux  
     Uvxx(0:Nx, 1:Ns) => Uxx
     dt =  (tf-t0)/Nt 
     Time = [ (t0 + dt*it, it=0, Nt ) ] 
     
     q = Nx 
    
     do i=1, Ns 
         dx =  ( x_interfaces(i+1) - x_interfaces(i) )/Nx 
         xv(:,i) = [ (x_interfaces(i) + dx*j, j=0, Nx ) ] 
         call Grid_Initialization( "nonuniform", grids(i), xv(:,i), q )
     end do  
   
     U(0,:) = 0 
     do it=0, Nt-1 
         
       U(it+1,:) = Eulerf( U(it,:), Time(it), Time(it+1), Space_discretization )        
       
     end do 
     
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
     
  
contains


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
       
end function 

function Values_at_boundaries(x) result(F) 
        real, intent(in) ::  x(:)
        real :: F( size(x) ) 
            
        integer :: i
        
        
     Uv( 0, 1)    = x(1)   
     Uv( Nx,  Ns) = x(2)  
     call Derivative( grids(1),   1, Uv(:,1),   Uvx(:,1),   0 ) 
     call Derivative( grids(Ns),  1, Uv(:,Ns),  Uvx(:,Ns), Nx  )  
     
      
     F(1) =  BC_heat_equation( x0, Time(it), Uv(0,1),   Uvx(0,1) ) 
     F(2) =  BC_heat_equation( xf, Time(it), Uv(Nx,Ns), Uvx(Nx,Ns) ) 
                         
              
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

     
end subroutine 

function Eulerf( U1, t1, t2, F ) result(U2) 
    real, intent(in) ::  U1(:), t1, t2
    real :: U2( size(U1) ) 
    
    interface
      function F(U, t)
        real, target :: U(0:), t 
        real :: F(0:size(U)-1) 
      end function 
    end interface  
    
   U2 = U1 + (t2-t1) * F(U1, t1)  

end function 



 
!*************************************************************************
subroutine Heat_equation_2D

      integer, parameter :: Nx = 20, Ny = 20, Nt = 500
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny) 
       
       real :: x0 = -1, xf = 1, y0 = -1, yf = 1, t0 = 0, tf =  0.5  
       integer :: i, j, k, q = 4  
       real :: levels(10);
       integer, parameter :: Nl = 4 
       character(len=10) :: legends(Nl) = [ "(a)", "(b)", "(c)", "(d)" ]
       character(len=100) :: path(Nl) = [                                   &    
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Heat2Da", &
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Heat2Db", &  
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Heat2Dc", &     
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Heat2Dd"  ]
       
    write (*, *) 'Time solution of the 2D heat equation'
    Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
    x(0) = x0;  x(Nx) = xf; y(0) = y0;  y(Ny) = yf;
    
     call Grid_Initialization( "nonuniform", "x", x, q )
     call Grid_Initialization( "nonuniform", "y", y, q )
   
     U(0, :, :) = 0 
     
!    Heat equation 2D    
     call Initial_Boundary_Value_Problem(                               & 
           Time_Domain = Time, x_nodes = x, y_nodes = y,                &
           Differential_operator =  Heat_equation2D,                    & 
           Boundary_conditions   =  Heat_BC2D,  Solution = U ) 
    
     do i=1, Nl 
      k = Nt / Nl * i    
      call plot_contour(x, y, U(k,:,:), "$x$", "$y$", levels,      &
                        legends(i), path(i), "isolines" ) 
     end do
        
contains


real function Heat_equation2D(x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy) result(F) 
           real,intent(in) :: x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy 
             
             F = Uxx + Uyy 
end function

!-------------------------------------------------------
real function Heat_BC2D( x, y, t, U, Ux, Uy ) result (BC) 
     real, intent(in) :: x, y, t, U, Ux, Uy 

        if (x==x0) then
                                                    BC = U - 1
        else if (x==xf .or. y==y0 .or. y==yf ) then
                                                    BC = U
        else
            write(*,*)  "Error in Heat_BC2D"; stop
        end if
end function

end subroutine


!********************************************************************************************
subroutine Advection_diffusion_1D

       integer, parameter :: Nx = 20, Nt = 1000
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = -1, xf = 1, t0 = 0, tf = 1.5
       integer :: i, k, q = 4 
       integer, parameter :: Nl = 5 
       character(len=10) :: legends(0:Nl) 
       character(len=100) :: path(2) = [   & 
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/AD1Da", &
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/AD1Db"  ]
       
       
     write (*, '(A50)') 'Time solution of the advection diffusion equation'
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ] 
     do i=0, Nl; legends(i) = "t="//float_to_str(time(Nt/Nl*i)) ; enddo 
     x(0) = x0; x(Nx) = xf 
     call Grid_Initialization( "nonuniform", "x", x, q )
  
     U(0, :)  =  exp( -25*x**2 )
     
!    Advection diffusion 1D       
     call Initial_Boundary_Value_Problem(                               & 
                 Time_Domain = Time, x_nodes = x,                       &
                 Differential_operator =  Advection_equation1D,         &
                 Boundary_conditions   =  Advection_BC1D,  Solution = U ) 
     
     call plot_parametrics(x, transpose(U(0:Nt/2-1:200,:)), legends(0:Nl/2), "$x$", "$u(x,t)$", "(a)", path(1)) 
     call plot_parametrics(x, transpose(U(Nt/2+100:Nt:100,:)), legends(Nl/2+1:Nl), "$x$", "$u(x,t)$", "(b)", path(2))
                             
contains 

real function Advection_equation1D( x, t, u, ux, uxx) result(F) 
               real, intent(in) ::  x, t, u, ux, uxx
            
            real :: nu = 0.02
            
            F = - ux + nu * uxx
end function 

!-------------------------------------------------------
real function Advection_BC1D(x, t, u, ux) result(BC) 
         real, intent(in) :: x, t, u, ux 
  
        if (x==x0 .or. x==xf) then
                            BC = u 
        else
             write(*,*)  "Error in Advection_BC1D"; stop
        endif
end function



end subroutine 


!*************************************************************************
subroutine Advection_diffusion_2D

      integer, parameter :: Nx = 20, Ny = 20, Nt = 1000
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny) 
       
       real :: x0 = -1, xf = 1
       real :: y0 = -1, yf = 1
       real :: t0 = 0, tf = 1.5
       integer :: i, j, Order = 8
       real :: levels(10);
       integer, parameter :: Nl = 4 
       character(len=10) :: legends(Nl) = [ "(a)", "(b)", "(c)", "(d)" ]
       character(len=100) :: path(Nl) = [                                   &    
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/AD2Da", &
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/AD2Db", &  
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/AD2Dc", &     
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/AD2Dd"  ]
       
     write (*, *) 'Time solution of the 2D advection diffusion equation'
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x(0)= x0; x(Nx) = xf; y(0)=y0; y(Ny) = yf 
     do i=1, 10; levels(i) = 0.05*i; end do; 
     
     call Grid_Initialization( "nonuniform", "x", x, Order )
     call Grid_Initialization( "nonuniform", "y", y, Order )
  
     U(0, :, :) = Tensor_product( exp(-25*x**2), exp(-25*y**2) ) 
 
!    Advection diffusion 2D      
     call Initial_Boundary_Value_Problem(                             & 
         Time_Domain = Time, x_nodes = x, y_nodes = y,                &
         Differential_operator =  Advection_equation2D,               & 
         Boundary_conditions   =  Advection_BC2D,  Solution = U ) 
     
    
     
     do i=1, 4 
      call plot_contour(x, y, U(300*(i-1),:,:), "$x$", "$y$", levels,      &
                        legends(i), path(i), "isolines" ) 
     end do
     
contains


!----------------------------------------------------
function Advection_equation2D(x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy) result(F) 
           real,intent(in) :: x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy
           real :: F 
        
        real :: nu = 0.02

        F = - Ux + nu * ( Uxx + Uyy )

end function


!-------------------------------------------------------
real function Advection_BC2D( x, y, t, U, Ux, Uy ) result (BC) 
          real, intent(in) :: x, y, t, U, Ux, Uy

        if (x==x0 .or. y==y0 .or. y==yf ) then
                               BC = U  
        elseif (x==xf) then
                       
                               BC =  FREE_BOUNDARY_CONDITION                        
        else
            Write(*,*)  "Error in Advection_BC2D"; stop 
        end if

end function

end subroutine  






!******************************************************************
subroutine Wave_equation_1D

       integer, parameter :: Nx = 40, Nt = 1000, Nv = 2 
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv), x(0:Nx)  
       real ::  x0 = -1, xf = 1, t0 = 0, tf = 4. 
       integer :: i, Order = 4, k;  
       integer, parameter :: Nl = 6 
       character(len=10) :: legends(0:Nl) 
       character(len=100) :: path(2) = [   & 
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Waves1Da", &
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Waves1Db"  ]
       
       
  write (*, '(A50)') 'Time solution of the wave equation'
      
  Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]
  do i=0, Nl; legends(i) = "t="//float_to_str(time(Nt/Nl*i)) ; enddo   
  x(0) = x0; x(Nx) = xf;  
  call Grid_Initialization( "nonuniform", "x", x, Order )
    
  U(0, :,1) = exp( - 15 * x**2 );   U(0, :,2) = 0    
  
!    Wave equation 1D   
     call Initial_Boundary_Value_Problem(                          &  
            Time_Domain = Time, x_nodes = x,                       & 
            Differential_operator =  Wave_equation1D,              &
            Boundary_conditions   =  Wave_BC1D,  Solution = U ) 
  
  
    
     call plot_parametrics(x, transpose(U(0:Nt/2:Nt/Nl,:,1)), legends(0:Nl/2), "$x$", "$u(x,t)$", "(a)", path(1)) 
     call plot_parametrics(x, transpose(U(Nt/2:Nt:Nt/Nl,:,1)), legends(Nl/2:Nl), "$x$", "$u(x,t)$", "(b)", path(2))
    
contains 


function Wave_equation1D( x, t, u, ux, uxx) result(F)
     real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
     real :: F(size(u))
            
        real :: v, vxx, w 
        v = u(1); vxx = uxx(1); 
        w = u(2);  
        
        F  = [w,  vxx]
                 
end function 
!-------------------------------------------------------
function Wave_BC1D(x, t, u, ux) result(BC) 
  real, intent(in) :: x, t, u(:), ux(:) 
  real :: BC( size(u) ) 
         
        if (x==x0 .or. x==xf) then
                                 BC = u
        else
             write(*,*)  "Error in Waves_BC1D"; stop
        endif
end function



end subroutine 
 


!**************************************************************************
subroutine Wave_equation_2D

       integer, parameter :: Nx = 20, Ny = 20, Nt = 300, Nv = 2 
       real ::  x(0:Nx), y(0:Ny)
       real :: Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, Nv)  
       
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1, t0 = 0, tf =  2.
       integer :: i, j, Order = 8 
       
       
       integer, parameter :: Nl = 4, Nlev=30  
       real :: levels(Nlev)
       character(len=10) :: legends(Nl) = [ "(a)", "(b)", "(c)", "(d)" ]
       character(len=100) :: path(Nl) = [                                   &    
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Waves2Da", &
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Waves2Db", &  
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Waves2Dc", &     
           "./doc/chapters/Initial_Boundary_Value_Problem/figures/Waves2Dd"  ]
       
     write (*, *) 'Time solution of the 2D wave equation'
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     x(0) = x0;  x(Nx) = xf; y(0) = y0;  y(Ny) = yf;
     do i=1, Nlev; levels(i) = -1 + 2*(i-1)/real(Nlev-1); end do; 
         
!    Wave equation 2D           
     call Grid_Initialization( "nonuniform", "x", x, Order )
     call Grid_Initialization( "nonuniform", "y", y, Order )
  
     U(0, :, :, 1) = Tensor_product( exp(-10*x**2) , exp(-10*y**2) ) 
     U(0, :, :, 2) = 0
  
     call Initial_Boundary_Value_Problem(                             & 
                       Time_Domain = Time, x_nodes = x, y_nodes = y,  & 
                       Differential_operator = Wave_equation2D,       & 
                       Boundary_conditions = Wave_BC2D, Solution = U )
     
    write(*,*) " CPU_time_BC = ", CPU_time_BC          
    
     do i=1, Nl 
      call plot_contour(x, y, U(Nt/(Nl-1)*(i-1),:,:,1), "$x$", "$y$", levels,      &
                        legends(i), path(i), "isolines" ) 
     end do
     
contains 



function Wave_equation2D( x, y, t, u, ux, uy, uxx, uyy, uxy ) result(L)
      real, intent(in) ::  x,y,t,u(:),ux(:),uy(:),uxx(:),uyy(:),uxy(:)
      real :: L(size(u))
            
            real :: v, vxx, vyy, w 
            
            v = u(1);  vxx = uxx(1); vyy = uyy(1)  
            w = u(2); 
        
            L(1)  = w 
            L(2)  = vxx +vyy
                 
end function 



!-------------------------------------------------------
function Wave_BC2D(x,y, t, u, ux, uy) result(BC) 
  real, intent(in) :: x, y,  t, u(:), ux(:), uy(:) 
  real :: BC( size(u) ) 
  
        real :: v, w 
        v = u(1) 
        w = u(2)
        
        if (x==x0 .or. x==xf .or. y==y0 .or. y==yf ) then
                            BC = [v, w] 
        else
             write(*,*)  "Error in BC2D_waves"; 
             write(*,'(a, 4f7.2)' ) "x0, xf, y0, yf =", x0, xf, y0, yf  
             write(*,'(a, 4f7.2)' ) "x,  y =", x, y  
             stop 
        endif
      
end function

end subroutine 


character(len=5) function float_to_str(x) result(c) 
      real, intent(in) :: x

    write(c,'(f5.1)')  x
   
end function










!*******************************************************************************************************************************************
!  
!*******************************************************************************************************************************************   
subroutine Heat_equation_1D_CN

       integer, parameter :: Nx = 10, Nt = 100
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx)  
       
       real ::  x0 = 0, xf = 1
       real :: t0 = 0, tf = 0.1  
       integer :: i, j, Order = 2 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     
     x(0) = x0; x(Nx) = xf  
     call Grid_Initialization( "uniform", "x", x, Order  )
     
     
     U(0,:) =  0 !  x**2   
       
     call Initial_Boundary_Value_Problem(                               & 
         
                       Time_Domain = Time, x_nodes = x,   & 
                       Differential_operator =  Heat_equation,           & 
                       Boundary_conditions   =  Heat_BC,                 & 
                       Solution = U, Scheme = Euler ) ! Scheme = Runge_Kutta4 ) ! ! ,   ! S ! 
                      
     
     !call scrmod("reverse")
     !call qplot(x, U(0,:),  Nx+1)
     !call qplot(x, U(Nt,:), Nx+1)  
  
   
     
    write(*,*) " x = ", x(5), " U = ", U(Nt, 5) 
    write(*,*) " x = ", x(6), " U = ", U(Nt,6) 
  
     
contains 

!----------------------------------------------------------
real function Heat_equation( x, t, u, ux, uxx) 
        real, intent(in) ::  x, t, u, ux, uxx
            
            Heat_equation =   uxx
                 
end function 
!-------------------------------------------------------
real function Heat_BC(x, t, u, ux)
  real, intent(in) :: x, t, u, ux 
  
        if (x==x0) then
                          Heat_BC = u - 1 
                      !       Heat_BC = ux

        else if (x==xf) then
                           Heat_BC = u
                       !     Heat_BC = ux

        else
             write(*,*)  "Error in Heat_BC"
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif

        
end function


end subroutine

subroutine Heat_equation_1D_system 

       integer, parameter :: Nx = 10, Nt = 100, Nv = 1
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv)  
       
       real ::  x0 = 0, xf = 1
       real :: t0 = 0, tf = 0.1  
       integer :: i, j, Order = 2 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ]  
     
     x(0) = x0; x(Nx) = xf  
     call Grid_Initialization( "uniform", "x", x, Order  )
     
     
     U(0,:, 1) =  0 !  x**2   
       
     call Initial_Boundary_Value_Problem(                               & 
         
                       Time_Domain = Time, x_nodes = x,   & 
                       Differential_operator =  Heat_equation,           & 
                       Boundary_conditions   =  Heat_BC,                 & 
                       Solution = U, Scheme = Euler ) ! Scheme = Runge_Kutta4 ) ! ! ,   ! S ! 
                      
     
     !call scrmod("reverse")
     !call qplot(x, U(0,:),  Nx+1)
     !call qplot(x, U(Nt,:), Nx+1)  
  
   
     
    write(*,*) " x = ", x(5), " U = ", U(Nt, 5, 1) 
    write(*,*) " x = ", x(6), " U = ", U(Nt,6, 1) 
    
     
contains 

!----------------------------------------------------------
function Heat_equation( x, t, u, ux, uxx) result(F) 
        real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
        real :: F(size(u)) 
            
            F(1) =   uxx(1) 
            
                 
end function 
!-------------------------------------------------------
function Heat_BC(x, t, u, ux)
  real, intent(in) :: x, t, u(:), ux(:) 
  real :: Heat_BC(size(u)) 
  
        if (x==x0) then
                          Heat_BC = u - 1 
                      !       Heat_BC = ux

        else if (x==xf) then
                           Heat_BC = u
                       !     Heat_BC = ux

        else
             write(*,*)  "Error in Heat_BC"
             write(*,*) " x0 =", x0, " xf =", xf 
             write(*,*) " x = ", x 
             stop 
        endif

end function



end subroutine 

end module 
