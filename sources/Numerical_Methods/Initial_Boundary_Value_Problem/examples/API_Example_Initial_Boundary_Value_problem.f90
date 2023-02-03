module API_Example_Initial_Boundary_Value_Problem 

use Linear_systems
use Initial_Boundary_Value_Problems
use Collocation_methods
use Temporal_Schemes
use Stability_regions
use Temporal_error 
use Stability
use plots 
use Non_Linear_Systems
use Utilities
use Numerical_Recipes
implicit none 


contains 
 
subroutine IBVP_examples 

      call Heat_equation_1D
      call Stability_and_Error_Heat1D
      call Stability_Heat_equation_1D
      
      call Heat_equation_1D_implicit
      call Advection_Diffusion_1D
      call Wave_equation_1D
     
      call Heat_equation_2D
      call Stability_and_Error_Heat2D
      call Advection_Diffusion_2D
      call Wave_equation_2D
end subroutine 
 

subroutine Heat_equation_1D

       integer, parameter :: Nx = 20, Nt = 1000, Nv = 1
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv)  
       
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
       
     U(0, :, 1)  =  exp(-25*x**2 )
     call Initial_Boundary_Value_Problem(                              & 
                       Time_Domain = Time, x_nodes = x,                & 
                       Differential_operator =  Heat_equation1D,       & 
                       Boundary_conditions   =  Heat_BC1D,             & 
                       Solution = U ) 
     
     call plot_parametrics(x, transpose(U(0:Nt/2-1:200,:,1)),           & 
                   legends(0:Nl/2), "$x$", "$u(x,t)$", "(a)", path(1)) 
     call plot_parametrics(x, transpose(U(Nt/2+100:Nt:100,:,1)),        &
                   legends(Nl/2+1:Nl), "$x$", "$u(x,t)$", "(b)", path(2)) 
contains           

function Heat_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
          real :: F( size(u) ) 
            
            F(1) =   uxx(1)
end function 

function Heat_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u(:), ux(:) 
    real :: BC( size(u) ) 
  
        if (x==x0 .or. x==xf) then
                            BC(1) = u(1)  
        else
             write(*,*)  "Error in Heat_BC1D"; stop 
        endif
end function

end subroutine 

 
subroutine Stability_and_Error_Heat1D

       integer, parameter :: Nx = 20,  Nv = 1
       real ::  x(0:Nx), U(0:Nx, Nv),  F(0:Nx, Nv), t   
       real ::  A( (Nx+1)*Nv, (Nx+1)*Nv), R(0:Nx, Nv)
       complex ::  lambda( (Nx+1)*Nv ) 
       real ::  x0 = -1, xf = 1
       integer :: i, q, qmax 
 
 !   Polynomial Order=2,4,6, ... Nx   
     qmax = 2 ! Nx   
     do q=2, qmax, 2   
         
 !    Spatial discretization         
      x(0) = x0; x(Nx) = xf; t = 0 
      call Grid_Initialization( "nonuniform", "x", x, q )
      U = Test_U(Nv, x) 
      F = Spatial_discretization( Heat_equation1D, Heat_BC1D, x, U, t ) 
      call plot(x, U(:,1), "Test condition Heat1D U = (1-x2) exp(-5 x2) ") 
      call plot(x, F(:,1), "Spatial discretization Heat1D F = Uxx") 
      
 !    Spatial Truncation Error 
      R = Spatial_Truncation_Error( Nv, Heat_equation1D, Heat_BC1D, x, q, Test_U  )
      call plot(x(1:Nx-1), R(1:Nx-1,1), "Estimated spatial truncation error R")
      
      R = Test_Uxx(Nv, x) - F
      call plot(x(1:Nx-1), R(1:Nx-1,1), "Exact spatial truncation error R")
       
 !    Matrix of the linear system      
      write(*,*) " Linear operator (system matrix) "  
      A = Linear_operator( Nv, x, q, Heat_equation1D, Heat_BC1D )  
      do i=1, Nx+1
         write(*,'(100f8.2)') A(i,:) 
      end do 
      write(*,*) 
      
      call Eigenvalues_QR(A, lambda) 
      write(*,*) " Eigenvalues of the Linear operator (system matrix) "  
      do i=1, (Nx+1)*Nv
          write(*,'(i3,2f8.2)')  i, lambda(i) 
       end do 
       write(*,*)  
       
     end do 
     
contains     

function Test_U( Nv, x) result(U)
      integer, intent(in) :: Nv 
      real, intent(in) ::  x(0:) 
      real :: U( 0:size(x)-1, Nv ) 
            
   U(:,1)  =  (1-x**2) * exp(-5*x**2) 
          
end function 

function Test_Uxx( Nv, x) result(Uxx)
      integer, intent(in) :: Nv 
      real, intent(in) ::  x(0:) 
      real :: Uxx( 0:size(x)-1, Nv ) 
            
 ! U(:,1)   =  (1-x**2) * exp(-5*x**2)
 ! Ux(:,1)  =  -2* x * exp(-5*x**2) + (1-x**2) *( -10*x * exp(-5*x**2) ) 
 !          =  (-12*x + 10* x**3 ) * exp(-5 *x**2) 
      
   Uxx(:,1) = (-12 + 30*x**2) * exp(-5 * x**2) +  (-12*x + 10* x**3 ) * (-10*x) * exp( -5 *x**2 ) 
          
end function 


function Heat_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
          real :: F( size(u) ) 
            
            F(1) =   uxx(1)
end function 

function Heat_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u(:), ux(:) 
    real :: BC( size(u) ) 
  
        if (x==x0 .or. x==xf) then
                            BC(1) =  u(1)  
        else
             write(*,*)  "Error in Heat_BC1D"; stop 
        endif
end function

end subroutine     



subroutine Heat_equation_1D_implicit

       integer, parameter :: Nx = 20, Nt = 20, Nv = 1
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv)  
       
       real ::  x0 = -1, xf = 1, t0 = 0, tf = 1 
       integer :: i, j, k, q = 2  
   
      write (*, '(A50)') 'Implicit 1D Heat equation: Inverse Euler + FD2 ' 
      write(*,*) "press enter"
      read(*,*) 
      
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ] 
     
   
!    Heat equation 1D    
     x(0) = x0; x(Nx) = xf
     call Grid_Initialization( "nonuniform", "x", x, q )
       
     U(0,:,1)  =  exp(-25*x**2 )
     call Initial_Boundary_Value_Problem( Time, x, Heat_equation1D, Heat_BC1D, U, Inverse_Euler ) 
     
     call plot(x, U(0,:,1),"Inverse Euler integration, Initial condition") 
     call plot(x, U(Nt,:,1),"Inverse Euler Solution at t=1 (Nx=20, Nt=20)") 
     
     call Initial_Boundary_Value_Problem( Time, x, Heat_equation1D, Heat_BC1D, U, CranK_Nicolson ) 
     call plot(x, U(Nt,:,1),"Crank-Nicolson Solution at t=1 (Nx=20, Nt=20)") 
     
contains           

function Heat_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
          real :: F( size(u) ) 
            
            F(1) =   uxx(1)
end function 

function Heat_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u(:), ux(:) 
    real :: BC( size(u) ) 
  
        if (x==x0 .or. x==xf) then
                            BC(1) = u(1)  
        else
             write(*,*)  "Error in Heat_BC1D"; stop 
        endif
end function

end subroutine 


!******************************************************************
!* Eigenvalues of the heat equation with FD and Chebyshev 
!******************************************************************
subroutine Stability_Heat_equation_1D

   integer, parameter :: Nx = 20
   real ::   x(0:Nx), A(0:Nx, 0:Nx), W(0:Nx, 0:Nx) 
   complex :: lambda(Nx, 0:Nx) ! order=2,4,6,.., index eigenvalue=1,2,3....Nx-1
   real :: rlambda(Nx, 0:Nx) 
   
   integer, parameter ::N = 50, N_levels = 9 
   real :: x_R(0:N), x_I(0:N), dx_R = 5./N, dx_I = 8./N , Region(0:N,0:N), levels(0:N_levels) 
          
   real ::  x0 = 0, xf = 1   
   integer :: i, j, k, q  
   
   x(0) = x0; x(Nx) = xf
  
    
  ! Eigenvalues of discrete operator with different orders 
    do q=2,  Nx
       call Grid_Initialization( "nonuniform", "x", x, q )
       A =  Linear_operator( 1, x, q, Heat_equation1D, Heat_BC1D )
                
      call Eigenvalues_QR( A, lambda(q,:)) 
      rlambda(q,:) = -Real( lambda(q,:) ) 
      call sortr1( rlambda(q,:), Nx+1, "A") 
      !write(*,*) lambda(q,:) 
      !read(*,*) 
      
    end do 
    
    x_R = [ ( -4 + dx_R * i  , i=0,N )]
    x_I = [ ( -4 + dx_I * j  , j=0,N )]    
    call Absolute_Stability_Region(Runge_Kutta4, x_R, x_I, Region) 
    
    call plot_eigenvalues
    call plot_stability_region 
  
           
contains 

function Heat_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
          real :: F( size(u) ) 
            
            F(1) =   uxx(1)
end function 

function Heat_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u(:), ux(:) 
    real :: BC( size(u) ) 
    
    real :: h = 1 
    
        if (x==x0) then
                            BC(1) = u(1) 
        else if (x==xf) then
                            BC(1) = u(1)  
        else
             write(*,*)  "Error in Heat_BC1D"; stop 
        endif
end function

subroutine plot_eigenvalues

   real :: xmin, xmax, ymin, ymax
   real :: xi(0:Nx), yi(0:Nx) 
   character(len=10) :: col(15) = ["red","green", "cyan", "orange", "white", "blue", "blue", "cyan", "orange", "white", "green", "blue", "cyan", "orange", "white" ]
   real ::  lambda_exact(0:Nx),  lambda_FD2(0:Nx)  
   
   real :: dx, PI = 4*atan(1.)
   
    dx = (xf-x0)/Nx 
   
    xmin = 0; xmax = Nx; ymin = 0; ymax = 5000;  
    call dislin_ini(xmax, xmin, ymax, ymin) 
   
    do k=0, Nx
      lambda_exact(k)  = (k*PI)**2 
      lambda_FD2(k)  = 4 /dx**2 * ( sin(k*PI*dx/2) )**2 
    end do 
       
    xi = [ (k, k=0,Nx) ]
    do j= 2,  6, 2    
       yi = rlambda(j,:)
       call color(col(j)); call incmrk(1);   call marker(21);  
       call curve( xi, yi, Nx+1 )
    end do 
    
    j=1; yi = rlambda(Nx,:) 
    call color(col(j)); call incmrk(1);   call marker(21);  
    call curve( xi, yi, Nx+1 )
         
    j= 5; yi = lambda_exact
    call color(col(j)); call incmrk(1);   call marker(21);  
    call curve( xi, yi, Nx+1 )
    
    call plot_title( ["Eigenvalues of d2u/dx2 = lambda u ", &
                      " u(0)=0, u(1)=0", " N = 20 " ] ) 
    call plot_legends( ["q=2","q=4", "q=6","q=N(Chebyshev)","Exact eigenvalues" ] )
    
    call disfin 
 

end subroutine 

subroutine plot_stability_region

   real :: xmin=-4, xmax=1, ymin=-4, ymax=4, dt = 0.001 
   real :: xi(0:Nx), yi(0:Nx) 
       
   
    levels = [ ( j/real(N_levels)  , j=0, N_levels )] 
    call dislin_ini(xmax, xmin, ymax, ymin) 
    
    do i=0, N_levels 
       call contur( x_R, size(x_R), x_I, size(x_I), Region, levels(i) ) 
    end do 
    
    xi = -dt * rlambda(2,:)    
    yi = 0 
    call color("red");  call marker(21);  call incmrk(-1) 
    call curve( xi, yi, Nx+1 )
      
    call plot_title( [" Eigenvalues (q=2) multiplied by the time step (dt=0.001) ", &
                      " inside the stability region of a fourth order Runge-Kutta" ] )
    call disfin 
       

end subroutine 
end subroutine 






subroutine Heat_equation_2D

      integer, parameter :: Nx = 20, Ny = 20,  Nv = 1, Nt = 500
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, 1:Nv) 
       
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
   
     U(0, :, :, 1) = 0 
     
!    Heat equation 2D    
     call Initial_Boundary_Value_Problem(                               & 
           Time_Domain = Time, x_nodes = x, y_nodes = y,                &
           Differential_operator =  Heat_equation2D,                    & 
           Boundary_conditions   =  Heat_BC2D,  Solution = U ) 
    
     do i=1, Nl 
      k = Nt / Nl * i    
      call plot_contour(x, y, U(k,:,:,1), "$x$", "$y$", levels,      &
                        legends(i), path(i), "isolines" ) 
     end do
        
contains


function Heat_equation2D(x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy) result(F) 
   real,intent(in) :: x, y, t, U(:), Ux(:), Uy(:), Uxx(:), Uyy(:), Uxy(:) 
   real F( size(U) ) 
                     F(1) = Uxx(1) + Uyy(1) 
end function


function Heat_BC2D( x, y, t, U, Ux, Uy ) result (BC) 
     real, intent(in) :: x, y, t, U(:), Ux(:), Uy(:) 
     real :: BC(size(U))
        if      (x==x0) then; BC = U - 1              
        else if (x==xf) then; BC = U              
        else if (y==y0) then; BC = Uy - 5 * U    
        else if (y==yf) then; BC = Uy + 5 * U                                              
        else
            write(*,*)  "Error in Heat_BC2D"; stop
        end if
end function

end subroutine
                                    

subroutine Stability_and_Error_Heat2D


       integer, parameter :: Nx = 10, Ny =10,  Nv = 1
       real ::  x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv),  F(0:Nx, 0:Ny, Nv), R(0:Nx, 0:Ny, Nv), t   
       real ::  A( (Nx+1)*(Ny+1)*Nv, (Nx+1)*(Ny+1)*Nv)
       real, allocatable :: v(:)
       complex ::  lambda( (Nx+1)*(Ny+1)*Nv ) 
       real ::  x0 = -1, xf = 1, y0 = -1, yf = 1
       integer :: i, q, qmax       
     
!    Heat equation 2D  
     qmax = 2   
     do q=2, qmax, 2   
         
 !    Spatial discretization 
      x(0) = x0; x(Nx) = xf; y(0) = y0; y(Ny) = yf; t = 0 
      call Grid_Initialization( "nonuniform", "x", x, q )
      call Grid_Initialization( "nonuniform", "y", y, q )
      
      U = Test_U(Nv, x, y) 
      call plot_contour(x, y, U(:,:,1), "x", "y", legend="Test condition Heat2D U") 
      
      F = Spatial_discretization( Heat_equation2D, Heat_BC2D, x, y, U, t ) 
      call plot_contour(x, y, F(:,:,1), "x", "y", legend="Spatial discretization F = Uxx + Uyy" ) 
       
 !    Spatial Truncation Error
      R = Spatial_Truncation_Error( Nv, Heat_equation2D, Heat_BC2D, x, y, q, Test_U  )
      call plot_contour(x, y, R(:,:,1), "x", "y", legend="Estimated spatial truncation error" ) 
    
       
!     Matrix of the linear system        
      A = Linear_operator( Nv, x, y, q, Heat_equation2D, Heat_BC2D )  
      write(*,*) " Linear operator (system matrix) "  
      do i=1, (Nx+1)*(Ny+1)*Nv
         v = A(i,:) 
         write(*,'(A3,i3, 200f8.1)') "i=", i, pack(v, v/=0)
      end do 
      write(*,*) 
      
      call Eigenvalues_QR(A, lambda) 
      write(*,*) " Eigenvalues of the Linear operator (system matrix) "  
      do i=1, (Nx+1)*(Ny+1)*Nv
         write(*,'(i3,2f10.1)')  i, lambda(i) 
      end do 
      write(*,*)
      
     end do 
     
contains         

function Heat_equation2D(x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy) result(F) 
   real,intent(in) :: x, y, t, U(:), Ux(:), Uy(:), Uxx(:), Uyy(:), Uxy(:) 
   real F( size(U) ) 
                     F(1) = Uxx(1) + Uyy(1) 
end function


function Heat_BC2D( x, y, t, U, Ux, Uy ) result (BC) 
     real, intent(in) :: x, y, t, U(:), Ux(:), Uy(:) 
     real :: BC(size(U))
        if (x==x0) then
                                                    BC(1) = U(1) 
        else if (x==xf .or. y==y0 .or. y==yf ) then
                                                    BC(1) = U(1)
        else
            write(*,*)  "Error in Heat_BC2D"; stop
        end if
end function

function Test_U( Nv, x, y) result(U)
      integer, intent(in) :: Nv 
      real, intent(in) ::  x(0:), y(0:) 
      real :: U( 0:size(x)-1, 0:size(y)-1, Nv ) 
            
   U(:,:, 1)  =  Tensor_product( (1-x**2)*exp(-5*x**2), (1-y**2)*exp(-5*y**2) )
          
end function 


end subroutine   







subroutine Advection_diffusion_1D

       integer, parameter :: Nx = 20, Nt = 1000, Nv = 1 
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv)  
       
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
  
     U(0, :, 1)  =  exp( -25*x**2 )
     
!    Advection diffusion 1D       
     call Initial_Boundary_Value_Problem(                               & 
                 Time_Domain = Time, x_nodes = x,                       &
                 Differential_operator =  Advection_equation1D,         &
                 Boundary_conditions   =  Advection_BC1D,  Solution = U ) 
     
     call plot_parametrics(x, transpose(U(0:Nt/2-1:200,:,1)), legends(0:Nl/2), "$x$", "$u(x,t)$", "(a)", path(1)) 
     call plot_parametrics(x, transpose(U(Nt/2+100:Nt:100,:,1)), legends(Nl/2+1:Nl), "$x$", "$u(x,t)$", "(b)", path(2))
                             
contains 

function Advection_equation1D( x, t, u, ux, uxx) result(F) 
               real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
               real :: F( size(u) ) 
            
            real :: nu = 0.02
            
            F(1) = - ux(1) + nu * uxx(1)
end function 

!-------------------------------------------------------
function Advection_BC1D(x, t, u, ux) result(BC) 
         real, intent(in) :: x, t, u(:), ux(:)
         real :: BC ( size(u) ) 
  
        if (x==x0 .or. x==xf) then
                            BC(1) = u(1) 
        else
             write(*,*)  "Error in Advection_BC1D"; stop
        endif
end function

end subroutine 






subroutine Advection_diffusion_2D

      integer, parameter :: Nx = 20, Ny = 20, Nv =1, Nt = 1000
      real :: x(0:Nx), y(0:Ny), Time(0:Nt), U(0:Nt, 0:Nx, 0:Ny, 1:Nv) 
       
       real :: x0 = -1, xf = 1, y0 = -1, yf = 1
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
  
     U(0, :, :, 1) = Tensor_product( exp(-25*x**2), exp(-25*y**2) ) 
 
!    Advection diffusion 2D      
     call Initial_Boundary_Value_Problem(                             & 
         Time_Domain = Time, x_nodes = x, y_nodes = y,                &
         Differential_operator =  Advection_equation2D,               & 
         Boundary_conditions   =  Advection_BC2D,  Solution = U ) 
     
     do i=1, Nl 
      j = Nt*(i-1)/(Nl-1)   
      call plot_contour(x, y, U(j,:,:,1), "$x$", "$y$", levels,  &
                        legends(i), path(i), "isolines" ) 
     end do
     
contains


function Advection_equation2D(x, y, t, U, Ux, Uy, Uxx, Uyy, Uxy) result(F) 
    real,intent(in) :: x, y, t, U(:), Ux(:), Uy(:), Uxx(:), Uyy(:), Uxy(:)
    real :: F( size(U) )
        
        real :: nu = 0.02

        F(1) = - Ux(1) + nu * ( Uxx(1) + Uyy(1) )
end function


function Advection_BC2D( x, y, t, U, Ux, Uy ) result (BC) 
          real, intent(in) :: x, y, t, U(:), Ux(:), Uy(:)
          real :: BC( size(U) ) 

        if (x==x0 .or. y==y0 .or. y==yf ) then
                               BC(1) = U(1)  
        elseif (x==xf) then
                       
                               BC(1) =  FREE_BOUNDARY_CONDITION                        
        else
            Write(*,*)  "Error in Advection_BC2D"; stop 
        end if

end function

end subroutine 







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

end module 

