module MUSE_orbits

       
        use Linear_systems
        use Non_Linear_Systems
        use Cauchy_Problem
        use Temporal_Schemes
        use Temporal_scheme_interface
        use Stability_regions
        use Stability
        use Temporal_error 
        use show_orbits  
        use plots   
        use ODE_Interface
        use Numerical_Recipes
        
        implicit none 

       real, parameter :: PI = 4 * atan(1d0) 
       real :: mu = 0.0121505856
       private 
       public :: Orbits_and_Numerical_Methods, kepler_orbits
       
    contains  
    
 subroutine Orbits_and_Numerical_Methods
 
  integer :: option 
    
  option = 1   
  do while (option>0) 
  write(*,*) " Enter milestone number to execute " 
  write(*,*) " 0. Exit  "
  write(*,*) " 2. Kepler orbit  " 
  write(*,*) " 3. Error in orbits. Convergence  "
  write(*,*) " 4. Regions of absolute stability with examples  "
  write(*,*) " 5. Lagrange points and stability   "
  write(*,*) " 6. Orbits L1, L2 and L3. Arenstorf orbit  "
  read(*,*) option 
 
  select case(option) 
       
     case(0) 
          exit 
 
     case(2) 
        call Kepler_orbits 
        
     case(3) 
        call Error_in_orbits 
        
     case(4) 
        call Regions_of_absolute_stability   
        
     case(5) 
        call Lagrange_points_and_stability
        
      case(6) 
        call Arenstorf_orbit_and_L123   
        
 
     case default 
         write(*,*) " Not implemented"
         
      end select
  end do 
  
 
 end subroutine 

  
 subroutine Error_in_orbits
     
    call Error_Kepler_orbit
    call Convergence_rate_Euler_RK4_Kepler
    call Convergence_rate_GBS_Kepler
  
 end subroutine 
 
 subroutine Regions_of_absolute_stability
 
     integer :: N(3) = [ 1000, 400, 200 ] 
     real :: U0(4) = [1., 0., 0., 0.7]
     real :: U1(2) = [1., 0.]
     
    
    
    call Stability_region_examples
    call Stability_region_examples_high_order
    
    call orbits_and_schemes( 7., Oscillator, Euler, 2, U1, [400, 1500, 2000]  )
     
    call orbits_and_schemes( 1.2, Oscillator, Leap_frog, 2, U1, [150, 200, 300]   )

    call orbits_and_schemes( 1.2, Oscillator, Crank_Nicolson, 2, U1, [150, 200, 300]   )
   
    call orbits_and_schemes( 1.2, Oscillator, Runge_kutta4, 2, U1, [150, 200, 300]  )
    
    
 end subroutine 
 
 subroutine Lagrange_points_and_stability
 
      call CR3BP_Lagrange_points_and_stability
     
 
 end subroutine 
 
 subroutine Arenstorf_orbit_and_L123
 
    call Arenstorf_orbit
    call Shampine_Gordon_orbit
    call Lyapunov_orbits_L1_L2_L3
 
 
 end subroutine 
 
 
 
!****************************************************************
!  Problems: 
!     1) Oscillator 
!     2) Kepler orbit 
!     3) Circular restricted three body problem: CR3BP, CR3BP_2D
!****************************************************************
  
function Oscillator(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
       
     F =  [ U(2), - U(1) ] 

end function  
 
function F_Kepler(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
     
     real :: r(2), drdt(2)
      
     r = U(1:2); drdt = U(3:4);   

     F =  [ drdt, - r / norm2(r)**3 ] 

end function 
   
function F_Kepler_polar(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
     
     real :: r, drdt
      
     r = U(1); drdt = U(2);   

     F =  [ drdt, 1/r**2 * ( 1/r - 1 ) ] 

end function 

function CR3BP(U, t) result(F) 
  real :: U(:), t
  real :: F ( size(U) ) 
    
     real :: x, y, z, vx, vy, vz, dvxdt, dvydt, dvzdt   
     real :: d, r 
   
    
         
     x  = U(1);   y = U(2);   z = U(3);
     vx = U(4);  vy = U(5);  vz = U(6);
     
     d = sqrt( (x+mu)**2 + y**2 + z**2 ) 
     r = sqrt( (x-1+mu)**2 + y**2 + z**2 ) 
     
     dvxdt = x + 2 * vy - (1-mu) * ( x + mu )/d**3 - mu*(x-1+mu)/r**3
     dvydt = y - 2 * vx - (1-mu) * y/d**3 - mu * y/r**3
     dvzdt = - (1-mu)*z/d**3 - mu*z/r**3
     
     F = [ vx, vy, vz, dvxdt, dvydt, dvzdt ] 
     
              
end function

function CR3BP_2D(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
       
     real :: x, y , vx, vy, dxdt, dydt, dvxdt, dvydt  
     real :: D1, D2 
      
     x = U(1); y = U(2); vx = U(3); vy = U(4)     
 
     D1 = sqrt( (x+mu)**2     + y**2 )**3
     D2 = sqrt( (x-(1-mu))**2 + y**2 )**3
     
     dxdt = vx  
     dydt = vy  
     dvxdt = x + 2 * vy - (1-mu)*( x + mu )/D1 - mu*(x-(1-mu))/D2
     dvydt = y - 2 * vx - (1-mu) * y/D1 - mu * y/D2
    
     F = [ dxdt, dydt, dvxdt, dvydt ]
  
end function

 
 
!***************************************** 
!   Kepler_orbit and errors
!*****************************************
 subroutine Kepler_orbits
   integer, parameter :: N = 10000
   real :: Time(0:N), U(0:N,4,4) 
   real :: t0 =0, tf=12*PI 
   integer :: i
   character (len=30) :: legends(4) = ["Euler method", "Inverse Euler method", & 
                                       "Crank_Nicolson method", "Runge_kutta4  " ] 
 
   Time = [ (t0 + (tf -t0 ) * i/real(N), i=0, N ) ]
   do i=1, 4; U(0,:,i) = [1., 0., 0., 1.]; end do
       
   call Cauchy_ProblemS( Time, F_Kepler, U(:,:,1), Euler )
   call Cauchy_ProblemS( Time, F_Kepler, U(:,:,2), Inverse_Euler )
  
   call Cauchy_ProblemS( Time, F_Kepler, U(:,:,3), Crank_Nicolson )
   call Cauchy_ProblemS( Time, F_Kepler, U(:,:,4), Runge_kutta4 )
           
 
   do i=1,4
       call plot(U(:,1,i), U(:,2,i), legends(i) ) 
   end do      
    
  
 
 end subroutine
 
 
subroutine Kepler_orbit_polar
   integer, parameter :: N = 10000
   real :: Time(0:N), U(0:N,2)
   real :: t0 =0, tf=12*PI 
   
   
   integer :: i
 
   Time = [ (t0 + (tf -t0 ) * i/real(N), i=0, N ) ]
   
     U(0,:) = [1., 0.01]
   
     call Cauchy_ProblemS( Time_Domain=Time,                       & 
                           Differential_operator = F_Kepler_polar, & 
                           Solution = U, Scheme = Euler )
     
    write (*, *) 'Polar coordinates: Kepler orbit  '     
    write(*,*) "press enter " ; read(*,*) 
   
    call plot_parametrics(U(:,1), U(:,2:2), ["r"], "r", "v") 
 
end subroutine
  
subroutine Error_Kepler_orbit
        
    real :: t0 = 0, tf =   12*PI           ! Time domain for Kepler problem 
    integer, parameter :: N = 10000         ! Time steps of the first grid 
    integer, parameter :: Nv = 4           ! Number of variables         
    real :: Time (0:N)                     ! Time domain 
    real :: U(0:N, Nv)                     ! Numerical solution 
    real :: Error(0:N, Nv), E(0:N,Nv)      ! Error
    real :: Ue(0:N,Nv)                     ! Exact solution 
    
    integer :: i
       
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    U(0,:) = [1, 0, 0, 1]   !Initial conditions for the Kepler equation 
    Ue(0,:) = U(0,:)
    
    call Cauchy_ProblemS( Time, F_Kepler, Ue )    
    call Error_Cauchy_Problem( Time, F_Kepler, & 
                             !  Euler, 1, U, Error )
                               Inverse_Euler, 1, U, Error )  
                            !  Crank_Nicolson, 2, U, Error ) 
    E = Ue - U 
    call plot_parametrics( Time, Error(:,1:2), ["Error_x", "Error_y"], "$ t $  ", "$ Richardson extrapolation Error $") 
    call plot_parametrics( Time, E(:,1:2), ["x_exact-x", "y_exact-y"], "$ t $  ", "$ Error $")
    call plot_orbits(N, Nv, Ue, U) 
   
end subroutine   
   


!********************************************************************** 
!  Convergence rate of numerical methods with the number of time steps  
!**********************************************************************
subroutine Convergence_rate_Euler_RK4_Kepler
        
 
    real :: t0 = 0, tf = 4*pi           ! Time domain 
    integer, parameter :: N = 500     ! Time steps of the first grid 
    integer, parameter :: Nv = 4      ! Number of variables 
    integer, parameter :: M = 10      ! number of time grids  
    real :: U(0:N,Nv)                 ! Solution
    real :: Time (0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 2      ! number of numerical schemes to be evaluated 
    real :: log_E(M, Np), log_N(M)    ! Error versus  time steps 
                                      ! first index: grid, second index: scheme   
    real :: order
    integer :: i 
    character (len=20) :: names(Np) = [ "Euler", "RK4"] 
   
   
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    U0=[1,0,0,1]             !Initial conditions 
    U(0,:) = U0 
    
    call Cauchy_ProblemS( Time, F_Kepler,  U, Euler ) 
    
    call Temporal_convergence_rate( Time, F_Kepler,  U0,      & 
                                    Euler, order, log_E(:,1), log_N )
    write(*,*) "Order Euler = ", order 
       
    call Temporal_convergence_rate( Time, F_Kepler,  U0,      &
                                    Runge_Kutta4, order, log_E(:,2), log_N )
    write(*,*) "Order Runge_Kutta4 = ", order 
       
    call plot_parametrics( log_N, log_E, names, "$\log N $", "$\log E$ ") 
    
  
end subroutine

!********************************************************************** 
!  Convergence rate of numerical methods with the number of time steps  
!**********************************************************************
subroutine Convergence_rate_GBS_oscillator
        
 
    real :: t0 = 0, tf = 20          ! Time domain 
    integer, parameter :: N = 10      ! Time steps of the first grid 
    integer, parameter :: Nv = 2      ! Number of variables 
    integer, parameter :: M = 10      ! number of time grids  
    real :: U(0:N,Nv)                 ! Solution
    real :: Time (0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 8      ! number of numerical schemes to be evaluated 
    real :: log_E(M, Np), log_N(M)    ! Error versus  time steps 
                                      ! first index: grid, second index: scheme   
   
    real :: order
    integer :: i, j  
    character (len=20) :: names(Np)
   
       
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    U0=[1,0]      
    
    do j=1, Np 
       names(j) = "Levels"
       call set_GBS_levels(j) 
       call Temporal_convergence_rate( Time, oscillator,  U0,      &
                                       GBS_scheme, order, log_E(:,j), log_N )
       write(*,*) "Order GBS = ", order 
    end do     
    call plot_parametrics( log_N, log_E, names, "$\log N $", "$\log E$ ") 
    
  
end subroutine

!********************************************************************** 
!  Convergence rate of numerical methods with the number of time steps  
!**********************************************************************
subroutine Convergence_rate_GBS_Kepler
        
    
    real :: t0 = 0, tf = 4*pi         ! Time domain 
    integer, parameter :: N = 10      ! Time steps of the first grid 
    integer, parameter :: Nv = 4      ! Number of variables 
    integer, parameter :: M = 10      ! number of time grids  
    real :: U(0:N,Nv)                 ! Solution
    real :: Time (0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 8      ! number of numerical schemes to be evaluated 
    real :: log_E(M, Np), log_N(M)    ! Error versus  time steps 
                                      ! first index: grid, second index: scheme   
   
    real :: order
    integer :: i, j  
    character (len=20) :: names(Np)
     
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    U0 = [1,0,0,1] 
    
    do j=1, Np 
       names(j) = "Levels"
       call set_GBS_levels(j) 
       call Temporal_convergence_rate( Time, F_Kepler,  U0,      &
                                       GBS_scheme, order, log_E(:,j), log_N )
       write(*,*) "Order GBS = ", order 
    end do     
    call plot_parametrics( log_N, log_E, names, "$\log N $", "$\log E$ ") 
    
  
end subroutine

subroutine Stability_region_examples
   integer, parameter ::N = 50
   integer :: i, j  
   real :: x(0:N), y(0:N), Region(0:N,0:N)
   real :: x0 = -4d0, xf = 4d0, dx
   real :: y0 = -4d0, yf = 4d0, dy
    
   integer, parameter :: N_levels = 9 
   real :: levels(0:N_levels)  
   
   dx = (xf-x0) / N; dy = (yf-y0) / N
   x = [ ( x0 + dx * i  , i=0,N )]
   y = [ ( y0 + dy * j  , j=0,N )]
    
   levels = [ ( j/real(N_levels)  , j=0, N_levels )]
   
   
   call Absolute_Stability_Region(Euler, x, y, Region) 
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines") 
   
   write (*, *) "Region of absolute stability: Inverse Euler"  
   call Absolute_Stability_Region(Inverse_Euler, x, y, Region) 
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines") 
   
    write (*, *) "Region of absolute stability: RK4"  
   call Absolute_Stability_Region(Runge_kutta4, x, y, Region) 
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines") 
   
   write (*, *) "Region of absolute stability: Crank_Nicolson"  
   call Absolute_Stability_Region(Crank_Nicolson, x, y, Region) 
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines") 
end subroutine 

!******************************************************** 
! Efficiency of high order methods. 
! Number of time steps as a function of error tolerance 
!********************************************************
subroutine Temporal_effort_with_tolerance_GBS_RK
        
  
    real :: t0 = 0, tf = 80*pi               ! Time domain for Van der Pol oscillator
    integer, parameter :: N = 100           ! Time steps of the first grid 
    integer, parameter :: Nv = 4            ! Number of variables
    integer, parameter :: M = 9             ! number of different tolerances           
    real :: Time (0:N)                      ! Time domain 
    real :: U0(Nv)                          ! Initial conditions  
    real :: U(0:N,Nv)                       ! Solution 
    integer, parameter :: Np = 2            ! number of numerical schemes to be evaluated 
    real :: log_mu(M), log_effort(M,Np)     ! log time steps versus  log ( 1/ tolerance)
                                            ! first index: grid, second index: scheme  
    
    
    character (len=20) :: names(Np), family(Np)   
    integer :: i, j
    
      
    family = ["GBS", "eRK"]
    names = ["GBS", "RK87"]
    U0 = [1,0,0,1] 
    U(0,:) = U0 
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    log_mu = [( real(i), i=1, M ) ]
    do j=1, Np 
     call set_solver(family_name = family(j), scheme_name = names(j) )
     call set_tolerance( 1d-5)
     call Cauchy_ProblemS( Time, F_Kepler, U ) 
     call plot_parametrics( U(:,1), U(:,2:2), names(2:2), "$x $", "$y$ ") 
  
     call Temporal_effort_with_tolerance( Time, F_kepler, U0, log_mu, log_effort(:,j)  )
    end do 
    stop 
    
    call plot_parametrics( log_mu, log_effort(:,1:Np), names, "$-\log \epsilon $", "$\log$ Effort ") 
   
end subroutine



subroutine Stability_region_examples_high_order
   integer, parameter ::N = 200
   integer :: i, j  
   real :: x(0:N), y(0:N), Region(0:N,0:N)
   real :: x0 = -8, xf = 8, dx
   real :: y0 = -8, yf = 8, dy
    
   integer, parameter :: N_levels = 9 
   real :: levels(0:N_levels)  
   
   dx = (xf-x0) / N; dy = (yf-y0) / N
   x = [ ( x0 + dx * i  , i=0,N )]
   y = [ ( y0 + dy * j  , j=0,N )]
    
   levels = [ ( j/real(N_levels)  , j=0, N_levels )]
   
   
   call set_GBS_levels(NL=1) 
   call Absolute_Stability_Region(GBS_Scheme, x, y, Region)
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines")
   
   call set_GBS_levels(NL=2) 
   call Absolute_Stability_Region(GBS_Scheme, x, y, Region)
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines")
   
   call set_GBS_levels(NL=4) 
   call Absolute_Stability_Region(GBS_Scheme, x, y, Region)
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines")
   
   call set_GBS_levels(NL=8) 
   call Absolute_Stability_Region(GBS_Scheme, x, y, Region)
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines")
   
   call set_GBS_levels(NL=10) 
   call Absolute_Stability_Region(GBS_Scheme, x, y, Region) 
   call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$",  & 
                     levels = levels, graph_type ="isolines")
   
end subroutine 


subroutine Eigenvalues_Kepler_polar

    integer, parameter :: N = 2 
    real :: t,  A(N, N), U(N)
    integer :: i 
    complex :: lambda(N) 

    U = [ 1, 0] 
    t = 0 

    A =  System_matrix( F_Kepler_polar, U, t ) 
    
    do i=1, N 
        write(*,'(A, 5f8.2)') "A =", A(i,:) 
    end do 
    
    call Eigenvalues_QR(A, lambda) 
    
    do i=1, N 
        write(*,'(A, 2f8.2)') "lambda =", lambda(i) 
    end do 
    
    
end subroutine 


subroutine Orbits_and_schemes(xmax, F, Scheme, Nv, U0, N ) 
    procedure (ODES) :: F
    procedure (Temporal_Scheme), optional :: Scheme
    integer, intent(in) :: Nv 
    real, intent(in) :: xmax, U0(:)
    integer, intent(in) :: N(:) 
    
  
   real, allocatable :: Time(:), U(:,:)
   real :: t0 =0, tf=12*PI 
   integer :: i, j
   integer :: Np 
   real :: xf, yf
   character(len=10) :: colors(3) = [ "red", "green", "blue" ] 
   
  
   Np = size(N) 
   call ini_dislin(xmax = xmax, xmin = -xmax, ymax = xmax, ymin = -xmax ) 
   do i=1, Np 
       
    allocate( Time(0:N(i)), U(0:N(i), Nv) )   
    Time = [ (t0 + (tf -t0 ) * j/real(N(i)), j=0, N(i) ) ]
    U(0,:) = U0(:)
   
    call Cauchy_ProblemS( Time, F, U, Scheme )
    
    call color(colors(i)); 
    call marker(-1); call incmrk(1)
   ! call chncrv('color')
    call curve( U(:,1), U(:,2), N(i)+1) 
    
   
    call marker(21); call incmrk(-1)
    call hsymbl(120)
    xf = U(N(i),1); yf = U(N(i),2) 
    call curve([xf, xf], [yf, yf], 2 )
    
    deallocate( Time, U )
    
   end do    
   call disfin 
   
  
end subroutine
  

!*********************************** 
!  Phase space Kepler 
!***********************************
subroutine Phase_space_Kepler_orbit_polar
   integer, parameter :: N = 1000
   real :: Time(0:N), U(0:N,2)
   real :: t0 =0, tf=12*PI 
   
   integer, parameter :: Np = 10
   real :: r(0:N,0:Np), v(0:N,0:Np)
   real :: v0(0:Np) 
   real :: vi = -0.8, vf = 0.8 
   character(len=10) :: legends(0:Np)
   
   integer :: i
 
   Time = [ (t0 + (tf -t0 ) * i/real(N), i=0, N ) ]
   v0   = [ (vi + (vf -vi ) * i/real(Np), i=0, Np ) ]
   legends= "v0"
   write(*,*) " v0 =", v0 
   
   do i=0, Np 
     U = 0   
     U(0,:) = [1., v0(i)]
   
     call Cauchy_ProblemS( Time_Domain=Time,                       & 
                           Differential_operator = F_Kepler_polar, & 
                           Solution = U, Scheme = Runge_kutta4 )
     
     r(:, i) = U(:,1) 
     v(:, i) = U(:,2) 
   end do   
   
   call plot_parametrics(r(:,0:Np), v(:,0:Np), legends,"r","drdt")
 
end subroutine

subroutine Phase_space_Kepler_orbit
   integer, parameter :: N = 5000
   real :: Time(0:N), U(0:N,4) 
   real :: t0 =0, tf=48*PI 
   integer :: i
   
   integer, parameter :: Np = 10
   real :: x(0:N,0:Np), y(0:N,0:Np)
   real :: dx(0:N,0:Np), dy(0:N,0:Np)
   real :: r(0:N,0:Np), v(0:N,0:Np)
   real :: v0(0:Np) 
   real :: vi = -0.3, vf = 0.3 
   character(len=10) :: legends(0:Np)
   
 
   Time = [ (t0 + (tf -t0 ) * i/real(N), i=0, N ) ]
   v0   = [ (vi + (vf -vi ) * i/real(Np), i=0, Np ) ]
   legends = "v0"
   
   do i=0, Np 
    U(0,:) = [1., 0., v0(i), 1.+v0(i)]
   
    call Cauchy_ProblemS( Time_Domain=Time, Differential_operator = F_Kepler, & 
                    !    Solution = U, Scheme = Euler )
                    !     Solution = U, Scheme = Crank_Nicolson )
                         Solution = U, Scheme = Runge_kutta4 )
    x(:, i)  = U(:,1) 
    y(:, i)  = U(:,2)
    dx(:, i) = U(:,3) 
    dy(:, i) = U(:,4)
    
    r(:,i) = sqrt( x(:,i)**2 + y(:,i)**2 ) 
    v(:,i) = ( x(:,i)*dx(:,i) + y(:,i)*dy(:,i) )/r(:,i)
   end do 
 
    call plot_parametrics(x(:,0:Np), y(:,0:Np), legends,"x","y")
 
    call plot_parametrics(r(:,0:Np), v(:,0:Np), legends,"r","v")
 
end subroutine
 


!*********************************************************************
! Restricted three body problem integrated with different numerical schemes     
!*********************************************************************    
subroutine CR3BP_Lagrange_points_and_stability


    integer, parameter :: N= 20000 ! time steps 
    integer, parameter :: M = 6   ! number of variables 
    real :: U(0:N, M)             ! state vector 
    real :: E(0:N, M)             ! error for different time steps 
    real :: Time(0:N)             ! time domain  
    integer, parameter :: NL = 5  ! number of Lagrange points 
    real :: U0(NL, M)             ! Lagrange points
    real :: eps(M)                ! disturbances around Lagrange points
    real :: A(M,M)                ! Jacobian 
    complex :: lambda(M)          ! eigenvalues 
           
   
   
    real :: t0 = 0
    real :: tf = 4*PI/0.3            ! final integration time 
    integer :: i, j                  ! index 
 !   real :: mu =  0.012150586550569  ! Earth-Moon 
               
   
    Time = [ (t0 + (tf -t0 )*i/N, i=0, N ) ]
    
    U0(1,:) = [ 0.8, 0.6, 0., 0., 0., 0.  ]
    U0(2,:) = [ 0.8, -0.6, 0., 0., 0., 0.  ]
    U0(3,:) = [ -0.1, 0.0, 0., 0., 0., 0.  ]
    U0(4,:) = [ 0.1, 0.0, 0., 0., 0., 0.  ]
    U0(5,:) = [ 1.1, 0.0, 0., 0., 0., 0.  ]
  
    
do i=1, NL   
  ! Lagrange points L1, L2, L3, L4, L5 
    
    call Newton( G, U0(i, 1:3) ) 
    
  ! Jacobian matrix at Lagrange points   
    A =  System_matrix( F = CR3BP, U0 = U0(i,:) , t =0.)
  
  ! stability of the Lagrange points   
    call Eigenvalues_QR( A, lambda )
    A = System_matrix( F = CR3BP, U0 = U0(i,:), t =0.)
        
    write(*,'(A,6f12.4)') " Lagrange point =", U0(i,:) 
    write(*,'(A,6E16.4)') " F in Lagrange point =", CR3BP( U0(i,:), 0.0 )
    
    
    do j=1, M 
         write(*,*) " Eigenvalue =", lambda(j) 
    end do  
    read(*,*)

    
  ! orbit around the Lagrange point  
    call random_number(eps)
    eps = 1d-2 * eps 
    U(0,:) = U0(i,:) + eps 
    
    call Cauchy_ProblemS( Time, CR3BP, U )
   
    call plotm(U, mu, N)  
    call qplot( U(:,1), U(:,2), N+1)
    
  ! Error of the orbit   
    call Error_Cauchy_Problem( Time,  CR3BP, Runge_kutta4, 4, U, E  ) 
    call qplot( Time, E(:,1), N+1)                           

end do     
    
contains  
!-------------------------------------------------
function G(Y) 
  real, intent(in) :: Y(:)
  real :: G(size(Y)) 
    
      real :: X(6), GX(6) 
     X(1:3) = Y
     X(4:6) = 0.
     
     GX = CR3BP(X, 0.)
     G = GX(4:6) 
     
         
end function


end subroutine
  

subroutine Lyapunov_orbits_L1_L2_L3  

    real :: t0 = 0
    real :: tf = 10 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 2000   ! Number of time steps
    integer, parameter :: Nv = 6     ! Number of variables 
    real :: Time(0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 3     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
                                     ! first index: time, second: component, third: scheme 
    real :: Error(0:N, Np)           ! Error associated to different schmes 
    integer :: i, j                  ! indexes 
    character (len=40) :: names(Np)  ! name of parametric curves 
    character (len=40) :: family(Np) ! name of numerical scheme family 
    real :: tolerance = 1d-8         ! error tolerance  
  
     mu = 0.0121505856
     write(*,*) "Lyapunov orbits L1 L2 L3 "   
       
     U(0,:,1) = [0.6089, 0., 0., 0., 0.85246, 0. ]
     U(0,:,2) = [1.3220, 0., 0., 0., -0.60859, 0. ]
     U(0,:,3) = [-1.9118, 0., 0., 0., 1.70633, 0. ]
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
     
     call set_solver( "eRK", "RK87" )
     call set_tolerance(tolerance)
     do j=1, Np 
      call Cauchy_ProblemS( Time, CR3BP, U(:, :, j) )
     end do  
     call plot_parametrics( U(:, 1, :), U(:, 2, :), ["L1", "L2", "L3"], "$x$", "$y$" ) 
     
end subroutine


subroutine Arenstorf_orbit  

    real :: t0 = 0
    real :: tf = 2*17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 3000   ! Number of time steps
    integer, parameter :: Nv = 6     ! Number of variables 
    real :: Time(0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 4     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
                                     ! first index: time, second: component, third: scheme 
    integer :: i, j                  ! indexes 
    character (len=40) :: names(Np)  ! name of parametric curves 
    character (len=40) :: family(Np) ! name of numerical scheme family 
    real :: tolerance = 1d-12         ! error tolerance 
    
    integer, parameter :: M = 9             ! number of different tolerances   
    real :: log_mu(M), log_effort(M,Np)     ! log time steps versus  log ( 1/ tolerance)
                                            ! first index: grid, second index: scheme
    integer :: N0, N_effort 

    log_mu = [( real(i), i=1, M ) ]
    family = ["GBS",  "wGBS", "weRK", "eRK"    ] 
    names = [ "GBS",  "WODEX", "WDOP853", "RK87"  ] 
    
     mu = 0.012277471 ! Earth-Moon 
     U0 = [0.994, 0., 0., 0., -2.0015851063790825, 0. ]
     do j=1, Np 
       U(0,:,j) = U0
     end do   
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
     do j=1, Np 
       call set_solver( family(j), names(j) )
       call set_tolerance(tolerance)
       N0  = get_effort() 
       call Cauchy_ProblemS( Time, CR3BP, U(:, :, j) )
       N_effort = get_effort() - N0 
       write(*,*) " N_effort =", N_effort 
       call Temporal_effort_with_tolerance( Time, CR3BP, U0, log_mu, log_effort(:,j)  )
     end do 
     
     write(*,*) "Arenstorf orbit (GBS, WDOP853, RK87, WODEX)"
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" ) 
     call plot_parametrics( log_mu, log_effort(:,1:Np), names, "Tolerance", "$\log$ Effort ")
    
end subroutine

subroutine Shampine_Gordon_orbit  

    real :: t0 = 0
    real :: tf = 17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 2000   ! Number of time steps
    integer, parameter :: Nv = 6     ! Number of variables 
    real :: Time(0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 3     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
    real :: Error(0:N, Np)           ! Error associated to different schmes 
                                     ! first index: time, second: component, third: scheme 
    integer :: i, j                  ! indexes 
    character (len=40) :: names(Np)  ! name of parametric curves 
    character (len=40) :: family(Np) ! name of numerical scheme family 
    real :: tolerance = 1d-5         ! error tolerance  
    
   
    write(*,*) "Shampine Gordon orbit "     
    read  (*,*)
    
    
     do j=1, Np 
       U(0,:,j) = [1.2, 0., 0., 0., -1.049357509830, 0. ]
     end do   
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ] 
     
     call Cauchy_ProblemS( Time, CR3BP, U(:, :, 1) )
     call plot_parametrics( U(:, 1, 1:1), U(:, 2, 1:1), ["RK4"], "$x$", "$y$" ) 
     
     
     names = ["WDOP853", "DOPRI54", "RK87" ] 
     family = ["weRK", "eRK", "eRK" ]
     do j=1, Np 
       call set_solver( family(j), names(j) )
       call set_tolerance(tolerance)
       call Cauchy_ProblemS( Time, CR3BP, U(:, :, j) )
     end do 
     
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" ) 
 
     
    
    
end subroutine


end module  
