module API_Example_Cauchy_Problem

    use Cauchy_Problem
    use Temporal_Schemes
    use Stability_regions
    use Temporal_error 
    use Stability
    use plots 
    use MUSE_orbits
      
    implicit none
    
    
    
    contains  
  
!*********************************** 
!  Cauchy problem solutions and tools
!***********************************
subroutine Cauchy_problem_examples

   
   call First_order_ODE
   call Linear_Spring
   call Lorenz_Attractor
   call Stability_regions_RK2_RK4
   call Error_solution
   call Convergence_rate_RK2_RK4
   call Variable_step_with_Predictor_Corrector
   
   
end subroutine   

!*********************************** 
!  Advanced numerical methods 
!***********************************
subroutine Advanced_Cauchy_problem_examples
  integer :: option
  
  option = 1
  do while (option>0) 
     write(*,*) "Advanced methods, select an option: " 
     write(*,*) " 0. Exit/quit  "
     write(*,*) " 1. Van del Pol system " 
     write(*,*) " 2. Henon Heiles system "
     write(*,*) " 3. Variable time step versus constant time step "
     write(*,*) " 4. Convergence rate of Runge Kutta wrappers  "
     write(*,*) " 5. Arenstorf orbit (Embedded Runke-Kutta) "
     write(*,*) " 6. Arenstorf orbit (GBS methods, Wrapper ODEX)"
     write(*,*) " 7. Arenstorf orbit (ABM methods, Wrapper ODE113)"
     write(*,*) " 8. Computational effort Runge-Kutta methods" 
     write(*,*) " 9. Orbits and numerical methods (Master MUSE)  "
     read(*,*) option 
     select case(option)
      case(1)  
             call Van_der_Pol_oscillator
      case(2)      
             call Henon_Heiles_system 
      case(3)         
             call Variable_step_simulation
      case(4)        
             call Convergence_rate_Runge_Kutta_wrappers
      case(5) 
             call Runge_Kutta_wrappers_versus_original_codes 
      case(6) 
             call GBS_and_wrapper_ODEX
             call Arenstorf_with_GBS
      case(7) 
             call ABM_and_wrapper_ODE113
      case(8) 
             call Temporal_effort_with_tolerance_eRK
      case(9) 
             call Orbits_and_Numerical_Methods         
      case default
     end select 
  end do
end subroutine 

!*********************************** 
!  First order equation 
!***********************************
subroutine First_order_ODE
        
    real :: t0 = 0, tf = 4
    integer :: i 
    integer, parameter :: N = 1000  !Time steps
    real :: Time(0:N), U(0:N,1)

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
    U(0,1) =  1
   
    call Cauchy_ProblemS( Time_Domain = Time,                      & 
                          Differential_operator = F1, Solution = U ) 
    
    write(*,*) "Solution of du/dt - 2 u"     
    call plot_parametrics(Time, U, ["Solution of du/dt-2u"], "time", "u") 
  
contains

function F1( U, t ) result(F) 
    real :: U(:), t 
    real :: F(size(U)) 
    
    F(1) = -2*U(1)
  
end function 


end subroutine


!*********************************** 
!  Second order equation  
!***********************************
subroutine Linear_Spring 
    integer :: i 
    integer, parameter :: N = 100   !Time steps
    real :: t0 = 0, tf = 4, Time(0:N), U(0:N, 2)

    Time = [ (t0 + (tf -t0 )*i/real(N), i=0, N ) ]
    U(0,:) = [ 5, 0] 
    call Cauchy_ProblemS( Time_Domain = Time ,                     & 
                          Differential_operator = F_spring,        & 
                          Solution = U, Scheme = Crank_Nicolson    )
   
    write (*, *) 'Solution of the Cauchy problem:  ' 
    write (*, *) ' d2u/dt2 = -3 t u, u(0) = 5, du(0)/dt = 0' 
    call plot_parametrics(Time, U, ["Sd2u/dt2 = -3 t  u"], "time", "u") 
  
contains

function F_spring( U, t )  result(F) 

    real :: U(:), t 
    real :: F(size(U)) 
            
    real, parameter :: a = 3.0 
    
    F(1) = U(2)
    F(2) = -a * t * U(1)
     
end function

end subroutine


!*********************************** 
!  Chaos behaviour 
!***********************************
subroutine Lorenz_Attractor
   integer, parameter :: N = 10000
   real :: Time(0:N), U(0:N,3) 
   real :: a=10., b=28., c=2.6666666666
   real :: t0 =0, tf=25 
   integer :: i
 
   Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
   
   U(0,:) = [12, 15, 30]
   
    call Cauchy_ProblemS( Time_Domain=Time, Differential_operator=F_L, & 
                          Solution = U, Scheme = Runge_Kutta4 )
    
    
    write (*, *) 'Solution of the Lorenz attractor  '   
    call plot_parametrics(U(:,1),U(:,2:2), ["Lorenz attractor"],"x","y")
       
contains

 function F_L(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
     
     real :: x, y , z
      
     x = U(1); y = U(2); z = U(3)   

     F(1) =   a * ( y - x ) 
     F(2) =   x * ( b - z ) - y 
     F(3) =   x * y - c * z   

 end function
 
end subroutine



!*********************************** 
!  Regions of absolute stability  
!***********************************
subroutine Stability_regions_RK2_RK4
   integer, parameter ::N = 50
   integer :: i, j  
   real :: x(0:N), y(0:N), Region(0:N,0:N)
   real :: x0 = -4d0, xf = 1d0, dx
   real :: y0 = -4d0, yf = 4d0, dy
   character(len=100) :: path(2) =  & 
                   [ "./doc/chapters/Cauchy_problem/figures/RK2a", & 
                     "./doc/chapters/Cauchy_problem/figures/RK4b"  ] 
   character(len=5) :: legends(2) = [ "(a)", "(b)" ] 
   
   integer, parameter :: N_levels = 9 
   real :: levels(0:N_levels)  
   
   dx = (xf-x0) / N; dy = (yf-y0) / N
   x = [ ( x0 + dx * i  , i=0,N )]
   y = [ ( y0 + dy * j  , j=0,N )]
   write (*, *) "Region of absolute stability: Runge kutta 2, 4"  
   
   levels = [ ( j/real(N_levels)  , j=0, N_levels )]
   do j=1, 2  
     if (j==1)      then 
         call Absolute_Stability_Region(Runge_Kutta2, x, y, Region) 
     else if (j==2) then 
         call Absolute_Stability_Region(Runge_Kutta4, x, y, Region) 
     end if 
     
     call plot_contour(x, y, Region, "$\Re(z)$","$\Im(z)$", levels, & 
                       legends(j), path(j), "isolines")  
   end do 
   

end subroutine 



!*********************************************************** 
!  Error estimation by means of Richardson extrapolation 
!***********************************************************
subroutine Error_solution
        
    real :: t0 = 0, tf = 30       ! Time domain for Van der Pol oscillator
    integer, parameter :: N = 800 ! Time steps of the first grid 
    integer, parameter :: Nv = 2  ! Number of variables van de pol
    integer, parameter :: M = 10  ! number of time grids           
    real :: Time (0:N)            ! Time domain 
    real :: Solution(0:N, Nv)     ! Solution 
    real :: Error(0:N, Nv)        ! Error
    integer, parameter :: Np = 2  ! number of graphs: solution and error 
    real :: y(0:N, Np) 
    
    character (len=20) :: names(Np) = [ "Solution", "Error"] 
    character (len=200) :: path(Np) = [ "./doc/chapters/Cauchy_problem/figures/VanderPol_errora", &
                                        "./doc/chapters/Cauchy_problem/figures/VanderPol_errorb"  ] 
    integer :: i
    
    write(*,*) "Error determination by Richardson extrapolation"   
    write(*,*) "Temporal scheme : Rk2 "  
    
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    Solution(0,:) = [3, 4]   !Initial conditions for VanDerpol equation 
    
    call Error_Cauchy_Problem( Time, VanDerPol_equation,        &
                               Runge_Kutta2, 2, Solution, Error ) 
                               
    y(:,1) = Solution(:,1)
    y(:,2) = Error(:,1)
    
    call plot_parametrics( Time, y(:,1:1), [" "], "$ t $  ", " $x$ ", "(a)", path(1) ) 
    call plot_parametrics( Time, y(:,2:2), [" "], "$ t $  ", "$ E $  ", "(b)", path(2) ) 
    
end subroutine





!********************************************************************** 
!  Convergence rate of numerical methods with the number of time steps  
!**********************************************************************
subroutine Convergence_rate_RK2_RK4
        
 
    real :: t0 = 0, tf = 30           ! Time domain for Van der Pol oscillator
    integer, parameter :: N = 1000    ! Time steps of the first grid 
    integer, parameter :: Nv = 2      ! Number of variables van de pol
    integer, parameter :: M = 10      ! number of time grids  
    real :: U(0:N,Nv)                 ! Solution
    real :: Time (0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 2      ! number of numerical schemes to be evaluated 
    real :: log_E(M, Np), log_N(M)    ! Error versus  time steps 
                                      ! first index: grid, second index: scheme 
    real :: order
   
    integer :: i 
    character (len=20) :: names(Np) = [ "RK2", "RK4"] 
    character (len=200) :: path(Np) = [ "./doc/chapters/Cauchy_problem/figures/Convergencea", &
                                        "./doc/chapters/Cauchy_problem/figures/Convergenceb"  ] 
    
    write(*,*) "Convergence rate: Error versus number of time steps"  
    write(*,*) "Temporal scheme : Rk2, Rk4 "   
    
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    U0=[3, 4]             !Initial conditions for van der pol
    U(0,:) = U0 
    
    call Cauchy_ProblemS( Time, VanDerPol_equation,  U, Runge_Kutta2 ) 
    
    call Temporal_convergence_rate( Time, VanDerPol_equation,  U0,      & 
                                    Runge_Kutta2, order, log_E(:,1), log_N )
    write(*,*) "Order Runge_Kutta2 = ", order 
    
    call Temporal_convergence_rate( Time, VanDerPol_equation,  U0,      &
                                    Runge_Kutta4, order, log_E(:,2), log_N )
     write(*,*) "Order Runge_Kutta4 = ", order 
       
    call plot_parametrics( U(:,1), U(:,2:2), [" "], "$x$", "$y$", "(a)", path(1)) 
    call plot_parametrics( log_N, log_E, names, "$\log N $", "$\log E$ ", "(b)", path(2)) 
    
  
end subroutine








subroutine Van_der_Pol_oscillator
        
    real :: t0 = 0, tf = 30
    integer, parameter :: N = 350, Nv = 2
    real :: Time (0:N), U(0:N, Nv, 2)
    integer :: i

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
    
    U(0,:,1) =  [3, 4]  
    call set_solver("eRK", "RK87")
    call set_tolerance(1d-8)
    call Cauchy_ProblemS( Time, VanDerPol_equation, U(:,:,1) ) 
    
    U(0,:,2) =  [0, 1] 
    call set_solver("eRK", "Fehlberg87")
    call set_tolerance(1d-8)
    call Cauchy_ProblemS( Time, VanDerPol_equation, U(:,:,2) ) 

    write (*, *) "VanDerPol oscillator with RK87 and Fehlberg87 "    
    write(*,*) "press enter ";  read(*,*)
    call plot_parametrics(time,U(:,1, 1:2),["RK87","Fehlberg87"],"t","x")

end subroutine  

  


subroutine Henon_Heiles_system
    integer, parameter :: N = 1000, Nv = 4 , M = 1 !Time steps    
    real, parameter :: dt = 0.1
    real :: t0 = 0, tf = dt * N
    real :: Time (0:N), U(0:N, Nv), H(0:N)
    integer :: i

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
            
    U(0,:) = [0., 0., 0.6,0.] 
    call set_solver("GBS")
    call set_tolerance(1d-2)
    call Cauchy_ProblemS( Time, Henon_equation, U )
    
    write (*, *) 'Henon Heiles system  '     
    write(*,*) "press enter " ; read(*,*) 
    call plot_parametrics(U(:,1),U(:,2:2), ["Henon Heiles"],"x","y")
    

end subroutine

!*****************************************************************************************
! plot a simple trayectory with Predictor-Corrector of variable step  
!*****************************************************************************************
subroutine Variable_step_with_Predictor_Corrector

        
   
    integer, parameter :: N = 200  !Time steps
    real :: Time(0:N), U(0:N,2)
    real :: t0 = 0, tf = 80
    integer :: i 

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
    U(0,1) =  1
    
    write(*,*) "Oscillator integrated with:"
    write(*,*) "  1) variable step Predictor-Corrector"  
    write(*,*) "  2) constant step Runge-Kutta 2"  
   
    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F,& 
                          Solution = U, Scheme = Predictor_Corrector1 )  
    call scrmod("reverse")
    call qplot(Time, U(:,1), N+1) 
    
    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F,& 
                          Solution = U, Scheme =  Runge_Kutta2 ) 
    call scrmod("reverse")
    call qplot(Time, U(:,1), N+1) 
  
contains

function F( U, t ) 
    real :: U(:), t 
    real :: F(size(U)) 
    
    F = [ U(2), - U(1) ]
    
  
end function 

end subroutine




!************************************************************ 
! Comparison of variable time step versus constant time step  
!************************************************************
subroutine Variable_step_simulation
        
    real :: t0 = 0, tf = 30           ! Time domain for Van der Pol oscillator
    integer, parameter :: N = 340     ! Time steps of the first grid 
    integer, parameter :: Nv = 2      ! Number of variables van de pol
    real :: Time (0:N)                ! Time domain 
    integer, parameter :: Np = 2      ! number of solutions to be evaluated 
    real :: U(0:N, Nv, Np)            ! Solutions for constant time step and variable time step 
    
    character (len=60) :: names(Np) = [ "const $\Delta t $ ", "var $\Delta t $"] 
    character (len=200) :: path(Np) = [ "./doc/chapters/Cauchy_problem/figures/ConstantDTa", &
                                        "./doc/chapters/Cauchy_problem/figures/ConstantDTb"  ] 
    integer :: i, j
    
    write(*,*) "Comparison constant step and variable step"   
    write(*,*) "Temporal scheme : Runge Kutta 2 order "   
    
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    do j=1, Np 
     U(0,:, j) = [3, 4]         !Initial conditions for VanDerpol equation 
    end do  
      
    call set_solver(family_name="eRK", scheme_name="HeunEuler21")
    call set_tolerance(1e10) 
    call Cauchy_ProblemS( Time, VanDerPol_equation, U(:,:,1) )
    
    call set_solver(family_name="eRK", scheme_name="HeunEuler21")
    call set_tolerance(1e-6) 
    call Cauchy_ProblemS(Time, VanDerPol_equation, U(:,:,2) )
    
    call plot_parametrics( Time, U(:,1,:),     names, "$t$", "$x$", "(a)", path(1) ) 
    call plot_parametrics( U(:,1,:), U(:,2,:), names, "$x$", "$y$", "(b)", path(2) ) 
    
    
end subroutine



!**************************************************************** 
!  It checks if Embedded Runge Kutta wrappers are properly built  
!****************************************************************
subroutine Convergence_rate_Runge_Kutta_wrappers 
        
   
    real :: t0 = 0, tf = 30          ! Time domain for Van der Pol oscillator
    integer, parameter :: N = 1000   ! Time steps of the first grid 
    integer, parameter :: Nv = 2     ! Number of variables van de pol
    integer, parameter :: M = 9      ! number of time grids           
    real :: Time (0:N), U0(Nv)       ! Time domain and initial conditions
    real ::  U(0:N, Nv)              ! Solution
    integer, parameter :: Np = 2     ! number of wrappers 
    real :: log_E(M, Np), log_N(M)   ! Error versus  time steps 
                                     ! first index: grid, second index: scheme 
    
    character (len=20) :: names(Np) = [ "WDOPRI5", "WDOP853"] 
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/ConvergenceWa", &
                                        "./doc/chapters/Cauchy_problem/figures/ConvergenceWb"  ] 
  
    integer :: i
    real :: order 
    
    write(*,*) "Convergence rate: Error versus number of time steps"   
    write(*,*) "Temporal scheme : DOPRI5, DOP853  "  
    
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    U0=[3, 4]       
    
    call set_solver(family_name="weRK",scheme_name="WDOPRI5")
    call set_tolerance(1e6)
    call Temporal_convergence_rate(                                      &
         Time_domain = Time, Differential_operator = VanDerPol_equation, & 
         U0 = U0,  order = order, log_E = log_E(:,1), log_N = log_N )
    write(*,*) "Order WDOPRI5 = ", order 
    
    call set_solver(family_name="weRK",scheme_name="WDOP853")
    call set_tolerance(1e6)
    call Temporal_convergence_rate(                                      &
         Time_domain = Time, Differential_operator = VanDerPol_equation, & 
         U0 = U0,  order = order, log_E = log_E(:,2), log_N = log_N )
    write(*,*) "Order WDOP853 = ", order 
    
    U(0,:) = U0 
    call Cauchy_ProblemS( Time, VanDerPol_equation, U, Runge_Kutta2 ) 
    call plot_parametrics( U(:,1), U(:,2:2), [" "], "$x$", "$y$", "(a)", path(1) ) 
    call plot_parametrics( log_N, log_E, names, " $\log N $", " $\log E $", "(b)", path(2) ) 
  

end subroutine









!******************************************************** 
!  It compares original codes with wrappers
!********************************************************
subroutine Runge_kutta_wrappers_versus_original_codes 

    real :: t0 = 0
    real :: tf = 17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 1500   ! Time steps of the first grid 
    integer, parameter :: Nv = 4     ! Number of variables 
    real :: Time(0:N), U0(Nv)        ! Time domain and initial conditions
    integer, parameter :: Np = 3     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
                                     ! first index: grid, second index: scheme 
    real :: x(0:N), y(0:N, 1) 

    character (len=40) :: names(Np) =    & 
                          [ "$\epsilon=10^{-3}$", "$\epsilon=10^{-4}$", "$\epsilon=10^{-5}$"]
    real :: tolerances(Np) = [ 1d-3, 1d-4, 1d-5 ]  
    integer :: No, i, j  
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/ArenstorfWRKa", &
                                        "./doc/chapters/Cauchy_problem/figures/ArenstorfRKb"  ] 
       
     write(*,*) "Arenstorf orbit "     
     write(*,*) "Comparison of original codes with wrappers(press enter)"   
     read  (*,*)
     
     call  dr_dopri5
     call read_file( "./results/orbit.plt", Time, x, y(:,1), No)
     call plot_parametrics( x(0:No-1), y(0:No-1, 1:1),["DOPRI5"],"x","y" ) 
 
     do j=1, Np 
       U(0,:,j) = [0.994, 0., 0., -2.0015851063790825 ]
     end do   
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ] 
     
     call set_solver(family_name="weRK", scheme_name="WDOPRI5")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names,           &
                            "$x$", "$y$", "(a)", path(1) ) 
    
     call set_solver(family_name="eRK", scheme_name="DOPRI54")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,          &
                            "$x$", "$y$", "(b)", path(2) )
    
end subroutine

!******************************************************** 
!  It compares original codes with wrappers
!********************************************************
subroutine GBS_and_wrapper_ODEX

    real :: t0 = 0
    real :: tf = 17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 1000   ! Time steps of the first grid 
    integer, parameter :: Nv = 4     ! Number of variables 
    real :: Time(0:N)                ! Time domain 
    integer, parameter :: Np = 3     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
                                     ! first index: grid, second index: scheme 
    real :: x(0:N), y(0:N, 1) 
    
    real :: tolerances(Np) =[ 1d-1, 1d-2, 1d-3 ]
    character (len=40) :: names(Np) =                                          & 
             [ "$\epsilon=10^{-1}$", "$\epsilon=10^{-2}$", "$\epsilon=10^{-3}$"]
    
    integer :: No, i, j  
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/GBSa", &
                                        "./doc/chapters/Cauchy_problem/figures/GBSb"  ] 
       
     write (*, *) 'Arenstorf orbit  '     
     write (*, *) 'GBS and ODEX  ' 
     
     do j=1, Np 
       U(0,:,j) = [0.994, 0., 0., -2.0015851063790825 ] 
     end do   
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
     
     call set_solver(family_name="wGBS")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names,            &      
                            "$x$", "$y$", "(a)", path(1) )
            
     call set_solver(family_name="GBS")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names,            &
                            "$x$", "$y$", "(b)", path(2) )
end subroutine


subroutine Arenstorf_with_GBS

    real :: t0 = 0
    real :: tf = 2*17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 50       ! Time steps of the first grid 
    integer, parameter :: Nv = 4       ! Number of variables 
    real :: Time(0:N)                  ! Time domain 
    real :: U(0:N, Nv)                 ! Solutions
    real :: eps 
    real :: t1 
    
    integer ::  i, j, N1 
    character(len=10) :: scheme(2) = [ "GBS", "wGBS"] 
    character(len=10) :: c(2) = [ "red", "white"]  
    real :: xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5  
   
     write (*, *) 'Arenstorf orbit with GBS, press enter '  
    
    
     U(0,:) = [0.994, 0., 0., -2.0015851063790825 ]  
     Time(0:N) = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
     
     eps = 1e-12
     call metafl("xwin")
     call scrmod("reverse")
     call disini  
     call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
   
     
     do j=1, 2 
       call set_solver( scheme(j) )
       call set_tolerance(eps)
     
       call Cauchy_ProblemS( Time, Arenstorf_equations, U )
       call color( c(j) ) 
       call curve( U(:,1), U(:,2), N+1) 
     end do 
     call disfin 
     
   
end subroutine



!******************************************************** 
!  It compares original codes with wrappers
!********************************************************
subroutine ABM_and_wrapper_ODE113

    real :: t0 = 0
    real :: tf = 17.0652165601579625 ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 8000   ! Time steps of the first grid 
    integer, parameter :: Nv = 4     ! Number of variables 
    real :: Time(0:N)                ! Time domain 
    integer, parameter :: Np = 3     ! number of graphs
    real :: U(0:N, Nv, Np)           ! Solutions
                                     ! first index: grid, second index: scheme 
    real :: x(0:N), y(0:N, 1) 

    real :: tolerances(Np) =[ 1d-3, 1d-4, 1d-5 ]
    character (len=40) :: names(Np) =                                          & 
             [ "$\epsilon=10^{-3}$", "$\epsilon=10^{-4}$", "$\epsilon=10^{-5}$"]
    
    integer :: No, i, j  
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/ABMa", &
                                        "./doc/chapters/Cauchy_problem/figures/ABMb"  ] 
       
     write (*, *) 'Arenstorf orbit '     
     write (*, *) 'ABM and ODE113  ' 
    
     do j=1, Np 
       U(0,:,j) = [0.994, 0., 0., -2.0015851063790825 ] 
     end do  
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
     
     call set_solver(family_name="wABM")
     do j=1, Np  
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,         &
                            "$x$", "$y$", "(a)", path(1) )
    
     call set_solver(family_name="ABM")
     do j=1, Np 
       call set_tolerance(tolerances(j))
       call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, j) )
     end do 
     call plot_parametrics( U(:, 1, :), U(:, 2, :) , names,        &
                               "$x$", "$y$", "(b)", path(2) )
     
end subroutine





!******************************************************** 
! Efficiency of high order methods. 
! Number of time steps as a function of error tolerance 
!********************************************************
subroutine Temporal_effort_with_tolerance_eRK
        
  
    real :: t0 = 0, tf = 30                 ! Time domain for Van der Pol oscillator
    integer, parameter :: N = 300           ! Time steps of the first grid 
    integer, parameter :: Nv = 2            ! Number of variables van de pol
    integer, parameter :: M = 9             ! number of different tolerances           
    real :: Time (0:N)                      ! Time domain 
    real :: U0(Nv)                          ! Initial conditions  
    real :: U(0:N,Nv)                       ! Solution 
    integer, parameter :: Np = 10            ! number of numerical schemes to be evaluated 
    real :: log_mu(M), log_effort(M,Np)      ! log time steps versus  log ( 1/ tolerance)
                                            ! first index: grid, second index: scheme  
    
    
    character (len=20) :: names(Np) = [ "HeunEuler21", "RK21", "BogackiShampine", "CashKarp", "RK87", "RK65", "DOPRI54", "Fehlberg54", "Fehlberg87", "Verner65"] 
    character (len=200) :: path(2) =  [ "./doc/chapters/Cauchy_problem/figures/Tstepsa", &
                                        "./doc/chapters/Cauchy_problem/figures/Tstepsb"  ] 
    integer :: i, j
    
    write(*,*) "Convergence rate: Error versus tolerance "   
    write(*,*) "Temporal schemes :  Embedded Runge Kutta "  
    
    U0 = [3, 4]
    Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ]
    log_mu = [( i, i=1, M ) ]
       
    do j=1, Np 
      call set_solver(family_name = "eRK", scheme_name = names(j) )
      call Temporal_effort_with_tolerance( Time, VanDerPol_equation, U0,& 
                                           log_mu, log_effort(:,j)  )
     end do 
    
     call plot_parametrics( log_mu, log_effort(:,1:3), names(1:3),      &
                         "$-\log \epsilon$", "$\log M$ ", "(a)", path(1))
     call plot_parametrics( log_mu, log_effort(:,4:7), names(4:7),      &
                         "$-\log \epsilon$", "$\log M$ ", "(b)", path(2)) 
end subroutine



















!***********************************************************************
!*  If gives the linearized operator F(U,t) in some point U0 
!***********************************************************************
subroutine test_System_matrix
          
          integer, parameter :: N = 4 
          real :: A(N,N),  U0(N)   
          integer :: j  
          real :: x, y, dx ,dy, t 
                  
          t = 0
          x = cos(t) 
          y = sin(t) 
          dx = -sin(t) 
          dy = cos(t) 
          
          U0  = [ x, y, dx, dy ]
          A = System_matrix( Kepler_F, U0, t) 
          
          do j=1, N 
             write(*,'(10f8.3)') A(j, :) 
          end do 
          
        
contains 
         
             
end subroutine  
 
!----------------------
! Kepler force 
!----------------------
function Kepler_F(U, t) result(F)
  real :: U(:), t 
  real :: F( size(U) ) 

  F = [ U(3), U(4), -U(1) / norm2( U(1:2) )**3 ,  -U(2) / norm2( U(1:2) )**3 ] 
   
end function  


   



!*************************************************************************
! read a tecplot file format x, y 
!*************************************************************************
 subroutine read_file(file_name, time, x, y, N  ) 
   character(len=*), intent(in) :: file_name
   real, intent(out) :: time(:), x(:), y(:)
   integer, intent(out) :: N 

      integer :: unit = 103 
      logical :: available 
      integer :: i, err  
   
     open(unit, file=file_name ) 
     read(unit, *) ! skip header 
     
     i = 1
     do while (1) 
     
          read(unit, iostat=err, fmt=*) time(i), x(i), y(i) 
                   
          if (err<0) exit 
          
          i = i + 1 
          
     end do 
     
     N = i - 1 
     
     close(unit)


 end subroutine
 






!******************************************************
!* Functions for Examples
!******************************************************
function Arenstorf_equations(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
     real :: mu = 0.012277471
     
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

function VanDerPol_equation( U, t ) result(F) 
                     real :: U(:), t 
                     real :: F(size(U)) 
    
    real :: mu = 5., x, v 
    
    x = U(1); v = U(2) 
    F = [ v, mu * (1 - x**2) * v  - x ] 
  
end function


function oscillator( U, t ) result(F) 
             real :: U(:), t 
             real :: F(size(U)) 
    
    F = [ -U(2), U(1) ] 
     
  
end function 

function Henon_equation( U, t )  result(F) 
                 real :: U(:), t 
                 real :: F(size(U)) 
                 
    real :: x, y, px, py, lambda = -1
  
    x = U(1) ; y = U(2); px = U(3); py = U(4) 
    
    F = [ px, py, -x-2* lambda*x*y, -y-lambda*(x**2 - y**2) ]
    
end function

real function Henon_Energy( q, p) result(H)
      real, intent(in) :: q(:) , p(:) 
      
      real :: lambda = -1 
      H = (norm2(q)**2 + norm2(p)**2) / 2  &
        + lambda * ( q(2)*q(1)**2 - q(2)**3 / 3) 
        
      
end function





end module



