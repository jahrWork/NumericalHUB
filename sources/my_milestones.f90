
module my_milestones

        use dislin 
        use Linear_systems
        use Cauchy_Problem
        use Temporal_Schemes
        use Fourier_Interpolation
        use Chebyshev_interpolation
        use plots 
        use Boundary_value_problems
        use Collocation_methods
        use Initial_Boundary_Value_Problems
        use Utilities
        use API_Example_Lagrange_Interpolation
        
        use Interpolation
        use Lagrange_interpolation
        use Legendre_points
        
        use API_Example_Cauchy_Problem
        use Numerical_Recipes
        implicit none 

       
    contains  
  


subroutine  milestone_examples


      call Milestone1_2A
      call Milestone2B
      call Milestone3
      call Milestone4
      call Milestone5 
      call Milestone6 
      call Milestone7
 
end subroutine






!***********************************************************
!* Milestone 1 y 2A
!***********************************************************
subroutine Milestone1_2A 
 
 integer, parameter :: N=25, M=400 
 real :: x(0:N), f(0:N)    ! N+1 given points
                           ! M+1 interpolated points 
 real :: I_N(0:N, 0:M), fe(0:2, 0:M), Error(0:2, 0:M) 
 real :: Lebesgue_N(-1:N, 0:M),  PI_N(0:N, 0:M)    
 real :: xp(0:M), a=-1, b=1, theta(0:N), alpha     
 integer :: i  
 
 x  = [ (a + (b-a)*i/N, i=0, N) ] ! N+1 points
 theta  = [ (PI*i/N, i=0, N) ]
 x = -cos( theta ) 
  
 xp = [ (a + (b-a)*i/M, i=0, M) ] ! M+1 points
 
 ! Runge function 1/(1+25x**2)  = 1/2 [ 1/(1+5xi) + 1/(1-5xi) ]
 ! d^nf/dx^n = n!  (5i)**n [ (-1)**n (1+5xi)**(-n)  + (1-5xi)**(-n) ] / 2 
 f = 1/( 1 + 25*x**2)  
 !f = sin ( PI * x ) 
 
 I_N = Interpolant(x, f, N, xp) 
 Lebesgue_N = Lebesgue_functions( x, xp ) 
 PI_N = PI_error_polynomial( x, xp ) 
 
 fe(0,:)  = 1/( 1 + 25*xp**2) 
 fe(1,:) =  -50 * xp / (1 + 25*xp**2)**2  
 fe(2,:) =  -50  / (1 + 25*xp**2)**2  +(50 * xp)**2 / (1 + 25*xp**2)**2  
 
 !fe(0,:) =  sin ( PI * xp ) 
 !fe(1,:) =  PI * cos ( PI * xp )
 !fe(2,:) =  -PI**2 * sin ( PI * xp ) 
 
 Error = fe - I_N(0:2, :)  
 
 write(*,*) " Plot results of milestone 1 y 2A " 
 write(*,*) " press enter "
 read(*,*) 
 call plot1
 
 
contains 


subroutine plot1

  real :: xmin, xmax, ymin, ymax
  xmin = a; xmax = b; ymin = -2; ymax = 2

  call metafl("xwin")
  call page(4000, 4000)
  call scrmod("reverse")
  
  call disini 
  call winfnt("Courier New Bold")
  call titlin( "Interpolated function", 3); 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("white");  call curve( xp, I_N(0,:), M+1)
  call incmrk(-1);  call marker(21);
  call color("red"); call curve( x, f, N+1)
  call color("white"); call height(80);call title 
  call disfin
  
  call disini 
  call winfnt("Courier New Bold")
  call titlin( "PI Error function", 3); 
  call setscl(PI_N(0:2,:), M+1, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("white");call curve( xp, PI_N(0,:), M+1)
  call color("red"); call curve( xp, PI_N(1,:), M+1)
  call color("blue");  call curve( xp, PI_N(2,:), M+1)
  call color("white"); call height(80);  call title 
  call disfin
 
  
  call disini 
  call winfnt("Courier New Bold")
  call titlin( "Lebesgue function", 3); 
  call axsscl("log","y"); call labels("log","y"); call labdig(-1, "y")
  call setscl(Lebesgue_N(0,:), M+1, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10 ) 
  call curve( xp, Lebesgue_N(0,:), M+1)
  call color("white"); call height(80);  call title 
  call disfin
  
  call disini 
  call winfnt("Courier New Bold")
  
  call titlin( "Interpolation Error", 3) 
  call setscl(Error, M+1, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10 ) 
  call curve( xp, Error(0,:), M+1)
  call color("red"); call curve( xp, Error(1,:), M+1)
  call color("blue");  call curve( xp, Error(2,:), M+1)
  call incmrk(-1);  call marker(21);
  call color("white"); call curve( x, f-f, N+1)
  call color("white"); call height(80);  call title 
  
  call disfin
  
 
end subroutine 
end subroutine 
 








!***********************************************************
!* Milestone 2B  
!***********************************************************
subroutine Milestone2B  
 
 integer, parameter :: N=4, M=400 
 real :: x(0:N), f(0:N)    ! N+1 given points
                           ! M+1 interpolated points 
 real :: I_N(0:N, 0:M), P_N(2, 0:M), fe(0:M), Error(2, 0:M) 
 real :: xp(0:M), a=-1, b=1, theta(0:N)
 real :: c1(0:N), c2(0:N), PI = 4 * atan(1d0) 
 integer :: i, k  
 
 theta  = [ (PI*i/N, i=0, N) ]
 x = -cos( theta ) 
 xp = [ (a + (b-a)*i/M, i=0, M) ] ! M+1 points
 
 fe = cos( PI * xp )  
 f = cos ( PI * x ) 
 I_N = Interpolant(x, f, N, xp) 
 
 c1 = Chebyshev_transform(f) 
 P_N(1, :) = Chebyshev_interpolant( c1, xp ) 
 
 c2 = 0
 do k=0, N, 2 
     c2(k) = c_hat(k/2) 
 end do 
 P_N(2, :) = Chebyshev_interpolant( c2, xp ) 
 
  Error(1,:) = fe - I_N(0,:)
  Error(2,:) = fe - P_N(2,:)
!  Error(1,:) = fe - P_N(1,:) 

 
 write(*,*) " Plot results of milestone 2B " 
 write(*,*) " press enter "
 read(*,*) 
 call plot_milestone2B
 
 
contains

!******************************************************************
!  Cehebyshev expansion. Truncated series. 
!       cos( pi x ) = sum c_m T_m(x) = sum c_m cos( m theta ) 
! 
!  Demonstration:  
!
!  exp( i pi x ) = cos( pi x ) + i sin( pi x ) 
!  Taylor expansion 
!  exp( i pi x ) = sum ( (i pi x)**n / n! ( odd plus even) 
!                = sum  ( (-1)**(2n) pi x )**(2n) / (2n)! + i ....
!
!  cos( pi x ) = sum (-1)**(n) ( pi x )**(2n) / (2n)!
!  Since x = cos( theta ) then, 
!
!  cos( pi x ) = sum (-1)**(n) (pi)**(2n)/(2n)! ( cos theta )**(2n) 
! 
!  Moivre  ( cos theta )**n = ( exp(i theta) + exp(-i theta) )**n / 2**n 
!
!   (x + y)**n = sum from k=0 to k=n  (n | k) x**(n-k) y**k 
!
!   ( cos theta )**(2n) = 1/2**(2n) sum ( 2n | k) exp( i theta( 2n-2k) ) 
!                       = 1/2**(2n-1) sum _0 ^{n} ( 2n | k ) cos ( (2n-2 k) theta ) 
!
! If m = 2n -2 k then, n = m/2 + k 
!
!    c_m = sum from n=m to infinity 2(-1)**n ( PI/2)**(2n) / ( (n-m)! (n+m)! ) 
!
!******************************************************************
function c_hat( m ) result(c) 
  integer, intent(in) :: m
  real :: c 
  
  integer :: n
  real :: f 
  
  c = 0 
  do n=m, 100 
              f =  gamma( real(n-m+1) ) * gamma( real(n+m+1) )
              c = c + 2 * (-1)**n * (PI/2)**(2*n) / f 
  end do 
  if (m==0) c = c/2
  
end function 


subroutine plot_milestone2B 

  real :: xmin, xmax, ymin, ymax
  xmin = a; xmax = b; ymin = -2; ymax = 2

  call metafl("xwin")
  call page(4000, 4000)
  call scrmod("reverse")
  
  call disini 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("white");  call curve( xp, I_N(0,:), M+1)
  call color("blue");   call curve( xp, P_N(1,:), M+1)
 ! call color("red");  call curve( xp, P_N(2,:), M+1)
 ! call color("orange");  call curve( xp, fe, M+1)
  
  call incmrk(-1);  call marker(21);
  call color("red"); call curve( x, f, N+1) 
  
  call plot_legends( [ "I_N", "P_N",  "cos(PI x)" ] ) 
  call plot_title( ["Comparison (N=4): Interpolation versus Truncated series", & 
                     "I_N: Chebyshev interpolant, P_N: Truncated series"] ) 
  call disfin
   
  call disini 
  call setscl(Error, M+1, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("blue");  call curve( xp, Error(1,:), M+1)
  call color("red");   call curve( xp, Error(2,:), M+1)
  
  call incmrk(-1);  call marker(21);
  call color("white"); call curve( x, f-f, N+1) 
  
  call plot_legends( ["E_I", "E_P", "Points"] ) 
  call plot_title( ["Error comparison: Interpolation versus Truncated series with  N = 4", & 
                    "E_I: Interpolation error, E_P: Truncated series error  " ] )
  call disfin
  
  
end subroutine 

end subroutine 






!***********************************************************
!* Milestone 3
!***********************************************************
subroutine Milestone3
 
 integer, parameter ::  M=1000
 real ::  xp(0:M)
 real, allocatable :: x(:), y(:), PI_N(:, :), Lebesgue_N(:,:) 
 real :: a=-1, b=1, PI = 4 * atan(1.) 
 integer :: i, k, N   
 
 
write(*,*) " Plot results of milestone 3 " 
write(*,*) " press enter "
read(*,*) 
 
xp = [ (a + (b-a)*i/M, i=0, M) ] ! M+1 points


do N=2, 8, 2 
      
   call dislin_ini(xmax = b, xmin= a, ymax = 15., ymin = 0.0) 
   
   call plot_title( ["Finite difference formulas. First and second derivative" ] )

   allocate( x(0:N), y(0:N), PI_N(0:N,0:M) ) 
   x  = [ (a + (b-a)*i/N, i=0, N) ] 
   y = 0 
 
   PI_N = PI_error_polynomial( x, xp ) 
   
   call color("red");  call marker(0); call incmrk(0) 
   call curve(xp, PI_N(1,:), M+1 )
   
   call color("blue"); call marker(0); call incmrk(0) 
   call curve(xp, PI_N(2,:), M+1 )
   
   
   call color("white");  call marker(21);  call incmrk(-1) 
   call curve(x, y, N+1 ) 
   
   deallocate(x, y, PI_N ) 
   call plot_legends( [ "d PI/dx", "d2 PI/dx",  "points" ] ) 
   call disfin
   
end do 

do N=2, 4, 2 
      
   call dislin_ini(xmax = b, xmin= a, ymax = 30., ymin = 0.0) 
   
   call plot_title( ["Finite difference formulas. First and second derivative" ] )

   allocate( x(0:N), y(0:N), Lebesgue_N(-1:N,0:M) ) 
   x  = [ (a + (b-a)*i/N, i=0, N) ] 
   y = 1 
 
   
   Lebesgue_N = Lebesgue_functions( x, xp ) 
   
   call color("white");  call marker(0); call incmrk(0) 
   call curve(xp, Lebesgue_N(0,:), M+1 )
   
   call color("red");  call marker(0); call incmrk(0) 
   call curve(xp, Lebesgue_N(1,:), M+1 )
   
   call color("blue"); call marker(0); call incmrk(0) 
   call curve(xp, Lebesgue_N(2,:), M+1 )
   
   
   call color("white");  call marker(21);  call incmrk(-1) 
   call curve(x, y, N+1 ) 
   
   deallocate(x, y, Lebesgue_N ) 
   call plot_legends( [ "Lebesgue_N", "d Lebesgue_N/dx", "d2 Lebesgue_N/dx",  "points" ] ) 
   call disfin
   
end do 
 
 
end subroutine 




!********************************************************************
!* Milestone4
!*****************************************************************
subroutine Milestone4
 
 integer :: q                      ! interpolant order 2, 4, 6, 8 
 integer :: Nq = 6                 ! max interpolant order 
 integer :: N                      ! # of nodes (piecewise pol. interpol.) 
 integer :: k = 2                  ! derivative order 
 integer :: p = 0                  ! where error is evaluated p=0, 1,...N
 integer, parameter :: M = 100     ! number of grids ( N = 10,... N=10**4) 
 real :: epsilon = 1d-12           ! order of the random perturbation 
 
 real :: PI = 4 * atan(1d0), logN  
 integer ::  j, l
 
 real, allocatable :: x(:), f(:), dfdx(:)       ! function to be interpolated 
 real, allocatable :: dIdx(:)                   ! derivative of the interpolant 
 real, allocatable :: log_Error(:,:),log_dx(:)  ! Error versus Dx for q=2, 4, 6, 8
 character(len=40) :: title 
 
 
 write(*,*) "Second derivative error versus spatial step for q=2,4,6,8 " 
 write(*,*) " Test function:  f(x) = cos pi x  " 
 write(*,*) " press enter "
 read(*,*) 
 
    
 
do k=1, 2 
 allocate(  log_Error(M,4), log_dx(M)   )  
 log_dx = 0; log_Error = 0   

 l = 0 
 do q=2, Nq, 2 
  l = l +1    
  do j=1, M 
     
   logN = 1 + 3.*(j-1)/(M-1)   
   N = 2*int(10**logN)
   
   allocate( x(0:N), f(0:N),  dfdx(0:N), dIdx(0:N)   ) 
   x(0) = -1; x(N) = 1 
 
   call Grid_Initialization( "uniform", "x", x, q ) 
   
   call random_number(f)
   f = cos ( PI * x ) + epsilon * f
   if (k==1) then 
       dfdx = - PI * sin ( PI * x ) 
       title = " First derivative "
   elseif (k==2) then 
       dfdx = - PI**2 * cos ( PI * x ) 
        title = " Second derivative "
   end if 
   
   call Derivative( "x", k, f, dIdx )
   
   log_dx(j) = log( x(1) - x(0) ) 
   log_Error(j, l) = log( abs(dIdx(p) - dfdx(p)) )    
 
   deallocate( x, f, dIdx, dfdx ) 
  end do 
 end do  
 
 call scrmod("reverse")
 call plot_parametrics( log_dx, log_Error, ["Error q=2", "Error q=4", "Error q=6"], & 
                       "log_dx","log_Error", title)
 
 deallocate( log_Error, log_dx ) 
end do 
 
end subroutine 



!********************************************************************
!* Milestone5
!*****************************************************************
subroutine Milestone5

    integer, parameter :: N = 40,  Nv = 1, Np = 3
    real :: x(0:N), U(0:N,Nv,Np)
    real :: x0 = -1 , xf = 1
    integer :: i, p 
    real :: pi = 4*atan(1.) 
   
    write (*, *) 'Milestone5: Solution of boundary value problems ' 
    write (*, *) ' yxx + exp(-x**2)  yx + - y = 100 sin( pi x) sin(5 pi x) ' 
    write (*, *) ' y(-1) = 0, dydx(+1) = 0 '
    write(*,*) " press enter "
    read(*,*) 
      
    
 x(0) = x0; x(N) = xf  
 do p=1, Np 
    call Grid_Initialization( grid_spacing = "nonuniform", &
                              direction = "x",   q = 2*p, nodes = x )
!   Legendre solution   
    call Boundary_Value_Problem( x = x,                                & 
                                 Differential_operator = ODES,     & 
                                 Boundary_conditions   = BCs, & 
                                 Solution = U(:,:,p) )
 end do 
 
call plot(x, transpose(U(:,1,:)),"Milestone5. BVP with order q=2,4,6"  )  
 
contains 


!****** Differential operator *********
function ODES(x, y, yx, yxx) result(L)
   real, intent(in) :: x, y(:), yx(:), yxx(:)   
   real :: L(size(y)) 
   
    integer :: n = 6
        
    L =  yxx + exp(-x**2) * yx - y - 100*sin(pi*x) * sin(5*pi*x)
       
end function 
    
!********* Boundary conditions *********
function BCs(x, y, yx) 
           real, intent(in) :: x, y(:), yx(:)   
           real :: BCs(size(y))

        if (x==x0) then
                           BCs = y
        else if (x==xf ) then
                           BCs = yx                    
        else 
            write(*,*) " Error BCs x=", x; stop  
        endif            
                 
end function  

end subroutine 


!********************************************************************
!* Milestone6
!*****************************************************************
subroutine Milestone6

   integer, parameter :: N = 1000, Nv=2, Np=2
   real :: Time(0:N), U(0:N, Nv, Np) 
   real :: a=10., b=28., c=2.6666666666
   real :: t0 =0, tf=25 
   integer :: i, p 
 
   Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
  
   do p=1, Np  
    U(0,:, p) = [0, 1]
    if (p==1) then    
        call Cauchy_ProblemS( Time_Domain=Time, Differential_operator=ODE, & 
                              Solution = U(:,:,p), Scheme = Euler )
    else if (p==2) then 
        call Cauchy_ProblemS( Time_Domain=Time, Differential_operator=ODE, & 
                              Solution = U(:,:,p), Scheme = Runge_Kutta4 )
    end if 
        
   end do 
    
   write (*, *) 'Solution of Milestone 6  '  
   write(*,*) " press enter "
   read(*,*) 
   call plot_parametrics(Time, U(:,1,:), ["Euler", "RK4"],"t","x", "Milestone6: Cauchy problem with different time schemes")
   
   call Stability_regions_Euler_RK4  
   
contains

function ODE(U, t) result(F)           
     real :: U(:),t
     real :: F(size(U))
     
     real :: x, dxdt
      
     x = U(1); dxdt = U(2);   

     F = [ dxdt, cos( 2*t ) - sin(x) ]
     

 end function
 
end subroutine

subroutine Stability_regions_Euler_RK4
   integer, parameter ::N = 50
   integer :: i, j  
   real :: x(0:N), y(0:N), Region(0:N,0:N)
   real :: x0 = -4d0, xf = 1d0, dx
   real :: y0 = -4d0, yf = 4d0, dy
   
   integer, parameter :: N_levels = 9 
   real :: levels(0:N_levels)  
   
   dx = (xf-x0) / N; dy = (yf-y0) / N
   x = [ ( x0 + dx * i  , i=0,N )]
   y = [ ( y0 + dy * j  , j=0,N )]
   
   levels = [ ( j/real(N_levels)  , j=0, N_levels )]
   do j=1, 2  
     if (j==1)      then 
         call Absolute_Stability_Region(Euler, x, y, Region) 
     else if (j==2) then 
         call Absolute_Stability_Region(Runge_Kutta4, x, y, Region) 
     end if 
     
     call plot_contour(x, y, Region, "Re(z)","Im(z)", levels, graph_type ="isolines"   )  
   end do 
   

end subroutine 



!********************************************************************
!* Milestone7
!*****************************************************************
subroutine Milestone7

       integer, parameter :: Nx = 20, Nt = 1000, Nv = 1
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv)  
       
       real ::  x0 = 0, xf = 1, t0 = 0, tf = 0.2 
       integer :: i, j, k, q = 2 
       integer, parameter :: Nl = 5 
       character(len=10) :: legends(0:Nl) 
     
       
     write (*, '(A50)') 'Milestone 7: Time solution of the 1D heat equation'
     write(*,*) " press enter "
     read(*,*) 
     
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ] 
     do i=0, Nl; legends(i) = "t="//float_to_str(time(Nt/Nl*i)) ; enddo 
         
     x(0) = x0; x(Nx) = xf
     
!    Heat equation 1D  
     call Grid_Initialization( "nonuniform", "x", x, q )
       
     U(0, :, 1)  =  exp(-25*(x-0.5)**2 )
     call Initial_Boundary_Value_Problem(                              & 
                       Time_Domain = Time, x_nodes = x,                & 
                       Differential_operator =  Heat_equation1D,       & 
                       Boundary_conditions   =  Heat_BC1D,             & 
                       Solution = U ) 
     
     call plot_parametrics(x, transpose(U(0:Nt:200,:,1)),           & 
                   legends, "$x$", "$u(x,t)$", "Milestone 7: Heat equation") 
     
     call Stability_Heat_equation_1Dm
     call Error_Heat1D
   

end subroutine 


function Heat_equation1D( x, t, u, ux, uxx) result(F)
          real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
          real :: F( size(u) ) 
            
            F(1) =   uxx(1)
end function 

function Heat_BC1D(x, t, u, ux) result(BC) 
    real, intent(in) :: x, t, u(:), ux(:) 
    real :: BC( size(u) ) 
    
    real ::  x0 = 0, xf = 1 
    
        if (x==x0) then
                            BC(1) = u(1) 
        else if (x==xf) then
                            BC(1) = u(1)  
        else
             write(*,*)  "Error in Heat_BC1D"; stop 
        endif
end function     
     


!******************************************************************
!   Eigenvalues of the heat equation with FD and Chebyshev 
!
!           d2u/dx2 = lambda u    with u(0)=0 and u(1)=0 
!
!    exact eigenfunctions u_k(x) = sin( k pi x ) with 
!          eigenvalues lambda_k = - ( k pi )**2   
!******************************************************************
subroutine Stability_Heat_equation_1Dm

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
    do q=2,  Nx, 2 
       call Grid_Initialization( "nonuniform", "x", x, q )
       A =  Linear_operator( 1, x, q, Heat_equation1D, Heat_BC1D )
                
      call Eigenvalues_QR( A, lambda(q,:)) 
      rlambda(q,:) = -Real( lambda(q,:) ) 
      call sortr1( rlambda(q,:), Nx+1, "A") 
      !write(*,*) "lambda =", lambda(q,:) 
      !read(*,*) 
      
    end do 
    
    x_R = [ ( -4 + dx_R * i  , i=0,N )]
    x_I = [ ( -4 + dx_I * j  , j=0,N )]    
    call Absolute_Stability_Region(Runge_Kutta4, x_R, x_I, Region) 
    
    call plot_eigenvalues
    call plot_stability_region 
    call plot_interpolation 
  
           
contains 


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
    call plot_legends( ["q=2","q=4", "q=6", "q=N(Chebyshev)","Exact eigenvalues" ] )
    
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



subroutine plot_interpolation

 integer, parameter :: N=20, M=500 
 real :: x(0:N), f(0:N)    ! N+1 given points
                           ! M+1 interpolated points 
 real :: xp(0:M), fe(0:M), a=0, b=1, y(0:M, 4)    
 integer :: i, j, k, q(3) = [2, 6, N]
 real, allocatable :: I_N(:,:) 
    
 x  = [ (a + (b-a)*i/N, i=0, N) ] ! N+1 points
 xp = [ (a + (b-a)*i/M, i=0, M) ] ! M+1 points
 
 k = 15
 fe  = sin ( k * PI * xp ) 
 
 do i=1, 3
 call Grid_Initialization( "nonuniform", "x", x, q(i)) 
 f = sin ( k * PI * x ) 
 
 allocate( I_N(0:q(i), 0:M) ) 
 I_N = Interpolant(x, f, q(i), xp) 
 y(:, i) = I_N(0,:) 
 deallocate(I_N) 
end do 
 
 y(:, 4) = fe 
 
 call plot1
 call plot_parametrics(xp, y, ["q=2", "q=6", "q=N", "sin(k pi x)" ], "x", "y", title = "Interpolation with N=20 and k=15")
 
 
 contains
 subroutine plot1

  real :: xmin, xmax, ymin, ymax
  xmin = a; xmax = b; ymin = -2; ymax = 2

  call metafl("xwin")
  call page(4000, 4000)
  call scrmod("reverse")
  
  call disini 
  call winfnt("Courier New Bold")
  call titlin( "Interpolation of sin(k PI x): Chebyshev(blue),  q=6(orange), N=20, k=15", 1) 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  call color("white");  call curve( xp, y(:,4), M+1) 
  call color("blue");  call curve( xp, y(:,3), M+1)
  call color("orange");  call curve( xp, y(:,2), M+1)
  
  call incmrk(-1);  call marker(21);
  call color("red"); call curve( x, f, N+1)
  
 
  call color("white"); call height(20);call title 
  call disfin
  
 
   
end subroutine 
 
end subroutine





subroutine Error_Heat1D

       integer, parameter :: Nx = 40,  Nv = 1
       real ::  x(0:Nx), U(0:Nx, Nv),  F(0:Nx, Nv), t   
       real ::  R(3, 0:Nx, Nv)
       real ::  x0 = 0, xf = 1
       integer :: i, q, qmax 
  
     do i=1, 2 
         
 !    Polynomial Order=2,4,6, ... Nx           
      q = 2*i 
      
 !    Spatial discretization         
      x(0) = x0; x(Nx) = xf; t = 0 
      call Grid_Initialization( "nonuniform", "x", x, q )
      U = Test_U(Nv, x) 
      F = Spatial_discretization( Heat_equation1D, Heat_BC1D, x, U, t ) 
      
 !    Spatial Truncation Error 
      R(i,:,:) = Spatial_Truncation_Error( Nv, Heat_equation1D, Heat_BC1D, x, q, Test_U  )   
      
     end do  
     R(3,:,:) = Test_Uxx(Nv, x) - F
   
!   call plot(x, F(:,1), "Spatial discretization Heat1D F = Uxx")     
    call plot(x, U(:,1), "Test condition Heat1D U ")   
    call plot(x(1:Nx-1), R(:, 1:Nx-1,1), "Spatial truncation error R with q=2,4 (Richardson & Exact)")
      
     
contains     

function Test_U( Nv, x) result(U)
      integer, intent(in) :: Nv 
      real, intent(in) ::  x(0:) 
      real :: U( 0:size(x)-1, Nv ) 
            
   U(:,1)  =   exp(-80*(x-0.5)**2) 
          
end function 

function Test_Uxx( Nv, x) result(Uxx)
      integer, intent(in) :: Nv 
      real, intent(in) ::  x(0:) 
      real :: Uxx( 0:size(x)-1, Nv ) 
            
 ! U(:,1)   =   exp(-80*(x-0.5)**2) 
 ! Ux(:,1)  =  -160* (x-0.5)  * exp(-80*(x-0.5)**2)  
   Uxx(:,1) = -160 * exp(-80*(x-0.5)**2) -160*(x-0.5) * exp(-80*(x-0.5)**2)  *(-160)* (x-0.5)
          
end function 

end subroutine 






end module  

