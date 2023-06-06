module API_Example_IBVP_Chebyshev

  
    use API_Example_Chebyshev_interpolation
    use Initial_Boundary_Value_Problems
    use Temporal_Schemes
    use Collocation_methods
    use Temporal_error
    use plots
    use Interpolation
    use Lagrange_interpolation
    use Stability
    use Stability_regions
    use Linear_systems
    use Initial_Boundary_Value_Problem1D
    use Utilities
    use Numerical_Recipes
    
implicit none

contains

subroutine Chebyshev_problems
integer :: option 

option = 1            
do while (option>0) 
    
     write(*,*) "Welcome to Advanced methods" 
     
     write(*,*) " select an option " 
     write(*,*) " 0. Exit/quit  "
     write(*,*) " 1. Chebyhev interpolations  "
     write(*,*) " 2. Chebyshev spectral derivative versus FD derivatives  "
     write(*,*) " 3. Heat equation and Wave equation with FD2 and Chebyshev collocation" 
     write(*,*) " 4. Advection diffusion equation with Chebyshev collocation"
     
     read(*,*) option 
     
     select case(option)
     case(1) 
         call Chebyshev_interpolation_examples
         
     case(2)  
         call Derivative_comparison_Chebyshev  
         
     case(3) 
         call Heat_and_wave_equations_with__high_order
           
     case(4)
         call Advection_diffusion_1D_Chebyshev
         
     case default
              
     end select 
     
end do
    
end subroutine     
    
    

subroutine Heat_and_wave_equations_with__high_order

       integer, parameter :: Nx = 10, Nt = 1000, M = 200, Nv = 1
       real ::  x(0:Nx)
       real :: Time(0:Nt), U(0:Nt,0:Nx, Nv)  
       real :: xp(0:M), yp(0:Nx,0:M), y_exact(0:M)  
       real, parameter :: PI = 4 * atan(1.) 
       real ::  x0 = -1, xf = 1, t0 = 0, tf = 2.0
       integer :: i, k, q = Nx, problem 
       integer, parameter :: Nl = 5 
       character(len=10) :: legends(0:Nl) 
       real :: w = 2 * PI 
       
     xp = [ (x0 + (xf-x0)*i/M, i=0, M ) ]   
     Time = [ (t0 + (tf-t0)*i/Nt, i=0, Nt ) ] 
     do i=0, Nl; legends(i) = "t="//float_to_str(time(Nt/Nl*i)) ; enddo 
     x(0) = x0; x(Nx) = xf 
     
      
     do problem=1, 2; do q=2, Nx, Nx-2
         
       call Grid_Initialization( "nonuniform", "x", x, q )
       U(0, :, 1)  =  0
     
       call Initial_Boundary_Value_Problem( Time, x, Problem_eq,  Problem_BC,  U ) 
     
       yp = Interpolant( x, U(Nt,:,1), Nx, xp)
       call plot  
       
     end do; end do 
     
     
contains
 
function Problem_eq( x, t, u, ux, uxx) result(F) 
               real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
               real :: F(size(u)) 
            
    real :: c, nu 
     
    if (problem==1) then 
          c = 0; nu = 1 
          
    elseif (problem==2) then 
         c = 1; nu = 0 
    end if 
            
   F = - c * ux + nu * uxx
   
end function 

!-------------------------------------------------------
function Problem_BC(x, t, u, ux) result(BC) 
         real, intent(in) :: x, t, u(:), ux(:) 
         real :: BC( size(u) ) 
  
        if (x==x0 ) then
                            BC = u - sin( w*t )
        else  if (x==xf) then
            
                     if (problem==1) then 
                         BC = u 
                     else if (problem==2) then 
                         BC = FREE_BOUNDARY_CONDITION 
                     end if  
        else
             write(*,*)  "Error in Advection_BC1D"; stop
        endif
end function
    
subroutine  plot

    real :: xmin, xmax, ymin, ymax
    integer :: j, N 
  
    xmin = minval(x); xmax = maxval(x) 
    if (problem == 1) then 
         ymax = +0.5; ymin = -0.5;
    else if (problem==2) then 
        ymax = +1; ymin =-1;
    end if 
    call dislin_ini(xmax, xmin, ymax, ymin) 
      
    call color("red"); call incmrk(1);   call marker(21);  
    call curve( x, U(Nt,:,1), Nx+1 ) 
    
    call color("blue"); call incmrk(1);   call marker(-1);  
    call curve( xp, yp(0,:), M+1 )
    
    
    if (problem==1) then 
        
     call plot_legends(["t=2 (nodal points)", "t=2 (interpolated)" ])   
     if (q==2) then 
        
        call plot_title( ["Heat equation with Finite differences 2th order",  & 
                          " du/dt = d2u/dx2 = 0 with N=10(sampling rule : N > 4 k)", & 
                          " u(-1,t) = sin ( 2 * pi t ), u(+1,t) = 0  "  ] ) 
      else if (q==Nx) then 
   
        
      call plot_title( ["Heat equation with Chebyshev interpolant",  & 
                       " du/dt = d2u/dx2 = 0 with N=10(sampling rule : N > 4 k)", & 
                       " u(-1,t) = sin ( 2 * pi t ), u(+1,t) = 0 "  ] ) 
      end if   
      
    else if (problem==2) then 
        
      call plot_legends(["t=2 (nodal points)", "t=2 (interpolated)", "t=2 (exact)"])     
     
      y_exact = sin( w*(tf-xp) )   
      call color("white"); call incmrk(1);   call marker(-1);  
      call curve( xp, y_exact, M+1 )  
      if (q==2) then 
        
        call plot_title( ["Wave equation with Finite differences 2th order",  & 
                          " du/dt + du/dx = 0 with N=10(sampling rule : N > 4 k)", & 
                          " u(-1,t) = sin ( 2 * pi t ), x=+1 outflow "  ] ) 
      else if (q==Nx) then 
        
      call plot_title( ["Wave equation with Chebyshev interpolant",  & 
                       " du/dt + du/dx = 0 with N=10(sampling rule : N > 4 k)", & 
                       " u(-1,t) = sin ( 2 * pi t ), x=+1 outflow"  ] ) 
      end if 
    end if
    
    call disfin 

end subroutine   
end subroutine 




!**************************************************************
!* Derivative comparison: Chebyshev versus FD
!**************************************************************
 subroutine Derivative_comparison_Chebyshev

    integer, parameter :: Nx = 20
    real :: x(0:Nx), y(0:Nx), yx(2, 0:Nx)  
    real ::  x0 = 0, xf = 1 
    integer :: i, Order
    real, parameter :: PI = 4 * atan(1.)
         
     x(0) = x0; x(Nx) = xf
     
     Order = Nx 
     call Grid_Initialization( "nonuniform", "x", x, Order )
     y = sin(PI*x) 
     call Derivative( "x", 1, y, yx(1, :) )  
     
     call Grid_Initialization( "Chebyshev", "x", x )
     call Derivative( "x", 1, y, yx(2, :) )  
     
     write(*,*) " maxval yx Cheb =", maxval(yx(2,:)-yx(1,:)) 
     call plot(x, yx, "First derivative: Nonuniform Order=Nx versus Chebyshev extrema") 
     
     call Grid_Initialization( "nonuniform", "x", x, Order )
     y = sin(PI*x) 
     call Derivative( "x", 2, y, yx(1, :) )  
     
     call Grid_Initialization( "Chebyshev", "x", x )
     call Derivative( "x", 2, y, yx(2, :) )  
     
     write(*,*) " maxval yx Cheb =", maxval(yx(2,:)-yx(1,:)) 
     call plot(x,  yx, "Second derivative: Nonuniform Order=Nx versus Chebyshev extrema") 
    
end subroutine 



!******************************************************************
! Chebyshev with collocation ponits and with Chebyshev expansions 
!******************************************************************
 subroutine Advection_diffusion_1D_Chebyshev

    integer, parameter :: Nx = 32, Nt = 100000, Nv = 1 
    real :: x(0:Nx), Time(0:Nt), U(0:Nt,0:Nx, Nv)  
    real ::  x0 = 0, xf = 1, t0 = 0, tf = 10.0
    integer :: i, Order = Nx 
    real, parameter :: PI = 4 * atan(1.), w = 16 * PI, nu = 0.0005
         
     Time = [ (t0 + (tf-t0)/Nt*i, i=0, Nt ) ] 
     x(0) = x0; x(Nx) = xf
     U(0, :, 1)  = 0
               
     x = [ ( (x0+xf)/2 + (x0-xf)/2 *cos( i*PI/Nx), i=0, Nx) ]  
     call Grid_Initialization( "unmodified", "x", x, Order )
     call Initial_Boundary_Value_Problem( Time, x, AD_eq, AD_BC, U ) 
     call plotN("Chebyshev")   
    
     call Grid_Initialization( "Chebyshev", "x", x )
     call Initial_Boundary_Value_Problem( Time, x, AD_eq, AD_BC, U ) 
     call plotN("Chebyshev") 
     
     Order = 12 
     call Grid_Initialization( "nonuniform", "x", x, q=Order )
     call Initial_Boundary_Value_Problem( Time, x, AD_eq, AD_BC, U ) 
     call plotN("FD12")      
        
contains 


function AD_eq( x, t, u, ux, uxx) result(L)
        real, intent(in) ::  x, t, u(:), ux(:), uxx(:)
        real :: L(size(u)) 
            
         L(1) = -ux(1) + nu * uxx(1) 
                 
end function 

function AD_BC(x, t, u, ux) result(BC) 
  real, intent(in) :: x, t, u(:), ux(:)  
  real :: BC(size(u)) 
  
        if (x==x0) then
                            BC(1) = u(1) - sin( w * t ) 

        else if (x==xf) then
                            BC(1) = FREE_BOUNDARY_CONDITION 
        else
             write(*,*)  "Error in AD_BC x=", x; stop 
             
        endif

end function

subroutine  plotN(tit)
    character(len=*), intent(in) :: tit 
    real :: xmin, xmax, ymin, ymax
    integer :: j, N 
    integer, parameter ::  M = 200 
    real ::  U_exact(0:M), xp(0:M), Up(0:Order, 0:M)
    real :: alpha 
    
   
    xmin = minval(x); xmax = maxval(x) 
    ymax = 1; ymin = -1
    call dislin_ini(xmax, xmin, ymax, ymin)  
  
    xp = [ (x0 + (xf-x0)*i/M, i=0, M ) ] 
    Up(0:Order,:) = Interpolant(x, U(Nt,:,1), Order, xp)  
    
    alpha = ( -1/nu + sqrt( 1/nu**2 + 4 * w**2 ) )/2
    U_exact = exp( - alpha * xp  )  * sin( w*(tf-xp) ) 
    
    call color("red"); call incmrk(1);   call marker(21);  
    call curve( x, U(Nt,:,1), Nx+1 )
    
    call color("blue"); call incmrk(1);   call marker(-1);  
    call curve( xp, Up(0,:), M+1 )
    
    call color("white"); call incmrk(1);   call marker(-1);  
    call curve( xp, U_exact, M+1 )
    
    call plot_title( ["Advection diffusion equation with collocation Chebyshev",  & 
                        " du/dt + du/dx = nu d2u/dx2 with N=32(sampling rule : N > 4 k)", & 
                        " u(0,t) = sin ( w t ), x=+1 outflow with w=16*PI, nu=0.0005 "  ] ) 
    
    call plot_legends( [tit//"(N=32)", tit // " interpolated",  "Exact solution " ] ) 
    call disfin 

end subroutine 
end subroutine 
    
    
 





end module
    