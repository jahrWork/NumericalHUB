module API_Example_Finite_Differences

   
    !use Finite_differences
    use Collocation_methods 
    use Non_Linear_Systems
    use Linear_systems
    use plots 
    
    implicit none
    
contains  
  

subroutine Finite_difference_examples

   
   call Derivative_function_x
   call Derivative_function_xy
   call Derivative_error   
  

end subroutine

subroutine Derivative_function_x

    integer, parameter :: Nx = 20, Order = 4
    real :: x(0:Nx)
    real :: x0 = -1, xf = 1
    integer :: i
    real :: pi = 4 * atan(1.) 
    real :: u(0:Nx), uxk(0:Nx, 2),  ErrorUxk(0:Nx, 2)
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]   
    
    call Grid_Initialization( "nonuniform", "x", x, Order )
         
    u = sin(pi * x) 
       
    call Derivative( 'x' , 1 , u , uxk(:,1) ) 
    call Derivative( 'x' , 2 , u , uxk(:,2) )  


    ErrorUxk(:,1) = uxk(:,1) - pi* cos(pi * x) 
    ErrorUxk(:,2) = uxk(:,2) + pi**2 * u
    
    write (*, *) 'Finite differences formulas: 4th order ' 
    write (*, *) 'First and second derivative of sin pi x  ' 
    write(*,*) "press enter ";   read(*,*)
   
    call plot_parametrics(x, Uxk, ["ux", "uxx"], "x", "y")  
     

end subroutine



!***********************************************************************
!
!***********************************************************************

subroutine Derivative_function_xy

    integer, parameter :: Nx = 20, Ny = 20, Order = 6
    real :: x(0:Nx), y(0:Ny)
    real :: x0 = -1, xf = 1, y0 = -1, yf = 1
    integer :: i, j
    real :: pi = 4 * atan(1.0) 
    real :: u(0:Nx,0:Ny), uxx(0:Nx,0:Ny), uy(0:Nx,0:Ny), uxy(0:Nx,0:Ny)
    real :: Erroruxx(0:Nx,0:Ny), Erroruxy(0:Nx,0:Ny)
    
    x = [ (x0 + (xf-x0)*i/Nx, i=0, Nx) ]
    y = [ (y0 + (yf-y0)*j/Ny, j=0, Ny) ]    
    
    call Grid_Initialization( "nonuniform", "x", x, Order )
    call Grid_Initialization( "nonuniform", "y", y, Order )
    
    u = Tensor_product( sin(pi*x), sin(pi*y) )
   
    call Derivative( ["x", "y"], 1, 2, u, uxx )  
    call Derivative( ["x", "y"], 2, 1, u, uy ) 
    call Derivative( ["x", "y"], 1, 1, uy, uxy ) 

    Erroruxx = uxx + pi**2 * u
    Erroruxy = uxy - pi**2 * Tensor_product( cos(pi*x), cos(pi*y) )
      
    write (*, *) '2D Finite differences formulas: 6th order ' 
    write (*, *) 'Second partial derivative with respect x' 
    write (*, *) 'of  u(x,y) = sin pi x sin pi y  '
    write(*,*) "press enter ";  read(*,*) 
   
    call plot_contour(x, y, uxx, "x", "y" )
    
end subroutine

    


!********************************************************************
!* 
!*****************************************************************
subroutine Derivative_error
 
 integer :: q                      ! interpolant order 2, 4, 6, 8 
 integer :: Nq = 8                 ! max interpolant order 
 integer :: N                      ! # of nodes (piecewise pol. interpol.) 
 integer :: k = 2                  ! derivative order 
 integer :: p = 0                  ! where error is evaluated p=0, 1,...N
 integer, parameter :: M = 100     ! number of grids ( N = 10,... N=10**4) 
 real :: log_Error(M,4),log_dx(M)  ! Error versus Dx for q=2, 4, 6, 8
 real :: epsilon = 1d-12           ! order of the random perturbation 
 
 real :: PI = 4 * atan(1d0), logN  
 integer ::  j, l=0
 
 real, allocatable :: x(:), f(:), dfdx(:)  ! function to be interpolated 
 real, allocatable :: dIdx(:)              ! derivative of the interpolant 
 
 do q=2, Nq, 2 
  l = l +1    
  do j=1, M 
     
   logN = 1 + 3.*(j-1)/(M-1)   
   N = 2*int(10**logN)
   
   allocate( x(0:N), f(0:N),  dfdx(0:N), dIdx(0:N) ) 
   x(0) = -1; x(N) = 1 
 
   call Grid_Initialization( "uniform", "x", x, q ) 
   
   call random_number(f)
   f = cos ( PI * x ) + epsilon * f
   dfdx = - PI**2 * cos ( PI * x ) 
   
   call Derivative( "x", k, f, dIdx )
   
   log_dx(j) = log( x(1) - x(0) ) 
   log_Error(j, l) = log( abs(dIdx(p) - dfdx(p)) ) 
 
   deallocate( x, f, dIdx, dfdx ) 
  end do 
 end do 
 
 call scrmod("reverse")
 write(*,*) "Second derivative error versus spatial step for q=2,4,6,8 " 
 write(*,*) " Test function:  f(x) = cos pi x  " 
 write(*,*) "press enter " ; read(*,*)
 
 call plot_parametrics( log_dx, log_Error, ["E2", "E4", "E6", "E8"], & 
                       "log_dx","log_Error")
 
end subroutine 




end module 
