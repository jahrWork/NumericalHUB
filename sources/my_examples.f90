
module my_examples

        use dislin 
        use Linear_systems
        use Cauchy_Problem
        use Temporal_Schemes
        use Fourier_Interpolation
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
        use API_Example_Chebyshev_Fourier_interpolation
        implicit none 

       
    contains  

subroutine  myExamples

      call myexampleA
      call myexampleB
      call myexampleC
      call myexampleD
      call myExampleE
      
      call Embedded_Runge_kutta_Arenstorf_orbit
      
      call Legendre_points_example      
       
      call Chebyshev_Fourier_interpolation_examples
      


end subroutine  
    


!*****************************************************************************************
! Plot a simple graph 
!*****************************************************************************************
subroutine myexampleA
 
 integer, parameter :: N=200
 real :: x(0:N), y(0:N)
 integer :: i 
 real, parameter :: PI = 4 * atan(1d0)
 real :: a = 0, b = 2 * PI 
 
   x  = [ (a + (b-a)*i/N, i=0, N) ] 

   y = sin( 2 * x ) 
 
   call scrmod("reverse")
   call qplot(x, y, N+1) 
 
end subroutine


!*****************************************************************************************
! plot a simple trayectory 
!*****************************************************************************************
subroutine myexampleB

        
   
    integer, parameter :: N = 200  !Time steps
    real :: Time(0:N), U(0:N,1)
    real :: t0 = 0, tf = 8
    integer :: i 

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
    U(0,1) =  1
   
    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F,& 
                          Solution = U, Scheme = Euler )
    
    call scrmod("reverse")
    call qplot(Time, U, N+1) 
  
contains

function F( U, t ) 
    real :: U(:), t 
    real :: F(size(U)) 
    
    F(1) = - U(1)
  
end function 

end subroutine
 
!*****************************************************************************************
! Plot a simple graph and generates latex file in a specific location 
!*****************************************************************************************
subroutine myexampleC 
    real, parameter :: PI = 4 * atan(1d0)
    integer, parameter :: N=200, Np = 3 
    real :: x(0:N), y(0:N, Np), a = 0, b = 2 * PI 
    integer :: i 
    character(len=100) :: path(4) =                    & 
    ["./results/myexampleCa", "./results/myexampleCb", & 
     "./results/myexampleCc", "./results/myexampleCd"  ]
    x  = [ (a + (b-a)*i/N, i=0, N) ] 
    y(:, 1)  = sin(x); y(:, 2)  = cos(x); y(:, 3)  = sin(2*x)

   call plot_parametrics( x, y, ["$\sin x$", "$\cos x$", "$\sin 2x$"], & 
                         "$x$", "$y$", "(a)", path(1) ) 
   call plot_parametrics( y(:,1), y(:,:), ["O1", "O2", "O3"],          & 
                         "$y_2$", "$y_1$", "(b)", path(2) ) 
   call plot_parametrics( y(:,1), y(:,2:2), ["O2"],  "$y_2$", "$y_1$", & 
                          "(c)", path(3) ) 
   call plot_parametrics( y(:,1), y(:,3:3), ["O3"], "$y_2$", "$y_1$",  & 
                          "(d)", path(4) ) 
end subroutine    
    
!*****************************************************************************************
! Plot a contour graph and generates latex file in a specific location 
!*****************************************************************************************
subroutine myexampleD 
 
    integer, parameter :: N=20, Nl = 29 
    real :: x(0:N), y(0:N), z(0:N, 0:N)
    real, parameter :: PI = 4 * atan(1d0)
    real :: levels(0:Nl), a = 0, b = 2 * PI  
    integer :: i 
    character(len=100) :: path(2) =  ["./results/myexampleDa", & 
                                      "./results/myexampleDb"  ] 
    x  = [ (a + (b-a)*i/N, i=0, N) ] 
    y  = [ (a + (b-a)*i/N, i=0, N) ]
    a = -1; b = 1 
    levels  = [ (a + (b-a)*i/Nl, i=0, Nl) ]
    z = Tensor_product( sin(x), sin(y) ) 
   
    call plot_contour(x, y, z, "x", "y", levels, "(a)",path(1),"color") 
    call plot_contour(x, y, z, "x", "y", levels, "(b)",path(2),"isolines") 
end subroutine       



!*****************************************************************************************
! Extraction of the matrix of a linear system of equations 
!*****************************************************************************************
subroutine myexampleE 

    integer :: N
    integer :: i 
    real, allocatable :: A(:,:)   
    
    N = 4 
    A = Linear_operatorF(Laplace1D, N)
    do i=1, N+1
        write(*,'(100f5.1)') A(i,:) 
    end do 
    

contains 

function Laplace1D(U) result(F) 
  real, intent(in) :: U(:) 
  real :: F(size(U)) 
  
  integer :: N, i 
  
  N = size(U) 
 
  do i=2, N-1  
      F(i) = U(i+1) - 2*U(i) - U(i-1) - s1D(i) 
  end do  
  
  F(1) = U(2) - 2*U(1) 
  F(N) =  - 2*U(N) - U(N-1) 
  
  
end function




function s1D(i) 
  integer, intent(in) :: i 
  real :: s1D 
  
   s1D = i**2 
   
end function 

end subroutine

function s2D(i,j) 
  integer, intent(in) :: i, j  
  real :: s2D 
  
   s2D = i**2 + j**2
   
end function 







subroutine Legendre_points_example
 
 
 real, allocatable :: x(:), y(:), alpha(:) 
 integer :: p, i 
 real, parameter :: PI = 4 * atan(1.) 
 
 
 
       
   p = 6
   allocate( x(0:p), alpha(0:p), y(0:p) ) 
   
 !  call Legendre_Gauss_Lobatto_points(p, x, alpha) 
   x(0) = 0 
   x(p) = 2
   call Legendre_Grid_Initialization("x", x, alpha )
      
   do i=0, p  
       write(*,*) " i, xi, wi = ", i+1, x(i), alpha(i)  
   end do 
   
  y = sin( PI*x ) 
  
  call plot(x, y, "alpha" ) 
   
   
 
end subroutine








!******************************************************** 
!  It compares original codes with wrappers
!********************************************************
subroutine Embedded_Runge_kutta_Arenstorf_orbit 

    real :: t0 = 0
    real :: tf = 5*17.0652165601579625  ! Time domain for Arenstorf orbit 
    integer, parameter :: N = 100000    ! Time steps of the first grid 
    integer, parameter :: Nv = 4        ! Number of variables 
    real :: Time(0:N), U0(Nv)           ! Time domain and initial conditions
    integer, parameter :: Np = 1        ! number of graphs
    real :: U(0:N, Nv, Np)              ! Solutions
                                        ! first index: grid, second index: scheme 
    real :: x(0:N), y(0:N, 1) 

    character (len=40) :: names(Np) =    [ "$\epsilon=10^{-5}$"]
    real :: tolerances(Np) = [ 1d-13 ]  
    integer :: i 
   
    
     U(0,:,1) = [0.994, 0., 0., -2.0015851063790825 ]
     Time = [ (t0 + (tf -t0 ) * i / real(N), i=0, N ) ] 
            
         
     call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, 1) )
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" )
          
     call set_solver(family_name="eRK", scheme_name="Fehlberg87")
     call set_tolerance(tolerances(1))
     call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, 1) )
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" )
     
     call set_solver(family_name="GBS")
     call set_tolerance(tolerances(1))
     call Cauchy_ProblemS( Time, Arenstorf_equations, U(:, :, 1) )
     call plot_parametrics( U(:, 1, :), U(:, 2, :), names, "$x$", "$y$" )
     
    
end subroutine


end module  

