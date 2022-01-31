module API_Example_Boundary_Value_Problem

    use Boundary_value_problems 
    !use Finite_differences
    use Collocation_methods 
    use Non_Linear_Systems
    use plots
    use Stability
    
    implicit none
    
           
interface  
    function FunctionRN_RN(x) result(F) 
      real, intent(in) :: x(:) 
      real :: F( size(x) ) 
     end function 
end interface  

contains  
  
    
subroutine BVP_examples 
   call Legendre_1D
   call Elastic_beam_1D
   
   call Poisson_2D
   call Elastic_Plate_2D
   call Elastic_Nonlinear_Plate_2D
end subroutine


subroutine Legendre_1D

    integer, parameter :: N = 40, q = 6  
    real :: x(0:N), U(0:N,1), Error(0:N,1) 
    real :: x0 = -1 , xf = 1
    integer :: i
    character(len=100) :: path(2) = [   & 
           "./doc/chapters/Boundary_Value_Problem/figures/Legendrea", &
           "./doc/chapters/Boundary_Value_Problem/figures/Legendreb"  ]

    write (*, *) 'Solution of boundary value problems ' 
    write (*, *) ' (1 - x**2) * yxx - 2 * x * yx + n * (n + 1) * y = 0 ' 
    write (*, *) ' y(-1) = +1, y(+1) = 1 '
      
    x(0) = x0; x(N) = xf  
    call Grid_Initialization( grid_spacing = "nonuniform", &
                              direction = "x",   q = q, nodes = x )
!   Legendre solution   
    call Boundary_Value_Problem( x_nodes = x,                         & 
                                 Differential_operator = Legendre,     & 
                                 Boundary_conditions   = Legendre_BCs, & 
                                 Solution = U(:,1) )
    
    Error(:,1) = U(:,1) - ( 231 * x**6 - 315 * x**4 + 105 * x**2 - 5 )/16.
    
    call plot_parametrics(x, U, ["P6"], "$x$", "$y$", "(a)", path(1))
    call plot_parametrics(x, Error, ["E"], "$x$", "$E$", "(b)", path(2))
contains 


!****** Differential operator *********
real function Legendre(x, y, yx, yxx) result(L)
   real, intent(in) :: x, y, yx, yxx   
   
    integer :: n = 6
        
    L = (1 - x**2) * yxx - 2 * x * yx + n * (n + 1) * y
       
end function 
    
!********* Boundary conditions *********
real function Legendre_BCs(x, y, yx) result(BCs)
           real, intent(in) :: x, y, yx            

        if (x==x0 .or. x==xf ) then
                           BCs = y - 1 
        else 
            write(*,*) " Error BCs x=", x; stop  
        endif            
                 
end function  

end subroutine 

subroutine Test_BVP1D

    integer, parameter :: N = 20, q = 6  
    real :: x(0:N), U(0:N,1), Error(0:N,1) 
    real :: x0 = -1 , xf = 1
    integer :: i
    
    x(0) = x0; x(N) = xf  
    call Grid_Initialization( grid_spacing = "nonuniform", &
                              direction = "x",   q = q, nodes = x )
!
    call Boundary_Value_Problem( x_nodes = x,                         & 
                                 Differential_operator = test,     & 
                                 Boundary_conditions   = test_BCs, & 
                                 Solution = U(:,1) )
    
    call plot_parametrics(x, U(:,1:1), ["y"], "$x$", "$y$" )
contains 


!****** Differential operator *********
real function test(x, y, yx, yxx) result(L)
   real, intent(in) :: x, y, yx, yxx   
   
       L =  yxx - 1
           
end function 
    
!********* Boundary conditions *********
real function test_BCs(x, y, yx) result(BCs)
           real, intent(in) :: x, y, yx            

        if (x==x0 ) then
                           BCs = yx 
        else  if (x==xf ) then
                           BCs = y - 1                   
        else 
            write(*,*) " Error BCs x=", x; stop  
        endif            
                 
end function  

end subroutine 




subroutine Elastic_beam_1D

    integer, parameter :: N = 40, q = 4, Nv = 2  
    real :: x(0:N), U(0:N, Nv)
    real :: x0 = -1 , xf = 1
    integer :: i
    character(len=100) :: path(2) = [   & 
           "./doc/chapters/Boundary_Value_Problem/figures/Beam1Da", &
           "./doc/chapters/Boundary_Value_Problem/figures/Beam1Db"  ]

    write (*, *) 'Solution of an elastic beam 1D' 
       
    x(0) = x0; x(N) = xf  
    call Grid_Initialization( grid_spacing = "nonuniform", &
                              direction = "x",   q = q, nodes = x )
    
!   Legendre solution  
    call Boundary_Value_Problem( x_nodes = x,                         & 
                                 Differential_operator = Beam_equations, & 
                                 Boundary_conditions   = Beam_BCs, & 
                                 Solution = U  ) !, Solver = Itera )
    
    

    call plot_parametrics(x, U(:,1:1), ["W"], "$x$", "$w$", "(a)", path(1))
    call plot_parametrics(x, U(:,2:2), ["M"], "$x$", "$E$", "(b)", path(2))
contains 


!****** Differential operator *********
function Beam_equations(x, U, Ux, Uxx) result(L)
   real, intent(in) :: x, U(:), Ux(:), Uxx(:)  
   real :: L(size(U)) 
   
    real :: W, Wxx, M, Mxx 
    W   = U(1);       M  = U(2) 
    Wxx = Uxx(1);    Mxx = Uxx(2) 
  
    L = [ Wxx - M,  Mxx - 100 * x ]
       
end function 

 
!********* Boundary conditions *********
function Beam_BCs(x, U, Ux) result(BCs)
           real, intent(in) :: x, U(:), Ux(:)
           real :: BCs(size(U)) 

   real :: W, Wx
   w = U(1); Wx = Ux(1) 
           
        if (x==x0) then
                           BCs = [ W, Wx ]  
        else if (x==xf) then
                           BCs = [ W, Wx ]   
        endif            
                 
end function  

end subroutine 
 



subroutine Poisson_2D

    integer, parameter :: Nx = 30, Ny = 30, q = 11 ! mesh grid and order
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny)        ! grid and solution 
    integer :: i, j 
    real :: a=0, b=1, Umax, Umin
    integer, parameter :: Nl = 19                  ! number of isolines  
    real ::  levels(0:Nl)  = 0                      ! value of isolines 
    character(len=100) :: path(2) = [   & 
           "./doc/chapters/Boundary_Value_Problem/figures/Poissona", &
           "./doc/chapters/Boundary_Value_Problem/figures/Poissonb"  ]

    write (*, *) '2D Poisson solver ' 
    x(0) = a; x(Nx) = b; y(0) = a; y(Nx) = b; 
    
    call Grid_Initialization( "nonuniform", "x",  x, q )
    call Grid_Initialization( "nonuniform", "y",  y, q )

!   Poisson equation     
    call Boundary_Value_Problem( x_nodes = x, y_nodes = y,               & 
                                 Differential_operator = Poisson,        & 
                                 Boundary_conditions = PBCs, Solution = U) 
    
    
    call plot_contour(x, y, U, "$x$", "$y$",              & 
                      levels, "(b)", path(2), "isolines" ) 
    do i=0, Nx; do j=0, Ny; 
        U(i,j) = source( x(i), y(j) ) 
    end do; end do 
    call plot_contour(x, y, U, "$x$", "$y$", levels, "(a)", path = path(1), graph_type = "isolines" ) 
    

contains
    

!********* Function *********
real function Poisson(x, y, u, ux, uy, uxx, uyy, uxy) result(L) 
  real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
        
         L =  uxx + uyy - source(x,y)
       
end function

!********* Boundary conditions *********
real function PBCs(x, y, u, ux, uy) result(BCs)
real, intent(in) :: x, y, u, ux, uy

        if ( x==a .or. x==b .or. y==a .or. y==b ) then
                                                 BCs = u
        else
            write(*,*) " Error BCs x=", x;stop 
        endif

end function

!************** source term *********
real function source(x, y)
    real, intent(in) :: x, y
        
        real :: r1, r2, a=100
        
         r1 = norm2( [x, y] - [ 0.2, 0.5] ) 
         r2 = norm2( [x, y] - [ 0.8, 0.5] ) 
        
         source   = a * exp(-a*r1**2) + a * exp(-a*r2**2) 
       
end function

end subroutine 


!*************************************
! Plates 2D 
!*************************************
subroutine Elastic_Plate_2D
    integer, parameter :: Nx = 20, Ny = 20, Nv = 2, q= 4  
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv), f(0:Nx, 0:Ny) 
    integer, parameter :: Nl = 19                  ! number of isolines  
    real ::  levels(0:Nl)  = 0                     ! value of isolines 
    real :: x0 = -1 , xf = 1 , y0 = -1  , yf = 1
    integer :: i, j 
    character(len=100) :: path(2) = [   & 
           "./doc/chapters/Boundary_Value_Problem/figures/Platea", &
           "./doc/chapters/Boundary_Value_Problem/figures/Plateb"  ]
    
    write (*, *) '2D Linear plate'
    
    x(0) = x0 ; x(Nx) =  xf; y(0) = y0; y(Ny) = yf;  U = 1  
    
!   Elastic linear plate     
    call Grid_Initialization( "nonuniform", "x", x, q )
    call Grid_Initialization( "nonuniform", "y", y, q ) 
    
    call Boundary_Value_Problem( x_nodes = x, y_nodes = y,             & 
                                 Differential_operator = Elastic_Plate,& 
                                 Boundary_conditions = Plate_BCs,      &
                                 Solution = U )
    
    do i=0, Nx; do j=0, Ny; f(i,j) = load(x(i), y(j)); enddo; enddo  
    call plot_contour( x, y, f, "$x$", "$y$", levels,                   &
                       "(a)", path(1), "isolines" )
    
    call plot_contour( x, y, U(:,:,1), "$x$", "$y$", levels,            & 
                       "(b)", path(2),"isolines")
    write(*,*) " w = ", maxval(U(:,:,1)), minval(U(:,:,1)) 
          
    
contains

function Elastic_Plate(x, y, u, ux, uy, uxx, uyy, uxy) result(L)
   real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
   real :: L(size(u)) 
   
   real :: w, wxx, wyy, v, vxx, vyy

        w = u(1); wxx = uxx(1); wyy = uyy(1) 
        v = u(2); vxx = uxx(2); vyy = uyy(2) 

        L(1) =  wxx + wyy - v   
        L(2) =  vxx + vyy - load(x,y) 
      
end function


        
!********* Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function Plate_BCs(x, y, u, ux, uy) result(BCs)
  real, intent(in) :: x, y, u(:), ux(:), uy(:)
  real :: BCs(size(u))
  
        real :: w, v 
        
        w = u(1) 
        v = u(2) 

        if (x==x0 .or. x==xf .or. y==y0 .or. y==yf ) then
            
            BCs(1) = w
            BCs(2) = v 
            
        else
            write(*,*) " Error BCs x=", x; stop 
        endif

end function

end subroutine 

real  function load(x, y) 
   real, intent(in) :: x, y
   
   load =  100*y 
      
end function



subroutine Elastic_Nonlinear_Plate_2D

    integer, parameter :: Nx = 20, Ny = 20, Nv = 4, q = 4 ! Order 
    real :: x(0:Nx), y(0:Ny), U(0:Nx, 0:Ny, Nv)
    real :: x0 = -1 , xf = 1 , y0 = -1 , yf = 1 
    real :: mu=100,  levels(20) = 0
    character(len=5):: legends(2)  = [ "(a)", "(b)" ] 
    character(len=100) :: path(2) = [   & 
           "./doc/chapters/Boundary_Value_Problem/figures/NLPlatea", &
           "./doc/chapters/Boundary_Value_Problem/figures/NLPlateb"  ]
    integer :: i 
        
     write (*, *) '2D Solution of a non linear plate' 
    
     x(0) = x0 ; x(Nx) =  xf; y(0) = y0; y(Ny) = yf 
     call Grid_Initialization( "nonuniform", "x", x, q )
     call Grid_Initialization( "nonuniform", "y", y, q )
   
!   Elastic nonlinear plate    
    call Boundary_Value_Problem( x_nodes = x, y_nodes = y,             & 
                                 Differential_operator = NL_Plate,     & 
                                 Boundary_conditions   = NL_Plate_BCs, & 
                                 Solution = U )
    
    call plot_contour( x, y, U(:,:,1), "$x$", "$y$", levels,    & 
                       legends(1), path(1), "isolines" )
    call plot_contour( x, y, U(:,:,3), "$x$", "$y$", levels,    & 
                       legends(2), path(2), "isolines" )
    write(*,*) " w = ", maxval(U(:,:,1)), minval(U(:,:,1)) 
        
contains


!****** Function *********
function NL_Plate(x, y, u, ux, uy, uxx, uyy, uxy) result(L)
  real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: L(size(u))
  
  real ::   w, wxx, wyy, wxy, v, vxx, vyy
  real :: phi, phixx, phiyy, phixy, F, Fxx, Fyy
  
  w   = u(1);   wxx = uxx(1);   wyy = uyy(1) ;   wxy = uxy(1)
  v   = u(2);   vxx = uxx(2);   vyy = uyy(2) 
  phi = u(3); phixx = uxx(3); phiyy = uyy(3) ; phixy = uxy(3)
  F   = u(4);   Fxx = uxx(4);   Fyy = uyy(4) 
  
  L(1) = wxx + wyy - v   
  L(2) = vxx + vyy - load(x,y)                              & 
         - mu * Lb(wxx, wyy, wxy, phixx, phiyy, phixy)                                 
  L(3) = phixx + phiyy - F
  L(4) = Fxx + Fyy  +  Lb(wxx, wyy, wxy, wxx, wyy, wxy)
  
end function



real function Lb( wxx, wyy, wxy, pxx, pyy, pxy) 
  real, intent(in) ::  wxx, wyy, wxy, pxx, pyy, pxy
 
  Lb = wxx * pyy + wyy * pxx - 2 * wxy * pxy                                 
  
end function

        
!**** Boundary conditions *********
!Simply supported plate Mij = 0 <=> nabla **2 U_z = 0
function NL_Plate_BCs(x, y, u, ux, uy) result(BCs)
  real, intent(in) :: x, y, u(:), ux(:), uy(:)
  real :: BCs(size(u))

        if (x==x0 .or. x==xf .or. y==y0 .or. y==yf) then
                            BCs = u
        else
            write(*,*) " Error BCs x, y=", x, y; stop 
        endif

    end function

end subroutine 



!***********************************************************************
! Assembling for a specific problem 
!***********************************************************************
subroutine BVP_FD

    integer, parameter :: Nx = 40 ! grid points 
    integer ::  Order = 6         ! finite differences order
    real :: x(0:Nx)               ! Grid distribution  
    real :: u(0:Nx)               ! Solution u(x)
    
!   Spatial domain             
    x(0) = -1; x(Nx) = +1  
    
 !  Grid points    
    call Grid_Initialization("nonuniform", "x",  x, Order)
    
!   Initial guess
    u = 1  
    
!   Newton solution    
    call Newton(Equations, u)
    
!   Graph 
    call qplot(x, u, Nx+1) 
contains 
 
! Difference equations 
 function Equations(u) result(F)
    real, intent (in) :: u(0:)
    real :: F(0:size(u)-1) 
    
    real :: uxx(0:Nx)
       
    call Derivative( "x", 2, u, uxx )
    
    F =  uxx + sin(6*u)   !  inner points
    
    F(0)  = u(0) - 1      ! B.C. at x = -1
    F(Nx) = u(Nx)         ! B.C. at x =  1
    
 end function  

end subroutine

subroutine BVP1D

    integer, parameter :: N = 40, q = 6  
    real :: x(0:N), U(0:N), x0 = -1 , xf = 1
   
    x(0) = x0; x(N) = xf  
    call Grid_Initialization( "nonuniform","x", x, q )
    call Boundary_Value_Problem( x, L, BCs, U )
    call qplot(x, U, N+1) 
contains 

real function L(x, u, ux, uxx) 
   real, intent(in) :: x, u, ux, uxx   
      
    L =  uxx  + sin(6*u) 
       
end function 

real function BCs(x, u, ux) 
  real, intent(in) :: x, u, ux            

        if (x==x0) then 
                    BCs = u - 1 
        else if (x==xf) then 
                    BCs = u 
        else 
            write(*,*) " Error BCs x=", x; stop  
        endif            
                 
end function  

end subroutine



!***********************************************************************
! Assembling for a specific problem 
!***********************************************************************
subroutine Matrix_BVP_FD

     
     
    integer, parameter :: Nx = 10 ! grid points 
    integer ::  Order = 6         ! finite differences order
    real :: x(0:Nx)               ! Grid distribution  
    real :: u(0:Nx)               ! Solution u(x)
    real :: A(0:Nx, 0:Nx) 
    real :: t 
    integer :: j 
    
!   Spatial domain             
    x(0) = -1; x(Nx) = +1  
    
 !  Grid points    
    call Grid_Initialization("nonuniform", "x",  x, Order)
    
!   Initial guess
    u = 0  
    
    A = System_matrix( Equations, u, t) 
          
    do j=0, Nx 
             write(*,'(20f8.3)') A(j, :) 
    end do 
 
contains 
 
! Difference equations 
function Equations(u, t) result(F)
    real :: u(0:), t 
    
    real :: F(0:size(u)-1) 
    
    real :: uxx(0:Nx)
       
    call Derivative( "x", 2, u, uxx )
    
    F =  uxx + sin(6*u)   !  inner points
    
    F(0)  = u(0)          ! B.C. at x = -1
    F(Nx) = u(Nx)         ! B.C. at x =  1
    
 end function  

end subroutine


end module 












