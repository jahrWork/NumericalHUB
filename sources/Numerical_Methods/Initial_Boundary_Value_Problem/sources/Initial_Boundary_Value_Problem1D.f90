
!***********************************************************************
! It integrates in time the Initial value boundary problem 1D.
! Given the differential operator and the boundary equations,
! the discrete solution is calculated.
! author: Juan A Hernandez, juanantonio.hernandez@upm.es 
!***********************************************************************
module Initial_Boundary_Value_Problem1D

use Cauchy_Problem
use Temporal_scheme_interface 
use Collocation_methods 
use Non_Linear_Systems
use Utilities
use Dependencies
use Temporal_error

implicit none   

private

public :: IBVP1D, Spatial_discretization1D

public :: Spatial_Error_IBVP1D
public :: Truncation_Spatial_Error_IBVP1D  
public :: Temporal_Error_IBVP1D 

public :: IBVP1D_system, Spatial_discretization1D_system 



abstract interface  

       real function DifferentialOperator1D(x, t, u, ux, uxx) 
                        real, intent(in) :: x, t, u, ux, uxx 
       end function  

       real function BC1D(x, t, u, ux) 
           real, intent(in) :: x, t, u, ux 
       end function  
       
       function DifferentialOperator1D_system(x, t, u, ux, uxx) 
                        real, intent(in) :: x, t, u(:), ux(:), uxx(:) 
                        real ::  DifferentialOperator1D_system( size(u) ) 
       end function  
       
       function IC1D(x) result(U)
           real, intent(in) :: x(:)
           real :: U(size(x))
       end function 
    
       
       
       function BC1D_system(x, t, u, ux) 
           real, intent(in) :: x, t, u(:), ux(:) 
           real :: BC1D_system( size(u) )  
       end function  


       function DifferentialOperatorODE(U, t) result(F) 
                        real, intent(in) :: U(:), t 
                        real :: F(size(U)) 
       end function 

end interface

type Boundary1D_system 
     
        integer :: i 
        real, allocatable :: value(:) 
        real, allocatable :: equation(:)  ! Nv equations at the boundary point 
        logical, allocatable :: impose(:) ! at least Nv boudary contitions to impose per point  
        
end type


contains
!**********************************************************************************************************************
! Initial Boundary Value Problem 1D (scalar case) 
!***********************************************************************************************************************
subroutine IBVP1D( Time_Domain, x_nodes,  Differential_operator,  & 
                    Boundary_conditions, Solution, Scheme) 

     real, intent(in) :: Time_Domain(:),  x_nodes(0:)
     procedure (DifferentialOperator1D) :: Differential_operator 
     procedure (BC1D) ::  Boundary_conditions
     real, intent(out) :: Solution(0:,0:)
     procedure (Temporal_Scheme), optional :: Scheme
     
     
     integer :: Nx
     Nx = size(x_nodes) - 1    
          
     call Cauchy_ProblemS( Time_Domain, Space_discretization, & 
                           Solution, Scheme                  )
      
contains

function Space_discretization( U, t ) result(F) 
                      real ::  U(0:), t, F(0:size(U)-1) 
           
  F = Spatial_discretization1D(   & 
      Differential_operator, Boundary_conditions, x_nodes, U, t )
         
      
end function     

end subroutine


function Spatial_discretization1D( Differential_operator, Boundary_conditions, & 
                                   x_nodes, U, t) result(F) 

     procedure (DifferentialOperator1D) :: Differential_operator 
     procedure (BC1D) ::  Boundary_conditions 
     real, intent(in) ::  x_nodes(0:), U(0:), t
     real :: F(0:size(U)-1)
       
     logical, save :: dU(2) ! matrix of dependencies( order )
     integer :: Nx
     real, allocatable :: Ux(:), Uxx(:)
     integer :: k, N1, N2 
     logical :: BC(0:1) 
     
     Nx = size(x_nodes) - 1
     allocate( Ux(0:Nx), Uxx(0:Nx) )
 
                 
!  ***  Check if Differential operator depends on Ux and Uxx        
        if (t==0) dU = IBVP1D_Dependencies( Differential_operator ) 
        
!  ***  It solves one or two equations at  boundaries 
        call Boundary_points( x_nodes, U, t, BC, Boundary_conditions) 
              
!  *** inner grid points 
        if (dU(1)) call Derivative( "x", 1, U, Ux)
        if (dU(2)) call Derivative( "x", 2, U, Uxx)  
        
        do k = 0, Nx
         if ( (k==0.and.BC(0)) .or. (k==Nx.and.BC(1)) ) then 
           F(k) = 0  
         else 
           F(k) = Differential_operator(x_nodes(k), t, U(k), Ux(k), Uxx(k))
         end if 
        enddo  

end function 

!----------------------------------------------------------------------
subroutine Boundary_points( x, U, t, BC, Boundary_conditions )
         real, target ::   x(0:), U(0:), t
         logical, intent(out) :: BC(0:1)  
         procedure (BC1D) ::  Boundary_conditions
    
      integer :: N
      real :: U1(1), Uc(2)
      
          
        N = size(x)-1
        if (method == "Fourier") then 
            BC = .false. 
      
        else if (Boundary_conditions( x(0),  t, 1.,  2.) ==       & 
                                 PERIODIC_BOUNDARY_CONDITION) then             
            U(0)  = U(N)
            BC = [ .true., .false. ]
            
            
       else if (Boundary_conditions( x(N),  t, 1.,  2.) ==  & 
                                     FREE_BOUNDARY_CONDITION) then 
           
            U1 = U(0) 
            call Newton( BCs1, U1 )
            BC = [ .true., .false. ]
            
       else 
            Uc = U(0:N:N)  
            call Newton( BCs2, Uc )
            BC =  .true.
            
       end if 
       
contains 

!-----------------------------------------------------------------------
function BCs2(Y) result(G) 
        real, intent (in) :: Y(:)
        real :: G(size(Y))

    real :: Wx(0:N)
   
    U(0) = Y(1) 
    U(N) = Y(2)
    
    call Derivative( "x", 1, U, Wx)
   
    G(1)= Boundary_conditions( x(0),  t, U(0),  Wx(0) )
    G(2)= Boundary_conditions( x(N), t, U(N), Wx(N) )
  
end function 

!-----------------------------------------------------------------------
function BCs1(Y) result(G) 
        real, intent (in) :: Y(:)
        real :: G(size(Y))

    real :: Wx(0:N)
  
    U(0) = Y(1) 
    call Derivative( "x", 1, U, Wx)
    G(1)= Boundary_conditions( x(0),  t, U(0),  Wx(0))
    
end function 

end subroutine



subroutine Temporal_Error_IBVP1D( Differential_operator, Boundary_conditions, Initial_conditions,   & 
                         Scheme, time, x_nodes, Order, U, Error)

     procedure (DifferentialOperator1D) :: Differential_operator
     procedure (BC1D) ::  Boundary_conditions
     procedure (IC1D) ::  Initial_conditions
     procedure (Temporal_Scheme), optional :: Scheme
     integer, intent(in) :: Order
     real, intent(in) ::  time(0:)   
     real, intent(inout) :: x_nodes(0:), U(0:,0:) ! 0:Nt, 0:Nx 
     real, intent(out) :: Error(0:, 0:)
     
    
    call Grid_Initialization( "nonuniform", "x",  x_nodes, Order)
    U(0,:) = Initial_conditions(x_nodes) 
    
    
    call Error_Cauchy_Problem( Time, F_Cauchy, Scheme, 2, U, Error ) 
         
         
contains 
function F_Cauchy( U, t ) result(F) 
                      real ::  U(0:), t, F(0:size(U)-1) 
           
  F = Spatial_discretization1D( & 
      Differential_operator, Boundary_conditions, x_nodes, U, t )
         
      
end function
                         
 
end subroutine  


   

!**********************************************************************************
! It determines the Local Truncation Spatial Error of the Solution 
! by means of Richardson extrapolation 
! INPUTS : 
!           Differential_operator 
!           Boundary_conditions  
!           Initial_conditions 
!           x_nodes or collocation points 
!           Order of the interpolation 
! OUTPUTS:
!           R : local truncation error of the spatial discretization 
!
! Author: juanantonio.hernandez@um.es (May 2021) 
!********************************************************************************** 
subroutine Truncation_Spatial_Error_IBVP1D( Differential_operator, Boundary_conditions,   &
                                            Initial_conditions, x_nodes, Order, R ) 
   
     procedure (DifferentialOperator1D) :: Differential_operator 
     procedure (BC1D) ::  Boundary_conditions 
     procedure (IC1D) ::  Initial_conditions
     integer, intent(in) :: Order 
     real, intent(inout) ::  x_nodes(0:)
     real, intent(out) :: R(0:size(x_nodes)-1)

   integer :: i, N 
   real, allocatable :: x1(:), U1(:), R1(:), x2(:), U2(:), R2(:)   
   real :: t = 0 
       
       N = size(x_nodes)-1;  
       allocate ( x1(0:N),  U1(0:N), R1(0:N), x2(0:2*N), U2(0:2*N), R2(0:2*N) ) 
       
       
       call Grid_Initialization( "nonuniform", "x", x_nodes, Order )
       
       x1 = x_nodes
       U1 = Initial_conditions(x1)
       R1 = Spatial_discretization1D( Differential_operator, Boundary_conditions, x1, U1, t)  
       
       do i=0, N-1 
           x2(2*i)   = x1(i) 
           x2(2*i+1) = ( x1(i) + x1(i+1) )/2
       end do 
       x2(2*N) = x1(N)
       
       call Grid_Initialization( "unmodified", "x", x2, Order )
       U2 = Initial_conditions(x2)
       R2 = Spatial_discretization1D( Differential_operator, Boundary_conditions, x2, U2, t)
       
       do i=0, N 
            R(i) = ( R2(2*i)- R1(i) )/( 1 - 1./2**Order ) 
       end do  
             
end subroutine 

!**********************************************************************************
! It determines the Spatial Error of the Solution by means of Richardson extrapolation 
! INPUTS : 
!           Differential_operator 
!           Boundary_conditions  
!           Initial_conditions 
!           time domain    
!           x_nodes or collocation points 
!           Order of the interpolation 
!               
! OUTPUTS:
!           U : solution ( 0:Nt, 0:Nx) ( time index, space index) 
!           Error of the spatial discretization   
!               
! Note: the default temporal scheme is a Runge-Kutta of 4th order. 
!       it is supposed that temporal error is much smaller than spatial error 
!       to asure that time step is chosen 1/10 of the original time step 
! Author: juanantonio.hernandez@um.es (May 2021) 
!**********************************************************************************               
subroutine Spatial_Error_IBVP1D( Differential_operator, Boundary_conditions,        & 
                                 Initial_conditions, time, x_nodes, Order, U, Error ) 
   
     procedure (DifferentialOperator1D) :: Differential_operator 
     procedure (BC1D) ::  Boundary_conditions 
     procedure (IC1D) ::  Initial_conditions
     integer, intent(in) :: Order 
     real, intent(in) ::  time(0:)
     real, intent(inout) ::  x_nodes(0:)
     real, intent(out) ::  U(0:, 0:), Error(0:, 0:)  ! 0:Nt, 0:Nx 

   integer :: i, j, Nt, Nx 
   real, allocatable :: x1(:), U1(:,:), x2(:), U2(:, :), t(:) 
   real :: t0, tf, dt 
       
       Nx = size(x_nodes)-1;   Nt = size(time)-1;
       allocate ( x1(0:Nx),  U1(0:10*Nt, 0:Nx), x2(0:2*Nx), U2(0:10*Nt, 0:2*Nx) ) 
       allocate ( t(0:10*Nt) )
       
       t0 = Time(0); tf = Time(Nt); dt = (tf-t0)/(10*Nt) 
       t = [ (t0 + dt*i, i=0, 10*Nt ) ]
       
       call Grid_Initialization( "nonuniform", "x", x_nodes, Order )
       
       x1 = x_nodes
       U1(0, :) = Initial_conditions(x1)
       call IBVP1D( t, x1, Differential_operator, Boundary_conditions, U1 ) 
       
       do i=0, Nx-1 
           x2(2*i)   = x1(i) 
           x2(2*i+1) = ( x1(i) + x1(i+1) )/2
       end do 
       x2(2*Nx) = x1(Nx)
       
       call Grid_Initialization( "unmodified", "x", x2, Order )
       U2(0, :) = Initial_conditions(x2)
       call IBVP1D( t, x2, Differential_operator, Boundary_conditions, U2 ) 
            
       do i =0, Nt
         do j=0, Nx 
            Error(i, j) = ( U2(10*i, 2*j)- U1(10*i, j) )/( 1 - 1./2**Order )  
            U(i, j) = U1(10*i, j)
         end do 
       end do   
      
       
end subroutine 

                                 
                     
                                 

                                 
                                 
  
                                 
!**********************************************************************************************************************
! Initial Boundary Value Problem 1D (vectorial case) 
!***********************************************************************************************************************
subroutine IBVP1D_system( Time_Domain, x_nodes,  Differential_operator,  Boundary_conditions, Solution, Scheme ) 

     real, intent(in) :: Time_Domain(:) 
     real, intent(inout) :: x_nodes(0:)
     procedure (DifferentialOperator1D_system) :: Differential_operator 
     procedure (BC1D_system) ::  Boundary_conditions
     real, intent(out) :: Solution(0:,0:, :) 
     procedure (Temporal_Scheme), optional :: Scheme
     
     real, pointer :: U_Cauchy(:, :) 
     integer :: Nt, Nx, Nv
     Nt = size(Time_Domain)-1;  Nx = size(x_nodes)-1; Nv = size(Solution, dim=3) 
    
    call my_reshape( Solution, Nt+1, Nv*(Nx+1), U_Cauchy ) 
    
    call Cauchy_ProblemS( Time_Domain, F_Cauchy, U_Cauchy, Scheme )
      
    
contains
function  F_cauchy( U, t ) result(F) 
          real ::  U(:), t, F(size(U))   
  
   F = Space_discretization( U, t )
           
end function 

function  Space_discretization( U, t ) result(F) 
          real, target ::  U(:), t, F(size(U))   
  
  real, pointer :: Uv(:,:), Fv(:,:)  
  
  Uv(0:Nx, 1:Nv) => U 
  Fv(0:Nx, 1:Nv) => F 
  
  Fv = Spatial_discretization1D_system(   & 
       Differential_operator, Boundary_conditions, x_nodes, Uv, t )
   
          
end function  
end subroutine

function Spatial_discretization1D_system( Differential_operator, Boundary_conditions, & 
                                          x,  U, t) result(F) 

     procedure (DifferentialOperator1D_system) :: Differential_operator 
     procedure (BC1D_system) ::  Boundary_conditions
     real, intent(in) ::  x(0:), t
     real ::   U(0:, :)
     real :: F(0:size(U)-1, size(U,dim=2))

        
    integer :: i,k, Nx, Nv
    real, allocatable :: Ux(:,:), Uxx(:,:), Uc(:), Fv(:)   
    logical, save, allocatable :: dU(:, :) 
    
  ! indexes of variables iU0):) ( at x0)  and iUN(:) ( at xN ) 
  ! where boundary conditions are imposed 
    integer,  allocatable :: iU0(:), iUN(:) 
    
  ! number of BCs at x0 and at xN   
    integer :: N0, N1
    
        Nx = size(x) - 1;  Nv = size(U, dim=2)
        allocate(  Ux(0:Nx, Nv), Uxx(0:Nx, Nv), Fv(Nv) ) 
        if (t==0) then 
                 if (allocated(dU)) deallocate(dU) 
                 allocate( dU(Nv, 2) ) 
                 dU = .true. 
        end if  
   
!  ***  solve a system of two equations to yield the values at the boundaries
        call Equations_at_boundary( x = x(0),  W = U(0,:),  Wx = Ux(0,:), iU = iU0 )
        call Equations_at_boundary( x = x(Nx), W = U(Nx,:), Wx = Ux(Nx,:), iU =  iUN )
        
        N0 = size(iU0); N1 = size(iUN) 
     
        allocate( Uc(N0+N1) ) 
        do k=1, N0;     Uc(k)    = U(0,  iU0(k) ); end do 
        do k=1, N1;     Uc(k+N0) = U(Nx, iUN(k) ); end do 
      
        call Newton( BCs, Uc )
   
!  *** inner grid points
        do i=1, Nv 
           if (dU(i,1)) call Derivative( "x", 1, U(0:,i), Ux(0:,i) )
           if (dU(i,2)) call Derivative( "x", 2, U(0:,i), Uxx(0:,i) )
        end do    
          
        do k=0, Nx
            Fv = Differential_operator( x(k), t, U(k,:), Ux(k,:), Uxx(k, :) )
            if (k==0) then 
                            do i=1, N0;  Fv(iU0(i)) = 0;  end do
             else if (k==Nx) then 
                             do i=1, N1;  Fv(iUN(i)) = 0;  end do
            else 
                F(k, :) = Fv 
            end if 
        enddo
  
contains 

subroutine Equations_at_boundary( x, W,  Wx, F, iU )
                real, intent(in) ::  x, W(:),  Wx(:)
                real, allocatable, optional, intent(out) ::  F(:)
                integer, allocatable, optional, intent(out) ::  iU(:) 

        integer :: i, k,  Nc           
        real :: F1( size(W) )   
         
        
        F1 = Boundary_conditions( x, t, W(:),  Wx(:) ) 
       
        Nc = count(F1 /= FREE_BOUNDARY_CONDITION ) 
    
        if (present(F) ) allocate( F(Nc) ) 
        if (present(iU)) allocate( iU(Nc) ) 
        
        k = 0 
        do i=1, Nv 
            if ( F1(i) /= FREE_BOUNDARY_CONDITION ) then 
              k = k + 1 
              if (present(F))   F(k) = F1(i) 
              if (present(iU))  iU(k) = i 
           end if
        end do  

end subroutine 
        

function BCs(Y) result(G) 
     real, intent (in) :: Y(:)
     real :: G(size(Y))
     
     real, allocatable :: G0(:), GN(:)   
     integer :: i,k    
     
     real :: Wx(0:Nx, Nv)

      
     do k=1, N0;       U(0, iU0(k) ) = Y(k);    end do
     do k=1, N1;       U(Nx, iUN(k) ) = Y(k+N0); end do 
     
     do i=1, Nv 
           call Derivative( "x", 1, U(0:, i), Wx(0:,i) )
     end do    
     
     call Equations_at_boundary( x(0),  U(0, :), Wx( 0,:), G0 )
     call Equations_at_boundary( x(Nx), U(Nx, :), Wx(Nx,:), GN )
      
     G = [ G0, GN ] 
    

end function 

end function

     
     
     
end module 
    
    