
!***********************************************************************
! Initial value boundary problem 1D.
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
!use Dependencies
use Dependencies_BC
use Temporal_error
use plots

implicit none   

private
public :: IBVP1DS,                    & 
          Spatial_discretization1DS,  & 
          Linear_operator1DS,         & 
          Spatial_Truncation_Error1DS


abstract interface  
  
       function DifferentialOperator1DS(x, t, u, ux, uxx) 
                        real, intent(in) :: x, t, u(:), ux(:), uxx(:) 
                        real ::  DifferentialOperator1DS( size(u) ) 
       end function  
   
       function BC1DS(x, t, u, ux) 
           real, intent(in) :: x, t, u(:), ux(:) 
           real :: BC1DS( size(u) )  
       end function  
       
end interface

! Boundaty point 1D with Nv variables per point
type BoundaryP1DS
     
        integer :: index                   ! grid index of the boundary condition 
        logical, allocatable :: impose(:)  ! impose at last Nv boudary contitions per point  
        integer, allocatable :: B_index(:) ! boundary index of variable to determine (Nb unknowns) 
                                           ! at last Nv unknowns per bounday point
                                           ! e.g. Nx = 40, Nv = 3, Nb = 5
                                           ! BC(1) % B_index = [ 1, 2, 3] 
                                           ! BC(2) % B_index = [ 4, 5]  
        
        logical, allocatable :: interfaces(:)  ! at last Nv interfaces per point 
        
end type


contains
                                 
!**********************************************************************
! Initial Boundary Value Problem 1D (vectorial case) 
!***********************************************************************
subroutine IBVP1DS( Time_Domain, x_nodes, Differential_operator,  & 
                    Boundary_conditions, Solution, Scheme            ) 

     real, intent(in) :: Time_Domain(:) 
     real, intent(inout) :: x_nodes(0:)
     procedure (DifferentialOperator1DS) :: Differential_operator 
     procedure (BC1DS) ::  Boundary_conditions
     real, intent(out) :: Solution(0:,0:, :) 
     procedure (Temporal_Scheme), optional :: Scheme
     
     real, pointer :: U_Cauchy(:, :) 
     integer :: Nt, Nx, Nv
     Nt=size(Time_Domain)-1; Nx=size(x_nodes)-1; Nv=size(Solution, dim=3) 
     call my_reshape( Solution, Nt+1, Nv*(Nx+1), U_Cauchy ) 
    
     call Cauchy_ProblemS( Time_Domain, F_Cauchy, U_Cauchy, Scheme )
      
contains
function  F_cauchy( U, t ) result(F) 
          real ::  U(:), t, F(size(U))   
  
   F = Space_discretization1D( Nx, Nv, Differential_operator,  & 
                               Boundary_conditions, x_nodes, U, t )
end function 
end subroutine

!*************************************************************************************************
! It gives the spatial discretization of differential operator together with boundary conditions
!*************************************************************************************************
function Space_discretization1D( Nx, Nv, Differential_operator, Boundary_conditions, x, U, t ) result(F) 
          integer, intent(in) :: Nx, Nv 
          procedure (DifferentialOperator1DS) :: Differential_operator 
          procedure (BC1DS) ::  Boundary_conditions
          real, intent(in) :: x(0:), t 
          real, target ::  U(:), F(size(U))   
  
  real, pointer :: Uv(:,:), Fv(:,:)  
  
  Uv(0:Nx, 1:Nv) => U 
  Fv(0:Nx, 1:Nv) => F 
  
  Fv = Spatial_discretization1DS(   & 
       Differential_operator, Boundary_conditions, x, Uv, t )
  !write(*,*) Fv 
  !read(*,*) 
  
end function 


function Spatial_discretization1DS(Differential_operator,      & 
                                   Boundary_conditions, x, U, t) result(F) 

     procedure (DifferentialOperator1DS) :: Differential_operator 
     procedure (BC1DS) ::  Boundary_conditions
     real, intent(in) ::  x(0:), t
     real :: U(0:, :), F(0:size(U)-1, size(U,dim=2))
        
    integer :: i, k, l, ib 
    integer :: Nx, Nv
    real, allocatable :: Ux(:,:), Uxx(:,:), Fi(:), Ub(:)  
    integer, save ::  Nb 
    logical, allocatable, save:: dU(:,:),  dBC(:,:) 
    type (BoundaryP1DS), save :: BC(2)
      
    Nx = size(x)-1; Nv = size(U, dim=2) 
 
! ***  Allocate internal variables and identifies the boundary unknowns 
       !if (t==0) then 
         call Initialize
         call Boundary_unknowns 
       !end if 
       allocate (  Ux(0:Nx, Nv), Uxx(0:Nx, Nv), Fi(Nv), Ub(Nb) )
        
! ***  It solves boundary equations to yield values at boundaries     
       call Determine_boundary_points( Ub, dBC, BC, Nx, Nv,  x, U, t, & 
                                       Boundary_conditions ) 
       
! ***  It calculates only derivatives which are used 
       call Involved_Derivatives(  Nx, Nv, dU, U, Ux, Uxx )  
       do i=0, Nx
           Fi = Differential_operator(x(i), t, U(i,:), Ux(i,:), Uxx(i,:)) 
            
      ! ** Boundary point
           k = Boundary_index(Nx, i) 
           if (k>0) then 
            do l=1, Nv
             if (BC(k) % impose(l) .or. BC(k) % interfaces(l)) & 
                                Fi(l) = IMPOSE_ZERO
            end do
           end if 
      ! ** inner point  
           F(i, :) = Fi(:)
       enddo 
contains 
    

subroutine Initialize

       integer :: j 
                      
       if (allocated(dU)) then 
            deallocate(dU, dBC) 
        end if      
        allocate(dU(Nv,2), dBC(Nv,2)) 
      
  
       dU  = .true. !IBVP1D_Dependencies_system( Nv, Differential_operator )
       dBC = .true. !BC_IBVP1D_Dependencies_system( Nv, Boundary_conditions, x(0), x(Nx) ) 
       
             
       do j=1, 2 
           
             if (  allocated(BC(j) % B_index )) then 
                 deallocate( BC(j) %  B_index, BC(j) % impose, BC(j) % interfaces ) 
             end if 
             
             allocate( BC(j) %  B_index(Nv) ) 
             allocate( BC(j) %  impose(Nv) ) 
             allocate( BC(j) % interfaces(Nv) )  
       end do
           
     
end subroutine 

subroutine Boundary_unknowns
     
 real :: u0(Nv), ux0(Nv), uy0(Nv), eq(Nv)   
 integer :: i_BC(2), k, i, l   
      
  i_BC = [0, Nx ]      
  u0 = 1.; ux0 = 2.; uy0 = 3  
  Nb = 0 
  
  
!  Fourier  
 if (method == "Fourier") then 
       do k=1, size(BC) 
           BC(k) % impose = .false.     
       end do   
       
 ! FD or Chebyshev       
 else   
  
  do k=1, size(BC)
      
    i = i_BC(k) 
    
    eq = Boundary_conditions( x(i),  t, u0,  ux0 ) 
    BC(k) % index = i 
    do l=1, Nv 
      BC(k) % interfaces(l) = eq(l) == INTERFACE_CONDITION
      BC(k) % impose(l) = eq(l) /= FREE_BOUNDARY_CONDITION  .and.  & 
                          eq(l) /= PERIODIC_BOUNDARY_CONDITION .and. & 
                          eq(l) /= INTERFACE_CONDITION 
      if (BC(k) % impose(l) ) then 
          Nb = Nb + 1  
          BC(k) % B_index(l) = Nb  
      end if 
    end do 
  end do
 
end if   
    
end subroutine       
end function 





subroutine Determine_boundary_points( Ub, dBC, BC, Nx, Nv, x, U, t, Boundary_conditions ) 
        real :: Ub(:) 
        logical, intent(in) :: dBC(:,:) 
        type (BoundaryP1DS) :: BC(:)
        integer, intent(in) :: Nx, Nv 
        real, intent(in) :: x(0:Nx), t
        real, intent(inout) :: U(0:Nx, Nv)
        procedure (BC1DS) ::  Boundary_conditions
        
        integer :: i, k, l, m
        
!  *** Ínitial condition for Ub        
       do k=1, size(BC); do l=1, Nv
           i = BC(k) % index
           m = BC(k) % B_index(l)  
          if( BC(k)% impose(l) ) Ub(m) = U(i, l)
       end do; end do 
        
 !  *** Solve boundary equations  
        call Newton( Boundary_equations, Ub )
          
contains 

function Boundary_equations(Z) result(G) 
        real, intent (in) :: Z(:)
        real :: G(size(Z))
    
       real :: Ux(0:Nx, Nv)
       integer :: i, k, m, l  
       real :: Gc(Nv), eq  
      
  
  !** update unknowns      
  do k=1, size(BC); do l=1, Nv
      i = BC(k) % index
      m = BC(k) % B_index(l)   
     
      if( BC(k)% impose(l) ) U(i, l) = Z(m)
      
  end do; end do 
  
  do l=1, Nv     
  !   write(*,*) " x_label =", x_label 
     call Derivative( "x", 1, U(:, l), Ux(:, l) )
  end do 
  
  !** Newton equations 
  do k=1, size(BC) 
       i = BC(k) % index 
       Gc = Boundary_conditions( x(i),  t, U(i, :),  Ux(i, :) )  
       
       do l=1, Nv  
         m = BC(k) % B_index(l) 
         if ( BC(k)% impose(l) ) G(m) = Gc(l)
       end do 
  end do
      
end function
end subroutine 




!****************************************************************
! vector of dependencies(derivatives) , ux ,uy, uxx uyy, uxy  
!****************************************************************
subroutine Involved_Derivatives( Nx, Nv, dU, U, Ux,  Uxx )
  integer, intent(in) :: Nx, Nv  
  logical, intent(in) :: dU(Nv, 2) 
  real, intent(in) ::   U(0:Nx, Nv)
  real, intent(out) :: Ux(0:Nx, Nv), Uxx(0:Nx, Nv)  
        
     integer :: k 
  
     !write(*,*) " Enter Involved derivatives Uxx ="
     !read(*,*) 
     
!  *** inner grid points
        do k=1, Nv 
                    
           if ( dU(k,1) ) call Derivative( "x", 1, U(:,k), Ux(:,k)  )
           if ( dU(k,2) ) call Derivative( "x", 2, U(:,k), Uxx(:,k) )
       
        end do  

end subroutine 
                                           
                                            
 
function Boundary_index(N, i) result(k) 
    integer, intent(in) :: N, i 
    integer :: k
    
    
       if (i==0) then  
                            k = 1
       else if (i==N) then 
                            k = 2 
       else 
                            k = -1 
       end if 
                
end function                                   
       









!**********************************************************************************
! It determines the Local Truncation Spatial Error of the Solution 
! by means of Richardson extrapolation 
! INPUTS : 
!           Differential_operator 
!           Boundary_conditions  
!           x_nodes or collocation points 
!           U solution to be discretized 
!           Order of the interpolation 
! OUTPUTS:
!           R : local truncation error of the spatial discretization 
!
! Author: juanantonio.hernandez@um.es (May 2022) 
!********************************************************************************** 
function Spatial_Truncation_Error1DS( Nv, Differential_operator,         &
                                      Boundary_conditions,               &
                                      x, Order, Test_function ) result(R) 
   
     procedure (DifferentialOperator1DS) :: Differential_operator 
     procedure (BC1DS) ::  Boundary_conditions 
     real, intent(inout) ::  x(0:) 
     interface 
        function Test_function(Nv,x) result(U) 
           integer, intent(in) :: Nv 
           real, intent(in) :: x(0:) 
           real :: U(0:size(x)-1, Nv) 
        end function 
     end interface 
     integer, intent(in) :: Nv, Order 
     real :: R(0:size(x)-1, Nv)

   integer :: i, Nx  
   real, allocatable :: x1(:), x2(:), U1(:,:), F1(:,:), U2(:,:), F2(:,:)   
   real :: t
   
       t = 0     
       Nx = size(x)-1;  
       allocate ( x1(0:Nx), x2(0:2*Nx), U1(0:Nx, Nv),F1(0:Nx, Nv), U2(0:2*Nx, Nv), F2(0:2*Nx, Nv) ) 
        
       x1 = x
       call Grid_Initialization( "nonuniform", "x", x1, Order )
       U1 = Test_function(Nv, x1) 
       F1 = Spatial_discretization1DS( Differential_operator, Boundary_conditions, x1, U1, t)  
         
       do i=0, Nx-1 
           x2(2*i)   = x1(i) 
           x2(2*i+1) = ( x1(i) + x1(i+1) )/2
       end do 
       x2(2*Nx) = x1(Nx)
       call Grid_Initialization( "unmodified", "x", x2, Order ) 
       
       U2 = Test_function(Nv, x2) 
       F2 = Spatial_discretization1DS( Differential_operator, Boundary_conditions, x2, U2, t)
       
       do i=0, Nx
            R(i,:) = ( F2(2*i,:) - F1(i,:) )/( 1 - 1./2**Order )  
       end do  
             
end function 


!*******************************************************************
! Given a vector function F: RN -> RN. 
! If the F (differential operator) is linear (F = A U + b), 
! it gives the system matrix A
!*******************************************************************
function Linear_operator1DS( Nv, x, Order, Differential_operator,  & 
                             Boundary_conditions ) result(A)

     integer, intent(in) :: Nv, Order
     real, intent(inout) :: x(0:) 
     procedure (DifferentialOperator1DS) :: Differential_operator 
     procedure (BC1DS) ::  Boundary_conditions
     real :: A(size(x)*Nv, size(x)*Nv) 
     
          
     real ::  U( size(x)*Nv ), b( size(x)*Nv ), t  
     integer :: j 
     
     call Grid_Initialization( "unmodified", "x", x, Order )
     
     t = 0 
     U = 0
     b = Space_discretization1D( size(x)-1, Nv, Differential_operator, Boundary_conditions, x, U, t ) 
     
     do j=1, size(x)*Nv 
 
         U = 0 
         U(j) = 1   
         A(:, j) = Space_discretization1D( size(x)-1, Nv, Differential_operator, Boundary_conditions, x, U, t ) - b
         
      end do 
      
             
end function  
  

                                 
end module 
    
    
     