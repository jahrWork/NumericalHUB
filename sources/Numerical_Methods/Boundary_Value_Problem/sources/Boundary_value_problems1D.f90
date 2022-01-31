!******************************************************************************
!
! Author: Juan A Hernandez (juanantonio.hernandez@upm.es) & Javier Escoto 
!******************************************************************************
module Boundary_value_problems1D 

use Linear_systems
use Non_Linear_Systems 
use Collocation_methods 

use Dependencies1D
use Linearity
use Stability
 
implicit none  

private
public :: Boundary_Value_Problem1D
public :: Boundary_Value_Problem1D_system

public :: Linear_Boundary_Value_Problem1D
public :: Linear_Boundary_Value_Problem1D_system

public :: Non_Linear_Boundary_Value_Problem1D
public :: Non_Linear_Boundary_Value_Problem1D_system


abstract interface  

       real function DifferentialOperator1D(x, u, ux, uxx) 
                        real, intent(in) :: x, u, ux, uxx 
       end function  

       real function      BC1D(x, u, ux) 
           real, intent(in) :: x, u, ux 
       end function  
       
       function DifferentialOperator1D_system(x, u, ux, uxx) 
                        real, intent(in) :: x, u(:), ux(:), uxx(:) 
                        real :: DifferentialOperator1D_system(size(u)) 
       end function  

       function BC1D_system(x, u, ux) 
           real, intent(in) :: x, u(:), ux(:) 
           real :: BC1D_system(size(u)) 
       end function  

       subroutine  NonLinearSolver(Function, x)
          use Jacobian_module
          procedure(FunctionRN_RN) :: Function
          real, intent(inout) :: x(:)
       end subroutine  


 end interface

  
contains 
    
!**************************************************************************************************************
! Boundary value problem 
!
!   Differential_operator(x, u, ux, uxx) = 0. Linear and non linear 
!   Boundary_conditions(x, u, ux) at x=0 and x=L 
!
!**************************************************************************************************************
subroutine Boundary_Value_Problem1D(                   & 
                x_nodes, Differential_operator,        & 
                Boundary_conditions, Solution, Solver )
       
     real, intent(in) :: x_nodes(0:)
     procedure (DifferentialOperator1D) :: Differential_operator
     procedure (BC1D) ::  Boundary_conditions
     real, intent(inout) :: Solution(0:) 
     procedure (NonLinearSolver), optional:: Solver
     
     
     logical :: linear1D 
     linear1D = Linearity_BVP1D( Differential_operator )
 
     
     if (linear1D) then 
         
       call Linear_Boundary_Value_Problem1D( x_nodes,                   &
            Differential_operator, Boundary_conditions, Solution)
     else 
         
       call Non_Linear_Boundary_Value_Problem1D( x_nodes,               &
            Differential_operator, Boundary_conditions, Solution, Solver)
     end if 

end subroutine 


!**************************************************************************************************************
! Linear Boundary Value Problem 1D
!
!       Differential_operator(x, u, ux, uxx)
!       Boundary_conditions(x, u, ux) at x=0 and x=L 
!
!**************************************************************************************************************
subroutine Linear_Boundary_Value_Problem1D( x_nodes,                     &
                                            Differential_operator,       &  
                                            Boundary_conditions, Solution)
     real, intent(in) :: x_nodes(0:)
     procedure (DifferentialOperator1D) :: Differential_operator
     procedure (BC1D) ::  Boundary_conditions
     real, intent(out) :: Solution(0:)
     
!  *** auxiliary variables
       integer ::  i, Nx 
       real, allocatable ::  b(:), U(:), A(:,:)
       logical :: dU(2) ! matrix of dependencies( order )     
   
!  *** Integration domain 
       Nx = size(x_nodes) - 1 
       allocate( b(0:Nx), U(0:Nx), A(0:Nx, 0:Nx) )
       dU =  BVP1D_dependencies( Differential_operator ) 

          
!  *** independent term  F = A U - b  ( U = inverse(A) b )         
       U = 0 
       b = -BVP_discretization(U) 
              
!  *** Kronecker delta to calculate the difference operator        
       do i=0, Nx
          U = 0 
          U(i) = 1 
          A(0:Nx, i) = BVP_discretization(U) + b
       enddo 
  
!  *** solve the linear system of equations  
       Solution = Gauss(A, b) 
contains 
!-----------------------------------------------------------------
function BVP_discretization(U) result(F)
           real, intent(in) ::  U(0:)
           real :: F(0:size(U)-1)  

            call FD_Equations( dU,  Nx, x_nodes, U, F,  & 
                          Differential_operator, Boundary_conditions  ) 
end function 

end subroutine 

!**************************************************************************************************************
! Non linear BVP 
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Non_Linear_Boundary_Value_Problem1D( x_nodes,             &
                        Differential_operator,  Boundary_conditions, & 
                        Solution, Solver)

     real, intent(in) :: x_nodes(0:)
     procedure (DifferentialOperator1D) :: Differential_operator
     procedure (BC1D) ::  Boundary_conditions
     real, intent(inout) :: Solution(0:) 
     procedure (NonLinearSolver), optional:: Solver
     
!  *** Nummber of grid points
       integer ::  Nx  
       logical :: dU(2) ! matrix of dependencies( order )     
       dU = BVP1D_dependencies( Differential_operator )
       Nx = size(x_nodes) - 1
  
!  *** Non linear solver
       if (present(Solver))  then
                   call Solver(BVP_discretization, Solution)
       else
                   call Newton(BVP_discretization, Solution)
       end if
contains

function BVP_discretization( U) result(F)
  real, intent (in) :: U(0:)
  real :: F(0:size(U)-1) 

  call FD_Equations( dU,  Nx, x_nodes, U, F,  & 
                          Differential_operator, Boundary_conditions  ) 
end function

end subroutine

    

subroutine FD_Equations( dU,  Nx, x, W, F,  & 
                          Differential_operator, Boundary_conditions  ) 
   logical :: dU(2) 
   integer :: Nx 
   real :: x(0:Nx), W(0:Nx), F(0:Nx)  
   procedure (DifferentialOperator1D) :: Differential_operator
   procedure (BC1D) ::  Boundary_conditions
   
     real :: Wx(0:Nx), Wxx(0:Nx) 
     integer :: i 
    
     if (dU(1)) call Derivative( "x", 1, W, Wx  )
     if (dU(2)) call Derivative( "x", 2, W, Wxx )
     call Derivative( "x", 1, W, Wx, 0  ) ! Derivative always at x=0  
     call Derivative( "x", 1, W, Wx, Nx ) ! Derivative always at x=Nx 
 
     F(0)  = Boundary_points( x(0),   W(0),   Wx(0), Wxx(0)  ) 
     F(Nx) = Boundary_points( x(Nx),  W(Nx),  Wx(Nx), Wxx(Nx) ) 
     do i=1, Nx-1
           F(i) = Differential_operator( x(i), W(i), Wx(i), Wxx(i) ) 
     enddo 
contains 
                           
function  Boundary_points(x, W, Wx, Wxx) result(F) 
             real, intent(in) :: x, W, Wx, Wxx  
             real :: F
    real :: C, D  
        
       C = Boundary_conditions( x, W, Wx ) 
       D = Differential_operator( x, W,  Wx, Wxx )   
     
      if (C == FREE_BOUNDARY_CONDITION) then 
                        F = D 
      else           
                        F = C 
      end if

end function  

end subroutine





















 
!**************************************************************************************************************
! Boundary value problem 1D system 
!
!       Differential_operator(x, u, ux, uxx)
!       Boundary_conditions(x, u, ux)  at boundary [a, b] 
!
!**************************************************************************************************************
subroutine Boundary_Value_Problem1D_system(x_nodes,  Differential_operator,   & 
                                           Boundary_conditions, Solution, Solver)
                                        
     real, intent(in) :: x_nodes(0:)
     procedure (DifferentialOperator1D_system) :: Differential_operator
     procedure (BC1D_system) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, :) 
     procedure (NonLinearSolver), optional:: Solver
    
     integer :: Nv
     logical :: linear1D 
     Nv = size(Solution, dim=2)
        
     linear1D = Linearity_BVP1D_system( Differential_operator, Nv )  
    
     if (linear1D) then 
         
              call Linear_Boundary_Value_Problem1D_system(     x_nodes, Differential_operator, & 
                                                               Boundary_conditions, Solution)
     else 
              call Non_Linear_Boundary_Value_Problem1D_system( x_nodes, Differential_operator, &
                                                               Boundary_conditions, Solution, Solver )         
     end if 
   
  
end subroutine   



!**************************************************************************************************************
!
!       Differential_operator(x, u, ux, uxx)
!       Boundary_conditions(x, u, ux)  at boundary [a, b]
!
!**************************************************************************************************************
subroutine Linear_Boundary_Value_Problem1D_system(x_nodes, Differential_operator,  & 
                                                  Boundary_conditions, Solution)
                                        
     real, intent(in) :: x_nodes(0:)
     procedure (DifferentialOperator1D_system) :: Differential_operator
     procedure (BC1D_system) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, :) 
     

!  *** variables specification
       integer ::   i, Nx, Nv, Ns  
       real, allocatable :: Difference_operator(:,:), V(:), b(:)
       logical ::  dU(size(Solution,dim=2), 2)  ! matrix of dependencies( variable, derivative )
       
!  ***  Integration domain 
        Nx = size(x_nodes) - 1 
        Nv = size(Solution, dim=2)  
        Ns = (Nx+1)*Nv 
       
        dU = BVP1D_system_dependencies( Differential_operator, Nv )

        allocate( Difference_operator(Ns, Ns), b(Ns), V(Ns) ) 
        V = 0 
      
!  ***  independent term  A U = b  ( U = inverse(A) b )    
        b = -Equations(V, 0.0)
       
!  *** Delta kronecker to calculate the difference operator  
       Difference_operator = System_matrix( Equations, V, 0.0)   
              
!  *** solve the linear system of equations  
       V = Gauss( Difference_operator, b )
   
       Solution = reshape( V, [ Nx+1, Nv ] ) 
     
 contains 
                                                  
function Equations(V, t) result(G)
    real :: V(:), t 
    real :: G(size(V)) 
    
    call FD_Equations_system( dU, Nv, Nx, x_nodes, V, G, & 
                        Differential_operator, Boundary_conditions )
    
end function      


 end subroutine 



!**************************************************************************************************************
!
!       Differential_operator_system(x, u, ux, uxx)
!       Boundary_conditions_system(x, u, ux)  at boundary [a, b]
!
!**************************************************************************************************************
subroutine Non_Linear_Boundary_Value_Problem1D_system( x_nodes, Differential_operator, &
                                                       Boundary_conditions, Solution, Solver )
                                                      
     real, intent(in) :: x_nodes(0:)
     procedure (DifferentialOperator1D_system) :: Differential_operator
     procedure (BC1D_system) ::  Boundary_conditions
     real, intent(inout) :: Solution(0:, :) 
     procedure (NonLinearSolver), optional:: Solver


!  *** variables specification
       integer ::  i,  k, ik, Nx, Nv, M
       real :: x0, xf
       real, allocatable ::  F(:), Us(:)
       logical, allocatable :: dU(:, :)  ! matrix of dependencies( variable, derivative )
 
!  *** Integration domain
       Nx = size(x_nodes) - 1
       Nv = size(Solution, dim=2)
       M = (Nx+1)*Nv  
       allocate(F(M), Us(M), dU(Nv,2))
       
       dU = .true. !check also BC ! BVP1D_system_dependencies( Differential_operator, Nv )


!  *** initial guess for iteration
       Us = reshape(Solution, [ M ])

!  *** Non linear solver
        if (present(Solver))  then
                   call Solver(System_BVP, Us)
        else
                   call Newton(System_BVP, Us)
        end if
        
        Solution = reshape( Us, [Nx+1, Nv] )

contains
  
function System_BVP(V) result(G)
    real, intent(in) :: V(:)
    real :: G(size(V)) 
    
    call FD_Equations_system( dU, Nv, Nx, x_nodes, V, G, & 
                       Differential_operator, Boundary_conditions )
    
end function      
 
end subroutine





subroutine  FD_Equations_system( du, Nv, Nx, x, W, F,  & 
                          Differential_operator, Boundary_conditions  ) 
   logical :: dU(Nv, 2) 
   integer :: Nv, Nx 
   real :: x(0:Nx), W(0:Nx, Nv), F(0:Nx, Nv)  
   procedure (DifferentialOperator1D_system) :: Differential_operator
   procedure (BC1D_system) ::  Boundary_conditions
   
     real :: Wx(0:Nx, Nv), Wxx(0:Nx, Nv) 
     integer :: i 
   
     do i=1, Nv 
       if (dU(i,1)) call Derivative( "x", 1, W(0:,i), Wx(0:,i)  )
       if (dU(i,2)) call Derivative( "x", 2, W(0:,i), Wxx(0:,i) )
       call Derivative( "x", 1, W(0:,i), Wx(0:,i), 0 )  ! Derivative always at x=0 (BCs) 
       call Derivative( "x", 1, W(0:,i), Wx(0:,i), Nx ) ! Derivative always at x=Nx(BCs) 
     end do 
 
     F(0,:)  = Boundary_points( x(0),   W(0,:),   Wx(0,:), Wxx(0,:)  ) 
     F(Nx,:) = Boundary_points( x(Nx),  W(Nx,:),  Wx(Nx,:), Wxx(Nx,:) ) 
     do i=1, Nx-1
           F(i, :) = Differential_operator( x(i), W(i,:), Wx(i,:), Wxx(i,:) ) 
     enddo 

contains 
                           
function  Boundary_points(x, W, Wx, Wxx) result(F) 
             real, intent(in) :: x, W(:), Wx(:), Wxx(:)  
             real :: F(size(W)) 
   
        real :: C(size(W)), D(size(W)) 
        integer :: l, Nv 
        Nv = size(W) 
        
       C = Boundary_conditions( x, W, Wx ) 
       D = Differential_operator( x, W,  Wx, Wxx )   
                                                                       
       do l=1, Nv 
         if (C(l) == FREE_BOUNDARY_CONDITION) then 
                        F(l) = D(l) 
         else           
                        F(l) = C(l) 
         end if
      end do 

end function  

end subroutine




    
end module
