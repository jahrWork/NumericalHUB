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
public :: Boundary_Value_Problem1DS
public :: Linear_Boundary_Value_Problem1DS
public :: Non_Linear_Boundary_Value_Problem1DS


abstract interface  

       
       function DifferentialOperator1DS(x, u, ux, uxx) 
                        real, intent(in) :: x, u(:), ux(:), uxx(:) 
                        real :: DifferentialOperator1DS(size(u)) 
       end function  

       function BC1DS(x, u, ux) 
           real, intent(in) :: x, u(:), ux(:) 
           real :: BC1DS(size(u)) 
       end function  

       subroutine  NonLinearSolver(Function, x)
          use Jacobian_module
          procedure(FunctionRN_RN) :: Function
          real, intent(inout) :: x(:)
       end subroutine  


 end interface

  
contains 
 
 
!*********************************************************************************
! Boundary value problem 1D system 
!
!       Differential_operator(x, u, ux, uxx) = 0 at inner points x_i
!       Boundary_conditions(x, u, ux)  = 0  at boundary points
!
!**********************************************************************************
subroutine Boundary_Value_Problem1DS( x, Differential_operator,        & 
                                 Boundary_conditions, Solution, Solver)   
                                        
     real, intent(in) :: x(0:)
     procedure (DifferentialOperator1DS) :: Differential_operator
     procedure (BC1DS) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, :) 
     procedure (NonLinearSolver), optional:: Solver
    
     integer :: Nv
     logical :: linear1D 
     Nv = size(Solution, dim=2)
        
     linear1D = Linearity_BVP1DS( Differential_operator, Nv )  
    
     if (linear1D) then 
         
        call Linear_Boundary_Value_Problem1DS & 
             (x, Differential_operator,  Boundary_conditions, Solution)
     else 
         
      call Non_Linear_Boundary_Value_Problem1DS( x, & 
         Differential_operator, Boundary_conditions, Solution, Solver)  
      
     end if 
     
end subroutine   



subroutine Linear_Boundary_Value_Problem1DS(x, Differential_operator,  & 
                                            Boundary_conditions, Solution)
                                        
     real, intent(in) :: x(0:)
     procedure (DifferentialOperator1DS) :: Differential_operator
     procedure (BC1DS) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, :) 
     

!  *** variables specification
       integer ::   j, Nx, Nv, M 
       real, allocatable :: A(:,:), U(:), b(:)
       logical ::  dU(size(Solution,dim=2), 2)  ! matrix of dependencies( variable, derivative )
       
!  ***  Integration domain 
        Nx = size(x) - 1 
        Nv = size(Solution, dim=2)  
        M = (Nx+1)*Nv 
       
        dU = BVP1D_system_dependencies( Differential_operator, Nv )

        allocate( A(M, M), b(M), U(M) ) 
         
      
!  ***  independent term  F = A U - b  ( U = inverse(A) b )   
        U = 0
        b = - v_Difference_equations(U)
       
!  ***  Delta kronecker to calculate the difference operator  
        do j=1, size(x)*Nv 
 
         U = 0 
         U(j) = 1   
         A(:, j) = v_Difference_equations(U) + b
         
       end do 
              
!  *** solve the linear system of equations  
       U = Gauss( A, b )
   
       Solution = reshape( U, [ Nx+1, Nv ] ) 
     
contains 

function v_Difference_equations(U) result(F)
     real, target, intent (in) :: U(:)
     real, target :: F(size(U))

    real, pointer :: Ui(:,:), Fi(:,:) 
       
    Ui(0:Nx, 1:Nv) => U
    Fi(0:Nx, 1:Nv) => F
  
    Fi = Difference_equations( dU, Nv, Nx, x, Ui,  & 
                               Differential_operator, Boundary_conditions )
end function         

end subroutine 

 
 
 
 
 
 
 

subroutine Non_Linear_Boundary_Value_Problem1DS( x,  &
          Differential_operator, Boundary_conditions, Solution, Solver )
     real, intent(in) :: x(0:)
     procedure (DifferentialOperator1DS) :: Differential_operator
     procedure (BC1DS) ::  Boundary_conditions
     real, intent(inout) :: Solution(0:, :) 
     procedure (NonLinearSolver), optional:: Solver
     
!  *** variables specification
       integer ::  Nx, Nv, M
       real, allocatable ::  F(:), Us(:)
       logical, allocatable :: dU(:, :)  
 
!  *** Integration domain
       Nx = size(x) - 1; Nv = size(Solution, dim=2); M = (Nx+1)*Nv
       allocate(F(M), Us(M), dU(Nv,2)); dU = .true. 

!  *** initial guess for iteration
       Us = reshape(Solution, [ M ])

!  *** Non linear solver
        if (present(Solver))  then
                   call Solver(Space_discretization, Us)
        else
                   call Newton(Space_discretization, Us)
        end if
        Solution = reshape( Us, [Nx+1, Nv] )
contains
                                               
function Space_discretization(V) result(F)
    real, intent(in) :: V(:)
    real :: F(size(V)) 
    
    F = v_Difference_equations( V )
    
end function                                                    
 
function v_Difference_equations(U) result(F)
     real, target, intent (in) :: U(:)
     real, target :: F(size(U))

    real, pointer :: Ui(:,:), Fi(:,:) 
       
    Ui(0:Nx, 1:Nv) => U
    Fi(0:Nx, 1:Nv) => F
  
    Fi = Difference_equations( dU, Nv, Nx, x, Ui,  & 
                               Differential_operator, Boundary_conditions )
 end function                                               
                                                 
     
 
end subroutine





function Difference_equations( du, Nv, Nx, x, W, & 
                               Differential_operator, Boundary_conditions  ) result(F)
   logical :: dU(Nv, 2) 
   integer :: Nv, Nx 
   real :: x(0:Nx), W(0:Nx, Nv)
   procedure (DifferentialOperator1DS) :: Differential_operator
   procedure (BC1DS) ::  Boundary_conditions
   real ::  F(0:Nx, Nv)  
   
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

end function 




    
end module
