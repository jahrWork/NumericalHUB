!******************************************************************************
!
! Authors: Juan A Hernandez (juanantonio.hernandez@upm.es)  
!          Javier Escoto (javier.escoto.lopez@alumnos.upm.es)
!******************************************************************************
module Boundary_value_problems2D 

use Linear_systems
use Non_Linear_Systems 
use Collocation_methods 

use Linearity
use Dependencies2D

 
implicit none  

private
public :: Boundary_Value_Problem2DS
public :: Linear_Boundary_Value_Problem2DS
public :: Non_Linear_Boundary_Value_Problem2DS


public :: linear2D

   
abstract interface  

     
    subroutine  NonLinearSolver(Function, x)
       use Jacobian_module
       procedure(FunctionRN_RN) :: Function
       real, intent(inout) :: x(:)
    end subroutine  

    function DifferentialOperator2DS(x, y, u, ux, uy, uxx, uyy, uxy) 
                     real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
                     real :: DifferentialOperator2DS(size(u))
    end function  

    function BC2DS(x, y, u, ux, uy) 
        real, intent(in) :: x, y, u(:), ux(:), uy(:) 
        real :: BC2DS(size(u)) 
    end function  

end interface

logical  :: linear2D
 
contains 
 
!***************************************************************************************
! Boundary value problem 2D system 
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy) = 0  at inner points (x_i, y_j ) 
!       Boundary_conditions(x, u, ux, uy) = 0   at boundary points 
!
!****************************************************************************************
subroutine Boundary_Value_Problem2DS(x, y,  Differential_operator,   & 
                                     Boundary_conditions, Solution, Solver)
                                        
     real, intent(in) :: x(0:), y(0:)
     procedure (DifferentialOperator2DS) :: Differential_operator
     procedure (BC2DS) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, 0:, :) 
     procedure (NonLinearSolver), optional:: Solver
    
     integer :: Nv
     Nv = size(Solution, dim=3)
        
     linear2D = Linearity_BVP2D_system( Differential_operator, Nv )    
       
     if (linear2D) then 
              
         call Linear_Boundary_Value_Problem2DS(   x, y, Differential_operator, & 
                                                  Boundary_conditions, Solution)
     else 
         call Non_Linear_Boundary_Value_Problem2DS( x, y, Differential_operator, &
                                                    Boundary_conditions, Solution, Solver )         
     end if 
   
  
end subroutine   


subroutine Linear_Boundary_Value_Problem2DS(x, y,  Differential_operator,  & 
                                            Boundary_conditions, Solution)
     real, intent(in) :: x(0:), y(0:)
     procedure (DifferentialOperator2DS) :: Differential_operator
     procedure (BC2DS) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, 0:, :) 
     

!  *** variables specification
       integer ::   Nx, Ny, Nv, M, i, j, k, ijk
       real, target, allocatable ::  A(:,:), F(:), b(:) 
       real, pointer ::  Fi(:,:,:), bi(:,:,:)
       logical ::  dUs(size(Solution,dim=3), 5)  ! matrix of dependencies( variable, derivative )
       real :: kappa 
       
!  *** Integration domain 
       Nx = size(x) - 1; Ny = size(y) - 1; Nv = size(Solution, dim=3)  
       M = (Nx+1)*(Ny+1)*Nv 
       allocate( A(M, M), F(M), b(M) )
       
       dUs = BVP2D_system_dependencies( Differential_operator, Nv ) 
       
       Fi(0:Nx, 0:Ny, 1:Nv) => F(1:M)
       bi(0:Nx, 0:Ny, 1:Nv) => b(1:M)
      
!  *** independent term  F = A U - b  
       Solution = 0
       bi =  Difference_equations( dUs, x, y,  Solution, Differential_operator, Boundary_conditions ) 
       
!  *** Delta kronecker to calculate the difference operator   
       do k=1, Nv; do i=0, Nx; do j=0, Ny   
       
              Solution = 0;   Solution(i,j,k) = 1.0 
      
              Fi =  Difference_equations( dUs, x, y,  Solution, Differential_operator, Boundary_conditions ) - bi 
              
              ijk = i+1 + j*(Nx+1) + (Nx+1)*(Ny+1)*(k-1) 
                          
              A(:,ijk) = F
       
       end do; end do; end do
  
!  *** solve the linear system of equations  
   !   ! kappa = Condition_number(Difference_operator) 
       F = Gauss(A, b)
          
       Solution = reshape( F, [Nx+1, Ny+1, Nv] )
  

end subroutine                                            
                                           
  

subroutine Non_Linear_Boundary_Value_Problem2DS( x, y,  Differential_operator, &
                                                 Boundary_conditions, Solution, Solver )
     real, intent(in) :: x(0:), y(0:)
     procedure (DifferentialOperator2DS) :: Differential_operator
     procedure (BC2DS) ::  Boundary_conditions
     real, intent(inout) :: Solution(0:,0:, :) 
     procedure (NonLinearSolver), optional:: Solver

!  *** variables specification
       integer :: Nx, Ny, Nv, M
       real, allocatable ::   U(:)
       logical ::  dUs(size(Solution,dim=3), 5)  ! matrix of dependencies( variable, derivative )
       
      
!  *** Integration domain
       Nx = size(x) - 1; Ny = size(y) - 1; Nv = size(Solution, dim=3)
       M = (Nx+1)*(Ny+1)* Nv 
       allocate(U(M))
       
       dUs = BVP2D_system_dependencies( Differential_operator, Nv )
   
!  *** initial guess for iteration
       U = reshape(Solution, [ M ]);

!  *** Non linear solver
        if (present(Solver))  then
                   call Solver(F, U)
        else

                   call Newton(F, U)
              
        end if

        Solution = reshape( U, [Nx+1, Ny+1, Nv] )

   
contains
 function F(U) 
      real,  intent (in) :: U(:)
      real :: F(size(U))

   F =  Difference_equations_v( U ) 
    
 end function                                                        

 function Difference_equations_v(U) result(F)
      real, target, intent (in) :: U(:)
     real, target :: F(size(U))

    real, pointer :: Ui(:,:,:), Fi(:,:,:) 
       
    Ui(0:Nx, 0:Ny, 1:Nv) => U
    Fi(0:Nx, 0:Ny, 1:Nv) => F
    
    Fi = Difference_equations( dUs, x, y, Ui, Differential_operator, Boundary_conditions ) 
    
 end function 
end subroutine

  






function Difference_equations( dU, x, y,  W, Differential_operator, Boundary_conditions ) result(F) 
     logical, intent(in) ::  dU(:,:) 
     real, intent(in) :: x(0:), y(0:), W(0:,0:, :)
     procedure (DifferentialOperator2DS) :: Differential_operator  
     procedure (BC2DS) ::  Boundary_conditions
     real ::  F(0:size(x)-1, 0:size(y)-1, size(W,dim=3) )  
         
    integer ::  i, j, k, Nx, Ny, Nv 
    real, allocatable :: Fc(:), Fi(:) 
    real, allocatable ::  Wx(:,:,:), Wxx(:,:,:), Wy(:,:,:), Wyy(:,:,:), Wxy(:,:,:)
    
    Nx = size(x) -1; Ny = size(y) -1; Nv = size(W, dim=3)
    allocate(   Wx(0:Nx, 0:Ny, Nv),  Wy(0:Nx, 0:Ny, Nv),    &
                Wxy(0:Nx, 0:Ny, Nv), Wxx(0:Nx, 0:Ny, Nv), Wyy(0:Nx, 0:Ny, Nv) )
    allocate( Fc(Nv), Fi(Nv) )
    
    do i = 1, Nv 
       if ( dU(i,1) ) call Derivative( ["x","y"], 1, 1, W(0:,0:,i), Wx(0:,0:,i)  )
       if ( dU(i,2) ) call Derivative( ["x","y"], 2, 1, W(0:,0:,i), Wy(0:,0:,i)  )
       if ( dU(i,3) ) call Derivative( ["x","y"], 1, 2, W(0:,0:,i), Wxx(0:,0:,i) )
       if ( dU(i,4) ) call Derivative( ["x","y"], 2, 2, W(0:,0:,i), Wyy(0:,0:,i) )
       
       if ( ( dU(i,5) ).and.( dU(i,1) ) ) then
           
           call Derivative( ["x","y"], 2, 1, Wx(0:,0:,i), Wxy(0:,0:,i) )
           
       elseif ( ( dU(i,5) ).and.( dU(i,2) ) ) then
           
           call Derivative( ["x","y"], 1, 1, Wy(0:,0:,i), Wxy(0:,0:,i) )
           
       elseif ( dU(i,5) ) then
           
           call Derivative( ["x","y"], 1, 1, W(0:,0:,i), Wx(0:,0:,i)  )
           call Derivative( ["x","y"], 2, 1, Wx(0:,0:,i), Wxy(0:,0:,i) )
           
       end if
    end do     
    
  
     
    do i=0, Nx; do j=0, Ny
        
           Fi =  Differential_operator( x(i), y(j), W(i,j, :), Wx(i,j, :), Wy(i,j, :), & 
                                        Wxx(i,j, :), Wyy(i,j, :), Wxy(i,j, :) ) 
           
           if (i==0 .or. i==Nx .or. j==0 .or. j==Ny) then 
               Fc = Boundary_conditions( x(i), y(j),  W(i,j,:), Wx(i,j,:), Wy(i,j,:) )
               do k=1, Nv 
                  if (Fc(k) /= FREE_BOUNDARY_CONDITION ) then 
                      Fi(k) = Fc(k) 
                  end if 
               end do 
           end if 
           
           F(i,j,:) = Fi 
        
    enddo; enddo 
    
    
end function 

    
end module

