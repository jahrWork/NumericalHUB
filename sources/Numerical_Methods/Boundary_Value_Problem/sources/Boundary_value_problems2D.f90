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
public :: Boundary_Value_Problem2D
public :: Boundary_Value_Problem2D_system

public :: Linear_Boundary_Value_Problem2D
public :: Linear_Boundary_Value_Problem2D_system

public :: Non_Linear_Boundary_Value_Problem2D
public :: Non_Linear_Boundary_Value_Problem2D_system


public :: linear2D

   
abstract interface  

     
    subroutine  NonLinearSolver(Function, x)
       use Jacobian_module
       procedure(FunctionRN_RN) :: Function
       real, intent(inout) :: x(:)
    end subroutine  


    real function DifferentialOperator2D(x, y, u, ux, uy, uxx, uyy, uxy) 
                     real, intent(in) :: x, y, u, ux, uy, uxx, uyy, uxy
    end function  

    real function      BC2D(x, y, u, ux, uy) 
        real, intent(in) :: x, y, u, ux, uy 
    end function  

    function DifferentialOperator2D_system(x, y, u, ux, uy, uxx, uyy, uxy) 
                     real, intent(in) :: x, y, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
                     real :: DifferentialOperator2D_system(size(u))
    end function  

    function BC2D_system(x, y, u, ux, uy) 
        real, intent(in) :: x, y, u(:), ux(:), uy(:) 
        real :: BC2D_system(size(u)) 
    end function  

       

end interface


logical :: optimization = .false.
logical  :: linear2D = .true. 
 
  
 
 
contains 
    
 

!**************************************************************************************************************
! Boundary value problem 2D
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Boundary_Value_Problem2D( x_nodes, y_nodes,  Differential_operator, Boundary_conditions, Solver, Solution)
       
     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2D) :: Differential_operator
     procedure (BC2D) ::  Boundary_conditions
     procedure (NonLinearSolver), optional:: Solver
     real, intent(inout) :: Solution(0:,0:)
     
  
   
     linear2D = Linearity_BVP2D( Differential_operator ) 
    
     if (linear2D) then 
        
              call Linear_Boundary_Value_Problem2D(    x_nodes, y_nodes,  Differential_operator,  & 
                                                       Boundary_conditions, Solution )         
     else 
             call Non_Linear_Boundary_Value_Problem2D( x_nodes, y_nodes,  Differential_operator, &
                                                       Boundary_conditions, Solver, Solution )
     end if 
   

end subroutine 

!**************************************************************************************************************
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Linear_Boundary_Value_Problem2D(x_nodes, y_nodes,  Differential_operator,  Boundary_conditions, Solution)
        
                                        
     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2D) :: Differential_operator
     procedure (BC2D) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, 0:) 
     

!  *** variables specification
       integer ::  i, j, ij, Nx, Ny, M  
       real, allocatable ::  bi(:,:), F(:,:), U(:,:), Ux(:,:), Uxx(:,:), Uy(:,:), Uyy(:,:), Uxy(:,:)
       real, allocatable :: Difference_operator(:,:), V(:), bm(:)

       logical ::  dUs(5) ! matrix of dependencies( variable, derivative )
       
       dUs = BVP2D_dependencies( Differential_operator ) 
   
       
!  *** Integration domain 
       Nx = size(x_nodes) - 1 
       Ny = size(y_nodes) - 1 
       M = (Nx+1)*(Ny+1)

       allocate( bi(0:Nx, 0:Ny), F(0:Nx, 0:Ny), U(0:Nx, 0:Ny), Ux(0:Nx, 0:Ny), Uy(0:Nx, 0:Ny), &
                    Uxy(0:Nx, 0:Ny), Uxx(0:Nx, 0:Ny), Uyy(0:Nx, 0:Ny) )
       allocate( Difference_operator(M, M), bm(M), V(M) ) 
       
     
!  *** independent term  A U = b  ( U = inverse(A) b )         
       U = 0 
       call Difference_equation(x_nodes, y_nodes, U, Ux, Uy, Uxx, Uyy, Uxy, bi) 
       
!  *** Delta kronecker to calculate the difference operator        
       do i=0, Nx 
           do j=0, Ny 
              U = 0 
              U(i,j) = 1.0 
      
              call Difference_equation(x_nodes, y_nodes, U, Ux, Uy, Uxx, Uyy, Uxy, F)
             
              F = F - bi
              ij = i + j * (Nx+1) + 1
              Difference_operator(1:M,ij) = reshape( F, [ M ])
           
          enddo 
       enddo 
       
  
!  *** solve the linear system of equations  
       call LU_Factorization(Difference_operator)
       
        bm = -reshape( bi, [ M ])
 
       V = Solve_LU( Difference_operator, bm ) 
       
       Solution = reshape( V, [Nx+1, Ny+1] ) 
   
       
      deallocate( U, Ux, Uxx, Uy, Uyy, Uxy, F, bi ) 
      deallocate( V, Difference_operator, bm ) 
      

contains 
!-----------------------------------------------------------------
subroutine Difference_equation(x, y, W, Wx, Wy, Wxx, Wyy, Wxy, Fxy) 
           real, intent(in) :: x(0:), y(0:), W(0:,0:)
           real, intent(out) :: Wx(0:,0:), Wy(0:,0:), Wxx(0:,0:), Wyy(0:,0:), Wxy(0:,0:), Fxy(0:, 0:)  
         

    integer ::  i, j, p  
    real :: C, D
    integer :: i0(2), iN(2), di(2), j0(2), jN(2), dj(2)    
   
     i0 = 0; iN = Nx; di = [ 1,  Nx ] 
     j0 = 0; jN = Ny; dj = [ Ny, 1 ]  
    
    
      if ( dUs(1) ) call Derivative( ["x","y"], 1, 1, W, Wx  )
      if ( dUs(2) ) call Derivative( ["x","y"], 2, 1, W, Wy  )
      if ( dUs(3) ) call Derivative( ["x","y"], 1, 2, W, Wxx )
      if ( dUs(4) ) call Derivative( ["x","y"], 2, 2, W, Wyy )
      
      if ( ( dUs(5) ).and.( dUs(1) ) ) then
           
           call Derivative( ["x","y"], 2, 1, Wx, Wxy )
           
       elseif ( ( dUs(5) ).and.( dUs(2) ) ) then
           
           call Derivative( ["x","y"], 1, 1, Wy, Wxy )
           
       elseif ( dUs(5) ) then
           
           call Derivative( ["x","y"], 1, 1,  W, Wx  )
           call Derivative( ["x","y"], 2, 1, Wx, Wxy )
           
       end if
        
    
!  ***  boundary conditions
    do p=1, 2 ! 4 edges for boundary conditions    
        do i = i0(p), iN(p), di(p)  ! z = a  and z = b 
          do j = j0(p), jN(p), dj(p)  
              C = Boundary_conditions(   x(i), y(j),  W(i,j),   Wx(i,j),   Wy(i,j) ) 
              D = Differential_operator( x(i), y(j),  W(i,j),  Wx(i,j),  Wy(i,j),    & 
                                                      Wxx(i,j), Wyy(i,j),  Wxy(i,j)  )   
               
              if (C == FREE_BOUNDARY_CONDITION) then 
                        Fxy(i, j) = D 
              else           
                        Fxy(i, j) = C 
              end if
              
        end do; end do; 
    end do 
           
!  *** inner grid points 
        do i=1, Nx-1
          do j=1, Ny-1
            Fxy(i,j) = Differential_operator( x(i), y(j), W(i,j), Wx(i,j), Wy(i,j), Wxx(i,j), Wyy(i,j), Wxy(i,j) ) 
          enddo 
        enddo 
           
end subroutine 

end subroutine 

!**************************************************************************************************************
! Non Linear Boundary Value Problem2D
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Non_Linear_Boundary_Value_Problem2D( x_nodes, y_nodes, Differential_operator, &
                                                          Boundary_conditions, Solver, Solution)
                                                      

     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2D) :: Differential_operator
     procedure (BC2D) ::  Boundary_conditions
     procedure (NonLinearSolver), optional:: Solver
     real, intent(inout) :: Solution(0:,0:)

!  *** variables specification
       integer ::  Nx, Ny, M
       real, allocatable ::  F(:), U_1(:)
       
       logical ::  dUs(5) ! matrix of dependencies( variable, derivative )
       dUs = BVP2D_dependencies( Differential_operator ) 

!  *** Integration domain
       Nx = size(x_nodes) - 1
       Ny = size(y_nodes) - 1
       M = (Nx+1)*(Ny+1)

       allocate(F(M), U_1(M))

!  *** initial guess for iteration
       U_1 = reshape(Solution, [ M ])

!  *** Non linear solver
        if (present(Solver))  then
                   call Solver(System_BVP, U_1)
        else

               !    Write(*,*)  "... running"

                   call Newton(System_BVP, U_1)
        end if

        Solution = reshape( U_1, [Nx+1, Ny+1] )

     deallocate(F, U_1)
    contains
    
!-----------------------------------------------------------------------
 function System_BVP(U) result(F)
                real, intent (in) :: U(:)
                real :: F(size(U))

                real :: UU(0:Nx,0:Ny), FF(0:Nx,0:Ny)

                UU = reshape( U, [Nx+1, Ny+1] )

                call Difference_equation(x_nodes, y_nodes, UU, FF)

                F = reshape(FF, [ M ]);

 end function 
 
!-----------------------------------------------------------------
 subroutine Difference_equation(x, y, W, Fxy)
            real, intent(in) :: x(0:), y(0:), W(0:,0:)
            real, intent(out) :: Fxy(0:, 0:)
   
    real :: Wx(0:Nx,0:Nx), Wy(0:Nx,0:Ny), Wxx(0:Nx,0:Ny), Wyy(0:Nx,0:Ny), Wxy(0:Nx,0:Ny)
    integer ::  i, j, p   
   
    real :: C, D
    integer :: i0(2), iN(2), di(2), j0(2), jN(2), dj(2)    
   
     i0 = 0; iN = Nx; di = [ 1,  Nx ] 
     j0 = 0; jN = Nx; dj = [ Ny, 1 ]  
    
    
     if (dUs(1))  call Derivative( ["x","y"], 1, 1, W, Wx  ) 
     if (dUs(2))  call Derivative( ["x","y"], 2, 1, W, Wy  )
     if (dUs(3))  call Derivative( ["x","y"], 1, 2, W, Wxx )
     if (dUs(4))  call Derivative( ["x","y"], 2, 2, W, Wyy )

     if ( ( dUs(5) ).and.( dUs(1) ) ) then
           
           call Derivative( ["x","y"], 2, 1, Wx, Wxy )
           
       elseif ( ( dUs(5) ).and.( dUs(2) ) ) then
           
           call Derivative( ["x","y"], 1, 1, Wy, Wxy )
           
       elseif ( dUs(5) ) then
           
           call Derivative( ["x","y"], 1, 1,  W, Wx  )
           call Derivative( ["x","y"], 2, 1, Wx, Wxy )
           
       end if   
    
!  ***  boundary conditions
    do p=1, 2 ! 4 edges for boundary conditions    
        do i = i0(p), iN(p), di(p)  ! z = a  and z = b 
          do j = j0(p), jN(p), dj(p)  
              C = Boundary_conditions(   x(i), y(j),  W(i,j),   Wx(i,j),   Wy(i,j) ) 
              D = Differential_operator( x(i), y(j),  W(i,j),  Wx(i,j),  Wy(i,j),    & 
                                                      Wxx(i,j), Wyy(i,j),  Wxy(i,j)  )   
               
              if (C == FREE_BOUNDARY_CONDITION) then 
                        Fxy(i, j) = D 
              else           
                        Fxy(i, j) = C 
              end if
              
        end do; end do; 
    end do 
           
!  *** inner grid points 
        do i=1, Nx-1
          do j=1, Ny-1
            Fxy(i,j) = Differential_operator( x(i), y(j), W(i,j), Wx(i,j), Wy(i,j), Wxx(i,j), Wyy(i,j), Wxy(i,j) ) 
          enddo 
        enddo 
           
end subroutine

end subroutine

























    
    
    
    
    
    
    
    
!**************************************************************************************************************
! Boundary value problem 2D system 
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Boundary_Value_Problem2D_system(x_nodes, y_nodes,  Differential_operator,   & 
                                           Boundary_conditions, Solution, Solver)
                                        
     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2D_system) :: Differential_operator
     procedure (BC2D_system) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, 0:, :) 
     procedure (NonLinearSolver), optional:: Solver
    
     integer :: Nv
     Nv = size(Solution, dim=3)
        
    linear2D = Linearity_BVP2D_system( Differential_operator, Nv )  
     
     if (linear2D) then 
         
              call Linear_Boundary_Value_Problem2D_system(     x_nodes, y_nodes, Differential_operator, & 
                                                               Boundary_conditions, Solution)
     else 
              call Non_Linear_Boundary_Value_Problem2D_system( x_nodes, y_nodes, Differential_operator, &
                                                               Boundary_conditions, Solution, Solver )         
     end if 
   
  
end subroutine   




!**************************************************************************************************************
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Linear_Boundary_Value_Problem2D_system(x_nodes, y_nodes,  Differential_operator,  & 
                                                  Boundary_conditions, Solution)
        
                                        
     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2D_system) :: Differential_operator
     procedure (BC2D_system) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, 0:, :) 
     

!  *** variables specification
       integer ::  i, j, ijk, k, Nx, Ny, Nv, Ns  
       real, allocatable ::  bi(:,:,:), F(:,:, :), U(:,:,:), Ux(:,:,:), Uxx(:,:,:), Uy(:,:,:), Uyy(:,:,:), Uxy(:,:,:)
       real, allocatable :: Difference_operator(:,:), V(:), bm(:)

       logical ::  dUs(size(Solution,dim=3), 5)  ! matrix of dependencies( variable, derivative )
       
!  *** Integration domain 
       Nx = size(x_nodes) - 1 
       Ny = size(y_nodes) - 1 
       Nv = size(Solution, dim=3)  
       Ns = (Nx+1)*(Ny+1)*Nv 
       
       dUs = BVP2D_system_dependencies( Differential_operator, Nv ) 

       allocate( bi(0:Nx, 0:Ny, Nv), F(0:Nx, 0:Ny, Nv), U(0:Nx, 0:Ny, Nv), Ux(0:Nx, 0:Ny, Nv), Uy(0:Nx, 0:Ny, Nv), &
                    Uxy(0:Nx, 0:Ny, Nv), Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv) )
       allocate( Difference_operator(Ns, Ns), bm(Ns), V(Ns) ) 
       
      
!  *** independent term  A U = b  ( U = inverse(A) b )         
       U = 0 
       call Difference_equation(x_nodes, y_nodes, U, Ux, Uy, Uxx, Uyy, Uxy, bi) 
       
!  *** Delta kronecker to calculate the difference operator   
       do k=1, Nv 
       do i=0, Nx 
           do j=0, Ny 
              U = 0 
              U(i,j,k) = 1.0 
      
             call Difference_equation(x_nodes, y_nodes, U, Ux, Uy, Uxx, Uyy, Uxy, F)
             
             F = F - bi
    !          ijk = i+1 + j*(Nx+1) + (Nx+1)*(Ny+1)*(k-1) 
   !          ijk = k + i*Nv + j * (Nx+1)*Nv
              ijk = v_index(i, j, k) 
     !        Difference_operator(1:M,ijk) = reshape( F, [ M ])
         !     Difference_operator(1:M,ijk) = reshape_to_vector( Nx, Ny, Nv, F )  
              Difference_operator(1:Ns,ijk) = reshape_to_vector( Nx, Ny, Nv, F )   
           
          end do 
       end do 
       end do 
       
  
!  *** solve the linear system of equations  
       call LU_Factorization(Difference_operator)
       
 !      bm = -reshape( bi, [ M ])
       bm = reshape_to_vector( Nx, Ny, Nv, bi ) 
 
       V = Solve_LU( Difference_operator, bm ) 
       
   !    Solution = reshape( V, [Nx+1, Ny+1, Nv] ) 
       call  reshape_to_matrix( Nx, Ny, Nv, V )
   
       
      deallocate( U, Ux, Uxx, Uy, Uyy, Uxy, F, bi ) 
      deallocate( V, Difference_operator, bm ) 
      

contains 
!-----------------------------------------------------------------
subroutine Difference_equation(x, y, W, Wx, Wy, Wxx, Wyy, Wxy, Fxy) 
           real, intent(in) :: x(0:), y(0:), W(0:,0:, :)
           real, intent(out) :: Wx(0:,0:,:), Wy(0:,0:,:), Wxx(0:,0:,:), Wyy(0:,0:,:), Wxy(0:,0:,:), Fxy(0:,0:,:)  
         

    integer ::  i, j, l, p 
    real :: C(Nv), D(Nv)
    integer :: i0(2), iN(2), di(2), j0(2), jN(2), dj(2)    
   
     i0 = 0; iN = Nx; di = [ 1,  Nx ] 
     j0 = 0; jN = Ny; dj = [ Ny, 1 ]  
    
     do i=1, Nv 
       if ( dUs(i,1) ) call Derivative( ["x","y"], 1, 1, W(0:,0:,i), Wx(0:,0:,i)  )
       if ( dUs(i,2) ) call Derivative( ["x","y"], 2, 1, W(0:,0:,i), Wy(0:,0:,i)  )
       if ( dUs(i,3) ) call Derivative( ["x","y"], 1, 2, W(0:,0:,i), Wxx(0:,0:,i) )
       if ( dUs(i,4) ) call Derivative( ["x","y"], 2, 2, W(0:,0:,i), Wyy(0:,0:,i) )
       
       if ( ( dUs(i,5) ).and.( dUs(i,1) ) ) then
           
           call Derivative( ["x","y"], 2, 1, Wx(0:,0:,i), Wxy(0:,0:,i) )
           
       elseif ( ( dUs(i,5) ).and.( dUs(i,2) ) ) then
           
           call Derivative( ["x","y"], 1, 1, Wy(0:,0:,i), Wxy(0:,0:,i) )
       else
           
           call Derivative( ["x","y"], 1, 1, W(0:,0:,i), Wx(0:,0:,i)  )
           call Derivative( ["x","y"], 2, 1, Wx(0:,0:,i), Wxy(0:,0:,i) )
           
       end if
       
     end do   
        
        
!  ***  boundary conditions
    do p=1, 2 ! 4 edges for boundary conditions    
        do i = i0(p), iN(p), di(p)  ! z = a  and z = b 
          do j = j0(p), jN(p), dj(p)  
              C = Boundary_conditions(   x(i), y(j),  W(i,j,:),   Wx(i,j,:),   Wy(i,j,:)    ) 
              D = Differential_operator( x(i), y(j),  W(i,j, :),  Wx(i,j, :),  Wy(i,j,:),    & 
                                                      Wxx(i,j,:), Wyy(i,j,:),  Wxy(i,j,:)  )   
                                                                       
              do l=1, Nv 
                if (C(l) == FREE_BOUNDARY_CONDITION) then 
                        Fxy(i, j, l) = D(l) 
                else           
                        Fxy(i, j, l) = C(l) 
                end if
              end do 
              
        end do; end do; 
    end do 
        
           
!  *** inner grid points 
    do i=1, Nx-1
       do j=1, Ny-1
           Fxy(i,j, :) = Differential_operator( x(i), y(j), W(i,j, :), Wx(i,j, :), Wy(i,j, :), Wxx(i,j, :), Wyy(i,j, :), Wxy(i,j, :) ) 
       enddo 
    enddo 
    
    end subroutine 
 
 subroutine  reshape_to_matrix( Nx, Ny, Nv, V ) 
    integer, intent(in) :: Nx, Ny, Nv 
    real, intent(in) :: V(Nv, 0:Nx, 0:Ny) 
    
     
     integer :: l, m, n
     
     do l=0, Nx
          do m=0, Ny 
              do n=1, Nv 
                 Solution(l, m, n) = V(n, l, m) 
              end do 
          end do 
      end do
     
    
 end subroutine 
 
 !------------------------------------------------------------
function  reshape_to_vector( Nx, Ny, Nv, A ) result(V) 
    integer, intent(in) :: Nx, Ny, Nv 
    real, intent(in) :: A(0:Nx, 0:Ny, Nv ) 
    real :: V( (Nx+1)*(Ny+1)*Nv )
   
    
      integer :: l, m, n,  lmn 
      
      do l=0, Nx
          do m=0, Ny 
              do n=1, Nv 
                  lmn =  v_index(l, m, n)
                  V(lmn) = A(l, m, n)
              end do 
          end do 
      end do
      
    
end function 
 
 !------------------------------------------------------------
integer function  v_index(i, j, k)
    integer, intent(in) :: i, j, k 
    
          
    v_index  = k + i * Nv + j * (Nx+1) * Nv 
    
 end function 

 

end subroutine 










!**************************************************************************************************************
!
!       Differential_operator_system(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions_system(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Non_Linear_Boundary_Value_Problem2D_system( x_nodes, y_nodes,  Differential_operator, &
                                                       Boundary_conditions, Solution, Solver )
                                                      
     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2D_system) :: Differential_operator
     procedure (BC2D_system) ::  Boundary_conditions
     real, intent(inout) :: Solution(0:,0:, :) 
     procedure (NonLinearSolver), optional:: Solver


!  *** variables specification
       integer ::  i, j, k, ijk, Nx, Ny, Nv, M
       real :: x0, xf, y0, yf
       real, allocatable ::  F(:), U_1(:)
       
       logical :: dUs(size(Solution,dim=3), 5)  ! matrix of dependencies( variable, derivative )
      
 
!  *** Integration domain
       Nx = size(x_nodes) - 1
       Ny = size(y_nodes) - 1
       Nv = size(Solution, dim=3)
       M = (Nx+1)*(Ny+1)* Nv 
   
       dUs = BVP2D_system_dependencies( Differential_operator, Nv ) 


       allocate(F(M), U_1(M))

!  *** initial guess for iteration
       U_1 = reshape(Solution, [ M ]);

!  *** Non linear solver
        if (present(Solver))  then
                   call Solver(System_BVP, U_1)
        else

                   call Newton(System_BVP, U_1)
              
        end if

        Solution = reshape( U_1, [Nx+1, Ny+1, Nv] )

     deallocate(F, U_1)

    contains
    
!-----------------------------------------------------------------------
 function System_BVP(U) result(F)
                real, intent (in) :: U(:)
                real :: F(size(U))

                real :: UU(0:Nx, 0:Ny, Nv), FF(0:Nx, 0:Ny, Nv) 
                integer:: i,j,k, ijk

                UU = reshape( U, [Nx+1, Ny+1, Nv] )
  
                call Difference_equation(x_nodes, y_nodes, UU, FF)

                F=reshape(FF, [ M ])
    
 end function 
!--------------------------------------------------------
 subroutine Difference_equation(x, y, W, Fxy)
            real, intent(in) :: x(0:), y(0:), W(0:,0:,:)
            real, intent(out) :: Fxy(0:, 0:,:)
            
            
     real ::   Wx(0:Nx, 0:Ny, Nv), Wy(0:Nx, 0:Ny, Nv)  
     real ::  Wxx(0:Nx,0:Ny, Nv), Wyy(0:Nx, 0:Ny, Nv), Wxy(0:Nx, 0:Ny, Nv)
   

    integer ::  i, j, l, p 
    real :: C(Nv), D(Nv)
    integer :: i0(2), iN(2), di(2), j0(2), jN(2), dj(2)    
   
     i0 = 0; iN = Nx; di = [ 1,  Nx ] 
     j0 = 0; jN = Ny; dj = [ Ny, 1 ]  
           
     do i = 1, Nv 
       if ( dUs(i,1) ) call Derivative( ["x","y"], 1, 1, W(0:,0:,i), Wx(0:,0:,i)  )
       if ( dUs(i,2) ) call Derivative( ["x","y"], 2, 1, W(0:,0:,i), Wy(0:,0:,i)  )
       if ( dUs(i,3) ) call Derivative( ["x","y"], 1, 2, W(0:,0:,i), Wxx(0:,0:,i) )
       if ( dUs(i,4) ) call Derivative( ["x","y"], 2, 2, W(0:,0:,i), Wyy(0:,0:,i) )
       
       if ( ( dUs(i,5) ).and.( dUs(i,1) ) ) then
           
           call Derivative( ["x","y"], 2, 1, Wx(0:,0:,i), Wxy(0:,0:,i) )
           
       elseif ( ( dUs(i,5) ).and.( dUs(i,2) ) ) then
           
           call Derivative( ["x","y"], 1, 1, Wy(0:,0:,i), Wxy(0:,0:,i) )
           
       elseif ( dUs(i,5) ) then
           
           call Derivative( ["x","y"], 1, 1, W(0:,0:,i), Wx(0:,0:,i)  )
           call Derivative( ["x","y"], 2, 1, Wx(0:,0:,i), Wxy(0:,0:,i) )
           
       end if
     end do     
     
        
!  ***  boundary conditions
    do p = 1, 2 ! 4 edges for boundary conditions    
        do i = i0(p), iN(p), di(p)  ! z = a  and z = b 
          do j = j0(p), jN(p), dj(p)  
              
              
              C = Boundary_conditions(   x(i), y(j),  W(i,j,:),   Wx(i,j,:),   Wy(i,j,:)    ) 
              D = Differential_operator( x(i), y(j),  W(i,j, :),  Wx(i,j, :),  Wy(i,j,:),    & 
                                                      Wxx(i,j,:), Wyy(i,j,:),  Wxy(i,j,:)  )   
                                                                       
              do l = 1, Nv 
                if (C(l) == FREE_BOUNDARY_CONDITION) then 
                        Fxy(i, j, l) = D(l) 
                        write(*,*) "FREE BCS"
                else           
                        Fxy(i, j, l) = C(l) 
                        
                end if
              end do 
              
        end do; end do; 
        end do 

           
!  *** inner grid points 
    do i = 1, Nx-1
       do j = 1, Ny-1
           Fxy(i,j, :) = Differential_operator( x(i), y(j), W(i,j, :), Wx(i,j, :), Wy(i,j, :), Wxx(i,j, :), Wyy(i,j, :), Wxy(i,j, :) ) 
       enddo 
    enddo 
           
end subroutine


end subroutine


    
end module


