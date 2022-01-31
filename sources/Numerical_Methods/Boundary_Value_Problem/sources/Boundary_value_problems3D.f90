!******************************************************************************
!
! Author: Juan A Hernandez (juanantonio.hernandez@upm.es) & Javier Escoto 
!******************************************************************************
module Boundary_value_problems3D 

use Linear_systems
use Non_Linear_Systems 
use Collocation_methods 

use Dependencies
 
implicit none  

private
public :: Boundary_Value_Problem3D_system 
public :: linear3D
   
abstract interface  

       subroutine  NonLinearSolver(Function, x)
          use Jacobian_module
          procedure(FunctionRN_RN) :: Function
          real, intent(inout) :: x(:)
       end subroutine  

       
       function DifferentialOperator3D_system(x, y, z, u, ux, uy, uz, uxx, uyy, uzz, uxy, uxz, uyz) 
                        real, intent(in) :: x, y, z, u(:), ux(:), uy(:), uz(:), uxx(:), uyy(:), uzz(:), uxy(:), uxz(:), uyz(:)  
                        real :: DifferentialOperator3D_system(size(u))
       end function  

       function      BC3D_system(x, y, z, u, ux, uy, uz) 
           real, intent(in) :: x, y, z, u(:), ux(:), uy(:), uz(:)  
           real :: BC3D_system(size(u)) 
       end function  
       

end interface


 
logical :: linear3D = .true. 
!logical, allocatable :: dUs(:,:,:) ! matrix of dependencies( direction, order, variable )

logical, allocatable :: dUs(:,:) ! matrix of dependencies( variable, derivative )
 
contains 
    
 
!**************************************************************************************************************
! Boundary value problem 3D
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Boundary_Value_Problem3D_system(x_nodes, y_nodes, z_nodes,  N_variables, Differential_operator,  Boundary_conditions, Solution)
        
                                        
     real, intent(inout) :: x_nodes(0:), y_nodes(0:), z_nodes(0:)
     integer, intent(in) ::  N_variables 
     procedure (DifferentialOperator3D_system) :: Differential_operator
     procedure (BC3D_system) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, 0:, 0:, :) 
     
     allocate( dUs(N_variables, 9) ) 
     dUs = .true. 
  !   linear2D = 
    
     if (linear3D) then 
        
              call Linear_Boundary_Value_Problem3D_system(x_nodes, y_nodes, z_nodes,  N_variables, Differential_operator,  Boundary_conditions, Solution)         
     else 
             
     end if 
     
     deallocate( dUs ) 
   

end subroutine 


!**************************************************************************************************************
!
!       Differential_operator(x, y, u, ux, uy, uxx, uyy, uxy)
!       Boundary_conditions(x, u, ux, uy)  at boundary [a, b]x[c,d] 
!
!**************************************************************************************************************
subroutine Linear_Boundary_Value_Problem3D_system(x_nodes, y_nodes, z_nodes,  N_variables, Differential_operator,  Boundary_conditions, Solution)
        
                                        
     real, intent(inout) :: x_nodes(0:), y_nodes(0:), z_nodes(0:)
     integer, intent(in) ::  N_variables 
     procedure (DifferentialOperator3D_system) :: Differential_operator
     procedure (BC3D_system) ::  Boundary_conditions
     real, intent(out) :: Solution(0:, 0:, 0:, :) 
     

!  *** variables specification
       integer ::  i, j, k, l, ijkl, Nx, Ny, Nz, Nv, M  
       real, allocatable ::  bi(:,:,:, :), F(:,:,:, :), U(:,:,:, :), Ux(:,:,:, :), Uxx(:,:,:, :), Uy(:,:,:, :), Uyy(:,:,:, :), Uxy(:,:,:, :)
       real, allocatable ::  Uz(:,:,:, :), Uzz(:,:,:, :), Uxz(:,:,:, :), Uyz(:,:,:, :)
       real, allocatable :: Difference_operator(:,:), V(:), bm(:)

!  *** Integration domain 
       Nx = size(x_nodes) - 1 
       Ny = size(y_nodes) - 1
       Nz = size(z_nodes) - 1
       Nv = N_variables 
       M = (Nx+1)*(Ny+1)*(Nz+1)*Nv 

       allocate( bi(0:Nx, 0:Ny, 0:Nz, Nv), F(0:Nx, 0:Ny, 0:Nz, Nv), U(0:Nx, 0:Ny, 0:Nz,  Nv) ) 
       allocate(  Ux(0:Nx, 0:Ny, 0:Nz, Nv),   Uy(0:Nx, 0:Ny, 0:Nz,  Nv),   Uz(0:Nx, 0:Ny, 0:Nz,  Nv) ) 
       allocate( Uxx(0:Nx, 0:Ny, 0:Nz, Nv),  Uyy(0:Nx, 0:Ny, 0:Nz, Nv),   Uzz(0:Nx, 0:Ny, 0:Nz, Nv) )
       allocate( Uxy(0:Nx, 0:Ny, 0:Nz,  Nv), Uxz(0:Nx, 0:Ny, 0:Nz, Nv),   Uyz(0:Nx, 0:Ny, 0:Nz, Nv) )
       allocate( Difference_operator(M, M), bm(M), V(M) ) 
       
  
      
!  *** independent term  A U = b  ( U = inverse(A) b )         
       U = 0 
       call Difference_equation(x_nodes, y_nodes, z_nodes, U, Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, bi) 
       
!  *** Delta kronecker to calculate the difference operator   
       do l=1, Nv
           do i=0, Nx 
               do j=0, Ny
                   do k=0, Nz     
       
                     U = 0 
                     U(i,j,k,l) = 1.0 
      
                     call Difference_equation(x_nodes, y_nodes, z_nodes, U, Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, F)
             
                     F = F - bi
                     ijkl = 1 + i + j*(Nx+1) + k*(Nx+1)*(Ny+1) + (l-1)*(Nx+1)*(Ny+1)*(Nz+1)
                     Difference_operator(1:M,ijkl) = reshape( F, [ M ])
         
                   end do 
               end do 
           end do 
       end do 
       
  
!  *** solve the linear system of equations  
       call LU_Factorization(Difference_operator)
       
       bm = - reshape( bi, [ M ])
 
       V = Solve_LU( Difference_operator, bm ) 
       
       Solution = reshape( V, [Nx+1, Ny+1, Nz+1, Nv] ) 
   
       
      deallocate( U, Ux, Uy, Uz, Uxx, Uyy, Uzz, Uxy, Uxz, Uyz, F, bi ) 
      deallocate( V, Difference_operator, bm ) 
      

contains 
!-----------------------------------------------------------------
subroutine Difference_equation(x, y, z, W, Wx, Wy, Wz, Wxx, Wyy, Wzz, Wxy, Wxz, Wyz, Fxyz) 
           real, intent(in) :: x(0:), y(0:), z(0:), W(0:,0:, 0:, :)
           real, intent(out) :: Wx(0:,0:, 0:, :),  Wy(0:,0:, 0:, :),  Wxx(0:,0:, 0:, :), Wyy(0:,0:, 0:, :),  Wxy(0:,0:, 0:, :)
           real, intent(out) :: Wz(0:,0:, 0:, :), Wzz(0:,0:, 0:, :),  Wxz(0:,0:, 0:, :), Wyz(0:,0:, 0:, :), Fxyz(0:,0:, 0:, :)  
         

    integer ::  i, j, k, l, p 
    real :: C(Nv), D(Nv)
    integer :: i0(3), iN(3), di(3), j0(3), jN(3), dj(3), k0(3), kN(3), dk(3)    
   
    i0 = 0; iN = Nx; di = [ 1,  1,  Nx ] 
    j0 = 0; jN = Nx; dj = [ 1,  Ny, 1  ] 
    k0 = 0; kN = Nx; dk = [ Nz, 1,  1  ] 
    
    
        
     do l=1, Nv 
         
        if ( dUs(l,1) ) call Derivative( ["x", "y", "z"], 1, 1, W(0:,0:, 0:, l), Wx(0:,0:, 0:, l) )
        if ( dUs(l,2) ) call Derivative( ["x", "y", "z"], 2, 1, W(0:,0:, 0:, l), Wy(0:,0:, 0:, l) )
        if ( dUs(l,3) ) call Derivative( ["x", "y", "z"], 3, 1, W(0:,0:, 0:, l), Wz(0:,0:, 0:, l) )
        
        if ( dUs(l,4) ) call Derivative( ["x", "y", "z"], 1, 2, W(0:,0:, 0:, l), Wxx(0:,0:, 0:, l) )
        if ( dUs(l,5) ) call Derivative( ["x", "y", "z"], 2, 2, W(0:,0:, 0:, l), Wyy(0:,0:, 0:, l) )
        if ( dUs(l,6) ) call Derivative( ["x", "y", "z"], 3, 2, W(0:,0:, 0:, l), Wzz(0:,0:, 0:, l) )
        
        if ( dUs(l,7) ) call Derivative( ["x", "y", "z"], 2, 1, Wx(0:,0:, 0:, l), Wxy(0:,0:, 0:, l) )
        if ( dUs(l,8) ) call Derivative( ["x", "y", "z"], 3, 1, Wx(0:,0:, 0:, l), Wxz(0:,0:, 0:, l) )
        if ( dUs(l,9) ) call Derivative( ["x", "y", "z"], 3, 1, Wy(0:,0:, 0:, l), Wyz(0:,0:, 0:, l) )
        
     end do    
        
!  ***  boundary conditions
    do p=1, 3 ! 6 planes for boundary conditions    
        do i = i0(p), iN(p), di(p)  ! z = a  and z = b 
          do j = j0(p), jN(p), dj(p)   
            do k = k0(p), kN(p), dk(p)  
              C = Boundary_conditions(   x(i), y(j), z(k), W(i,j,k,:),   Wx(i,j,k,:),   Wy(i,j,k,:),  Wz(i,j,k,:)    ) 
              D = Differential_operator( x(i), y(j), z(k), W(i,j,k, :),  Wx(i,j,k, :),  Wy(i,j,k,:),  Wz(i,j,k,:),   & 
                                                           Wxx(i,j,k,:), Wyy(i,j,k, :), Wzz(i,j,k, :),               & 
                                                           Wxy(i,j,k, :),Wxz(i,j,k, :), Wyz(i,j,k, :)                )  
              do l=1, Nv 
                if (C(i) == FREE_BOUNDARY_CONDITION) then 
                        Fxyz(i, j, k, l) = D(l) 
                else           
                        Fxyz(i, j, k, l) = C(l) 
                end if
              end do 
              
        end do; end do; end do  
    end do 
    
           
!  *** inner grid points 
        do i=1, Nx-1
          do j=1, Ny-1
            do k=1, Nz-1 
            Fxyz(i, j, k, :) = Differential_operator( x(i), y(j), z(k), W(i,j, k, :),              &
                                                      Wx(i,j,k, :), Wy(i,j,k, :), Wz(i,j,k, :),   & 
                                                      Wxx(i,j,k, :), Wyy(i,j,k, :), Wzz(i,j,k, :),& 
                                                      Wxy(i,j,k, :), Wxz(i,j,k, :), Wyz(i,j,k, :)     ) 
         
            end do 
          end do 
        end do 
           
end subroutine 

end subroutine 



 
    
    
    
    
end module
