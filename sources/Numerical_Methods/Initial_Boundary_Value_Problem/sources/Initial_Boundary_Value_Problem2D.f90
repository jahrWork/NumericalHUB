
!***********************************************************************
! It integrates in time the Initial value boundary problem 2D.
! Given the differential operator and the boundary equations,
! the discrete solution is calculated.
! author: Juan A Hernandez, juanantonio.hernandez@upm.es 
!***********************************************************************
module Initial_Boundary_Value_Problem2D

use Cauchy_Problem
use ODE_Interface 
use Temporal_scheme_interface
use Collocation_methods 
use Non_Linear_Systems
use Utilities
use Dependencies_IBVP2D
use Dependencies_BC
use plots

implicit none   

private
public :: IBVP2DS,                    & 
          Spatial_discretization2DS,  &  
          Linear_operator2DS,         & 
          Spatial_truncation_error2DS 


   
abstract interface  

 function DifferentialOperator2DS(x, y, t, u, ux, uy, uxx, uyy, uxy) 
  real, intent(in) :: x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real :: DifferentialOperator2DS (size(u)) 
  end function  

  function BC2DS(x, y, t, u, ux, uy) 
      real, intent(in) :: x, y, t, u(:), ux(:), uy(:) 
      real :: BC2DS(size(u)) 
  end function  

end interface

 
type BoundaryP2DS 
     
        integer :: i, j                    ! grid index (i,j)
        logical, allocatable :: impose(:)  ! at least Nv boudary contitions to impose per point  
        integer, allocatable :: B_index(:) ! boundary index of variable to determine (Nb unknowns) 
                                           ! at last Nv unknowns per bounday point
        
        logical, allocatable :: interfaces(:)  ! at last Nv interfaces per point 
        
end type


contains



!******************************************************************************
!*  Initial value boundary problem 2D. Vectorial case
!******************************************************************************
subroutine IBVP2DS( Time_Domain, x_nodes, y_nodes, Differential_operator, & 
                    Boundary_conditions, Solution, Scheme ) 

     real, intent(in) :: Time_Domain(0:) 
     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2DS) :: Differential_operator 
     procedure (BC2DS) ::  Boundary_conditions
     procedure (Temporal_Scheme), optional :: Scheme
     real, intent(out) :: Solution(0:, 0:, 0:, :) 
     
     real, pointer :: U_Cauchy(:, :) 
     integer :: Nx, Ny, Nt, Nv
     Nx = size(x_nodes)-1;  Ny = size(y_nodes)-1; Nt = size(Time_Domain)-1
     Nv = size(Solution, dim=4) 
         
    call my_reshape( Solution, Nt+1, (Nx+1)*(Ny+1)*Nv, U_Cauchy ) 
    
    CPU_time_BC = 0 
    call Cauchy_ProblemS( Time_Domain, F_Cauchy, U_Cauchy, Scheme )
    write(*, '("Boundary conditions, CPU Time=",f6.3," seconds.")') CPU_time_BC 
    
contains

function  F_Cauchy( U, t ) result(F) 
          real ::  U(:), t         
          real :: F(size(U)) 
          
    F = Space_discretization2D( Nx, Ny, Nv, Differential_operator, Boundary_conditions,  & 
                                x_nodes, y_nodes, U, t )
           
end function 
end subroutine

function  Space_discretization2D( Nx, Ny, Nv, Differential_operator, Boundary_conditions, & 
                                  x, y, U, t ) result(F) 

      integer, intent(in) :: Nx, Ny, Nv
      procedure (DifferentialOperator2DS) :: Differential_operator 
      procedure (BC2DS) ::  Boundary_conditions
      real, intent(in) :: x(0:), y(0:)
      real, target ::  U(:), t, F(size(U))   
  
  real, pointer :: Uv(:, :, :), Fv(:, :, :)  
  
  Uv(0:Nx, 0:Ny, 1:Nv) => U 
  Fv(0:Nx, 0:Ny, 1:Nv) => F 
  
  Fv = Spatial_discretization2DS(  Differential_operator,  Boundary_conditions, & 
                                   x, y, Uv, t )
 
  call write_residuals(t,  Fv) 
 
end function  


function Spatial_discretization2DS( Differential_operator, &     
                                    Boundary_conditions,   & 
                                    x, y, U, t) result(F) 

  procedure (DifferentialOperator2DS) :: Differential_operator 
  procedure (BC2DS) :: Boundary_conditions
  real, intent(in) :: x(0:), y(0:), t
  real, intent(inout) :: U(0:,0:,:) 
  real :: F(0:size(U, dim=1)-1, 0:size(U,dim=2)-1, size(U,dim=3))
  
    real, save, allocatable :: Ux(:,:,:), Uy(:,:,:), Uxx(:,:,:), Uxy(:,:,:),  Uyy(:,:,:)
    real, save, allocatable :: Fi(:), Ub(:)  
    integer :: Nx, Ny, Nv, i, j, k, l 
    type (BoundaryP2DS), save, allocatable :: BC(:)
    logical, save, allocatable :: dU(:,:)  !dependencies Operator(variables, derivative)
    integer, save :: Nb 
   
    
    Nx = size(x)-1; Ny = size(y)-1; Nv = size(U, dim=3)  
 
!  ***  Allocate internal variables and identifies the boundary unknowns    
        if (t==0) then 
          call Initialize_block
          call Boundary_unknowns 
          call Pointers_and_allocation
        end if 
       
    
!  ***  It solves boundary equations to yield values at boundaries  
        if (Nb>0) then 
         call Determine_boundary_points( Nb,  BC, x, y, Nv, U, Ux, Uy, t, Boundary_conditions ) 
        end if 
         
!  ***  It calculates only derivatives which are used 
        call Involved_DerivativesS( dU, U, Ux, Uy, Uxx, Uyy, Uxy )
   
        do i=0, Nx; do j=0, Ny
             
             Fi = Differential_operator( x(i), y(j), t, U(i,j,:), Ux(i,j,:), & 
                                         Uy(i,j,:), Uxx(i,j,:), Uyy(i,j,:), Uxy(i,j,:) ) 
             if (Nb>0) then 
          ! ** Boundary point when k>0
               k = Boundary_index(Nx, Ny, i, j)
               if (k>0) then     
                  do l=1, Nv
                    if (BC(k) % impose(l) .or. BC(k) % interfaces(l)) Fi(l) = IMPOSE_ZERO
                  end do 
               end if 
             end if 
             
        ! ** inner point  
             F(i, j, :) = Fi(:)
         
        enddo; enddo
contains 
  

subroutine Initialize_block                                      
                                     
       if (allocated(dU)) deallocate(dU)  
       allocate( dU(Nv,5) )
       dU  = IBVP2D_Dependencies_system( Nv, Differential_operator )
        
       if (allocated(BC)) deallocate(BC)
       allocate( BC( 2*(Ny+1) + 2*(Nx-1) ) )
       do k=1, size(BC) 
           allocate( BC(k) % impose(Nv), BC(k) % interfaces(Nv), BC(k) % B_index(Nv) ) 
       end do
       
end subroutine  

subroutine  Pointers_and_allocation

   if (allocated(Ux)) deallocate(Ux, Uy, Uxx, Uxy, Uyy, Fi, Ub) 
   allocate(  Ux(0:Nx, 0:Ny, Nv),  Uy(0:Nx, 0:Ny, Nv), & 
              Uxx(0:Nx, 0:Ny, Nv), Uxy(0:Nx, 0:Ny, Nv),  Uyy(0:Nx, 0:Ny, Nv) ) 
   allocate( Fi(Nv), Ub(Nb) ) 


end subroutine

subroutine Boundary_unknowns
   
       integer :: i, j, k, l, m 
       real :: u0(Nv), ux0(Nv), uy0(Nv), eq(Nv)   
     
  
  u0 = 1.; ux0 = 2.; uy0 = 3  
  
  Nb = 0   
!  Fourier  
 if (method == "Fourier") then 
       do k=1, size(BC) 
           BC(k) % impose = .false.     
       end do   
       
 ! FD or Chebyshev       
 else   
  
  Nb = 0 
  do k=1, size(BC)
      
    call  ij_index( Nx, Ny, k, i, j) 
      
    eq = Boundary_conditions (x(i),  y(j),  t, u0,  ux0,  uy0 ) 
    BC(k) % i = i 
    BC(k) % j = j 
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
   
 endif  
 write(*,*) " Number of unknown boundary points:", Nb 
  
end subroutine 
 
    
  
end function   
       


subroutine Determine_boundary_points( Nb, BC, x, y, Nv, U, Ux, Uy, t, Boundary_conditions ) 
        integer, intent(in) :: Nb  
        type (BoundaryP2DS) :: BC(:)
        real, intent(in) :: x(0:), y(0:), t
        integer, intent(in) :: Nv 
        real, intent(inout) :: U(0:, 0:, :),   Ux(0:, 0:, :),  Uy(0:, 0:, :) 
        procedure (BC2DS) ::  Boundary_conditions
 
      integer ::  i, j, k, m, l, r, rmax = 1000 
      real :: Gc(Nv), t0, tf  
      real :: dZ(Nb), dZ0(Nb), G(Nb), G0(Nb), Error, eps = 1d-7
     
  call cpu_time(t0)  
  dZ = 0; r = 1; Error = 1 
 
  do while (Error > eps .and. r < rmax) 
  
    do l=1, Nv           
       call Derivative( [ "x", "y" ], 1, 1, U(:, :, l), Ux(:, :, l) )
       call Derivative( [ "x", "y" ], 2, 1, U(:, :, l), Uy(:, :, l) )
    end do 
 
!** Boundary equations 
    do k=1, size(BC) 
       i = BC(k) % i 
       j = BC(k) % j 
       Gc = Boundary_conditions (x(i), y(j), t, U(i, j, :), Ux(i, j, :), Uy(i, j, :) )  
      
       do l=1, Nv  
           
         if ( BC(k)% impose(l) ) then   
             
             m = BC(k) % B_index(l) 
             G(m) = Gc(l) 
             
             ! Secant method
             if (r>1) then 
                 dZ(m) = - G(m) / (G(m)-G0(m)) * dZ0(m) 
             else 
                 dZ(m) = 1d-3
             end if 
             U(i, j, l) = U(i, j, l) + dZ(m)
             
         end if 
       end do 
    end do
    
    G0 = G; dZ0 = dZ 
    Error = norm2(dZ) 
    r = r + 1 
    !write(*,*) " Error =", Error; read(*,*) 
    
  end do 
  
  call cpu_time(tf) 
  CPU_time_BC = CPU_time_BC + tf-t0 
 
 

end subroutine 




!*********************************************************************************
! vector of dependencies(derivatives) , ux ,uy, uxx uyy, uxy  
!*********************************************************************************
subroutine Involved_DerivativesS( dU, U, Ux, Uy, Uxx, Uyy, Uxy )
  logical, intent(in) :: dU(:,:) 
  real, intent(in) :: U(:,:,:)
  real, intent(out) :: Ux(:,:,:), Uy(:,:,:), Uxx(:,:,:), Uyy(:,:,:), Uxy(:,:,:)  
        
     integer :: k 
  
!  *** inner grid points
        do k=1, size(DU, dim=1)
                    
           if ( dU(k,1) ) call Derivative( [ "x", "y" ], 1, 1, U(:,:,k), Ux(:,:,k)  )
           if ( dU(k,2) ) call Derivative( [ "x", "y" ], 2, 1, U(:,:,k), Uy(:,:,k)  )
           if ( dU(k,3) ) call Derivative( [ "x", "y" ], 1, 2, U(:,:,k), Uxx(:,:,k) )
           if ( dU(k,4) ) call Derivative( [ "x", "y" ], 2, 2, U(:,:,k), Uyy(:,:,k) )
           if ( dU(k,5) ) call Derivative( [ "x", "y" ], 2, 1, Ux(:,:,k),Uxy(:,:,k) )
           
        end do  
       
  

end subroutine 

!************************************************************************
!* Relation between i,j boundary points and k index associated to the boundary 
!    
!  j=5 | 6 17 18 19 20 12 
!  j=4 | 5             11  
!  j=3 | 4             10 
!  j=2 | 3  k values    9
!  j=1 | 2              8
!  j=0 | 1 13 14 15 16  7
!  ______________________  
!  i=   0  1  2  3  4  5 
!     
!************************************************************************
! It gives i,j from k index 
!************************************************************************    
subroutine ij_index(Nx, Ny, k, i, j) 
   integer, intent(in) :: Nx, Ny, k
   integer, intent(out) :: i, j 
   
    if (k<=Ny+1) then 
                       i = 0; j = k-1 
                       
    else if (k<=2*Ny+2) then 
        
                       i = Nx; j = k - Ny - 2 
                       
    else if (k<=2*Ny+Nx+1) then 
        
                       j = 0; i = k - 2 * Ny - 2   
    else 
                       j = Ny; i = k - 2*Ny - Nx - 1  
    endif 

end subroutine 


!**************************************************************
! It gives k from i,j indexes 
!**************************************************************
integer function Boundary_index(Nx, Ny, i,j) result(k) 
         integer, intent(in) :: Nx, Ny, i, j 

    if (i==0) then 
        
               k = j + 1  
               
    else if (i==Nx) then 
        
               k = Ny + 2 + j 
               
    else if (j==0) then 
        
               k =  2*Ny + 2  + i 
               
    else if (j==Ny) then  
        
               k =  2*Ny + Nx + 1  + i 
    else 
               k = -1 
    end if 
    

end function     
    


subroutine write_residuals(t,  F)
 real, intent(in) :: t, F(:, :, :) 
 
  integer :: i, Nv 
  real, save :: t1
  logical :: write_R 
    
  Nv = size(F, dim=3) 

  if (t==0) then 
      t1 = 0 
      write_R = .true. 
      
  else if (int(10*(t-t1)) == 1) then 
      write_R = .true.
      
  end if 
  
  if (write_R) then 
   write(*,'(A, f5.2, A, 10e15.7)') " t=", t, "   norm2(Residuals) =",  (norm2( F(:,:,i) ), i=1, Nv)  
   t1 = t 
  end if 

end subroutine 






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
! Author: juanantonio.hernandez@um.es (May 2021) 
!********************************************************************************** 
function Spatial_Truncation_Error2DS( Nv, Differential_operator,       &
                                      Boundary_conditions,  x, y,      &
                                      Order, Test_function ) result(R) 

     procedure (DifferentialOperator2DS) :: Differential_operator 
     procedure (BC2DS) :: Boundary_conditions 
     integer, intent(in) :: Nv, Order   
     real, intent(inout) :: x(0:), y(0:)
     interface 
        function Test_function(Nv, x, y) result(U) 
           integer, intent(in) :: Nv 
           real, intent(in) :: x(0:), y(0:) 
           real :: U(0:size(x)-1, 0:size(y)-1, Nv) 
        end function 
     end interface 
     real :: R( 0:size(x)-1, 0:size(y)-1, Nv ) 

     integer :: i, j, Nx, Ny
     real, allocatable :: x1(:), y1(:), U1(:,:,:), F1(:,:,:)
     real, allocatable :: x2(:), y2(:), U2(:,:,:), F2(:,:,:)
     real :: t 
       
       t = 0 
       Nx = size(x)-1;  Ny = size(y)-1
       allocate( x1(0:Nx), y1(0:Ny), U1(0:Nx, 0:Ny, Nv), F1(0:Nx, 0:Ny, Nv) ) 
       allocate( x2(0:2*Nx), y2(0:2*Ny), U2(0:2*Nx, 0:2*Ny, Nv), F2(0:2*Nx, 0:2*Ny, Nv) )
         
       x1 = x; y1 = y
       call Grid_Initialization( "nonuniform", "x", x1, Order )
       call Grid_Initialization( "nonuniform", "y", y1, Order )
       
       U1 = Test_function(Nv, x1, y1) 
       F1 = Spatial_discretization2DS( Differential_operator, Boundary_conditions, x1, y1, U1, t)  
          
       do i=0, Nx-1 
           x2(2*i)   = x1(i) 
           x2(2*i+1) = ( x1(i) + x1(i+1) )/2
       end do 
       x2(2*Nx) = x1(Nx)
       
       do j=0, Ny-1 
           y2(2*j)   = y1(j) 
           y2(2*j+1) = ( y1(j) + y1(j+1) )/2
       end do 
       y2(2*Ny) = y1(Ny)
       call Grid_Initialization( "unmodified", "x", x2, Order )
       call Grid_Initialization( "unmodified", "y", y2, Order )
       
       U2 = Test_function(Nv, x2, y2) 
       F2 = Spatial_discretization2DS( Differential_operator, Boundary_conditions, x2, y2, U2, t)
     
       do i=0, Nx; do j=0, Ny
            R(i, j, :) = ( F2(2*i, 2*j, :) - F1(i, j, :) )/( 1 - 1./2**Order )  
       end do; end do   
                
end function 

!*******************************************************************
! Given a vector function F: RN -> RN. 
! If the F (differential operator) is linear (F = A U + b), 
! it gives the system matrix A
!*******************************************************************
function Linear_operator2DS( Nv, x, y, Order, Differential_operator,  & 
                             Boundary_conditions ) result(A)

     integer, intent(in) :: Nv, Order
     real, intent(inout) :: x(0:), y(0:)  
     procedure (DifferentialOperator2DS) :: Differential_operator 
     procedure (BC2DS) ::  Boundary_conditions
     real :: A(size(x)*size(y)*Nv, size(x)*size(y)*Nv) 
     
          
     real ::  U(size(x)*size(y)*Nv), b(size(x)*size(y)*Nv), t  
     integer :: j
        
     call Grid_Initialization( "unmodified", "x", x, Order )
     call Grid_Initialization( "unmodified", "y", y, Order )
     t = 0 
     U = 0
     b = Space_discretization2D( size(x)-1, size(y)-1, Nv, Differential_operator, Boundary_conditions, x, y, U, t ) 
        
     do j=1, size(x)*size(y)*Nv
         
         U = 0 
         U(j) = 1  
         A(:, j) = Space_discretization2D( size(x)-1, size(y)-1, Nv, Differential_operator, Boundary_conditions, x, y, U, t ) - b
         
     end do 
             
end function  




end module 
    
    