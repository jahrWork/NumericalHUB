
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

implicit none   

private

public :: IBVP2D, IBVP2D_system
public :: Spatial_discretization2D
public :: Spatial_discretization2D_system
   
abstract interface  

             
  real function DifferentialOperator2D(x, y, t, u, ux, uy, uxx, uyy, uxy) 
                   real, intent(in) :: x, y, t, u, ux, uy, uxx, uyy, uxy
  end function  

  real function      BC2D(x, y, t, u, ux, uy) 
      real, intent(in) :: x, y, t, u, ux, uy 
  end function 
  
 function DifferentialOperator2D_system(x, y, t, u, ux, uy, uxx, uyy, uxy) 
                   real, intent(in) :: x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
                   real :: DifferentialOperator2D_system (size(u)) 
  end function  

  function      BC2D_system(x, y, t, u, ux, uy) 
      real, intent(in) :: x, y, t, u(:), ux(:), uy(:) 
      real :: BC2D_system(size(u)) ! maximum number of boundary conditions at each point
  end function  


  function DifferentialOperatorODE(U, t) result(F) 
                   real, intent(in) :: U(:), t 
                   real :: F(size(U)) 
  end function 

end interface


 type Boundary2D 
     
        integer :: i, j     ! associated i,j grid points for thid boundary point 
        real :: value       ! value of the function at this boundary point 
        real :: equation    ! value of the equation in the boundary 
        logical :: impose   ! true if there is boundary condition
        
 end type
 
type Boundary2D_system 
     
        integer :: i, j 
        real, allocatable :: value(:) 
        real, allocatable :: equation(:)  ! Nv equations at the boundary point 
        logical, allocatable :: impose(:) ! at least Nv boudary contitions to impose per point  
        
end type
 

contains



!*******************************************************************
!* Initial value boundary problem 2D. Escalar case 
!*******************************************************************
subroutine IBVP2D( Time_Domain, x_nodes, y_nodes, Differential_operator,  Boundary_conditions, Solution, Scheme ) 

     real, intent(in) :: Time_Domain(:) 
     real, intent(inout) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2D) :: Differential_operator 
     procedure (BC2D) ::  Boundary_conditions  
     real, intent(out) :: Solution(0:, 0:, 0:) 
     procedure (Temporal_Scheme), optional :: Scheme
   

     integer :: Nt, Nx, Ny, M 
     real, pointer :: U_Cauchy(:, :) 
     Nt = size(Time_Domain)-1;  Nx = size(x_nodes)-1; Ny = size(y_nodes)-1
     M = (Nx+1) * (Ny+1)   
       
     call my_reshape( Solution, Nt+1, M, U_Cauchy ) 
    
     call Cauchy_ProblemS( Time_Domain, F_Cauchy, U_Cauchy, Scheme )
   
contains

function  F_Cauchy( U, t ) result(F) 
          real ::  U(:), t         
          real :: F(size(U)) 
          
 F = Space_discretization( U, t )
           
end function 

function  Space_discretization( U, t ) result(F) 
          real, target ::  U(:), t, F(size(U))   
  
  real, pointer :: Uv(:,:), Fv(:,:)  
  
  Uv(0:Nx, 0:Ny) => U 
  Fv(0:Nx, 0:Ny) => F 
  
  Fv = Spatial_discretization2D(  Differential_operator,  & 
                                  Boundary_conditions,    & 
                                  x_nodes, y_nodes, Uv, t )
   
end function  

end subroutine

function Spatial_discretization2D( Differential_operator, Boundary_conditions, & 
                                    x, y, U, t) result(F) 

     procedure (DifferentialOperator2D) :: Differential_operator 
     procedure (BC2D) ::  Boundary_conditions
     real, intent(in) ::  x(0:), y(0:), t
     real ::   U(0:, 0:)
     real :: F(0:size(U, dim=1)-1, 0:size(U,dim=2)-1)
  
    real, allocatable :: Ux(:,:), Uy(:,:), Uxx(:,:), Uxy(:,:),  Uyy(:,:)
    integer :: Nx, Ny
    logical ::  NO_BC
    integer :: i, j, k
    type (Boundary2D), save, allocatable :: BC(:)
    integer :: Nb 
    
    logical, save :: dU(5)  ! Dependencies of the differential operator
                            ! (derivatives) , ux ,uy, uxx uyy, uxy
    logical, save :: dBC(2) ! Dependencies of the BC operator
                            ! (derivatives) , ux ,uy 
    
    Nx = size(x)-1; Ny = size(y)-1      
    if (t==0) then 
       dU  = IBVP2D_Dependencies( Differential_operator ) 
       dBC = BC_IBVP2D_Dependencies( Boundary_conditions, & 
                                     x(0), x(Nx), y(0), y(Ny) ) 
       if (allocated(BC)) deallocate(BC)
       allocate( BC( 2*(Ny+1) + 2*(Nx-1) ) )  
    end if 
    
   
    allocate(  Ux(0:Nx, 0:Ny),  Uy(0:Nx, 0:Ny), & 
               Uxx(0:Nx, 0:Ny), Uxy(0:Nx, 0:Ny),  Uyy(0:Nx, 0:Ny) ) 
   
!  ***  It solves boundary  equations to yield values at  boundaries 
        call Boundary_points( dBC, BC, Nx, Ny, x, y, U, t, Boundary_conditions ) 
        
!  ***  It calculates only derivatives which are used 
        call involved_Derivatives( Nx, Ny, dU, U, Ux, Uy, Uxx, Uyy, Uxy ) 
      
!  *** inner grid points
        do i=0, Nx
             do j=0, Ny
                 k = boundary_index(Nx, Ny, i, j) 
                 if (k>0) then 
                    NO_BC = .not. BC(k) % impose
                 else
                    NO_BC = .true.
                 end if 
                 
                 if (NO_BC) then 
                       F(i,j) = Differential_operator( x(i), y(j), t, & 
                                U(i,j), Ux(i,j), Uy(i,j), Uxx(i,j), Uyy(i,j), Uxy(i,j) )
                 else 
                      F(i,j) = 0
                 end if   
             enddo
        enddo

end function 


subroutine Boundary_points( dBC, BC, Nx, Ny, x, y, U, t, Boundary_conditions ) 
        logical, intent(in) :: dBC(:)
        type (Boundary2D) :: BC(:)
        integer, intent(in) :: Nx, Ny 
        real, intent(in) :: x(0:Nx), y(0:Ny), t
        real, intent(inout) :: U(0:Nx, 0:Ny)
        procedure (BC2D) ::  Boundary_conditions
      
    real, allocatable :: Ub(:)
    integer ::  k, m, i, j  
     
  

 if (method == "Fourier") then 
       
           BC % impose = .false.     
 else      
   if ( all( dBC ==.false. ) ) then 
       
       do k=1, size(BC)
            call  ij_index( Nx, Ny, k, i, j) 
            U(i,j) = -Boundary_conditions ( x(i),  y(j),  t, 0.,  0.,  0. )  
        end do 
   else 
            
        call Boundary_unknowns() 
       
        call Newton( BCs, Ub )
        
        m = 1 
        do k=1, size(BC)
            call  ij_index( Nx, Ny, k, i, j) 
            if ( BC(k) % impose ) then 
               U(i,j) =  Ub(m) 
               m = m+ 1 
           endif 
        end do
          
   end if      
 end if 
       
          
contains 


subroutine Boundary_unknowns
   
       integer :: i, j, k
       real :: u0 = 1., ux0= 2., uy0 = 3. 
       integer :: Nb 
 
  do k=1, size(BC)
      
    call  ij_index( Nx, Ny, k, i, j) 
    
    BC(k) % equation = Boundary_conditions ( x(i),  y(j),  t, u0,  ux0,  uy0 ) 
    BC(k) % i = i 
    BC(k) % j = j 
    BC(k) % value = U(i, j) 
    BC(k) % impose = BC(k) % equation /= FREE_BOUNDARY_CONDITION .and. BC(k) % equation /= PERIODIC_BOUNDARY_CONDITION
   
  end do
  
  Nb = count( BC % impose ) 
   
  allocate ( Ub(Nb) )
  Ub = pack( BC(:) % value, BC % impose ) 
  
end subroutine 


function BCs(Z) result(G) 
    real, intent(in) :: Z(:)
    real :: G(size(Z))

 integer :: i, j, k, m  
 real :: Ux(0:Nx,0:Ny), Uy(0:Nx,0:Ny)
  
  m = 1           
  do k=1, size(BC) 
      
      if( BC(k)% impose) then 
            i = BC(k) % i
            j = BC(k) % j
         !   Solution(it,i,j) = Y(m)
            U(i,j) = Z(m)
            m = m + 1 
      end if 
      
  end do 
            
  call Derivative( [ "x", "y" ], 1, 1, U, Ux)
  call Derivative( [ "x", "y" ], 2, 1, U, Uy)

  m = 1 
  do k=1, size(BC) 
       if (BC(k)% impose) then 
           
           i = BC(k) % i 
           j = BC(k) % j 
           G(m) = Boundary_conditions (x(i),  y(j), t, U(i, j),  Ux(i, j),  Uy(i, j) )  
           m = m + 1 
       end if 
  end do

end function 

end subroutine

!-------------------------------------------------------------
! vector of dependencies(derivatives) , ux ,uy, uxx uyy, uxy  
!-------------------------------------------------------------
subroutine involved_Derivatives( Nx, Ny, dU, U, Ux, Uy, Uxx, Uyy, Uxy )
  integer, intent(in) :: Nx, Ny 
  logical, intent(in) :: dU(5) 
  real, intent(in) ::   U(0:Nx, 0:Ny)
  real, intent(out) :: Ux(0:Nx, 0:Ny),  Uy(0:Nx, 0:Ny), & 
                      Uxx(0:Nx, 0:Ny), Uyy(0:Nx, 0:Ny), Uxy(0:Nx, 0:Ny)  
        
        if ( dU(1) ) call Derivative( [ "x", "y" ], 1, 1, U, Ux  )
        if ( dU(2) ) call Derivative( [ "x", "y" ], 2, 1, U, Uy  )
        if ( dU(3) ) call Derivative( [ "x", "y" ], 1, 2, U, Uxx )
        if ( dU(4) ) call Derivative( [ "x", "y" ], 2, 2, U, Uyy )
        
        if ( dU(5).and.dU(1) ) then
            
                      call Derivative( ["x","y"], 2, 1, Ux, Uxy )
                      
        elseif ( dU(5).and.dU(2) ) then
            
                      call Derivative( ["x","y"], 1, 1, Uy, Uxy )
                      
        elseif ( dU(5) ) then
            
                     call Derivative( ["x","y"], 1, 1,  U, Ux  )
                     call Derivative( ["x","y"], 2, 1, Ux, Uxy )
        end if

end subroutine 





!-------------------------------------------------------------
! vector of dependencies(derivatives) , ux ,uy, uxx uyy, uxy  
!-------------------------------------------------------------
subroutine involved_Derivatives_system( Nx, Ny, Nv, dU, U, Ux, Uy, Uxx, Uyy, Uxy )
  integer, intent(in) :: Nx, Ny, Nv  
  logical, intent(in) :: dU(Nv, 5) 
  real, intent(in) ::   U(0:Nx, 0:Ny, Nv)
  real, intent(out) :: Ux(0:Nx, 0:Ny, Nv),  Uy(0:Nx, 0:Ny, Nv), & 
                      Uxx(0:Nx, 0:Ny, Nv), Uyy(0:Nx, 0:Ny, Nv), Uxy(0:Nx, 0:Ny, Nv)  
        
     integer :: k 
  
!  *** inner grid points
        do k=1, Nv 
                    
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
integer function boundary_index(Nx, Ny, i,j) result(k) 
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
    
    
    


































!********************************************************************************************************************************
!*  Initial value boundary problem 2D. Vectorial case
!********************************************************************************************************************************
subroutine IBVP2D_system( Time_Domain, x_nodes, y_nodes, Differential_operator, Boundary_conditions, Solution, Scheme ) 

     real, intent(in) :: Time_Domain(:) 
     real, intent(in) :: x_nodes(0:), y_nodes(0:)
     procedure (DifferentialOperator2D_system) :: Differential_operator 
     procedure (BC2D_system) ::  Boundary_conditions
     procedure (Temporal_Scheme), optional :: Scheme
     real, intent(out) :: Solution(0:, 0:, 0:, :) 
     
     real, pointer :: U_Cauchy(:, :) 
     integer :: Nx, Ny, Nt, Nv
     Nx = size(x_nodes)-1;  Ny = size(y_nodes)-1; Nt = size(Time_Domain)-1
     Nv = size(Solution, dim=4) 
         
    call my_reshape( Solution, Nt+1, (Nx+1)*(Ny+1)*Nv, U_Cauchy ) 
    
    call Cauchy_ProblemS( Time_Domain, F_Cauchy, U_Cauchy, Scheme )
    
contains

function  F_Cauchy( U, t ) result(F) 
          real ::  U(:), t         
          real :: F(size(U)) 
          
 F = Space_discretization( U, t )
 write(*,'(A6, f5.3)') " t = ", t   
           
end function 

function  Space_discretization( U, t ) result(F) 
          real, target ::  U(:), t, F(size(U))   
  
  real, pointer :: Uv(:, :, :), Fv(:, :, :)  
  
  Uv(0:Nx, 0:Ny, 1:Nv) => U 
  Fv(0:Nx, 0:Ny, 1:Nv) => F 
  
  Fv = Spatial_discretization2D_system(  Differential_operator,  & 
                                         Boundary_conditions,    & 
                                         x_nodes, y_nodes, Uv, t )
   
end function  

end subroutine






function Spatial_discretization2D_system(                        & 
                                          Differential_operator, &
                                          Boundary_conditions,   & 
                                           x, y, U, t) result(F) 

     procedure (DifferentialOperator2D_system) :: Differential_operator 
     procedure (BC2D_system) ::  Boundary_conditions
     real, intent(in) ::  x(0:), y(0:), t
     real ::   U(0:, 0:, :)
     real :: F(0:size(U, dim=1)-1, 0:size(U,dim=2)-1, size(U,dim=3))
  
    real, allocatable :: Ux(:,:,:), Uy(:,:,:), Uxx(:,:,:), Uxy(:,:,:),  Uyy(:,:,:)
    real, allocatable :: Fv(:) 
    integer :: Nx, Ny, Nv 
    logical ::  NO_BC
    integer :: i, j, k, l
    type (Boundary2D_system), save, allocatable :: BC(:)
    logical, save, allocatable :: dU(:,:)  !dependencies Operator(variables, derivative)
    logical, save, allocatable :: dBC(:,:) !dependencies BC(variables, derivative)
     
    integer :: Nb 
    
    Nx = size(x)-1; Ny = size(y)-1; Nv = size(U, dim=3)  
    
    if (t==0) then 
       if (allocated(dU)) deallocate(dU, dBC)  
       allocate( dU(Nv,5), dBC(Nv,2)  )
       dU  = IBVP2D_Dependencies_system( Nv, Differential_operator )
       dBC = BC_IBVP2D_Dependencies_system( Nv, Boundary_conditions, & 
                                            x(0), x(Nx), y(0), y(Ny) )   
       if (allocated(BC)) deallocate(BC)
       allocate( BC( 2*(Ny+1) + 2*(Nx-1) ) )
       do k=1, size(BC) 
           allocate( BC(k) % value(Nv), BC(k) % equation(Nv), BC(k) % impose(Nv)  ) 
       end do
    end if 
    
    allocate(  Ux(0:Nx, 0:Ny, Nv),  Uy(0:Nx, 0:Ny, Nv), & 
               Uxx(0:Nx, 0:Ny, Nv), Uxy(0:Nx, 0:Ny, Nv),  Uyy(0:Nx, 0:Ny, Nv) ) 
    allocate( Fv(Nv) ) 
    
!  ***  It solves boundary equations to yield values at boundaries     
        call Boundary_points_system( dBC, BC, Nx, Ny, Nv, x, y, U, t, Boundary_conditions ) 
         
!  ***  It calculates only derivatives which are used 
        call involved_Derivatives_system( Nx, Ny, Nv, dU, U, Ux, Uy, Uxx, Uyy, Uxy )  
   
        F = ZERO
        do i=0, Nx
             do j=0, Ny
                 Fv = Differential_operator( x(i), y(j), t, U(i,j,:), Ux(i,j,:), & 
                                             Uy(i,j,:), Uxx(i,j,:), Uyy(i,j,:), Uxy(i,j,:) ) 
                 k = boundary_index(Nx, Ny, i, j)                   
            ! ** Boundary point when k>0
                 if (k>0) then  
                    do l=1, Nv
                      if ( .not. BC(k) % impose(l) ) then
                              F(i, j, l) = Fv(l) 
                      end if   
                    end do                    
            ! ** inner point    
                 else 
                              F(i, j, :) = Fv(:)
                 end if 
             enddo
        enddo
        
end function 

!------------------------------------------------------------------------------
subroutine Boundary_points_system( dBC, BC, Nx, Ny, Nv, x, y, U, t, Boundary_conditions ) 
        logical, intent(in) :: dBC(:, :) 
        type (Boundary2D_system) :: BC(:)
        integer, intent(in) :: Nx, Ny, Nv 
        real, intent(in) :: x(0:Nx), y(0:Ny), t
        real, intent(inout) :: U(0:Nx, 0:Ny, Nv)
        procedure (BC2D_system) ::  Boundary_conditions
      
    real, allocatable :: Ub(:)
    integer ::  i, j, k, m, l 
    real :: t1, t2
    real :: zero(Nv) 
    
  zero = 0    
  call CPU_time(t1)    
  

 if (method == "Fourier") then 
       do k=1, size(BC) 
           BC(k) % impose = .false.     
      end do      
 else   
   if ( all( dBC ==.false. ) ) then 
       
       do k=1, size(BC)
            call  ij_index( Nx, Ny, k, i, j) 
            U(i,j, :) = -Boundary_conditions ( x(i),  y(j),  t, zero,  zero, zero )  
        end do 
   else 
     
!  ***  It identifies the boundary value unknowns       
        call Boundary_unknowns() 
        
         
 !  *** Solve boundary equations  
        call Newton( BCs, Ub )
        
         m = 1 
         do k=1, size(BC)
            call  ij_index( Nx, Ny, k, i, j) 
            do l=1, Nv 
              if ( BC(k) % impose(l) ) then 
                U(i,j,l) =  Ub(m) 
                m = m + 1 
              end if    
            end do   
         end do
   end if 
   
 end if 
 
 call CPU_time(t2)  
 CPU_time_BC = CPU_time_BC + t2-t1
       
          
contains 
!-----------------------------------------------------------------------
subroutine Boundary_unknowns() 
   
       integer :: i, j, k, l, m, Nb  
       real :: u0(Nv), ux0(Nv), uy0(Nv), eq  
     
  
  u0 = 1.; ux0 = 2.; uy0 = 3  
  Nb = 0 
  do k=1, size(BC)
      
    call  ij_index( Nx, Ny, k, i, j) 
    
    BC(k) % equation(:) = Boundary_conditions (x(i),  y(j),  t, u0,  ux0,  uy0 ) 
    BC(k) % i = i 
    BC(k) % j = j 
    BC(k) % value(:) = U(i, j, :)
    BC(k) % impose = BC(k)% equation /= FREE_BOUNDARY_CONDITION  .and. BC(k) % equation /= PERIODIC_BOUNDARY_CONDITION
    Nb = Nb + count( BC(k)% impose ) 
    
  end do
  
  allocate ( Ub(Nb) ) 
  m = 1 
  do k=1, size(BC)
     do l=1, Nv 
         eq = BC(k)% equation(l)
         if ( eq == PERIODIC_BOUNDARY_CONDITION ) then
             i = BC(k) % i
             j = BC(k) % j
             if (i==0)     U(0, :, l) = U(Nx, :, l)
             if (j==0)     U(:, 0, l) = U(:, Ny, l)
             
         else if (eq /= FREE_BOUNDARY_CONDITION ) then 
             Ub(m) = BC(k)% value(l)
             m = m + 1 
             
         end if 
     end do  
  end do
  
  
end subroutine 

!-----------------------------------------------------------------------
function BCs(Z) result(G) 
        real, intent (in) :: Z(:)
        real :: G(size(Z))
    
       real :: Ux(0:Nx, 0:Ny, Nv), Uy(0:Nx, 0:Ny, Nv)
       integer :: i, j, k, m, l  
       real :: Gv(Nv), eq  
      
  m = 1             
  do k=1, size(BC) 
     do l=1, Nv  
                 if( BC(k)% equation(l) /= FREE_BOUNDARY_CONDITION) then 
                    i = BC(k)% i
                    j = BC(k)% j
                    U(i, j, l) = Z(m)
                    m = m + 1 
                 end if 
     end do 
  end do 
  
  do l=1, Nv           
     call Derivative( [ "x", "y" ], 1, 1, U(:, :, l), Ux(:, :, l) )
     call Derivative( [ "x", "y" ], 2, 1, U(:, :, l), Uy(:, :, l) )
  end do 
 
 
  m = 1 
  do k=1, size(BC) 
       i = BC(k) % i 
       j = BC(k) % j 
       Gv = Boundary_conditions (x(i),  y(j),  t, U(i, j, :),  Ux(i, j, :),  Uy(i, j, :) )  
       
       do l=1, Nv 
                   eq = BC(k)% equation(l)
                   if ( eq /= FREE_BOUNDARY_CONDITION .and. eq /= PERIODIC_BOUNDARY_CONDITION) then 
                      G(m) = Gv(l) 
                      m = m + 1 
                   end if 
       end do 
  end do
      
end function

end subroutine 
























!********************************************************************************************************************************
!*  Initial value boundary problem 2D. Vectorial case
!********************************************************************************************************************************
!subroutine IBVP2D_system_( Time_Domain, x_nodes, y_nodes, Differential_operator, Boundary_conditions, Solution, Scheme ) 
!
!     real, intent(in) :: Time_Domain(:) 
!     real, intent(in) :: x_nodes(0:), y_nodes(0:)
!     procedure (DifferentialOperator2D_system) :: Differential_operator 
!     procedure (BC2D_system) ::  Boundary_conditions
!     procedure (Temporal_Scheme), optional :: Scheme
!     real, intent(out) :: Solution(0:, 0:, 0:, :) 
!   
!     
!     real, pointer :: U_Cauchy(:, :) 
!     type (Boundary2D_system), allocatable :: BC(:)
!     integer :: Nb 
!     logical, allocatable :: dU(:,:) ! dependencies Operator(variables, derivative)
!     logical, allocatable :: dBC(:,:) !dependencies BC(variables, derivative)
!     
!     integer :: Nx, Ny, Nt, Nv
!     Nx = size(x_nodes) - 1;  Ny = size(y_nodes) - 1; Nt = size(Time_Domain) - 1
!     Nv = size(Solution, dim=4) 
!     allocate( dU(Nv,5), dBC(Nv,2)  )
! 
!    dU  = IBVP2D_Dependencies_system( Nv, Differential_operator )
!    dBC = BC_IBVP2D_Dependencies_system( Nv, Boundary_conditions, & 
!                                         x_nodes(0), x_nodes(Nx), y_nodes(0), y_nodes(Ny) ) 
!         
!    call my_reshape( Solution, Nt+1, (Nx+1)*(Ny+1)*Nv, U_Cauchy ) 
!    
!    call Cauchy_ProblemS( Time_Domain, Space_discretization, U_Cauchy, Scheme )
!    
!contains
!
!!----------------------------------------------------------------------
!function  Space_discretization( U, t ) result(F) 
!          real ::  U(:), t         
!          real :: F(size(U))   
!
!         call Space_discretization_2D_system( U, t, F )
!         write(*,'(A6, f5.3)') " t = ", t   
!         
!end function    
!!-----------------------------------------------------------------------
!subroutine  Space_discretization_2D_system( U, t, F )
!          real :: U(0:Nx,0:Ny, Nv), t 
!          real :: F(0:Nx,0:Ny, Nv)
!              
!    integer :: i, j, k, l
!    real :: Ux(0:Nx,0:Ny, Nv), Uxx(0:Nx,0:Ny, Nv), Uxy(0:Nx,0:Ny, Nv)
!    real :: Uy(0:Nx,0:Ny, Nv), Uyy(0:Nx,0:Ny, Nv)
!    real :: Fv(Nv)
!    
!!  ***  It solves boundary equations to yield values at boundaries     
!        call Boundary_points_system_( dBC, BC, Nx, Ny, Nv, x_nodes, y_nodes, U, t, Boundary_conditions ) 
!         
!!  ***  It calculates only derivatives which are used 
!        call involved_Derivatives_system( Nx, Ny, Nv, dU, U, Ux, Uy, Uxx, Uyy, Uxy )  
!   
!        F = ZERO
!        do i=0, Nx
!             do j=0, Ny
!                 Fv = Differential_operator( x_nodes(i), y_nodes(j), t, U(i,j,:), Ux(i,j,:), & 
!                                             Uy(i,j,:), Uxx(i,j,:), Uyy(i,j,:), Uxy(i,j,:) ) 
!                 k = boundary_index(Nx, Ny, i, j)                   
!            ! ** Boundary point when k>0
!                 if (k>0) then  
!                    do l=1, Nv
!                      if ( .not. BC(k) % impose(l) ) then
!                              F(i, j, l) = Fv(l) 
!                      end if   
!                    end do                    
!            ! ** inner point    
!                 else 
!                              F(i, j, :) = Fv(:)
!                 end if 
!             enddo
!        enddo
!        
!end subroutine
!end subroutine
!
!!------------------------------------------------------------------------------
!subroutine Boundary_points_system_( dBC, BC, Nx, Ny, Nv, x, y, U, t, Boundary_conditions ) 
!        logical, intent(in) :: dBC(:, :) 
!        type (Boundary2D_system), allocatable :: BC(:)
!        integer, intent(in) :: Nx, Ny, Nv 
!        real, intent(in) :: x(0:Nx), y(0:Ny), t
!        real, intent(inout) :: U(0:Nx, 0:Ny, Nv)
!        procedure (BC2D_system) ::  Boundary_conditions
!      
!    real, allocatable :: Ub(:)
!    integer ::  i, j, k, m, l 
!    real :: t1, t2
!    real :: zero(Nv) 
!    
!  zero = 0    
!  call CPU_time(t1)    
!  if (.not.allocated(BC)) then 
!         allocate( BC( 2*(Ny+1) + 2*(Nx-1) ) ) 
!         do k=1, size(BC) 
!           allocate( BC(k) % value(Nv), BC(k) % equation(Nv), BC(k) % impose(Nv)  ) 
!         end do
!  end if 
!
! if (method == "Fourier") then 
!       do k=1, size(BC) 
!           BC(k) % impose = .false.     
!      end do      
! else   
!   if ( all( dBC ==.false. ) ) then 
!       
!       do k=1, size(BC)
!            call  ij_index( Nx, Ny, k, i, j) 
!            U(i,j, :) = -Boundary_conditions ( x(i),  y(j),  t, zero,  zero, zero )  
!        end do 
!   else 
!     
!!  ***  It identifies the boundary value unknowns       
!        call Boundary_unknowns() 
!        
!         
! !  *** Solve boundary equations  
!        call Newton( BCs, Ub )
!        
!         m = 1 
!         do k=1, size(BC)
!            call  ij_index( Nx, Ny, k, i, j) 
!            do l=1, Nv 
!              if ( BC(k) % impose(l) ) then 
!                U(i,j,l) =  Ub(m) 
!                m = m + 1 
!              end if    
!            end do   
!         end do
!   end if 
!   
! end if 
! 
! call CPU_time(t2)  
! CPU_time_BC = CPU_time_BC + t2-t1
!       
!          
!contains 
!!-----------------------------------------------------------------------
!subroutine Boundary_unknowns() 
!   
!       integer :: i, j, k, l, m, Nb  
!       real :: u0(Nv), ux0(Nv), uy0(Nv), eq  
!     
!  
!  u0 = 1.; ux0 = 2.; uy0 = 3  
!  Nb = 0 
!  do k=1, size(BC)
!      
!    call  ij_index( Nx, Ny, k, i, j) 
!    
!    BC(k) % equation(:) = Boundary_conditions (x(i),  y(j),  t, u0,  ux0,  uy0 ) 
!    BC(k) % i = i 
!    BC(k) % j = j 
!    BC(k) % value(:) = U(i, j, :)
!    BC(k) % impose = BC(k)% equation /= FREE_BOUNDARY_CONDITION  .and. BC(k) % equation /= PERIODIC_BOUNDARY_CONDITION
!    Nb = Nb + count( BC(k)% impose ) 
!    
!  end do
!  
!  allocate ( Ub(Nb) ) 
!  m = 1 
!  do k=1, size(BC)
!     do l=1, Nv 
!         eq = BC(k)% equation(l)
!         if ( eq == PERIODIC_BOUNDARY_CONDITION ) then
!             i = BC(k) % i
!             j = BC(k) % j
!             if (i==0)     U(0, :, l) = U(Nx, :, l)
!             if (j==0)     U(:, 0, l) = U(:, Ny, l)
!             
!         else if (eq /= FREE_BOUNDARY_CONDITION ) then 
!             Ub(m) = BC(k)% value(l)
!             m = m + 1 
!             
!         end if 
!     end do  
!  end do
!  
!  
!end subroutine 
!
!!-----------------------------------------------------------------------
!function BCs(Z) result(G) 
!        real, intent (in) :: Z(:)
!        real :: G(size(Z))
!    
!       real :: Ux(0:Nx, 0:Ny, Nv), Uy(0:Nx, 0:Ny, Nv)
!       integer :: i, j, k, m, l  
!       real :: Gv(Nv), eq  
!      
!  m = 1             
!  do k=1, size(BC) 
!     do l=1, Nv  
!                 if( BC(k)% equation(l) /= FREE_BOUNDARY_CONDITION) then 
!                    i = BC(k)% i
!                    j = BC(k)% j
!                    U(i, j, l) = Z(m)
!                    m = m + 1 
!                 end if 
!     end do 
!  end do 
!  
!  do l=1, Nv           
!     call Derivative( [ "x", "y" ], 1, 1, U(:, :, l), Ux(:, :, l) )
!     call Derivative( [ "x", "y" ], 2, 1, U(:, :, l), Uy(:, :, l) )
!  end do 
! 
! 
!  m = 1 
!  do k=1, size(BC) 
!       i = BC(k) % i 
!       j = BC(k) % j 
!       Gv = Boundary_conditions (x(i),  y(j),  t, U(i, j, :),  Ux(i, j, :),  Uy(i, j, :) )  
!       
!       do l=1, Nv 
!                   eq = BC(k)% equation(l)
!                   if ( eq /= FREE_BOUNDARY_CONDITION .and. eq /= PERIODIC_BOUNDARY_CONDITION) then 
!                      G(m) = Gv(l) 
!                      m = m + 1 
!                   end if 
!       end do 
!  end do
!      
!end function
!
!end subroutine 
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!!*******************************************************************
!!* Initial value boundary problem 2D. Escalar case 
!!*******************************************************************
!subroutine IBVP2D_old( Time_Domain, x_nodes, y_nodes, Differential_operator,  Boundary_conditions, Solution, Scheme ) 
!
!     real, intent(in) :: Time_Domain(:) 
!     real, intent(inout) :: x_nodes(0:), y_nodes(0:)
!     procedure (DifferentialOperator2D) :: Differential_operator 
!     procedure (BC2D) ::  Boundary_conditions  
!     real, intent(out) :: Solution(0:, 0:, 0:) 
!     procedure (Temporal_Scheme), optional :: Scheme
!   
!
!     integer :: Nt, Nx, Ny, M 
!     real, pointer :: U_Cauchy(:, :) 
!     real :: t_BC 
!     integer ::  it  
!     type (Boundary2D), allocatable :: BC(:)
!     integer :: Nb 
!     logical :: dU(5)  ! Dependencies of the differential operator
!                       ! (derivatives) , ux ,uy, uxx uyy, uxy
!     logical :: dBC(2) ! Dependencies of the BC operator
!                       ! (derivatives) , ux ,uy
!     
!     Nt = size(Time_Domain) -1;  Nx = size(x_nodes) - 1; Ny = size(y_nodes) - 1
!     M = (Nx+1) * (Ny+1)   
!   
!     dU  = IBVP2D_Dependencies( Differential_operator ) 
!     dBC = BC_IBVP2D_Dependencies( Boundary_conditions, & 
!                                   x_nodes(0), x_nodes(Nx), y_nodes(0), y_nodes(Ny) ) 
!          
!     call my_reshape( Solution, Nt+1, M, U_Cauchy ) 
!    
!     call Cauchy_ProblemS( Time_Domain, Space_discretization, U_Cauchy, Scheme )
!   
!contains
!!----------------------------------------------------------------------
!function  Space_discretization( Uc, t ) result(F) 
!          real ::  Uc(:), t         
!          real :: F(size(Uc))    
!
!         call Space_discretization2D( Uc, t, F) 
!         write(*,*) "IBVP2D Space discretization t = ", t 
!        
!         
!end function 
!
!subroutine Space_discretization2D( U, t, F )
!          real :: U(0:Nx, 0:Ny), t, F(0:Nx, 0:Ny) 
!    
!    real :: Ux(0:Nx, 0:Ny),  Uy(0:Nx, 0:Ny), & 
!           Uxx(0:Nx, 0:Ny), Uxy(0:Nx, 0:Ny),  Uyy(0:Nx, 0:Ny)
!    logical ::  NO_BC
!    integer :: i, j, k
!   
!!  ***  It solves boundary  equations to yield values at  boundaries 
!        call Boundary_points_old( dBC, BC, Nx, Ny, x_nodes, y_nodes, U, t, Boundary_conditions ) 
!        
!!  ***  It calculates only derivatives which are used 
!        call involved_Derivatives_old( Nx, Ny, dU, U, Ux, Uy, Uxx, Uyy, Uxy ) 
!      
!!  *** inner grid points
!       ! F = ZERO
!        do i=0, Nx
!             do j=0, Ny
!                 
!                 k = boundary_index(Nx, Ny, i, j) 
!                 if (k>0) then 
!                    NO_BC = .not. BC(k) % impose
!                 else
!                    NO_BC = .true.
!                 end if 
!                 
!                 if (NO_BC) then 
!                       F(i,j) = Differential_operator( x_nodes(i), y_nodes(j), t, & 
!                                U(i,j), Ux(i,j), Uy(i,j), Uxx(i,j), Uyy(i,j), Uxy(i,j) )
!                 else 
!                      F(i,j) = 0
!                 end if   
!             enddo
!        enddo
!
!end subroutine 
!
!end subroutine 
!
!!------------------------------------------------------------------------------
!subroutine Boundary_points_old( dBC, BC, Nx, Ny, x, y, U, t, Boundary_conditions ) 
!        logical, intent(in) :: dBC(:)
!        type (Boundary2D), allocatable :: BC(:)
!        integer, intent(in) :: Nx, Ny 
!        real, intent(in) :: x(0:Nx), y(0:Ny), t
!        real, intent(inout) :: U(0:Nx, 0:Ny)
!        procedure (BC2D) ::  Boundary_conditions
!      
!    real, allocatable :: Ub(:)
!    integer ::  k, m, i, j  
!     
!  if (.not.allocated(BC)) then 
!         allocate( BC( 2*(Ny+1) + 2*(Nx-1) ) )  
!  end if 
!
! if (method == "Fourier") then 
!       
!           BC % impose = .false.     
! else      
!   if ( all( dBC ==.false. ) ) then 
!       
!       do k=1, size(BC)
!            call  ij_index( Nx, Ny, k, i, j) 
!            U(i,j) = -Boundary_conditions ( x(i),  y(j),  t, 0.,  0.,  0. )  
!        end do 
!   else 
!            
!            
!        call Boundary_unknowns() 
!       
!        call Newton( BCs, Ub )
!        
!        m = 1 
!        do k=1, size(BC)
!            call  ij_index( Nx, Ny, k, i, j) 
!            if ( BC(k) % impose ) then 
!               U(i,j) =  Ub(m) 
!               m = m+ 1 
!           endif 
!        end do
!          
!   end if      
! end if 
!       
!          
!contains 
!
!!-----------------------------------------------------------------------
!subroutine Boundary_unknowns
!   
!       integer :: i, j, k
!       real :: u0 = 1., ux0= 2., uy0 = 3. 
!       integer :: Nb 
! 
!  do k=1, size(BC)
!      
!    call  ij_index( Nx, Ny, k, i, j) 
!    
!    BC(k) % equation = Boundary_conditions ( x(i),  y(j),  t, u0,  ux0,  uy0 ) 
!    BC(k) % i = i 
!    BC(k) % j = j 
!    BC(k) % value = U(i, j) 
!    BC(k) % impose = BC(k) % equation /= FREE_BOUNDARY_CONDITION .and. BC(k) % equation /= PERIODIC_BOUNDARY_CONDITION
!   
!  end do
!  
!  Nb = count( BC % impose ) 
!   
!  allocate ( Ub(Nb) )
!  Ub = pack( BC(:) % value, BC % impose ) 
!  
!end subroutine 
!
!!-----------------------------------------------------------------------
!function BCs(Z) result(G) 
!    real, intent(in) :: Z(:)
!    real :: G(size(Z))
!
! integer :: i, j, k, m  
! real :: Ux(0:Nx,0:Ny), Uy(0:Nx,0:Ny)
!  
!  m = 1           
!  do k=1, size(BC) 
!      
!      if( BC(k)% impose) then 
!            i = BC(k) % i
!            j = BC(k) % j
!         !   Solution(it,i,j) = Y(m)
!            U(i,j) = Z(m)
!            m = m + 1 
!      end if 
!      
!  end do 
!            
!  call Derivative( [ "x", "y" ], 1, 1, U, Ux)
!  call Derivative( [ "x", "y" ], 2, 1, U, Uy)
!
!  m = 1 
!  do k=1, size(BC) 
!       if (BC(k)% impose) then 
!           
!           i = BC(k) % i 
!           j = BC(k) % j 
!           G(m)   = Boundary_conditions (x(i),  y(j), t, U(i, j),  Ux(i, j),  Uy(i, j) )  
!           m = m + 1 
!           
!       end if 
!       
!  end do
!
!end function 
!
!end subroutine
!
!!-------------------------------------------------------------
!! vector of dependencies(derivatives) , ux ,uy, uxx uyy, uxy  
!!-------------------------------------------------------------
!subroutine involved_Derivatives_old( Nx, Ny, dU, U, Ux, Uy, Uxx, Uyy, Uxy )
!  integer, intent(in) :: Nx, Ny 
!  logical, intent(in) :: dU(5) 
!  real, intent(in) ::   U(0:Nx, 0:Ny)
!  real, intent(out) :: Ux(0:Nx, 0:Ny),  Uy(0:Nx, 0:Ny), & 
!                      Uxx(0:Nx, 0:Ny), Uyy(0:Nx, 0:Ny), Uxy(0:Nx, 0:Ny)  
!        
!        if ( dU(1) ) call Derivative( [ "x", "y" ], 1, 1, U, Ux  )
!        if ( dU(2) ) call Derivative( [ "x", "y" ], 2, 1, U, Uy  )
!        if ( dU(3) ) call Derivative( [ "x", "y" ], 1, 2, U, Uxx )
!        if ( dU(4) ) call Derivative( [ "x", "y" ], 2, 2, U, Uyy )
!        
!        if ( ( dU(5) ).and.( dU(1) ) ) then
!            
!                      call Derivative( ["x","y"], 2, 1, Ux, Uxy )
!                      
!        elseif ( ( dU(5) ).and.( dU(2) ) ) then
!            
!                      call Derivative( ["x","y"], 1, 1, Uy, Uxy )
!                      
!        elseif ( dU(5) ) then
!            
!                     call Derivative( ["x","y"], 1, 1,  U, Ux  )
!                     call Derivative( ["x","y"], 2, 1, Ux, Uxy )
!        end if
!
!end subroutine 





end module 
    
    
 
  