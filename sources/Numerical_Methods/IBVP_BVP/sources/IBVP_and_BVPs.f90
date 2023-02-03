module IBVPs_and_BVPs

use Cauchy_Problem
use Temporal_scheme_interface
use Collocation_methods 
use Linear_Systems
use Non_Linear_Systems
use Boundary_value_problems
use Utilities

implicit none   
private
public :: IBVP_and_BVP
   
abstract interface  
      
 function L_uv(       x, y, t, u,    ux,    uy,    uxx,    uyy,    uxy, &
                               v,    vx,    vy,    vxx,    vyy,    vxy  ) 
  real, intent(in) :: x, y, t, u(:), ux(:), uy(:), uxx(:), uyy(:), uxy(:)
  real, intent(in) ::          v(:), vx(:), vy(:), vxx(:), vyy(:), vxy(:)
  real :: L_uv (size(u)) 
 end function  

 function BC2DS(x, y, t, u, ux, uy) 
       real, intent(in) :: x, y, t, u(:), ux(:), uy(:) 
       real :: BC2DS(size(u)) ! maximum number of BCs at each point
 end function  
 
subroutine B_data(u, ux, uy, uxx, uyy) 
       real, intent(in) :: u(:,:,:), ux(:,:,:),  uy(:,:,:), uxx(:,:,:), uyy(:,:,:) 
end subroutine  
 
 
end interface


contains

    
subroutine IBVP_and_BVP(  Time, x, y, L_u, L_v, BC_u, BC_v, &
                          Ut, Vt, Scheme, Boundary_data ) 

   real, intent(in) :: Time(0:), x(0:), y(0:)
   procedure (L_uv) ::   L_u, L_v
   procedure (BC2DS) :: BC_u, BC_v
   real, intent(out) :: Ut(0:, 0:, 0:, :),  Vt(0:, 0:, 0:, :)
   procedure (Temporal_Scheme), optional :: Scheme
   procedure (B_data), optional :: Boundary_data   
  
   real, pointer :: U_Cauchy(:, :) 
   real :: t_BC 
   integer ::  it, Nx, Ny, Nt,  Nu, Nv, M1, M2, M3, M4 
   real, allocatable :: Ux(:,:,:), Uxx(:,:,:), Uxy(:,:,:),            & 
                        Uy(:,:,:), Uyy(:,:,:) 
   
   Nx = size(x) - 1 ; Ny = size(y) - 1; Nt = size(Time) - 1
   Nu = size(Ut, dim=4); Nv = size(Vt, dim=4) 
   allocate( Ux(0:Nx,0:Ny, Nu), Uxx(0:Nx,0:Ny, Nu), Uxy(0:Nx,0:Ny, Nu), & 
             Uy(0:Nx,0:Ny, Nu), Uyy(0:Nx,0:Ny, Nu) )        
   M1 = Nu*(Ny-1); M2 = Nu*(Ny-1); M3 = Nu*(Nx-1); M4 = Nu*(Nx-1) 
   
   
   call my_reshape( Ut, Nt+1, Nu*(Nx+1)*(Ny+1), U_Cauchy ) 
   
   call Cauchy_ProblemS( Time, BVP_and_IBVP_discretization, U_Cauchy, Scheme )
   
 
     
contains

!----------------------------------------------------------------------
function BVP_and_IBVP_discretization( U, t ) result(F) 
      real ::  U(:), t         
      real :: F(size(U))   

         call BVP_and_IBVP_discretization_2D( U, t, F )          
         
end function    

!-----------------------------------------------------------------------
subroutine BVP_and_IBVP_discretization_2D( U, t, F_u )
          real :: U(0:Nx,0:Ny, Nu), t, F_u(0:Nx,0:Ny, Nu) 
          
    integer :: i, j, k 
    real :: Vx(0:Nx,0:Ny, Nv), Vxx(0:Nx,0:Ny, Nv), Vxy(0:Nx,0:Ny, Nv)
    real :: Vy(0:Nx,0:Ny, Nv), Vyy(0:Nx,0:Ny, Nv), Uc(M1+M2+M3+M4) 
    
        t_BC = t;  call Binary_search(t_BC, Time, it)
 !  *** Step 1. BVP for V        
       call Boundary_Value_Problem( x, y, L_v_R, BC_v_R, Vt(it,0:,0:,:))
       
!  *** Derivatives for V
       do k=1, Nv 
        call Derivative( ["x","y"], 1, 1, Vt(it, 0:,0:, k), Vx (0:,0:,k) )
        call Derivative( ["x","y"], 2, 1, Vt(it, 0:,0:, k), Vy (0:,0:,k) )
        call Derivative( ["x","y"], 1, 2, Vt(it, 0:,0:, k), Vxx(0:,0:,k) )
        call Derivative( ["x","y"], 2, 2, Vt(it, 0:,0:, k), Vyy(0:,0:,k) )
        call Derivative( ["x","y"], 2, 1, Vx(0:    ,0:, k), Vxy(0:,0:,k) )
       end do   
       if (present(Boundary_data))  then 
        call Boundary_data( Vt(it, 0:,0:, :), Vx, Vy, Vxx, Vyy)  
       end if        
        
 !  *** initial boundary value :  Uc 
        call Asign_BV2s( U(0, 1:Ny-1, 1:Nu ), U( Nx, 1:Ny-1, 1:Nu ),  & 
                         U(1:Nx-1, 0, 1:Nu ), U( 1:Nx-1, Ny, 1:Nu ), Uc) 
 !  *** Step 2. Boundary points Uc from inner points U 
        call Newton( BCs, Uc )
!   *** asign boundary points Uc  to U        
        call Asign_BVs(Uc, U( 0, 1:Ny-1, 1:Nu ), U( Nx, 1:Ny-1, 1:Nu  ), &
                           U( 1:Nx-1, 0, 1:Nu ), U( 1:Nx-1, Ny, 1:Nu )  )         
!  *** Derivatives of U for inner grid points
       do k=1, Nu 
          call Derivative( ["x","y"], 1, 1, U(0:,0:, k),  Ux (0:,0:,k) )
          call Derivative( ["x","y"], 1, 2, U(0:,0:, k),  Uxx(0:,0:,k) )
          call Derivative( ["x","y"], 2, 1, U(0:,0:, k),  Uy (0:,0:,k) )
          call Derivative( ["x","y"], 2, 2, U(0:,0:, k),  Uyy(0:,0:,k) )
          call Derivative( ["x","y"], 2, 1, Ux(0:,0:,k),  Uxy(0:,0:,k) )
       end do  
!  *** Step 3. Differential operator L_u(U,V) at inner grid points
 F_u=0   
 do i=1, Nx-1; do j=1, Ny-1
    F_u(i, j, :) = L_u(  x(i), y(j), t, U(i, j, :),                   &
                         Ux(i, j, :), Uy(i, j, :), Uxx(i, j, :),      &
                         Uyy(i, j, :), Uxy(i, j, :), Vt(it, i, j, :), &
                         Vx(i, j, :), Vy(i, j, :), Vxx(i, j, :),      & 
                         Vyy(i, j, :), Vxy(i, j, :)  )
    
 end  do; end do 
end subroutine


!-----------------------------------------------------------------------
function BCs(Y) result(G) 
        real, intent (in) :: Y(:)
        real :: G(size(Y))
    
     real :: G1(M1), G2(M2), G3(M3), G4(M4) 

! ** Asign Newton's iteration Y to Solution     
     call Asign_BVs(Y, Ut(it, 0 ,1:Ny-1, 1:Nu), Ut(it, Nx, 1:Ny-1, 1:Nu),&
                       Ut(it, 1:Nx-1, 0, 1:Nu), Ut(it, 1:Nx-1, Ny, 1:Nu) ) 
! ** Calculate boundary conditions G      
     call Asign_BCs( G1,  G2, G3,  G4 ) 
  
     G = [ G1, G2, G3, G4 ] 

end function

!------------------------------------------------------------------------------------
subroutine Asign_BV2s(   U1,  U2,  U3,  U4, Uc  )
   real, intent(in) :: U1(M1), U2(M2), U3(M3), U4(M4) 
   real, intent(out) :: Uc(M1+M2+M3+M4)
   
     Uc = [ U1, U2, U3, U4 ] 
   
end subroutine

!---------------------------------------------------------------------------------
subroutine Asign_BVs( Y, U1,  U2,  U3, U4) 
  real, intent(in) :: Y(M1+M2+M3+M4)
  real, intent(out) :: U1(M1), U2(M2), U3(M3), U4(M4)  
   
    integer :: i1, i2, i3, i4 
  
     i1 = 1  + M1; i2 = i1 + M2; i3 = i2 + M3; i4 = i3 + M4
     
     U1 = Y(1  : i1-1)
     U2 = Y(i1 : i2-1) 
     U3 = Y(i2 : i3-1) 
     U4 = Y(i3 : i4-1) 
end subroutine


subroutine Asign_BCs(  G1, G2,  G3,  G4  ) 
   real, intent(out) :: G1(1:Ny-1,Nu), G2(1:Ny-1,Nu),                  &
                        G3(1:Nx-1,Nu), G4(1:Nx-1,Nu)  
  
   real :: Wx(0:Nx, 0:Ny, Nu), Wy(0:Nx, 0:Ny, Nu)
   integer :: i, j, k 
    
   do k=1, Nu 
    call Derivative( ["x","y"], 1, 1, Ut(it, 0:, 0:, k), Wx(0:, 0:, k) )
    call Derivative( ["x","y"], 2, 1, Ut(it, 0:, 0:, k), Wy(0:, 0:, k) )
   end do 
        
   do j = 1, Ny-1
     G1(j,:) = BC_u( x(0), y(j), t_BC,   &
                     Ut(it, 0, j, : ), Wx(0, j,:), Wy(0, j,:)) 
     G2(j,:) = BC_u( x(Nx), y(j), t_BC,   &
                     Ut(it, Nx, j, : ), Wx(Nx, j, :), Wy(Nx, j, :)) 
   end do
     
   do i = 1, Nx-1
     G3(i,:) = BC_u( x(i), y(0), t_BC,   &
                     Ut(it, i, 0,:), Wx(i, 0, :), Wy(i, 0,:)) 
     G4(i,:) = BC_u( x(i), y(Ny), t_BC,  &
                     Ut(it, i, Ny, : ), Wx(i, Ny, :), Wy(i, Ny, :)) 
   end do
   
end subroutine


!-----------------------------------------------------------------------
function L_v_R(xr, yr, V, Vx, Vy, Vxx, Vyy, Vxy) result(Fv)
     real, intent(in) :: xr, yr, V(:), Vx(:), Vy(:), Vxx(:), Vyy(:), Vxy(:)
     real :: Fv(size(V))      
 
     integer :: ix, iy 
     
     call Binary_search(xr, x, ix)  
     call Binary_search(yr, y, iy)
     
     Fv = L_v( xr, yr, t_BC, V(:),  Vx(:),  Vy(:),  &
               Vxx(:), Vyy(:), Vxy(:),              & 
               Ut(it, ix, iy, :),  Ux(ix, iy, :),   &
               Uy(ix, iy, :), Uxx(ix, iy, :),       &
               Uyy(ix, iy, :), Uxy(ix, iy, :)    ) 
      
end function


!-----------------------------------------------------------------------
function BC_v_R(xr, yr, V, Vx, Vy) result (G_v)
     real, intent(in) :: xr, yr, V(:), Vx(:), Vy(:)
     real :: G_v(size(V))      
 
    G_v = BC_v( xr, yr, t_BC, V, Vx, Vy )
    
        
end function


end subroutine   

end module 
    
    
  
