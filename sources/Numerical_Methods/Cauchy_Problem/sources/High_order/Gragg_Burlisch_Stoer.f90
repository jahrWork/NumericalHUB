module Gragg_Burlisch_Stoer 
    
    use ODE_Interface
    use Stability
    
    implicit none
   
    real, save :: Tolerance = 1e-8  
    integer, save :: NL_fixed = 0  
    integer, save :: N_GBS_effort = 0 
    
    contains     
    
subroutine set_GBS_tolerance( eps)
      real, intent(in) :: eps
     
      Tolerance = eps
      NL_fixed  = 0 
           
end subroutine

subroutine set_GBS_levels( NL )
    integer, intent(in) :: NL
     
      NL_fixed = NL 
      
end subroutine


integer function get_GBS_effort() result(N) 
     
      N = N_GBS_effort
      
end function 


!*******************************************************************************
!    Gragg-Bulirsch-Stoer (GBS) scheme
!    2019. Francisco Javier Escoto Lopez, javier.escoto.lopez@gmail.com 
!    2020. Juan A. Hernandez Ramos, juanantonio.hernandez@upm.es 
!    2021. GBS improvement. Variable step (jahr)
!******************************************************************************* 
subroutine GBS_Scheme( F, t1, t2, U1, U2, ierr) 
      procedure (ODES) ::   F 
      real, intent(in) ::   t1, t2, U1(:)
      real, intent(out) ::  U2(:)
      integer, intent(out) :: ierr
    
     real  ::  dt, t1s, t2s, U1s( size(U1) )            
     complex :: lambda( size(U1) ) 
     integer :: i, N_steps
    
     
     if (NL_fixed>0) then
         
         call GBS_solution_NL( F, t1, t2, U1, U2, NL_fixed) 
     else 
         
         N_steps = 0 
         U1s = U1; t1s = t1; t2s = t1   
         
         do while( t2s < t2 ) 
              lambda = Eigenvalues_Jacobian( F, U1s, t1s)
              dt  = 0.1 / maxval( abs(lambda) ) 
              
              if (t1s + dt > t2) then 
                  t2s = t2
             else
                  t2s = t1s + dt 
              end if 
              
              call GBS_Solution( F, t1s, t2s, U1s, U2, ierr) 
              U1s = U2; t1s = t2s  
              N_steps = N_steps + 1 
         end do 
       ! write(*,*) " N_steps = ", N_steps 
            
     end if 
     
     ierr = 0 
  
end subroutine


!*******************************************************************************
!    Gragg-Bulirsch-Stoer (GBS) scheme with NL constant number of levels 
!    2019. Francisco Javier Escoto Lopez, javier.escoto.lopez@gmail.com 
!    2020. Juan A. Hernandez Ramos, juanantonio.hernandez@upm.es 
!******************************************************************************* 
subroutine GBS_solution_NL( F, t1, t2, U1, U2, NL) 
      procedure (ODES) ::   F 
      real, intent(in) ::   t1, t2, U1(:)
      real, intent(out) ::  U2(:)
      integer, intent(in) :: NL
      
     real, allocatable :: U(:,:) 
     integer, allocatable :: N(:)
     integer :: i, Nv     
   
     Nv = size(U1) 
     if (NL<1) then 
         write(*,*) " ERROR: NL must be greater or equal than 1, NL =", NL 
         stop 
     else 
         allocate( U(NL, Nv), N(NL) )
     end if   
          
   ! *** Partition sequence definition
         N = mesh_refinement(NL) 
       
   ! *** Modified midpoint scheme for each level
         do i = 1, NL
           call Modified_midpoint_scheme( F, t1, t2, U1, U(i,:), N(i) )
         end do
   
   ! *** Richardson extraplation       
         U2 = Corrected_solution_Richardson( N, U ) 
    
  
end subroutine


subroutine GBS_Solution( F, t1, t2, U1, U2, ierr) 
      procedure (ODES) ::   F 
      real, intent(in) ::   t1, t2, U1(:)
      real, intent(out) ::  U2(:)
      integer, intent(out) :: ierr
    
     integer :: NL                   ! number of levels of mesh refinement  
     integer, parameter :: NLmax = 8 ! max number of allowed levels 
     real  :: Ucs( size(U1) )        ! corrected solution with NL+1 levels 
     real  :: Uc( size(U1) )         ! corrected solution with NL levels
     real  :: UL( NLmax+1,  size(U1))! solutions from l=1.. NLmax levels 
     real  :: Error                  ! estimated error from solutions 
                                     ! with NL and NL+1 levels 
     
     complex :: lambda( size(U1) ) 
     integer :: i 
     logical :: next_level 
     
         !write(*,*) "GBS_solution"
         NL = 1;  Error = 1 
         do while (Error > Tolerance .and. NL < NLmax ) 
              NL = NL + 1 
              call GBS_solutionL( F, t1, t2, U1, UL, Uc, Ucs, NL) 
              Error = norm2( Uc - Ucs ) 
              !write(*,*) " Error = ", Error
              !write(*,'(a10, e13.7,a10, e13.7,i3)') " Error =", Error, "Tolerance =", Tolerance, NL
              
         end do  
         U2 = Ucs    
     
     ierr = 0 
  
end subroutine


!*******************************************************************************
! It calculated two solutions U2 with NL levels  and Us with NL+1 levels
! with Gragg-Bulirsch-Stoer (GBS) scheme. 
! 2019. Francisco Javier Escoto Lopez, javier.escoto.lopez@gmail.com 
! 2020. Juan A. Hernandez Ramos, juanantonio.hernandez@upm.es 
!******************************************************************************* 
subroutine GBS_solutionL( F, t1, t2, U1, UL, Uc, Ucs, NL) 
      procedure (ODES) ::   F 
      real, intent(in) ::   t1, t2, U1(:)
      real, intent(inout) ::  UL(:, :)
      real, intent(out) :: Uc(:), Ucs(:) 
      integer, intent(in) :: NL
     
     integer :: i, N(NL+1)
   
 ! *** Partition sequence definition
       N = mesh_refinement(NL+1) 
        
 ! *** Modified midpoint scheme for each level
       if (NL==2) then  
         call Modified_midpoint_scheme( F, t1, t2, U1, UL(1,:), N(1) )
         call Modified_midpoint_scheme( F, t1, t2, U1, UL(2,:), N(2) )
       endif    
         
      call Modified_midpoint_scheme( F, t1, t2, U1, UL(NL+1,:), N(NL+1) )
         
! *** Corrected solution with NL levels NL and with NL+1 levels 
      Uc  = Corrected_solution_Richardson( N(1:NL),   UL(1:NL,:)   ) 
      Ucs = Corrected_solution_Richardson( N(1:NL+1), UL(1:NL+1,:) ) 
     
  
end subroutine


function mesh_refinement(Levels) result(N) 
       integer, intent(in) :: Levels 
       integer :: N(Levels) 
       
 integer :: N_Romberg(10) = [ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 ] 
 integer :: N_Burlirsch(10) = [ 1, 2, 3, 4, 6, 8, 12, 16, 24, 32 ] 
 integer :: N_Harmonic(10) = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 ]
 
 !N = N_Romberg(1:Levels)
 !N = N_Burlirsch(1:Levels) 
  N = N_Harmonic(1:Levels) 
       
end function        
 
  
!*******************************************************************************
!    Modified Midpoint scheme
!    given an interval [t1, t2] divided in 2n = number of steps in the interval 
!    and an initial value U1,
!    Gives back the value U2 
!
!    Francisco Javier Escoto Lopez (2019)
!    Juan A Hernandez Ramos (2020)
!*******************************************************************************  
subroutine Modified_midpoint_scheme( F, t1, t2, U1, U2, N )   
       procedure (ODES) ::    F
       real, intent(in) ::    t1, t2 , U1(:)
       real, intent(out) ::   U2(:)
       integer, intent(in) :: N
   
       real :: ti, h, U( size(U1), 0:2*N+1 ) ! number of steps is even 
       integer ::  i 

        h = (t2 - t1) / ( 2*N )

        U(:,0) = U1 
        U(:,1) = U(:,0) + h * F( U(:,0), t1 ) 

      ! Leap Frog goes to t2 + h = t_(2N+1) 
        do i=1, 2*N
            ti = t1 + i*h 
            U(:, i+1)  = U(:, i-1) + 2 * h * F( U(:,i), ti ) 
        end do
      ! average value at t2 = t1 + 2*N h       
        U2 = ( U(:, 2*N+1) + 2 *  U(:, 2*N) + U(:, 2*N-1)  )/4. 
        
        N_GBS_effort = N_GBS_effort + 2*N + 1 
       
end subroutine

!*************************************************************************
!  Richardson extrapolation coefficients
!  Solution as a polynomial function of h    
!      U(h) = a0 + a1 h**2 + a2 h**4 + a3 h**6 + ... 
!      U1 = a0 + a1 h1**2 + a2 h1**4 + a3 h1**6 + ... 
!      U2 = a0 + a1 h2**2 + a2 h2**4 + a3 h2**6 + ...
!
!  Lagrange interpolation formula with 
!      x = h**2 and Li(x) = Lagrange_i(x) 
!      Li(x) = product( (xj - x)/(xj-xi), with i/=j ) 
!      
!      U(x) = L1(x) U1 + L2(x) U2 + L3(x) U3 +.....
!      Uc (corrected solution) = L1(0) U1 + L2(0) U2 + L3(0) U3  + ... 
!
!  Barycentric formula 
!   Lj(x) = product(xj - x) / (xj-x) / product( (xj - x)/(xj-xi), with i/=j ) 
!
!   If wj =  1/ product( (xj - x)/(xj-xi), with i/=j ) 
!   and L(x) = product(xj - x), then 
!      
!   I(x) = sum ( Lj(x) Uj ) = sum ( wj/(x-xj) Uj ) / L(x) 
!   Since 1 = sum ( wj/(x-xj) ) / L(x), then L(x) = sum ( wj/(x-xj) ), 
!   
!     I(x) = sum ( wj/(x-xj) Uj ) / sum ( wj/(x-xj) )
!
!  2020, Juan A Hernandez Ramos juanantonio.hernandez@upm.es 
!*************************************************************************
function Corrected_solution_Richardson( N, U ) result (Uc) 
       integer, intent(in) :: N(:) 
       real, intent(in) :: U(:,:) 
       real :: Uc( size(U, dim=2) ) 
         
     integer :: j,  NL  ! number of levels
     real, allocatable ::  Lagrange(:), h(:), x(:), w(:)
     NL = size(N) 
     allocate( Lagrange(NL), h(NL), x(NL), W(NL) ) 
     
     h = 1. / (2*N) ! Leap Frog 
     x = h**2       ! even power of h (time step) 
     
     if (NL==1) then 
          Lagrange = 1 
          w = 1 
     else 
        do j=1, NL 
          Lagrange(j) =  product( x / ( x - x(j) ), x /= x(j) ) 
          w(j) = 1 / product( x(j) - x, x /= x(j) ) 
        end do 
     end if 
     
  !  Barycentric formula  
     Uc = matmul( w/x, U ) / sum( w/x ) 
  !  Lagrange formula   
  !  Uc = matmul( Lagrange, U)
end function 



end module 