

    
!********************************************************************************
!*  Temporal scheme for the solution of the Cauchy problem 
!*
!       U^{n+1} = G( U^n... U^{n-1+p}, F^n... F^{n-1+p}, dt ) 
!*
!*        Inputs: 
!*                F(U) vector valued function of the system of ordinary differential equations 
!*                t1 : initil time 
!*                t2 : final time  
!*                U1 :  vector for the initial condition 
!*
!*        Outputs:
!*                U2   : vector solution for the final state 
!*                ierr : integer variable to inform of internal errors 
!* 
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 2016
!********************************************************************************
module Temporal_Schemes

use Embedded_RKs
use Gragg_Burlisch_Stoer
use Adams_Bashforth_Moulton
use Wrappers 

use ODE_Interface
use Non_Linear_Systems 

implicit none

  private 
  public ::  &  
       Euler,               & ! U(n+1) <- U(n) + Dt F(U(n))
       Inverse_Euler,       & ! U(n+1) <- U(n) + Dt F(U(n+1))
       Crank_Nicolson,      & ! U(n+1) <- U(n) + Dt/2 ( F(n+1) + F(n) )
       Leap_Frog,           & ! U(n+1) <- U(n-1) + Dt/2 F(n) 
       Runge_Kutta2,        & ! U(n+1) <- U(n) + Dt/2 ( F(n)+F(U_Euler) ) 
       Runge_Kutta4,        & ! Runge Kutta method of order 4
       Adams_Bashforth2,    & ! U(n+1) <- U(n) + Dt/2 ( 3 F(n)-F(U(n-1) ) 
       Adams_Bashforth3,    & ! Adams Bashforth method of Order 3
       Predictor_Corrector1,& ! Variable step methods
  
       set_tolerance,  & 
       set_solver,     & 
       family,         & 
       get_effort,     & 
      
       ERK_scheme,     & 
       GBS_scheme,     & 
       set_GBS_levels, &
       PC_ABM, & 
       WERK,           & 
       WODEX,          & 
       WODE113
      
 character(len=50), save:: family = " "    

contains 
  

!*******************************************************************************
! Explicit Euler  
!*******************************************************************************  
 subroutine Euler(F, t1, t2, U1, U2, ierr )      
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
  
       U2 = U1 + (t2 - t1) * F(U1, t1) 
       ierr = 0    
 
 end subroutine 
!*******************************************************************************
! Implicit Inverse Euler    
!*******************************************************************************  
subroutine Inverse_Euler(F, t1, t2, U1, U2, ierr )   
     procedure (ODES) :: F
     real, intent(in) :: t1, t2, U1(:) 
     real, intent(out) ::  U2(:)
     integer, intent(out) :: ierr 

      real :: dt

      dt = t2-t1
      U2 = U1 
       
    ! Try to find a zero of the residual of the inverse Euler  
      call Newtonc( F = Residual_IE, x0 = U2 )
                
      ierr = 0
contains 

function Residual_IE(X) result(G) 
         real, target :: X(:)
         real :: G(size(X))
      
      G = X - U1 - dt * F(X, t2) 
      where (F(X, t2)==IMPOSE_ZERO) G = 0 
      
end function 
end subroutine 
  

  
!*******************************************************************************
! Crank Nicolson scheme    
!*******************************************************************************  
subroutine Crank_Nicolson(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 

      real :: dt
      real ::  a(size(U1))

      ierr = 0 
      dt = t2-t1
    
      a = U1  +  dt/2 * F( U1, t1)   
      U2 = U1 
        
      call Newtonc( F = Residual_CN, x0 = U2 )
     

  contains 
!---------------------------------------------------------------------------
!  Residual of the Crank Nicolson scheme 
!---------------------------------------------------------------------------
function Residual_CN(X) result(G)
         real, target :: X(:)
         real :: G(size(X)) 
       
         real :: Fc(size(X)) 
         
         Fc = F(X, t2) 
         
         G = X - dt/2 *  Fc - a 
         
         where (Fc==IMPOSE_ZERO) G = 0 
   
end function 


end subroutine   
  
  
  
  
 
   
!*******************************************************************************
! Runge Kutta 2 stages 
!*******************************************************************************  
subroutine Runge_Kutta2(F, t1, t2, U1, U2, ierr )
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
   
       real :: t, dt  
       real,save :: t_old
       integer :: N 
       real, save, allocatable :: k1(:), k2(:) 

       if (.not.allocated(k1)) then 
                         N = size(U1) 
                         allocate( k1(N), k2(N) )
       elseif (t1 < t_old) then
                         N = size(U1) 
                         deallocate( k1, k2 )
                         allocate( k1(N), k2(N) )
       endif 


       dt = t2-t1;   t = t1 

       k1 = F( U1, t)
       k2 = F( U1 + dt * k1, t + dt )
   
 
       U2 = U1 + dt * ( k1 + k2 )/2


      ierr = 0 
      t_old = t2

    end subroutine 

!*******************************************************************************
! Runge Kutta 4 stages 
!*******************************************************************************  
subroutine Runge_Kutta4( F, t1, t2, U1, U2, ierr )
     procedure (ODES) :: F
     real, intent(in) :: t1, t2, U1(:) 
     real, intent(out) ::  U2(:)
     integer, intent(out) :: ierr 
       
       real :: t, dt 
       real :: k1(size(U1)), k2(size(U1)), k3(size(U1)), k4(size(U1))
      
       dt = t2-t1;   t = t1 

       k1 = F( U1, t)
       k2 = F( U1 + dt * k1/2, t + dt/2 )
       k3 = F( U1 + dt * k2/2, t + dt/2 )
       k4 = F( U1 + dt * k3,   t + dt   )
 
       U2 = U1 + dt * ( k1 + 2*k2 + 2*k3 + k4 )/6

       ierr = 0 

    end subroutine 


!*******************************************************************************
!   Leap Frog 
!*******************************************************************************  
subroutine Leap_Frog(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
   

       real :: t, dt
       real, save :: t_old = 0;  
       integer ::  N 
       real, save, allocatable ::  U0(:)
       

       dt = t2 - t1;   t = t1
       
       N = size(U1)  
       
       if ( t1 < t_old .or. t_old == 0 ) then 
           
              if ( allocated(U0) ) then 
                               deallocate ( U0 )
              end if 
              allocate(  U0(N)  )
              U0 = U1 
              U2 = U1 + dt * F( U1, t) 
                   
       else 
          
            U2 = U0 + 2 * dt * F( U1, t) 
             
            U0 = U1      
                  
                               
       endif 

       ierr = 0
       t_old = t2
 
    end subroutine


!*******************************************************************************
! Second order Adams Bashforth    
!*******************************************************************************  
subroutine Adams_Bashforth2(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
   


       real :: t, dt
       real, save :: t_old = 0; 
       integer :: N
       real, save, allocatable ::   F0(:)


       dt = t2-t1;   t = t1 
       N = size(U1) 

       if (.not.allocated(F0)) then 
                         allocate( F0(N)  ) 
                         F0 = F( U1, t)  
                         U2 = U1 + dt * F0
       
       elseif (t1 < t_old ) then
                         deallocate(F0)
                         allocate( F0(N) ) 
                         F0 = F( U1, t)
                         U2 = U1 + dt * F0
                        
       else                                                      
                  U2 = U1 + dt/2 * ( 3*F( U1,t ) - F0 ) 
                  F0 = F( U1, t); 

                               
       endif 

       ierr = 0 
       t_old = t2
 
    end subroutine


!*******************************************************************************
! Third order Adams-Bashforth    
!*******************************************************************************  
subroutine Adams_Bashforth3(F, t1, t2, U1, U2, ierr )   
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 

       real :: t, dt
       real, save :: t_old = 0; 
       integer :: N
       real, save, allocatable ::   F1(:), F0(:), Fm1(:) 
       integer, save :: ipass; 


       dt = t2-t1;   t = t1 
       N = size(U1) 
              
       if (.not.allocated(F1)) then 
                         N = size(U1) 
                         allocate( F1(N), F0(N),Fm1(N) )
                         ipass = - 1 
                         
       elseif (t1 < t_old) then
                         N = size(U1) 
                         deallocate( F1, F0, Fm1 )
                         allocate( F1(N), F0(N),Fm1(N) )
                         ipass = - 1; 
       endif 
       
       ipass = ipass + 1 
       
       
      if (ipass ==0) then
                         F1 = F( U1, t)
                         U2 = U1 + dt * F1
                         F0 = F1 
                         
       elseif (ipass==1) then 
       
                        F1 = F( U1, t)
                        U2 = U1 + dt/2 * ( 3*F1 - F0) 
                        Fm1 = F0
                        F0 = F1 
                    
       else 
                        F1 = F( U1, t)
                                                                         
                        U2 = U1 + dt/12 * ( 23*F1 - 16*F0 + 5*Fm1) 
                        Fm1 = F0
                        F0 = F1 
                               
       endif 

       ierr = 0 
       t_old = t2
 
    end subroutine




!*******************************************************************************
! Predictor- corrector with variable time step 
!*******************************************************************************  
subroutine Predictor_Corrector1(F, t1, t2, U1, U2, ierr )
            procedure (ODES) :: F 
            real, intent(in) :: t1, t2, U1(:)
            real, intent(out) ::  U2(:)
            integer, intent(out) :: ierr 
            
       real  :: t, dt, dt0, Tolerance = 1d-5
       real :: Residual, b1, b2 
       real :: Up( size(U1) ), U( size(U1) ), F0( size(U1) ), F1( size(U1) )
       integer :: M 
    
  ! ** First step: Euler with dt small enough to achieve Tolerance 
  !    Since Euler local error  Error =  dt**2     
       dt = sqrt( Tolerance )   
       F0 = F( U1, t)
       U = U1 + dt * F0 
       t = t1 + dt  
       dt0 = dt 
       M = 0 
       
 ! ** Next steps with PC variable step      
      do while( t < t2 ) 
                                 
               ! ** Predictor: Adams-Bashforth orden 2  (variable time step) 
               !    page 83 Calculo Numerico en ecuaciones difenciales ordinarias (JA Hernandez)     
                    F1 = F(U, t); b1 = 2 + dt0/dt;  b2 = - dt0/dt
                    
                    Up = U + dt/2 * ( b1 * F1 + b2 * F0 )  
                    
               ! ** Corrector: Adams-Moulton orden 2  (variable time step)     
                    U  = U + dt/2 * ( F1 + F(Up, t) ) 
                    
                    dt0 = dt 
                    F0 = F1
                    
               ! ** New time step (page 93-94 same book) 
                    Residual  =  norm2(Up-U) 
                    dt =  dt * sqrt( Tolerance / Residual ) 
                    
                    dt = min( dt, t2-t) 
                    M = M + 1 
                !  write(*,*) " dt = ", dt, " Residual = ", Residual 
                    t = t + dt  
       enddo 
    
       U2 = U 
       ierr = 0 
       write(*,*) " Number of steps in t2-t1 =", M, " Residual = ",  Residual 
 
end subroutine


  

subroutine set_tolerance( tolerance ) 
      real, intent(in) :: tolerance
   
         if (family == " ") then 
   
            write(*,*) "ERROR: temporal scheme is not selected "
            stop 
   
         else 
   
             select case(family) 
             case("eRK") 
                         call set_eRK_tolerance(tolerance)
             case("GBS") 
                         call set_GBS_tolerance(tolerance)
             case("weRK") 
                         call set_weRK_tolerance(tolerance)
             case("wGBS") 
                         call set_wGBS_tolerance(tolerance) 
             case("wABM") 
                         call set_wABM_tolerance(tolerance)  
              case("ABM") 
                         call set_ABM_tolerance(tolerance)             
             case default 
                    write(*,*) " ERROR set_tolerance " 
                    stop 
               end select 
   
         end if 
        write(*,*) "set_Tolerance =", tolerance  

end subroutine 

subroutine set_solver( family_name, scheme_name) 
      character(len=*), intent(in) :: family_name
      character(len=*), intent(in), optional :: scheme_name
      
      character(len=30):: scheme
  
     family = family_name
     
     if(present(scheme_name))  then
       scheme = scheme_name
     else
       scheme ="Dopri54"
     end if
     
     if(family=="eRK") then
         call set_eRK_name(trim(scheme))
         
     elseif (family=="weRK") then
         call set_weRK_name(trim(scheme))
      
     elseif (family=="GBS") then
         
     elseif (family=="wGBS") then    
     
     elseif (family=="wABM") then 
         
     elseif (family=="ABM") then     
         
     else 
          write(*,*) " ERROR set_solver: not known family ", family_name  
         stop 
     end if
      
end subroutine
  
integer function get_effort(  ) result(N) 
     
         if (family == " ") then 
   
            write(*,*) "ERROR: temporal scheme is not selected "
            stop 
   
         else 
   
             select case(family) 
             case("eRK") 
                         N = get_eRK_effort()
             case("weRK") 
                         N = get_weRK_effort()
             case("GBS") 
                         N = get_GBS_effort()
             case("wGBS") 
                         N = get_wGBS_effort()
             case("ABM") 
                         N = get_ABM_effort()   
             case("wABM") 
                         N = get_wABM_effort() 
                        
             case default 
                    write(*,*) " ERROR set_tolerance " 
                    stop 
               end select 
   
         end if 

end function 


end module 





 
 
  
