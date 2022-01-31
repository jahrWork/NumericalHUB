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
!*          Manú J. Soto-Aranaz (manu.soto-aranaz.gonzalez@alumnos.upm.es) 2019
!********************************************************************************
module Wrappers

use ODE_Interface

implicit none

   public :: set_weRK_name, set_weRK_tolerance, set_wGBS_tolerance
   
   real, save:: wRK_Tolerance = 0
   real, save:: wGBS_Tolerance = 0 
   real, save:: wABM_Tolerance = 0 
   character (len=20), save :: WRK_Method = " " 
   
   integer, save:: N_weRK_effort = 0
   integer, save:: N_wGBS_effort = 0 
   integer, save:: N_wABM_effort = 0 
 
contains 

!*********************************************************
! Runge Kutta wrapper selection, tolerance and time steps 
!*********************************************************    
subroutine set_weRK_name( name)
      character(len=*), intent(in) :: name 
      
      wRK_Method = name
      
end subroutine

subroutine set_weRK_tolerance( eps)
      real, intent(in) :: eps
      
     wRK_Tolerance = eps
    
end subroutine    

integer function get_weRK_effort() result(N) 
     
      N = N_weRK_effort 
      
end function 
    

!*********************************************
! GBS wrapper tolerance and time steps 
!********************************************* 
subroutine set_wGBS_tolerance( eps)
      real, intent(in) :: eps
     
      wGBS_Tolerance = eps            
      
end subroutine    

integer function get_wGBS_effort() result(N) 
     
      N = N_wGBS_effort 
      
end function 



!********************************************
! ABM wrapper tolerance and time steps 
!****************************************** 
subroutine set_wABM_tolerance( eps)
      real, intent(in) :: eps
      
      wABM_Tolerance = eps            
      
end subroutine  

integer function get_wABM_effort() result(N) 
     
      N = N_wABM_effort 
      
end function 



!*******************************************************************************
! Wrapper of dopri5
!*******************************************************************************  
subroutine WERK(F, t1, t2, U1, U2, ierr )
        procedure (ODES) :: F
        real, intent(in) :: t1, t2 
        real, intent(in) ::  U1(:)
        real, intent(out) ::  U2(:)
        integer, intent(out) :: ierr  !intent out changed for test
                    
        integer :: N, NRDENS, LWORK, LIWORK, iout, itol, idid
        integer, allocatable :: IWORK(:), IPAR(:)
        real, allocatable ::  WORK(:), RPAR(:), RTOL(:), ATOL(:)
        external dopri5
        external dop853
        
        if (WRK_Method == " ") WRK_Method = "WDOPRI5"
        
        if (WRK_Tolerance == 0) WRK_Tolerance = 1d-4
       
        
!       DIMENSION OF THE SYSTEM
        N = size(U1) 
        NRDENS = 2
        IF (WRK_Method=="WDOPRI5") THEN
        LWORK = 8*N + 5*NRDENS + 21
        ELSEIF (WRK_Method=="WDOP853") THEN
        LWORK = 11*N+8*NRDENS+21
        END IF   
    
        LIWORK = NRDENS + 21
        allocate( IWORK(LIWORK), IPAR(2) ) 
        allocate( WORK(LWORK), RPAR(2), RTOL(N),ATOL(N) )
        
!       OUTPUT ROUTINE (AND DENSE OUTPUT) USED DURING INTEGRATION
        IOUT=2
!       DEFAULT VALUES FOR PARAMETERS
        IWORK = 0; WORK = 0. 
!       DENSE OUTPUT IS USED FOR THE TWO POSITION COORDINATES 1 AND 2
!       IWORK(5)=NRDENS; IWORK(21)=1; IWORK(22)=2  
        
        WORK(7)=t2-t1

!       REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL = 0; RTOL = WRK_Tolerance;  ATOL = RTOL
        U2 = U1
        IF (WRK_Method=="WDOPRI5") THEN
            CALL DOPRI5(N, FUNC, t1, U2, t2,              &
                        RTOL, ATOL, ITOL,                  &
                        SOLOUT, IOUT,                      &
                        WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)
        ELSEIF (WRK_Method=="WDOP853") THEN
            CALL DOP853(N, FUNC, t1, U2, t2,              &
                        RTOL, ATOL, ITOL,                  &
                        SOLOUT, IOUT,                      &
                        WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)
        Else 
            write(*,*) "Invalid Method selection use WDOPRI5 or WDOP853"
            stop 
        END IF 
             
       ierr =0
contains 

    subroutine FUNC(N, x, y, Fv, RPAR, IPAR)
       integer :: N, ipar 
       real :: x, y(N), Fv(N), RPAR(2)
           
          Fv = F( y, x ) 
          N_weRK_effort = N_weRK_effort + 1 
      
    end subroutine

    subroutine SOLOUT(NR,XOLD,X,Y,N,CON,ICOMP,ND,RPAR,IPAR,IRTRN)
        integer :: NR,N,ICOMP,ND,IPAR,IRTRN
        real :: XOLD,X,Y,CON,RPAR
    end subroutine     

end subroutine 
!*******************************************************************************
! Wrapper of ODEX
!*******************************************************************************  
subroutine WODEX(F, t1, t2, U1, U2, ierr )
       procedure (ODES) :: F
       real, intent(in) :: t1, t2, U1(:)  
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr 
    
        integer :: N, NRDENS, LWORK, LIWORK, iout, itol, idid, KM
        integer, allocatable :: IWORK(:)
        real, allocatable ::  WORK(:), RTOL(:), ATOL(:), RPAR(:) 
        
        integer, allocatable :: IPAR(:) 
        REAL:: H
        external ODEX
        
        if (WGBS_Tolerance == 0) WGBS_Tolerance = 0.1
           
         
        
!       DIMENSION OF THE SYSTEM
        N = size(U1) 
        allocate( RTOL(N), ATOL(N), RPAR(2), IPAR(2) ) 
        NRDENS = 2
        KM = 9   !KM=9 IF IWORK(2)=0 / KM=IWORK(2) IF IWORK(2).GT.0
        LWORK = N*(KM+5)+5*KM+20+(2*KM*(KM+2)+5)*NRDENS
        LIWORK = 2*KM+21+NRDENS
        allocate( IWORK(LIWORK),  WORK(LWORK) ) 
       ! allocate(IPAR(2), RPAR(2), RTOL(N),ATOL(N) )
        H=0.01D0
!       OUTPUT ROUTINE (AND DENSE OUTPUT) USED DURING INTEGRATION
        IOUT=2
!       DEFAULT VALUES FOR PARAMETERS
        IWORK = 0; WORK = 0. 
!       DENSE OUTPUT IS USED FOR THE TWO POSITION COORDINATES 1 AND 2
        IWORK(8)=NRDENS; IWORK(21)=1; IWORK(22)=2      

!       REQUIRED (RELATIVE AND ABSOLUTE) TOLERANCE
        ITOL = 0; RTOL = WGBS_Tolerance;  ATOL = RTOL
        U2 = U1 
        CALL ODEX(N, FVDPOL, t1, U2, t2, H,              &
                  RTOL, ATOL, ITOL,  SOLOUT, IOUT,       &
                  WORK, LWORK, IWORK, LIWORK, RPAR, IPAR, IDID)
     
       ierr =0     
contains 

    subroutine FVDPOL(N, X, Y, FV, RPAR, IPAR)
        integer :: N, ipar 
        real :: x, y(N), Fv(N), RPAR(2) 
        
            Fv = F( y, x ) 
            N_wGBS_effort = N_wGBS_effort + 1 
            
    end subroutine

    subroutine SOLOUT(NR,XOLD,X,Y,N,CON,NCON,ICOMP,ND,RPAR,IPAR,IRTRN)
        real :: NR,N,ND,IPAR,IRTRN, NCON
        real :: XOLD,X,Y(N),CON(NCON),ICOMP(ND),RPAR
    end subroutine     

end subroutine

!*******************************************************************************
! Wrapper of ODE113
!*******************************************************************************

subroutine WODE113(F, t1, t2,  U1, U2, ierr)
      
       procedure (ODES) :: F
 
       real, intent(in) :: t1, t2 
       real, intent(in) ::  U1(:)
       
       real, intent(out) ::  U2(:)
       integer, intent(out) :: ierr
       
       ! **** Variables for the wrapping
       
    real :: WORK(226),RPAR(1),RTOL(1),ATOL(1), H
    real :: relerr,abserr
    integer:: i, N, LWORK, IDID, IWORK(39), LIWORK,ITOL,IPAR(1), IOUT, iflag!LWORK87!LIWORK21
    integer :: GL
    
    real :: V(size(U1)) 
    
    if (wABM_Tolerance == 0) wABM_Tolerance = 0.1
    
    WORK=0.
    !RPAR=0.
   
    !H = 0.
    !LWORK=226
    IWORK=0
    !LIWORK=39
    !ITOL=0
    !IPAR=0
    !IOUT=0
    
    abserr = wABM_tolerance 
    relerr = wABM_tolerance 
    iflag = 1
    
    GL = size(U1) 
    U2 = U1
    call ode(F_call, GL, U2, t1, t2, relerr, abserr, iflag, work, iwork)
    
    
    ierr = 0
contains 
subroutine F_call(t, y, yp) !Adaptado para las subrutina ODE
    real t, y(GL), yp(GL)
        
        yp = F(y,t) 
        N_wABM_effort = N_wABM_effort + 1 
    
    end subroutine


end subroutine


end module 





 
 
  
