
!***********************************************************************
!   It integrates the Cauchy problem.   
!***********************************************************************    
module Cauchy_Problem

  use ODE_Interface
  use Temporal_scheme_interface
  use Temporal_Schemes
  implicit none   
   
private 
public ::           & 
  Cauchy_ProblemS,  & ! It calculates the solution of a Cauchy problem 
  set_tolerance,    & ! It sets the error tolerance of the integration 
  set_solver,       & ! It defines the family solver and the name solver
  get_effort,       & ! # function evaluations (ODES) after integration 
  set_GBS_levels      ! It fixes the number of levels of GBS schemes

contains   

!************************************************************************************************
! It integrates the following Cauchy problem 
!
!       dU/dt = F(U, t),      U(0) = U^0 (CI)   
!
!     Inputs: 
!            Time_Domain(:)        : time discretization 
!            Differential_operator : vector values function F(U, t)
!            Scheme                : Optional Temporal numerical scheme ( default: Runge_Kutta4)
!            Solution(0,:)         : Initial condition 
!
!     Outputs: 
!            Solution(:,:)         : first index represents time 
!                                    second index  represents the i-component of solution
!     
!*  Author: Juan A Hernandez (juanantonio.hernandez@upm.es) 
!***************************************************************************************************
subroutine Cauchy_ProblemS( Time_Domain, Differential_operator,  &       
                            Solution, Scheme ) 
     real, intent(in) :: Time_Domain(:) 
     procedure (ODES) :: Differential_operator
     real, intent(out) :: Solution(:,:) 
     procedure (Temporal_Scheme), optional :: Scheme
     
!  *** Initial and final time
       real :: start, finish, t1, t2  
       integer ::  i, N_steps, ierr
       
!  *** loop for temporal integration
       call cpu_time(start)  
       N_steps = size(Time_Domain) - 1; 
       do i=1, N_steps 
           t1 = Time_Domain(i) ; t2 = Time_Domain(i+1);
         
           if (present(Scheme))  then 
               call Scheme( Differential_operator, t1, t2,            & 
                            Solution(i,:), Solution(i+1,:), ierr )   
           
           else if (family/=" ") then
               call Adavanced_Scheme  
             
           else 
               call Runge_Kutta4( Differential_operator, t1, t2,      & 
                                  Solution(i,:), Solution(i+1,:), ierr ) 
           endif 
           if (ierr>0) exit 
       enddo 
      call cpu_time(finish)
 write(*, '("Cauchy_Problem, CPU Time=",f6.3," seconds.")') finish-start 
 write(*, *) 
contains 
      
subroutine Adavanced_Scheme  
                            
                            
          if (family=="eRK") then
              
             call ERK_Scheme( Differential_operator, t1, t2,        & 
                              Solution(i,:), Solution(i+1,:), ierr ) 
             
          elseif (family=="weRK") then 
              
                call WERK( Differential_operator, t1, t2,           & 
                           Solution(i,:), Solution(i+1,:), ierr )
                
           elseif (family== "GBS") then
               
                call GBS_Scheme( Differential_operator, t1, t2,      & 
                                 Solution(i,:), Solution(i+1,:), ierr )
                
            elseif (family== "wGBS") then
                
                call WODEX( Differential_operator, t1, t2,           & 
                            Solution(i,:), Solution(i+1,:), ierr )   
                
           elseif (family== "wABM") then
               
                call WODE113( Differential_operator, t1, t2,          & 
                              Solution(i,:), Solution(i+1,:), ierr )  
                
           elseif (family== "ABM") then
               
                call PC_ABM( Differential_operator, t1, t2,           & 
                            Solution(i,:), Solution(i+1,:), ierr )  
           else 
                write(*,*) " Advanced schem not known", family
                stop 
                
           end if      
 
end subroutine 

end subroutine                 

end module 
