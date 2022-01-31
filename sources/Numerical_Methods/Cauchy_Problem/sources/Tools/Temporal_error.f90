
module Temporal_error 


    use Cauchy_Problem
    use Temporal_scheme_interface 
    use ODE_Interface 
    use Linear_systems
    implicit none

private 
public ::                      & 
  Error_Cauchy_Problem,        & ! Richardson extrapolation 
  Temporal_convergence_rate,   & ! log Error versus log time steps 
  Temporal_effort_with_tolerance  ! log time steps versus log (1/tolerance) 

    
contains
    
!*******************************************************************
!* Given the time partition t1(1:N), it returns partition t2(1:2N) 
!* by halving each time step    
!*******************************************************************    
 function refine_mesh(t1) result(t2)   
    real, intent(in) :: t1(:) 
    real :: t2( 2*size(t1)-1 ) 
    
    integer :: i, N  
    
   N = size(t1) 
   
   t2(1:2*N-1:2) =  t1; 
   t2(2:2*N-1:2) = ( t1(1:N-1)  + t1(2:N) )/ 2 
   
    
 end function       
    
!**********************************************************************************
! It determines the Error of the Solution by means of Richardson extrapolation 
! INPUTS : 
!          Time_Domain           : initial partition time to integrate 
!          Differential_operator : physical problem 
!          Solution(0,:)         : Initial condition 
!          Scheme                : temporal scheme
!          order                 : order of Error of the temporal scheme 
! OUTPUTS:
!          Solution(:,:)         : first index represents time 
!                                  second index  represents the i-component of solution
!          Error(:,:)            : error of the solution  
! Author: juanantonio.hernandez@um.es
!********************************************************************************** 
subroutine Error_Cauchy_Problem( Time_Domain, Differential_operator, & 
                                 Scheme, order, Solution, Error ) 
             real, intent(in) :: Time_Domain(0:) 
             procedure (ODES) :: Differential_operator
             procedure (Temporal_Scheme) :: Scheme
             integer, intent(in) :: order
             real, intent(out) ::  Solution(0:,:), Error(0:,:) 
             
   integer :: i, N, Nv 
   real, allocatable :: t1(:), U1(:,:), t2(:), U2(:,:)  
       
       N = size(Time_Domain)-1;  Nv = size(Solution, dim=2) 
       allocate ( t1(0:N),  U1(0:N, Nv), t2(0:2*N), U2(0:2*N, Nv) ) 
       t1 = Time_Domain
       t2 = refine_mesh(t1)  
       
       U1(0,:) = Solution(0,:); U2(0,:) = Solution(0,:) 
       
       call  Cauchy_ProblemS(t1, Differential_operator, U1, Scheme) 
       call  Cauchy_ProblemS(t2, Differential_operator, U2, Scheme)    
       
       do i=0, N 
            Error(i,:) = ( U2(2*i, :)- U1(i, :) )/( 1 - 1./2**order ) 
       end do  
       Solution = U1 + Error 
       
       deallocate ( t1,  U1, t2, U2 ) 
             
end subroutine                  
    
  
!**********************************************************************************
! It determines log Error versus log N step for a solutions of ODES
! INPUTS : 
!          Time_Domain:             initial partition time to integrate 
!          Differential_operator:   physical problem 
!          U0:                      initial condition for the solution  
!          Scheme:                  temporal scheme
!          order:                   order of Error of the temporal scheme 
! OUTPUTS:
!          log_N:       vector for the different number of time partitions 
!          log_E:       error for each time partition     
!**********************************************************************************    
subroutine Temporal_convergence_rate( Time_Domain, Differential_operator,& 
                                      U0, Scheme, order, log_E, log_N)
                  real, intent(in) :: Time_Domain(:), U0(:)  
                  procedure (ODES) :: Differential_operator
                  procedure (Temporal_Scheme), optional :: Scheme
                  real, intent(out) :: order, log_E(:), log_N(:)
     
     real :: error, coef(2) 
     real, allocatable :: t1(:), t2(:),  U1(:, :), U2(:, :)
     integer :: i, j, m, N, Nv 
     
     N = size( Time_Domain ) - 1; Nv = size( U0 );  m = size(log_N) 
     allocate ( t1(0:N), U1(0:N, Nv) )
     
     U1(0,:) = U0(:); t1 = Time_Domain 
     call Cauchy_ProblemS( t1, Differential_operator, U1, Scheme ) 
    
     do i = 1, m ! simulations in different grids 
        N =  2 * N; allocate ( t2(0:N), U2(0:N, Nv) )
        t2 = refine_mesh(t1); U2(0,:) = U0(:)
                   
        call Cauchy_ProblemS( t2, Differential_operator, U2, Scheme )  
        
        error = norm2( U2(N, :) - U1(N/2, :) ) !/ (1 - 1./2**order) 
        log_E(i) = log10( error );  log_N(i) = log10( real(N) ) 
        deallocate( t1, U1 );  allocate ( t1(0:N), U1(0:N, Nv) )
        t1 = t2;  U1 = U2;     deallocate( t2, U2)
     end do
     
     do j=1, m; if (abs(log_E(j)) > 12 ) exit;  end do
     j = min(j, m); coef = linear_regression( log_N(1:j), log_E(1:j) )
     order = abs(coef(2)); log_E = log_E - log10( 1 - 1./2**order) 
end subroutine   
 
!**********************************************************************************
! It determines log Error versus log (1/tolerance) for a solutions of ODES
! INPUTS : 
!          Time_Domain           : initial partition time to integrate 
!          Differential_operator : physical problem 
!          U0                    : initial conditions                                
!                                      
! OUTPUTS:
!          log_mu      : vector for the different tolerances 
!          log_effort   : log number of time steps       
!**********************************************************************************                                 
subroutine Temporal_effort_with_tolerance( Time_Domain, Differential_operator,  U0,  & 
                                           log_mu, log_effort)

    real, intent(in) :: Time_Domain(:),  U0(:), log_mu(:)  
    procedure (ODES) :: Differential_operator
    real, intent(out) :: log_effort(:)
                                               
    integer :: Nv, N, M, i, j  
    real, allocatable :: U(:,:)
    real :: tolerance
    integer :: N_effort, N0  
                  
    Nv = size( U0 ) 
    N = size( Time_Domain ) - 1 
    M = size( log_mu ) 
    allocate( U(0:N, Nv)  )
   
    U(0,:) = U0 
                                          
    do j = 1, M  ! simulations for different tolerances 
        
        tolerance = 10** ( -log_mu(j) )
        call set_tolerance( tolerance )
        
        N0 = get_effort() 
        call Cauchy_ProblemS( Time_Domain, Differential_operator, U) 
        N_effort =  get_effort() - N0 
        !write(*,*) "N_effort =", N_effort
        !read(*,*)
        log_effort(j) = log10( real(N_effort)  )
            
    end do 
    
end subroutine                                          
                                          
                                          

endmodule
