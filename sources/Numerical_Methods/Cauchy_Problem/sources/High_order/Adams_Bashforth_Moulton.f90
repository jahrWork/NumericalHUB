module Adams_Bashforth_Moulton
    
    
    use ODE_Interface
    use Temporal_scheme_interface
    use Lagrange_interpolation
    implicit none
    
  
    private 
    public ::                   & 
              PC_ABM,           & ! Predictor Corrector ABM 
              set_ABM_tolerance,& ! set tolerance error 
              get_ABM_effort      ! # of function evaluation (ODES) 
    
    real, save, allocatable :: Y1(:,:) , Y2(:,:) , B(:,:), r(:)
    integer, save :: Steps
    
    real, save :: ABM_tolerance  = 0 
    integer, save :: N_ABM_effort = 0 
    
    logical, save :: IC_set = .false. 
    
    
contains    
   

!******************************************************************
!*  TODO: error tolerance allow to determine Steps and h (time step) 
!******************************************************************  
subroutine set_ABM_tolerance(eps)
      real, intent(in) :: eps
      
  ABM_Tolerance = eps
  Steps = 4
  IC_set  = .false. 
    
end subroutine

integer function get_ABM_effort() result(N) 
     
      N = N_ABM_effort
      
end function 


!******************************************************************
!*  Predictor Corrector based on Adams-Bashforth and Adams-Moulton  
!******************************************************************
subroutine PC_ABM( F, t1, t2, U1, U2, ierr )   
       procedure (ODES)     ::  F
       real, intent(in)     ::  t1, t2,  U1(:) 
       real, intent(out)    ::  U2(:)
       integer, intent(out) ::  ierr 
        
       real :: Up( size(U1) ) 
       
       if (.not.IC_set)  then 
           call set_ABM_IC( F, U1, t1, t2-t1)
           IC_set  = .true. 
       end if 
       
       
       call Predictor_AB( F, t1, t2, U1, Up, ierr )  
       call Corrector_AM( F, t1, t2, U1, Up, U2, ierr )  
       
       ierr = 0 
       
end subroutine

!********************************************************
!* Predictor Adams Bashforth 
!********************************************************
subroutine Predictor_AB( F, t1, t2, U1, U2, ierr )   
       procedure (ODES)     :: F
       real, intent(in)     :: t1, t2, U1(:) 
       real, intent(out)    :: U2(:)
       integer, intent(out) :: ierr 
   
       real :: h, alpha( size(U1) )
       integer :: i
       
        h = t2 - t1 
        Y1(1,:) = h * F(U1, t1) 
        
        Y2(0,:)  = matmul( B(0,:), Y1 ) 
        alpha = h * ( F(Y2(0, :), t2) ) - matmul( B(1,:), Y1 ) 
                 
        do i = 1, steps 
           Y2(i,:) = matmul( B(i,:), Y1 ) + alpha * r(i)
        end do
        
        U2 = Y2(0,:)
        ierr = 0 
        
        N_ABM_effort = N_ABM_effort + 2 
              
end subroutine
 
!********************************************************
!* Corrector Adams Moulton 
!********************************************************
subroutine Corrector_AM( F, t1, t2, U1, Up, U2, ierr )   
       procedure (ODES) ::    F
       real, intent(in) ::    t1, t2, U1(:), Up(:)   
       real, intent(out) :: U2(:)
       integer, intent(out) :: ierr
         
       real :: h, alpha( size(U1) ) 
       integer :: i
 
       h = t2 - t1
       Y1(1,:) = h * F(U1, t1) 
      
       alpha = h * F(Up, t2) - matmul( B(1,:), Y1 ) 
                 
        do i = 0, Steps 
           Y2(i,:) = matmul( B(i,:), Y1 ) + alpha * r(i)
        end do
        
        U2 = Y2(0,:)
        Y1 = Y2 
        ierr = 0
        
        N_ABM_effort = N_ABM_effort + 2 
        
end subroutine




!*******************************************************************************
!    Computes the mabtix B and coefficients r_i for a multivalue method
!    Author:  
!       Francisco Javier Escoto Lopez 2019 
!       Juan A Hernandez Ramos        2019 
!******************************************************************************* 
subroutine set_ABM_IC( F, U0, t0, h) 
        procedure(ODES) :: F
        real, intent (in) :: U0(:), t0, h
        
        integer :: i, j, Nv 
        real :: x(0:Steps-1), xp
        real Polynomials(-1:Steps-1, 0:Steps-1)
      
       
          Nv = size(U0)
          if (.not.allocated(Y1)) then 
            allocate ( Y1(0:Steps, Nv) , Y2(0:Steps, Nv) , B(0:Steps, 0:Steps), r(0:Steps) ) 
          end if   
        
    ! *** Grid definition for Lagrange interpolant
          x = [( i, i=0, Steps-1 )]
          xp = x(Steps-1)
          Polynomials = Lagrange_polynomials( x, xp ) 
       
          do i = 0, Steps
             r(i) = Polynomials(i-1, Steps-1) / Factorial(i) 
          end do
            
    ! *** Matrix for extrapolation computation
          B = 0
          do i = 0, Steps 
            do j = i, Steps 
                 B(i,j) = Factorial(j) / ( Factorial(i) * Factorial(j-i) )
            end do
          end do
            
    ! *** Initial conditions of multivalue methods   
          call Initial_conditions( F, U0, t0, h )
            
    
end subroutine
    
integer function Factorial(n)
    integer, intent(in) :: n
      
      Factorial = gamma( real(n+1) ) 
         
end function

!*******************************************************************************
!    Computes the initial conditions for Y1 which are derivatives 
!******************************************************************************* 
subroutine Initial_conditions( F, U0, t0, h )   
       procedure (ODES) ::    F 
       real, intent(in) :: U0(:), t0
       real, intent(in) :: h
    
       real :: Time(0:Steps-1), U( 0:Steps-1, size(U0))
       real :: Polynomials(-1:Steps-1,0:Steps-1 )  , FF( 0:Steps-1,size(U0))
       integer :: i, s 
       
       ! *** Grid definition for Lagrange interpolant
       !     and first "s" steps  
             s = Steps 
             Time =  [( t0 + i*h, i=0,s-1 )] 
       
       ! *** Lagrange Polynomials at t=t0
             Polynomials = Lagrange_polynomials( Time, t0 ) 
       
       ! *** Explicit Euler for first "s" steps
             U(0,:) = U0
             FF(0,:) = F( U0,Time(0) )
             do i = 0, s-2
                U(i+1,:) = U(i,:) + h * FF(i,:)
                FF(i+1,:) = F( U(i+1,:),Time(i+1) )
             end do
       
       ! *** Initial condition for derivatives from interpolant
             Y1(0,:) = U0
             Y1(1,:) = h* FF(0,:)
             do i = 2, s
               Y1(i,:) = matmul( Polynomials(i-1,:), FF ) * h**i / Factorial(i)
             end do 
       
 end subroutine  


                            
end module 




