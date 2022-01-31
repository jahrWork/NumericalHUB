
module my_examples

        use dislin 
        use Linear_systems
        use Cauchy_Problem
        use Temporal_Schemes
        use Fourier_Interpolation
        use plots 
        use Boundary_value_problems
        use Collocation_methods
        implicit none 

       real, parameter :: PI = 4 * atan(1d0) 
       
    contains  
  


!*****************************************************************************************
! Plot a simple graph 
!*****************************************************************************************
subroutine myexampleA
 
 integer, parameter :: N=200
 real :: x(0:N), y(0:N)
 integer :: i 
 real :: a = 0, b = 2 * PI 
 
   x  = [ (a + (b-a)*i/N, i=0, N) ] 

   y = sin( 2 * x ) 
 
   call scrmod("reverse")
   call qplot(x, y, N+1) 
 
end subroutine


!*****************************************************************************************
! plot a simple trayectory 
!*****************************************************************************************
subroutine myexampleB

        
   
    integer, parameter :: N = 200  !Time steps
    real :: Time(0:N), U(0:N,1)
    real :: t0 = 0, tf = 8
    integer :: i 

    Time = [ (t0 + (tf -t0 ) * i / (1d0 * N), i=0, N ) ]
    U(0,1) =  1
   
    call Cauchy_ProblemS( Time_Domain = Time,  Differential_operator = F,& 
                          Solution = U, Scheme = Euler )
    
    call scrmod("reverse")
    call qplot(Time, U, N+1) 
  
contains

function F( U, t ) 
    real :: U(:), t 
    real :: F(size(U)) 
    
    F(1) = - U(1)
  
end function 

end subroutine
 
!*****************************************************************************************
! Plot a simple graph and generates latex file in a specific location 
!*****************************************************************************************
subroutine myexampleC 
 
    integer, parameter :: N=200, Np = 3 
    real :: x(0:N), y(0:N, Np), a = 0, b = 2 * PI 
    integer :: i 
    character(len=100) :: path(4) =                    & 
    ["./results/myexampleCa", "./results/myexampleCb", & 
     "./results/myexampleCc", "./results/myexampleCd"  ]
    
    x  = [ (a + (b-a)*i/N, i=0, N) ] 
    y(:, 1)  = sin(x); y(:, 2)  = cos(x); y(:, 3)  = sin(2*x)

   call plot_parametrics( x, y, ["$\sin x$", "$\cos x$", "$\sin 2x$"], & 
                         "$x$", "$y$", "(a)", path(1) ) 
   call plot_parametrics( y(:,1), y(:,:), ["O1", "O2", "O3"],          & 
                         "$y_2$", "$y_1$", "(b)", path(2) ) 
   call plot_parametrics( y(:,1), y(:,2:2), ["O2"],  "$y_2$", "$y_1$", & 
                          "(c)", path(3) ) 
   call plot_parametrics( y(:,1), y(:,3:3), ["O3"], "$y_2$", "$y_1$",  & 
                          "(d)", path(4) ) 
end subroutine    
    
!*****************************************************************************************
! Plot a contour graph and generates latex file in a specific location 
!*****************************************************************************************
subroutine myexampleD 
 
    integer, parameter :: N=20, Nl = 29 
    real :: x(0:N), y(0:N), z(0:N, 0:N)
    real :: levels(0:Nl), a = 0, b = 2 * PI  
    integer :: i 
    character(len=100) :: path(2) =  ["./results/myexampleDa", & 
                                      "./results/myexampleDb"  ] 
    x  = [ (a + (b-a)*i/N, i=0, N) ] 
    y  = [ (a + (b-a)*i/N, i=0, N) ]
    a = -1; b = 1 
    levels  = [ (a + (b-a)*i/Nl, i=0, Nl) ]
    z = Tensor_product( sin(x), sin(y) ) 
   
    call plot_contour(x, y, z, "x", "y", levels, "(a)",path(1),"color") 
    call plot_contour(x, y, z, "x", "y", levels, "(b)",path(2),"isolines") 
end subroutine       



!*****************************************************************************************
! Fast Fourier Transform of a real 1D function 
!*****************************************************************************************
subroutine  myExampleE 


  integer, parameter :: N = 256 
  complex :: u(0:N-1) 
  

  integer :: i
  real :: k(0:N/2-1), ck(0:N/2-1)  
  real :: x(0:N-1) 
  
  x = [ ( 2*PI*i/N, i=0, N-1 ) ]  
  u =  cos(x) ! + sin (4 * x ) 
    
   call FFT( N, u )
  
   do i=0, N-1 
    write(*,*) " k = ", u(i)
   end do 
  
   
   u =  conjg(u)  / N    ! (Duhamel et al., 1988)
  
  call FFT(N, u)  
  u = conjg( u)  
  
  do i=0, N-1 
    write(*,*) " i = ", i,  u(i)
   end do 
  
  call scrmod("reverse")
  write (*, '(A40)') 'Reconstructed FFT from its harmonics ' 
  write(*,*) "press enter " 
  read(*,*)
  call qplot( x, real(u), N )
    
  
end subroutine



subroutine  allExamples


       !call    Arenstorf_orbit
       !call Temporal_effort_with_tolerance_GBS_RK 
  
       !call Test_Linear_regression  
       !call Test_Linear_regression  
     
       !call Heat_equation_1D_CN
       !call Heat_equation_1D_system
       !call myexampleA
       !call myexampleB
       !call myexampleC 
       !call myexampleD
       !call myexampleE
     
       !call Error_Kepler_orbit
      
       !call Regions_of_absolute_stability
       
       
       !call Advanced_Heat_equation_1D
       !call BVP_1D_2021(p1 =0.4839, p2 = 0.1002) 
       !stop 
       !call Error_Heat_equation_1D
          
       !call Two_solids_different_conductivity
       !call Heat_equation_1D
       !call Heat_equation2_1D
       !call Advection_Diffusion_1D
      
       !call Wave_equation_1D
       !call Heat_equation_2D      


end subroutine

 
subroutine BVP_1D_2021(p1, p2) 
   real, intent(in) :: p1, p2 

    integer, parameter :: N = 10 
    real :: x(0:N), U(0:N), x0 = -1 , xf = 1
    integer :: i, Order = 2 
    
    integer, parameter :: Na = 100 
    real ::  xa(0:Na), Ua(0:Na)
    
    x(0) = x0; x(N) = xf 
    xa(0) = x0; xa(Na) = xf  
    call Grid_Initialization( "nonuniform", "x",  xa, 4 )
    write(*,*) " xa =", xa
    call Boundary_Value_Problem( xa, L, BCs, Ua )
    
    x(0) = x0; x(N) = xf  
    call Grid_Initialization( "uniform", "x",  x, Order )
    call Boundary_Value_Problem( x, L, BCs, U )
    
    call scrmod("reverse") 
    call qplot(x, U, N+1)
    
    call scrmod("reverse") 
    call qplot(xa, Ua, Na+1)
    
   
    write(*,*) " U(N/2)= ", U(N/2)
    write(*,*) " Ua(Na/2)= ", Ua(Na/2)
      
  

contains 
 real function L(x, y, yx, yxx) 
    
        real, intent(in) :: x, y, yx, yxx   
           
    
        L =  p1 * yxx +  yx +   sin( 10 * p2*x) * y - cos(20 * p1*x)
       
           
    end function 
    real function BCs(x, y, yx) 
    
        real, intent(in) :: x, y, yx            

        if (x==x0) then
                           BCs = yx
        elseif (x==xf) then
                           BCs = y - p1
        else 
            write(*,*) " Error BCs x=", x  
            write(*,*) " a, b=", x0, xf
            stop  
        endif            
                 
    end function 

end subroutine 




end module  

