module API_Example_Fourier_series

    use Fourier_interpolation
    use Collocation_methods 
    use Linear_systems
    use plots 
    use API_Example_IBVP_Fourier
    use Utilities
    
    implicit none
    real :: PI = 4 * atan(1d0)
    
    
    
    contains
    
   
 subroutine Fourier_problems 
 
  integer :: option
    
  option = 1  
  do while (option>0) 
  write(*,*) " Enter example to execute " 
  write(*,*) " 0. Exit  "
  write(*,*) " 1. Fourier interpolation examples" 
  write(*,*) " 2. Fourier 2D derivative  " 
  write(*,*) " 3. Fourier IBVP examples  "
  write(*,*) " 4. Identify period of a signal with FFT  "
  read(*,*) option 
 
  select case(option) 
       
     case(0) 
          exit 
 
      case(1) 
          call Fourier_interpolation_examples 
        
      case(2)  
          call Test_Fourier_Spectral_Derivative2D 
          call Test_Fourier_Derivative2D
        
     case(3) 
        call Fourier_IBVP_examples
        
     case(4) 
         call Identify_period_with_FFT
   
 
     case default 
         write(*,*) " Not implemented"
         
      end select
  end do 
  
 
 end subroutine   
    
    
    
subroutine Fourier_interpolation_examples

     call Test_Fast_Fourier_Transform1D
     call Test_direct_and_inverse_FFT1D
     call Test_Fourier_Spectral_Derivative1D
     call Test_Fourier_Derivative1D
     call Fourier_Derivative_error1D
   
     
end subroutine  


!*****************************************
subroutine  Test_Fast_Fourier_Transform1D 
!*****************************************

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
   
   
   u = conjg( u ) / N    ! (Duhamel et al., 1988)
 
  call FFT(N, u)  
  u = conjg( u)  
  
  do i=0, N-1 
    write(*,*) " i = ", i,  u(i)
   end do 
  
  call scrmod("reverse")
  write (*, '(A40)') 'It plots cos x + sin 4 x  ' 
  write (*, '(A40)') 'Reconstructed FFT from its harmonics ' 
  write(*,*) "press enter " 
  read(*,*)
  call qplot( x, real(u), N )
    
  
end subroutine  

!*****************************************
subroutine  Test_direct_and_inverse_FFT1D
!*****************************************

  integer, parameter :: N = 256 
  complex :: u(0:N-1) 
  complex :: C(-N/2:N/2-1) 

  integer :: i, k 
 
  real :: x(0:N-1) 
  
  x = [ ( 2*PI*i/N, i=0, N-1 ) ]  
  U =  cos(x) ! + sin (4 * x ) 
    
  C =  FFT_1D( N, U )
      
  U =  iFFT_1D(N, C) 
  
  call scrmod("reverse")
  call qplot( x, real(U), N )
    
  
end subroutine    
  

!**********************************************
subroutine  Test_Fourier_Spectral_Derivative1D
!**********************************************

  integer, parameter :: N =  256 
  complex :: U(0:N-1), c(-N/2:N/2-1)  
  real :: Vxi(0:N-1), V(0:N-1) 
  
  integer :: i, j, k 
  real :: x(0:N-1) 
  
   x = [ ( 2*PI*i/N, i=0, N-1 ) ]  
  
   U =    Triangle( x ) 
   call scrmod("reverse")
   call qplot( x, real(U), N ) 
    
   C = FFT_1D(N, U )
   C = Spectral_derivative1D(N, C, 1)
   U = iFFT_1D(N, C)  
   
   write (*, '(A40)') 'Reconstructed derivative of u(x) from its harmonics ' 
   call qplot( x, real(U), N )
  
  
end subroutine 


!*****************************************
subroutine  Test_Fourier_Derivative1D 
!*****************************************

  integer, parameter :: N =  256 
  complex :: U(0:N-1), c(-N/2:N/2-1)  
  real :: Vxi(0:N-1), V(0:N-1) 
  
  integer :: i, j, k 
  real :: x(0:N-1) 
  
   x(0) = 0; x(N-1) = 1; 
   call Grid_Initialization( "Fourier", "x",  x )
  
   V = exp( -200*(x-0.5)**2 ) 
   call scrmod("reverse")
   call qplot( x, V, N )
  
   
   do k=1, 10 
     call Derivative( "x",  k, V, Vxi ) 
     call qplot( x, Vxi, N )
   end do 
    
  
end subroutine 



subroutine  plot_harmonics(N, u) 
     integer, intent(in) :: N 
     complex :: u(0:) 
     
      real :: k(0:N/2-1), ck(0:N/2-1), ak(0:N/2), bk(0:N/2)   
      integer :: i 
      complex :: II= (0., 1.) 
  
   
   do i=1, N/2-1 
    k(i) = i 
    ck(i) = 2 * abs(u(i)) / N 
    ak(i) = ( u(i) + u(N-i) )/ 2 
    bk(i) = ( u(i) - u(N-i) )/ (2*II)  
    ck(i) = 2 * sqrt( ak(i)**2 + bk(i)**2 ) / N 
    
   end do 
   
   call scrmod("reverse")
   call qplot( k, ck, N/2 ) 

end subroutine 
   
elemental real function unit_pulse(x) result(f) 
     real, intent(in) :: x 
     
     if ( x < PI/2 ) then 
                                  f = 0 
     else if ( x < 3*PI/2 ) then 
                                  f = 1 
     else 
                                  f = 0 
     end if 
     
     
end function  

elemental real function sawtooth(x) result(f) 
     real, intent(in) :: x 
     
     if ( x < PI/2 ) then 
                                  f = 0 
     else if ( x < 3*PI/2 ) then 
                                  f = x - PI/2 
     else 
                                  f = 0 
     end if 
     
     
end function 

elemental real function Triangle(x) result(f) 
     real, intent(in) :: x 
     
     if ( x < PI ) then 
                                  f = x
     else 
                                  f = 2*PI - x  
     end if 
     
     
end function  



!********************************************************************
!* 
!*****************************************************************
subroutine Fourier_Derivative_error1D
 
 integer :: q                      ! interpolant order 2, 4, 6, 8 
 integer :: N                      ! # of nodes (piecewise pol. interpol.) 
 integer :: k = 1                  ! derivative order 
 integer :: p = 5                  ! where error is evaluated p=0, 1,...N
 integer, parameter :: M = 10      ! number of grids ( N = 8 * 2**j) 
 real :: log_Error(M,4),log_dx(M)  ! Error versus Dx for q=2, 4, 6, 8
 real :: epsilon = 1d-12           ! order of the random perturbation 
 
 real :: PI = 4 * atan(1d0), logN  
 integer ::  j, l=0
 
 real, allocatable :: x(:), f(:), dfdx(:)  ! function to be interpolated 
 real, allocatable :: dIdx(:)              ! derivative of the interpolant 
 real :: beta
 
 
 do l=1, 4   ! l=1 Fourier and FD l=2,3,4
     
  do j=1, M  ! discretization points N=N(j) 
     
   N = 8 * 2**j 
   
   allocate( x(0:N-1), f(0:N-1),  dfdx(0:N-1), dIdx(0:N-1) ) 
   x(0) = 0; x(N-1) = 2*PI *(N-1)/N 
   
   q = 2*l
   if (l==1) then 
      call Grid_Initialization( "Fourier", "x",  x )
   else 
       call Grid_Initialization( "uniform", "x", x, q )
   end if 
  
   call random_number(f)
   beta = 60 
   epsilon = 1d-10 
   f = exp(- beta*(x-PI)**2 ) + epsilon * f
   dfdx = - beta * 2* (x-PI) * exp(- beta*(x-PI)**2 ) 
   call Derivative( "x", k, f, dIdx )
   
   if (j==2) then 
       write(*,*) " N = ", N 
       write(*,*) " blue: exact, red: numeric (press enter) "
     
       call metafl("xwin")
       CALL PAGE (4000, 4000)
       call scrmod("reverse")
       call disini 
        call graf(0., 2*PI, 0., 1., -10., 10., -10., 2.) 
        call color("blue"); call curve( x, dfdx, N) 
        call color("red");  call curve( x, dIdx, N)
       call disfin
   end if 
    
   p = N/2-1 
     
   log_dx(j) = log( x(1) - x(0) ) 
   log_Error(j, l) = log( abs(dIdx(p) - dfdx(p)) ) 
 
   deallocate( x, f, dIdx, dfdx ) 
  end do 
 end do  

 call scrmod("reverse")
 write(*,*) "Second derivative error versus spatial step for q=2,4,6,8 " 
 write(*,*) " Test function:  f(x) = cos x  " 
 write(*,*) "press enter " ; read(*,*)
 
 call plot_parametrics( log_dx, log_Error, ["Fourier","FD", "FD", "FD"], & 
                       "log_dx","log_Error")
 
end subroutine 





!**********************************************
subroutine  Test_Fourier_Spectral_Derivative2D 
!**********************************************

  integer, parameter :: N = 32, M = 32  
  real  :: U(0:N-1, 0:M-1)
  complex :: Uc(0:N-1, 0:M-1), c(-N/2:N/2-1, -M/2:M/2-1 )  
  real :: Vxi(0:N-1, 0:M-1), V(0:N-1) 
  
  integer :: i, j, k, l  
  real :: x(0:N-1), y(0:M-1)  
  integer, parameter :: N_levels = 9 
  real :: levels(0:N_levels) 
  real :: z0, zf 
  
   x = [ ( 2*PI*i/real(N), i=0, N-1 ) ]  
   y = [ ( 2*PI*j/real(M), j=0, M-1 ) ] 
   levels = [ ( z0 + (zf-z0) * j/real(N_levels)  , j=0, N_levels )]
   
   call Grid_Initialization( "Fourier", "x",  x )
   call Grid_Initialization( "Fourier", "y",  y )
   
   U =  Tensor_product( sin(x), sin(y) )  
   
   call scrmod("reverse")
   call plot_contour(x, y, U, "x","y", levels, graph_type ="isolines") 
    
   Uc = U 
   C = FFT_2D(N, M, Uc )
   
   C = Spectral_derivative2D(N, M, C, 1, 1)
   Uc = iFFT_2D(N, M, C )
   U = real(Uc)
   call plot_contour(x, y, U, "x","y", levels, graph_type ="isolines")
   
   C = Spectral_derivative2D(N, M, C, 1, 1)
   Uc = iFFT_2D(N, M, C )
   U = real(Uc)
   call plot_contour(x, y, U, "x","y", levels, graph_type ="isolines")
   
   C = Spectral_derivative2D(N, M, C, 1, 2)
   Uc = iFFT_2D(N, M, C )
   U = real(Uc)
   call plot_contour(x, y, U, "x","y", levels, graph_type ="isolines")
   
end subroutine  



!*****************************************
subroutine  Test_Fourier_Derivative2D 
!*****************************************

  integer, parameter :: N = 32, M = 32  
  real  :: U(0:N-1, 0:M-1), Ux(0:N-1, 0:M-1)
  real :: x(0:N-1), y(0:M-1)  
  integer, parameter :: N_levels = 9 
  real :: levels(0:N_levels) 
  real :: z0, zf 
  integer :: i, j 
  
   x = [ ( 2*PI*i/N, i=0, N-1 ) ]  
   y = [ ( 2*PI*j/N, j=0, M-1 ) ] 
   levels = [ ( z0 + (zf-z0) * j/real(N_levels)  , j=0, N_levels )]
   call Grid_Initialization( "Fourier", "x",  x )
   call Grid_Initialization( "Fourier", "y",  y )
  
  
   U =  Tensor_product( sin(x), sin(y) )  
   call scrmod("reverse")
   call plot_contour(x, y, U, "x","y", levels, graph_type ="isolines") 
   
   ! partial u / partial x 
   call Derivative( [ "x", "y" ], 1, 1, U, Ux ) 
   call plot_contour(x, y, Ux, "x","y", levels, graph_type ="isolines")
   
   ! partial u / partial y 
   call Derivative( [ "x", "y" ], 2, 1, U, Ux ) 
   call plot_contour(x, y, Ux, "x","y", levels, graph_type ="isolines")
  
   
end subroutine 



!*****************************************************************************************
! Identify the period of a periodic signal by considering the minimum number of samples
!*****************************************************************************************
subroutine  Identify_period_with_FFT

  
  integer :: N, Ns, Nf 
  
  complex, allocatable :: U(:), dU(:), dUr(:), C(:), dC(:)  
  complex, allocatable :: Us(:), dUs(:)   
  real, allocatable :: xs(:), x(:), xr(:), frequency(:), Amplitude(:) 
  real, allocatable :: Error(:) 
  real :: xmax, xmin, ymax, ymin 

  integer :: i, k
  real :: f 
 
! Total number of samples 
  Ns = 1000
! Fraction of samples to identify the period 
  N =  64
  
  allocate(  x(0:N-1),  U(0:N-1), dU(0:N-1),  dUr(0:N-1), Error(0:N-1) ) 
  allocate(  C(-N/2:N/2-1), dC(-N/2:N/2-1), frequency(0:N/2-1), Amplitude(0:N/2-1) ) 
  allocate( xs(0:Ns-1), Us(0:Ns-1), dUs(0:Ns-1) ) 
  
  xs  = [ (2*PI*i/real(Ns), i=0, Ns-1 ) ]
  Us  = sin(16*xs) + sin(24*xs) + cos(32*xs) 
  dUs = derivative_a( xs, real(Us) )
  
! First N samples of Us  
  x = xs(0:N-1) 
  U = Us(0:N-1) 
  C =   FFT_1D( N, U )
  
! Spectral Derivative and inverse result dUr    
  dC = Ns/N * Spectral_derivative1D(N, C, 1)
  dUr = iFFT_1D(N, dC) 
  
! Error between derivative of samling and reconstruction by means of  FFT
! Increase N until the Error is small   
  Error = real(dUs(0:N-1)) - real(dUr) 
  
  call scrmod("reverse")
  call qplot( x(2:N-3), Error(2:N-3), N-4)
  
  call disini  
  xmax = maxval(xs); xmin = 0.0; ymax = 10*2*PI; ymin = -10*2*PI 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10 , ymin, ymax, ymin, (ymax-ymin)/10)   
  call color("blue");  
  call curve( xs, real(dUs), Ns) 
  call color("red");  
  call curve( x, real(dUr), N)
  call disfin 
    
  do k=0, N/2-1 
    frequency(k) = k 
    Amplitude(k) = sqrt( abs(C(-k)**2) + abs(C(k)**2) ) 
   ! write(*,*) " k = ", k, C(k) 
  end do 
  Nf = N/2
  call disini  
  xmax = Nf; xmin = 0.0; ymax = 1.5; ymin = -0.1 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10 , ymin, ymax, ymin, (ymax-ymin)/10)    
  call marker(21); call incmrk(-1)
  call color("red");  
  call curve( frequency, Amplitude, Nf+1) 
  call disfin 
     
end subroutine

function derivative_a( x,  U) result(dU) 
   real, intent(in) :: x(0:), U(0:) 
   real dU( 0:size(U)-1 ) 
   
   integer :: i, N 
   
   N = size(x) - 1 
   
   dU = 0 
   do i=2, N-3      
     dU(i) = ( U(i-2) - 8 * U(i-1) + 8 * U(i+1) - U(i+2) )/( 6*(x(i+1) - x(i-1)) )  
   end do  
   
   
end function 




!*****************************************************************************************
! Fast Fourier Transform of a real 1D function 
!*****************************************************************************************
subroutine  FFT_example


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



end module