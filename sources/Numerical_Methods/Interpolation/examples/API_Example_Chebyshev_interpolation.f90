module API_Example_Chebyshev_interpolation

    use Chebyshev_interpolation
    use Collocation_methods 
    use Interpolation
    use Lagrange_interpolation
    use plots 
    
    implicit none
   
    
    contains
   

    
subroutine Chebyshev_interpolation_examples

    write (*,*) "**********************************************************" 
    write (*,*) " Chebyshev interpolation examples" 
    write (*,*) "**********************************************************"  
   
    call Truncation_error_cos_pikx 
    call Error_exact_versus_estimated 
    call Error_second_derivative_Chebyshev_versus_FD
   
    
    call Fast_chebyshev_transform_example 
      
    
end subroutine      
    
 
   

subroutine Truncation_error_cos_pikx
 
 integer, parameter :: N=10, M=200
 real :: x(0:N), y(0:N), xp(0:M), yp(0:M), y2p(0:M)  
 real :: Lebesgue_N(-1:N, 0:M),  PI_N(0:N, 0:M), I_N(0:N, 0:M)  
 real :: a=-1, b=1, PI = 4 * atan(1.) 
 integer :: i, k, degree   
 
   x  = [ (-cos(PI*i/N), i=0, N) ]
   xp = [ (a + (b-a)*i/M, i=0, M) ]
   
!  sampling rule:
!  Only 4 points are needed to resolve a wavelength of wavenumber k   
   k =  N/4  
   y = cos( k * PI * x ) 
   yp = cos( k * PI * xp ) 
   y2p = -(k*PI)**2 * cos( k * PI * xp ) 
 
   degree = N 
   Lebesgue_N = Lebesgue_functions( x, xp ) 
   PI_N = PI_error_polynomial( x, xp ) 
   I_N = Interpolant( x, y, degree, xp )
     
   call plot_interpolants 
   call plot_PI
   
contains 
subroutine  plot_interpolants

    real :: xmin, xmax, ymin, ymax
    integer :: j, N 
    character(len=10) :: col(4) = ["red","blue", "white", "orange" ]
   
    xmin = minval(x); xmax = maxval(x) 
    ymax = maxval(I_N(2,:) ) ; ymin =minval(I_N(2,:) );
    call scrmod("reverse") 
    call metafl("xwin")
   
    call page(4000, 4000)
    call disini 
    call graf(xmin, xmax, xmin, (xmax-xmin)/5, ymin, ymax, ymin, (ymax-ymin)/5 )
   
         
    j = 1 
    call color(col(j));! call incmrk(1);   call marker(21);  
    call curve( xp, yp, M+1 )
    
    j = 2 
    call color(col(j)); !call incmrk(1);   call marker(21);  
    call curve( xp, I_N(0,:), M+1 )
    
     j = 3 
    call color(col(j)); call incmrk(1);   call marker(21);  
    call curve( xp, y2p, M+1 )
    
    j = 4 
    call color(col(j)); call incmrk(1);   call marker(21);  
    call curve( xp, I_N(2,:), M+1 )
    
    call plot_legends(["y(x)", "I(x)","d2y/dx2", "d2I/dx2" ])
    call plot_title( ["Chebyshev interpolant and second derivative",  & 
                     " y(x)=cos(k pi x) with N=10(sampling rule : N > 4 k)" ] ) 
  
    call disfin 

end subroutine 
subroutine  plot_PI

    real :: xmin, xmax, ymin, ymax
    integer :: j, N 
    character(len=10) :: col(5) = ["red", "green", "cyan", "orange", "white" ]
  
    xmin = minval(x); xmax = maxval(x) 
    ymax = maxval( PI_N(2,:) ); ymin = 0
    call scrmod("reverse") 
    call metafl("xwin")
   
    call page(4000, 4000)
    call disini 
    call graf(xmin, xmax, xmin, (xmax-xmin)/5, ymin, ymax, ymin, (ymax-ymin)/5 )
    
    j = 1 
    call color(col(j)); call incmrk(1);   call marker(21);  
    call curve( xp, PI_N(2,:), M+1 )
    
    j = 5 
    call color(col(j)); call incmrk(1);   call marker(21);  
    call curve( xp, PI_N(0,:), M+1 )
    
    call plot_legends(["d2PI_N/dx2(x)", "PI_N(x)" ])
    call plot_title(["PI_N(x) Error function","Extrema Chebyshev points" ]) 
    call disfin 

end subroutine 
end subroutine 
   
    
    
    
    
    
    
    
    

!***********************************************************
!* Error_exact_versus_estimated
!***********************************************************
subroutine Error_exact_versus_estimated  
 
 
 real, allocatable :: x(:), f(:), theta(:)   
 real, allocatable :: I_N(:, :), I2_N(:), c2_N(:), fe(:, :), Error(:, :) 
 integer ::  N, M = 50,  i
 integer, parameter :: Np = 30
 real :: xp(Np), yp(2, Np), a = 1  
 real :: PI = 4 * atan(1d0)
 
 yp = 0 
 
 do i=1, Np  
 
  N = 10+ i   
  
  allocate(theta(0:N), x(0:N), f(0:N), I_N(0:N, 0:N), fe(0:2, 0:N), Error(2, 0:N)  ) 
  allocate ( I2_N(0:N), c2_N(0:N)  ) 
  theta  = [ (PI*i/N, i=0, N) ]
  x = -cos( theta ) 
  
  f = sin ( a * PI * x ) 
  I_N = Interpolant(x, f, N, x) 
    
! Error  of the second derivative by definition 
  fe(2,:) =  -(a * PI)**2 * sin ( a * PI * x ) 
  Error(1,:) = fe(2,:) - I_N(2, :) ! Function, first and second derivative 
  
! Error estimate of the second derivative 
! by substracting two interpolants
  c2_N = Chebyshev_transform(  I_N(2, :) ) 
  I2_N = Chebyshev_interpolant( c2_N(0:N-4), x ) 
  Error(2,:) = I_N(2, :) - I2_N 
    
  xp(i) = N 
  yp(1, i) = norm2( Error(1,:) )  
  yp(2, i) = norm2( Error(2,:) )
  
  deallocate(theta, x, f, I_N, fe, Error, c2_N, I2_N) 
  
 end do 
 
 call plot_error
 
contains 


subroutine plot_error
  real :: xmin, xmax, ymin, ymax
  xmin = xp(1); xmax = xp(Np); 

  call metafl("xwin")
  call page(4000, 4000)
  call scrmod("reverse")
  
  call disini   
 
  call axsscl("log","y"); call labels("log","y"); call labdig(-1, "y"); call labdig(0, "x")
  call setscl(yp, Np, "y") 
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  
  call incmrk(1); call marker(21);
  call color("white");   call curve( xp, yp(1,:), Np)
  
  call incmrk(1); call marker(21);
  call color("red");   call curve( xp, yp(2,:), Np)
  
  call plot_legends(["E1", "E2"])  
  call plot_title([ "Error of the second derivative cos(pi x)  versus N", & 
                    "Comparison of exact error E1 with estimated error E2"] ) 
  call disfin
  
 
end subroutine 
end subroutine 







!***********************************************************
!* Error_second_derivative_Chebyshev_versus_FD
!***********************************************************
subroutine Error_second_derivative_Chebyshev_versus_FD 
 
 
 real, allocatable :: x(:), f(:), theta(:)   
 real, allocatable :: I_N(:, :), I2_N(:, :), c2_N(:), fe(:, :), Error(:, :) 
 integer ::  N, M = 50,  i, q 
 integer, parameter :: Np = 15
 real :: xp(Np), yp(4, Np), a = 1,  PI = 4 * atan(1d0)  
 
 yp = 0 
 do i=1, Np  
 
  N = 10+ 2*i   
  
  allocate(theta(0:N), x(0:N), f(0:N), I_N(0:N, 0:N), fe(0:2, 0:N), Error(4, 0:N)  ) 
  allocate ( I2_N(0:N, 0:N), c2_N(0:N)  ) 
  
! Cbebyshec. Error  of the second derivative by definition 
  theta  = [ (PI*i/N, i=0, N) ]
  x = -cos( theta ) 
  
  f = sin ( a * PI * x ) 
  I_N = Interpolant(x, f, N, x)
  fe(2,:) =  -(a * PI)**2 * sin ( a * PI * x ) 
  Error(1,:) = fe(2,:) - I_N(2, :) ! Function, first and second derivative 
  
  ! FD. Error  of the second derivative by definition 
  call Grid_Initialization( "nonuniform", "x",  x, 2 )
  
  f = sin ( a * PI * x )
  I2_N = Interpolant(x, f, N, x)
  fe(2,:) =  -(a * PI)**2 * sin ( a * PI * x ) 
  Error(2,:) = fe(2,:) - I2_N(2, :) 
  
  
! FD. Error  of the second derivative by definition 
  call Grid_Initialization( "nonuniform", "x",  x, 4 )
  
  f = sin ( a * PI * x )
  I2_N = Interpolant(x, f, N, x)
  fe(2,:) =  -(a * PI)**2 * sin ( a * PI * x ) 
  Error(3,:) = fe(2,:) - I2_N(2, :) 
  
  call Grid_Initialization( "nonuniform", "x",  x, 6 )
  
  f = sin ( a * PI * x )
  I2_N = Interpolant(x, f, N, x)
  fe(2,:) =  -(a * PI)**2 * sin ( a * PI * x ) 
  Error(4,:) = fe(2,:) - I2_N(2, :) 
  
    
  xp(i) = N 
  yp(1, i) = maxval( Error(1,:) )  
  yp(2, i) = maxval( Error(2,:) ) 
  yp(3, i) = maxval( Error(3,:) ) 
  yp(4, i) = maxval( Error(4,:) )  
  
  write(*,*) " N = ", N, " Error =", yp(:, i)   
  
  deallocate(theta, x, f, I_N, fe, Error, c2_N, I2_N) 
  
 end do 
 
 call plot_error
 
contains 


subroutine plot_error
  real :: xmin, xmax, ymin, ymax
  xmin = xp(1); xmax = xp(Np); 

  call metafl("xwin")
  call page(4000, 4000)
  call scrmod("reverse")
  
  call disini   
 
  call axsscl("log","y"); call labels("log","y"); call labdig(-1, "y"); call labdig(0, "x")
 ! call setscl(yp, Np, "y") 
  ymax = -3; ymin = -15  
  call graf(xmin, xmax, xmin, (xmax-xmin)/10, ymin, ymax, ymin, (ymax-ymin)/10) 
  
  call incmrk(1); call marker(21);
  call color("white");   call curve( xp, yp(1,:), Np)
  
  call incmrk(1); call marker(21);
  call color("red");   call curve( xp, yp(2,:), Np)
  
  call incmrk(1); call marker(21);
  call color("blue");   call curve( xp, yp(3,:), Np)
  
  call incmrk(1); call marker(21);
  call color("orange");   call curve( xp, yp(4,:), Np)
   
  
  call plot_legends( ["Chebyshev q = N", "FD q =2", "FD q =4", "FD q =6"] ) 
  call plot_title( [ "Error of the second derivative sin(pi x)",  & 
                     " Chebyshev versus Finite Differences" ]    ) 

  
  call disfin
  
  
 
 
end subroutine 

end subroutine







subroutine Fast_chebyshev_transform_example
 
 integer, parameter :: N=8
 real :: x(0:N), y(0:N), yr(0:N), c(0:N) 
 real ::  PI = 4 * atan(1.) 
 integer :: i, k  
 
 
   x  = [ (-cos(PI*i/N), i=0, N) ]
  
 
   y =  1 + x + x**2  + x**3 + x**4 + x**5 + x**6 + x**7 + x**8
  ! y = x**2
   
   call plot(x, y, "y(x)" )
   
   c = Chebyshev_transform( y )
   do k=0,  N 
     write(*,*) k, "c_k =", c(k) 
   end do
   
  write(*,*)
  write(*,*) "Fast Chebyshev transform " 
  c = Fast_Chebyshev_transform(N, y) 
   do k=0,  N 
     write(*,*) k, "c_k =", c(k) 
   end do
    
   yr = iChebyshev_transform( c )
   call plot(x, yr, "Inverse Chebyshev y(x)" )
   write(*,*) " maxval y-yr =", maxval( abs(y - yr) )  
   
   yr = iFast_Chebyshev_transform( N, c )
   call plot(x, yr, "Inverse Fast Chebyshev y(x)" )
   write(*,*) " maxval y-yr =", maxval( abs(y - yr) ) 
   
   
end subroutine 


!
!
!!*****************************************
!subroutine  test_dct2
!!*****************************************
!
!  integer, parameter :: N = 64
!  complex :: U(0:N) 
!  
!
!  integer :: i
!  real :: k(0:N)
!  complex ::  ck(0:N), C1(0:N)   
!  real :: x(0:N), y(0:N),  PI = 4 * atan(1d0)  
!  real :: a= -1, b= 1
!  
!  
!   x = [ ( cos(PI*i/N), i=N, 0, -1) ]  
!   U = cos( PI * X ) 
!!  call qplot( x, real(U), N+1 )
! 
!   call DCT2(N, U )
!
!   ck(0) =  U(0) / N 
!   ck(1:N) = 2 * U(1:N)  / N
!!  call qplot([(real(i),i=0,N)], real(ck), N+1)
!  
!   U = 2 * U / N 
!  
!   call DCT2(N, U)  
!
!  C1 =  Chebyshev_Spectral_Derivative1D(Ck)
!  call DCT2(N, C1) 
! 
!  
!  
!  call scrmod("reverse")
!  call metafl("xwin")
!  write (*, *) 'It plots cos pi x and its dervative by doing:  ' 
!  write (*, *) '1. Cosine transform DTC2 of cos( pi x  ) to obtain Chebyshev coefficients' 
!  write (*, *) '2. Derivative in the spectral plane ' 
!  write (*, *) '3. Inverse cosine transform DTC2 to obtain physical values of the derivative ' 
!  write(*,*) "press enter " 
!  read(*,*)
!  
!  call disini  
!   call graf(-1., 1., -1., 0.2, -4., 4., -4., 1.); 
!   call curve( x, real(u), N+1 )
!   call curve( x, real(C1), N+1 )
!  call disfin 
!  
! 
!  
!  
!  
!  
!end subroutine      
    

end module 
    
    
    
    