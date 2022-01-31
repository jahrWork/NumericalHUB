module API_Example_Chebyshev_series

    use Chebyshev_interpolation
    use Collocation_methods 
    use plots 
    
    implicit none
    real :: PI = 4 * atan(1d0)
    
    
    contains
    
    
subroutine Chebyshev_examples

 !   call test_fft2
 !   call test_dct2
   
  !  all Test_fft2
       !call Test_FFTD
       !stop 
       !call Test_Fourier
       !stop 
       !call Fourier_Derivative_error
       !stop 
  !     call Fourier_IBVP_examples
  !     stop 
    
    
end subroutine  



!*****************************************
subroutine  test_dct2
!*****************************************

  integer, parameter :: N = 64
  complex :: U(0:N) 
  

  integer :: i
  real :: k(0:N)
  complex ::  ck(0:N), C1(0:N)   
  real :: x(0:N), y(0:N)  
  real :: a= -1, b= 1
  
  
   x = [ ( cos(PI*i/N), i=N, 0, -1) ]  
   U = cos( PI * X ) 
!  call qplot( x, real(U), N+1 )
 
   call DCT2(N, U )

   ck(0) =  U(0) / N 
   ck(1:N) = 2 * U(1:N)  / N
!  call qplot([(real(i),i=0,N)], real(ck), N+1)
  
   U = 2 * U / N 
  
   call DCT2(N, U)  

  call Chebyshev_Derivative(Ck, C1)
  call DCT2(N, C1) 
 
  
  
  call scrmod("reverse")
  call metafl("xwin")
  write (*, *) 'It plots cos pi x and its dervative by doing:  ' 
  write (*, *) '1. Cosine transform DTC2 of cos( pi x  ) to obtain Chebyshev coefficients' 
  write (*, *) '2. Derivative in the spectral plane ' 
  write (*, *) '3. Inverse cosine transform DTC2 to obtain physical values of the derivative ' 
  write(*,*) "press enter " 
  read(*,*)
  
  call disini  
   call graf(-1., 1., -1., 0.2, -4., 4., -4., 1.); 
   call curve( x, real(u), N+1 )
   call curve( x, real(C1), N+1 )
  call disfin 
  
 
  
  
  
  
end subroutine      
    

end module 
    
    
    
    