module Chebyshev_interpolation

 use Fourier_interpolation 
 implicit none 

 private 
 public :: DCT2, Chebyshev_Derivative 
contains 


  
!***************************************************************
! Discrete cosine Fourier Transform 
!***************************************************************
subroutine DCT2(N, U)
     integer, intent(in) :: N 
     complex, intent(inout)  :: U(0:N)
    
    complex   :: Uc(0:2*N-1)
     
    Uc(0:N) = U(0:N) 
    Uc(N+1:2*N-1) = U(N-1:1:-1) 
    
    call FFT(2*N, Uc)
  
    U(0:N) = Uc(0:N) / 2
 
  end subroutine 
  
!***************************************************************
! Derivative of Chebyshev expansion 
!***************************************************************
subroutine Chebyshev_Derivative(C, C1)
     complex, intent(in)  :: C(0:)
     complex, intent(out)  ::  C1(0:) 
    
    integer   :: k, N 
    N = size(C) - 1 
   
    C1(N) = 0 
    C1(N-1) = - 2 * N * C(N) 
    
    do k = N-1, 1, -1
       C1(k-1) = C1(k+1) - 2 * k * C(k)  !! instead of  + 2 * k * C(k) because x_0 = -1 
    end do
  
 
  end subroutine 


end module 
