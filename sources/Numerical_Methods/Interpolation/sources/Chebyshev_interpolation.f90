module Chebyshev_interpolation

 use Fourier_interpolation 
 use plots
 implicit none 

 !private 
 !public :: DCT2
 
 private 
 public ::                  & 
   Chebyshev_transform,             & ! Chebyshev coefficients c_k
   iChebyshev_transform,            & ! Inverse Chebyshev  f(x_j)  
   Chebyshev_Spectral_Derivative1D, & ! Derivative 1D in the spectral plane 
   Chebyshev_interpolant,           & ! I(x) =  sum( c_k T_k(x) ) 
   Chebyshev_Grid_Initialization,   & 
   Chebyshev_Derivative1D,          & 
   Fast_Chebyshev_transform,        & 
   iFast_Chebyshev_transform,       & 
   Chebyshev_extrema
 
 
  real, save :: x0, xf, y0, yf, z0, zf ! independent domain 
  real, save :: fx = 1, fy = 1, fz = 1 ! dimension factor for derivatives 
 
 
 
 
contains 


 
!***********************************************************************************
!*  It computes dimensional factors for Chebyshev derivatives 
!*
! Authors : Juan A Hernandez (juanantonio.hernandez@upm.es)  
!***********************************************************************************
subroutine Chebyshev_Grid_Initialization( direction, nodes ) 
      character(len=*),  intent(in) ::  direction 
      real, intent(inout) :: nodes(0:)
        
  integer :: N, i
  real, parameter :: PI = 4 * atan(1.) 

   N = size(nodes)-1 
   
   if (direction == "x") then 
       x0 = nodes(0); xf = nodes(N) 
       fx = 2  / ( xf - x0 ) 
       
       do i = 0, N
          nodes(i) = (x0+xf)/2 + (x0-xf)/2 *cos( i*PI/N) 
       enddo
       
    elseif (direction == "y") then 
       y0 = nodes(0); yf = nodes(N) 
       fy = 2  / ( yf - y0 ) 
       
       do i = 0, N
          nodes(i) = (y0+yf)/2 + (y0-yf)/2 *cos( i*PI/N) 
       enddo
       
    elseif (direction == "z") then 
       z0 = nodes(0); zf = nodes(N)  
       fz = 2  / ( zf - z0 ) 
       
       do i = 0, N
          nodes(i) = (z0+zf)/2 + (z0-zf)/2 *cos( i*PI/N) 
       enddo
       
    else 
        write(*,*) "ERROR: Chebyshev_Grid_Initialization "
        stop 
    end if 
    
  
    
   
    
end subroutine    
    

function Chebyshev_extrema( N, x0, xf ) result(x) 
     integer,  intent(in) ::  N
     real,  intent(in) ::  x0, xf 
     real :: x(0:N)
     
   real, parameter :: PI = 4 * atan(1.) 
   integer :: i 

  
    do i = 0, N
          x(i) = (x0+xf)/2 + (x0-xf)/2 *cos( i*PI/N) 
    enddo
   
    
end function    




!******************************************************************************
! Derivative1D by means of a Chebyshev expansion 
! 
!  I(x) = sum_m( c_k T_k(x) ) 
!
! Algorithm : 
!            1) given a set of colloction points    U = [ u_0, .... u_{n-1} ]
!            2) obtain its Chebyshev coefficients  c_k   k=0... N 
!            3) derive in the spectral plane: new coefficients: c1_k
!            4) transform to the physical space 
!   
!   Author: Juan A. Hernandez, Nov., 2022, juanantonio.hernandez@upm.es
!******************************************************************************
subroutine Chebyshev_Derivative1D( direction, derivative_order, U, Uxi )
   character(len=*), intent(in) :: direction
   integer, intent(in) :: derivative_order 
   real, intent(in) ::   U(0:)
   real, intent(out)::   Uxi(0:) 
   
     integer :: i, N 
     real ::  c(size(U))
     real, parameter :: PI = 4 * atan(1.) 
     N = size(U)-1 
          
   c = Chebyshev_transform( U )
 !  c = Fast_Chebyshev_transform(N, U) 
         
   do i=1, derivative_order
     c = Chebyshev_Spectral_Derivative1D( c )
     c = fx * c 
   end do 
  ! write(*,*) " c =", c 
   
   Uxi = iChebyshev_transform( c )  
!  Uxi = iFast_Chebyshev_transform(N, c) 
  
 !  write(*,*) " Uxi =", Uxi 
   
end subroutine 


!***************************************************************
! Derivative of Chebyshev expansion 
! Author: juanantonio.hernandez@upm.es, April, 2022.
!***************************************************************
function Chebyshev_Spectral_Derivative1D_old(C) result(C1) 
     real, intent(in)  :: C(0:)
     real  ::  C1(0:size(C)-1) 
    
    integer   :: k, N 
    N = size(C) - 1 
   
    C1(N) = 0 
    C1(N-1) = 2 * N * C(N)
    
     do k = N-1, 2, -1
       C1(k-1) = C1(k+1) + 2 * k * C(k)  
    end do
    C1(0) = C1(2) / 2 + C(1) 
    
end function 

!***************************************************************
! Derivative of Chebyshev expansion 
! Author: juanantonio.hernandez@upm.es, April, 2022.
!***************************************************************
function Chebyshev_Spectral_Derivative1D(C) result(C1) 
     real, intent(in)  :: C(0:)
     real  ::  C1(0:size(C)-1) 
    
    integer   :: k, N 
    N = size(C) - 1 
   
    C1(N) = 0 
    C1(N-1) = 2 * N * C(N)
    
     do k = N-1, 2, -1
       C1(k-1) = C1(k+1) + 2 * k * C(k)  
    end do
    C1(0) =  C(1) + C1(2) / 2 
    
end function 



!****************************************************
! Chebyshev expansion 
!   f(x) = sum _0 ^N  c_k T_k(x)  
!   with T_k(x) = cos( k theta ) and x = cos( theta ) 
!****************************************************
function iChebyshev_transform( c ) result(P) 
  real, intent(in) :: c(0:)
  real :: P(0:size(c)-1) 
  
  integer :: k, j, N
  real, allocatable :: x(:)
  real, parameter :: PI = 4 * atan(1.)
  
  N = size(c)-1 
  allocate( x(0:N) ) 
  x = [(-cos( PI *j / N), j=0, N) ]
  
  P = 0
  do j=0, N
    do k=0, N 
         P(j) = P(j) +  c(k) * Chebyshev(k, x(j) ) 
    end do 
  end do 
  
end function 



!*********************************************************
! Inverse Chebyshev expansion 
!   c_k = 1/gamma_k sum _0 ^N  alpha_j f(x_j)   T_k(x_j)  
!   with T_k(x) = cos( k theta ) and x = cos( theta ) 
!*********************************************************
function Chebyshev_transform( f ) result(c) 
  real, intent(in) ::  f(0:)
  real :: c(0:size(f)-1) 
  
  real, allocatable :: x(:), alpha(:), gamm(:)
  real, parameter :: PI = 4 * atan(1.) 
  integer :: k, j, N 
  N = size(f) -1 
    
  allocate( x(0:N), alpha(0:N), gamm(0:N) ) 
  
  alpha(0) = PI/(2*N);  alpha(1:N-1) = PI/N; alpha(N) = PI/(2*N)  
  gamm(0) = PI ;gamm(1:N-1) = PI/2;   gamm(N) = PI
  x = [(-cos( PI *j / N), j=0, N) ]  
    
  c = 0 
  do k=0, N 
    do j=0, N
       c(k) = c(k) + alpha(j) * f(j) * Chebyshev(k, x(j))  
   end do  
   c(k) = c(k) / gamm(k)
  end do 
  
end function 


real function Chebyshev( k, x ) result(Tk) 
  integer, intent(in) :: k 
  real, intent(in) :: x
  
    Tk = cos( k * acos(x) ) 
  
end function 



!****************************************************
! Chebyshev expansion 
!   f(x) = sum _0 ^N  c_k T_k(x)  
!   with T_k(x) = cos( k theta ) and x = cos( theta ) 
!****************************************************
function Chebyshev_interpolant( c, x ) result(P) 
  real, intent(in) :: c(0:), x(0:) 
  real :: P(0:size(x)-1) 
  
  integer :: k, j, N, M 
  real, parameter :: PI = 4 * atan(1.) 
  
  N = size(c) -1 ; M = size(x) - 1 
  P = 0
  do j=0, M
    do k=0, N 
         P(j) = P(j) +  c(k) * Chebyshev(k, x(j) ) 
    end do 
  end do 
  
end function 







!***************************************************************
! Discrete cosine Fourier Transform 
!***************************************************************
function Fast_Chebyshev_transform(N, U) result(c) 
     integer, intent(in) :: N 
     real, intent(in)  :: U(0:N)
     real  :: c(0:N)
    
    complex   :: Uf(0:2*N), cf(-N:N)
    integer :: j, k  
!    real :: x(0:2*N)
    real, parameter :: PI = 4 * atan(1.) 
     
    Uf(0:N) = U(0:N) 
    
    
!    x(0:N) = [ (- cos(PI*j/N), j=0, N ) ] 
    
    do j=0, N
        Uf(2*N-j) = Uf(j) 
 !       x(2*N-j) = 2 + x(j) 
    end do 
    
    !do j=0, 2*N 
    !    write(*,*) j, " Uf =", Uf(j)  
    !end do
    
  !  call plot(x, real(Uf), "extended")
    
    ! Fourier transform 
    !Cf(-N:N-1) = DFT(2*N, Uf) 
    !cf(N) = cf(-N) 
    
    ! Fast Fourier transform 
    call FFT(2*N, Uf)
    Uf(2*N) = Uf(0)  
    
    ! do j=0, 2*N 
    !    write(*,*) j, " coef Uf =", Uf(j) /(2*N)   
    !end do
    
    
    
    cf(0:N) = Uf(0:N) 
   ! cf(-N:-1) = Uf(N:2*N) 
  !  cf(-N+1:-1) = Uf(N:2*N-1) 
    
    
    do k=-1, -N, -1
        cf(k) = Uf(2*N+k) 
    end do 
    
    
       
   
    
    cf = cf / (2*N) 
    
    !do k=-N, N 
    !    write(*,*) k, " cf =", cf(k) 
    !end do
    
    !stop 
    !
    
    !WARNING 
    c(0) =  real( cf(0) )
    do k=1, N
      c(k) =  (-1)**k * ( cf(k) + conjg( cf(-k) ) )
    end do 
    c(N) = c(N)/2 
    
    
       
end function


!***************************************************************
!  U_j physical values from C_k  Chebyshev coefficients
!***************************************************************
function iFast_Chebyshev_transform(N, C) result(U) 
     integer, intent(in) :: N 
     real  :: c(0:N)
     real  :: U(0:N)
    
    complex   ::  Cf(0:2*N-1)
    integer :: k, j 
    
 !   do j=0, N
 !       Uf(2*N-j) = Uf(j) 
 !!       x(2*N-j) = 2 + x(j) 
 !   end do 
    
    
    Cf(0:N-1) = c(0:N-1)
    
     do k=N, 2*N-1
        cf(k) = c(2*N-k) 
     end do 
    
    !do k=0, 2*N-1 
    !    write(*,*) k, " cf =", cf(k) 
    !end do
    
    
     
   ! Cf(N+1:2*N-1) = c(N-1:1:-1)
     
    
    Cf = conjg( Cf ) /2
    Cf(0) = 2 * cf(0) 
    Cf(N) = 2 * cf(N) 
    
    
    call FFT(2*N, Cf) 
    
    Cf = conjg( Cf )
    
    !do j=0, 2*N-1 
    !    write(*,*) j, " j--- cf =", cf(j) 
    !end do
    
    do j=N, 0, -1
        U(N-j) = cf(j)
    end do
  !  
 !   U = real( Cf(0:N) ) 
       
end function



end module 
