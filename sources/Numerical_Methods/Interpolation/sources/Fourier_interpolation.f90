module Fourier_interpolation

 implicit none 
 real :: PI = 4 * atan(1d0)
 complex   :: II = (0, 1)
 
 
 private 
 public ::                  & 
     FFT,                   & ! Fast Fourier transform  
     FFT2D,                 & ! Fast Fourier transform  2D 
     FFT_1D,                & ! direct Fast Fourier Transform 1D
     FFT_2D,                & ! direct Fast Fourier Transform 2D
     iFFT_1D,               & ! inverse Fast Fourier Transform 1D
     iFFT_2D,               & ! inverse Fast Fourier Transform 2D
     Spectral_derivative1D, & ! Derivative 1D in the spectral plane 
     Spectral_derivative2D, & ! Derivative 1D in the spectral plane 
     Fourier_derivative1D,  & ! Fourier derivative 1D in the physical plane
     Fourier_derivative2D,  & ! Fourier derivative 1D in the physical plane
     Fourier_Grid_Initialization 
 
 
  real, save :: x0, xf, y0, yf, z0, zf ! independent domain 
  real, save :: fx = 1, fy = 1, fz = 1 ! dimension factor for derivatives 
 
    contains 
    
    
!*************************************************************************
! Cooley-Tukey FFT (Fast Fourier Transform)
!    
!        
!    1) Interpolant I(x) = sum_k( c_k exp(  I k x  ) )  
!    
!       from k= -N/2 to k = N/2 -1 forall x in [0, 2 PI ] 
!     
!    2) Conditions to obtain the interpolant     
!       u_j = sum_m( c_m exp(  2 PI m j/N I ) )  j=0,..N-1   
!    
!    3) Coefficients determination     
!       c_k = sum_j( u_j exp( -2 PI k j/N I ) )  sum from j=0,.. N-1 
!                                                forall k=-N/2.... N/2-1  
!    4) Cooley-Tukey algorithm     
!       c_k = Even_j + Odd_j     
!   
!    It can also be used for the inverse transform because
!    
!    conjugate(u_j) =  sum ( conjugate(c_m) exp( - 2 * PI m j/N I ) ) 
!    
!    which is the same expression used to calculate c_m 
!    by changing the input and the output with their conjugates.  
!        
!
!    Inputs:
!              - N        : number of collocation points(must be N=2^r)
!              - U(0:N-1) : collocation points 
!                           u_0, u_1,.. u_N-1
!
!    Outputs:
!              - U(0:N-1) : discrete fourier coefficients 
!                           N *(c_0, c_1, ... c_N/2, c_N/2-1, c_{-N/2} ..,c_{-3},  c_{-2}, c_{-1} )
!
!  If u_i is real, then c_k    = a_k + I b_k 
!                       c_{-k} = a_k - I b_k 
!    
!  If u_i is complex (u = u1 + I u2 ), then 
!                       c1_k =   1/2  ( c_k + conjg(c_{-k} ) )
!                       c2_k = - I/2  ( c_k - conjg(c_{-k} ) )    
!
!  Author: Juan A. Hernandez, Feb, 2018, Dec. 2020, May 2021    
!*************************************************************************
recursive subroutine FFT(N, U)
     integer, intent(in) :: N 
     complex, intent(inout)  :: U(0:N-1)
    
    complex   :: w, Even(0:N/2-1), Odd( 0:N/2-1 )
    integer   :: k 
     
    if(N <= 1) return
    
    Even = U(0:N-1:2)
    Odd  = U(1:N-1:2)
    
    call FFT( N/2, Even)
    call FFT( N/2, Odd)
   
    do k = 0, N/2-1
        
       w = exp( -2*PI*II * k/N )
       U(k)     = Even(k) + w * Odd(k)
       U(k+N/2) = Even(k) - w * Odd(k)  
       
     ! Note that  U(N/2)   ... U(N-1) 
     ! are:       C_{-N/2} ....C_{-1} 
     ! using periodic property of C_{k-N/2} 
    end do
 
end subroutine  

  

!****************************************************************
!   Fourier C_k coefficients from U_j physical values 
!****************************************************************
function FFT_1D(N, U) result(C) 
     integer, intent(in) :: N 
     complex, intent(in)  :: U(0:N-1)
     complex :: C(-N/2:N/2-1)
    
    
    complex :: Us(0:N-1) 
    
    Us = U 
    call FFT(N, Us) 
    
    C(0:N/2-1) = Us(0:N/2-1) 
    C(-N/2:-1) = Us(N/2:N-1) 
    C = C / N 
       
 end function 
 
 
!****************************************************************
!   U_j physical values from C_k  Fourier coefficients
!****************************************************************  
function iFFT_1D(N, C) result(U) 
     integer, intent(in) :: N  
     complex, intent(in) :: C(-N/2:N/2-1)
     complex  :: U(0:N-1)
    
     
    complex :: Cs(0:N-1) 
       
    Cs(0:N/2-1) = C(0:N/2-1) 
    Cs(N/2:N-1) = C(-N/2:-1)   
    
    Cs = conjg( Cs ) 
    call FFT(N, Cs) 
    Cs = conjg( Cs) 
    
    U = Cs
    
 
end function 



!*********************************************************************
! Derivative of Fourier expansion 
!  Inputs: 
!          c_m coefficients 
!          k   order of the derivative 
! Output: 
!         Dc_m coefficients of the kth derivative   
!  
!  I(x) = sum_m( c_m exp(  I m x  ) )   
!  
!  dkI/dxk (x) =  sum_k ( (I m)**k c_k exp( I m x ) )  
!  
!  sum from m= -N/2 to m = N/2 -1 forall x in [0, 2 PI ] 
!   
!   Author: Juan A. Hernandez, Nov., 2020, juanantonio.hernandez@upm.es    
!***********************************************************************
function Spectral_derivative1D(N, C, k) result(DC)
     integer, intent(in) :: N 
     complex, intent(in)  :: C(-N/2:N/2-1)
     integer, intent(in) :: k 
     complex  ::  DC(-N/2:N/2-1) 
    
    integer   :: m
   
    do m = -N/2, N/2-1 
       DC(m) =    ( II * m * fx )**k * C(m) 
    end do
    
end function 

!******************************************************************************
! Derivative1D by means of a Fourier expansions 
! 
!  I(x) = sum_m( c_m exp(  I m x  ) ) 
!
! Algorithm : 
!            1) given a set of colloction points    U = [ u_0, .... u_{n-1} ]
!            2) obtain its Fourier cofficients  C_k k=-n/2... N/2-1  
!            3) derive in the spectral plane: new coefficients: (I m)**k c_m 
!            4) transform to the physical space 
!   
!   Author: Juan A. Hernandez, Nov., 2020, juanantonio.hernandez@upm.es
!******************************************************************************
subroutine Fourier_Derivative1D( direction, derivative_order, U, Uxi )
   character(len=*), intent(in) :: direction
   integer, intent(in) :: derivative_order 
   real, intent(in) ::   U(0:)
   real, intent(out)::   Uxi(0:) 
   
     integer :: N
     complex :: Uc( size(U) ), C( size(U) ) 
     
   N = size(U) ! warning j=0,... N-1 (Fourier transform with N points= 2**r) 
   Uc = U 
   C = FFT_1D(N, Uc)
   C = Spectral_derivative1D(N, C, derivative_order)
   Uc = iFFT_1D(N, C)  
   
   Uxi  =  real(Uc) 
  
     
end subroutine 
 







   
!*************************************************************************
! FFT (Fast Fourier Transform) 2D 
!    
!        
!    1) Interpolant I(x,y) = sum_km( c_kl exp(  I k x + I m y  ) )  
!     
!    2) Conditions to obtain the interpolant     
!       u_ij = sum_km( c_km exp(  2 PI k i/N I ) exp(  2 PI m j/N I ))  i=0,...N-1, j=0,..M-1   
!    
!    3) Coefficients determination     
!
!       c_km = 1/(N M) sum_j  exp( -2 PI k j/N I ) [  sum_i( u_ij exp( -2 PI k i/N I ) ]  
!
!  Author: Juan A. Hernandez, Dec., 2020   
!*************************************************************************
recursive subroutine FFT2D(N, M, U)
     integer, intent(in) :: N, M  
     complex, intent(inout)  :: U(0:N-1, 0:M-1)
    
     integer :: i, j 
     
    do j=0, M-1
        call FFT(N, U(:,j) ) 
    end do 
    
    do i=0, N-1
        call FFT(M, U(i,:) ) 
    end do  
 
end subroutine      

!****************************************************************
!   Fourier C_km coefficients from U_ij physical values 
!****************************************************************
function FFT_2D(N, M, U) result(C) 
     integer, intent(in) :: N, M  
     complex, intent(in)  :: U(0:N-1, 0:M-1)
     complex :: C(-N/2:N/2-1, -M/2:M/2-1)
    
    
    complex :: Us(0:N-1, 0:M-1) 
    
    !write(*,*) " N = ", N , " M = ", M 
    Us = U
    call FFT2D(N, M, Us)
    
    C(0:N/2-1, 0:M/2-1) = Us(0:N/2-1, 0:M/2-1)
    C(-N/2:-1, 0:M/2-1) = Us(N/2:N-1, 0:M/2-1) 
    
    C(0:N/2-1, -M/2:-1) = Us(0:N/2-1, M/2:M-1) 
    C(-N/2:-1, -M/2:-1) = Us(N/2:N-1, M/2:M-1) 
    
    
    C = C / ( N * M ) 
    
       
end function 

!****************************************************************
!   U_j physical values from C_k  Fourier coefficients
!****************************************************************  
function iFFT_2D(N, M, C) result(U) 
     integer, intent(in) :: N, M   
     complex, target, intent(in) :: C(-N/2:N/2-1, -M/2:M/2-1)
     complex  :: U(0:N-1, 0:M-1)
    
     
    complex :: Cs(0:N-1, 0:M-1)
    
     Cs(0:N/2-1, 0:M/2-1) = C(0:N/2-1, 0:M/2-1)
     Cs(N/2:N-1, 0:M/2-1) = C(-N/2:-1, 0:M/2-1)
     
     Cs(0:N/2-1, M/2:M-1) = C(0:N/2-1, -M/2:-1)
     Cs(N/2:N-1, M/2:M-1) = C(-N/2:-1, -M/2:-1)  
    
    
    Cs = conjg( Cs )
    call FFT2D(N, M, Cs) 
    Cs = conjg(Cs) 
    
    U = Cs
    
 
end function 


function Spectral_derivative2D(N, M, C, order, coordinate) result(DC)
     integer, intent(in) :: N, M 
     complex, intent(in)  :: C(-N/2:N/2-1, -M/2:M/2-1)
     integer, intent(in) :: order, coordinate  
     complex  ::  DC(-N/2:N/2-1, -M/2:M/2-1) 
    
    integer   :: k, l
   
    if (coordinate == 1) then 
        
      do k = -N/2, N/2-1 
         DC(k, :) =  ( II * k * fx )**order * C(k, :) 
      end do
      
    else if (coordinate == 2) then 
        
      do l = -M/2, M/2-1 
         DC(:, l) =  ( II * l * fy )**order * C(:, l) 
      end do 
      
    else 
        write(*,*) " Error in Spectral_derivative2D "
        write(*,*) " wrong coordinate = ", coordinate 
        stop 
    end if 
    
    
end function  
!******************************************************************************
! Derivative2D by means of a Fourier expansions 
! 
!  I(x,y)          = sum_km( c_{km} exp(  I k x ) exp(  I m y ) ) 
!  dI/dx(x_i, x_j) = sum_km( I k c_{km} exp(  I k x_i ) exp(  I m y_j ) ) 
!  
! Algorithm : 
!            1) given a set of colloction points    U_{ij} = [ u_00, .... u_{n-1,m-1} ]
!            2) obtain its Fourier cofficients  C_km 
!            3) derive in the spectral plane: new coefficients: (I k)**order c_km 
!            4) transform to the physical space 
!   
!   Author: Juan A. Hernandez, Dec., 2020, juanantonio.hernandez@upm.es
!******************************************************************************
 subroutine Fourier_Derivative2D( direction, coordinate, derivative_order, U, Ux ) 

   character(len=*), intent(in) :: direction(1:2)
   integer, intent(in) :: coordinate, derivative_order
   real, intent(in) ::   U(0:, 0:)
   real, intent(out)::   Ux(0:, 0:) 
       
     integer :: N, M 
     complex, allocatable :: Uc(:, :), C(:, :) 
     
   ! Fourier transform with i=0,..N-1, j=0...M-1 with N=2**r1, M=2**r2
     N = size(U, dim=1);  M = size(U, dim=2); 
    
     allocate(  Uc(0:N-1, 0:M-1), C(-N/2:N/2-1, -M/2:M/2-1) ) 
    
     Uc = U 
     C = FFT_2D(N, M, Uc) 
     C = Spectral_derivative2D(N, M, C, derivative_order, coordinate)
     Uc = iFFT_2D(N, M, C) 
    
     Ux  =  real(Uc)   
   
 end subroutine


!***********************************************************************************
!*  It computes dimensional factors for Fourier derivatives 
!*
! Authors : Juan A Hernandez (juanantonio.hernandez@upm.es)  
!***********************************************************************************
subroutine Fourier_Grid_Initialization( direction, nodes ) 
      character(len=*),  intent(in) ::  direction 
      real, intent(inout) :: nodes(0:)
       
            
       
  real :: a , b;  
  integer :: N, i

   N = size(nodes) 
   
   if (direction == "x") then 
       x0 = nodes(0); xf = nodes(N-1) * N / real(N-1)
       fx = 2 * PI / ( xf - x0 ) 
       
    elseif (direction == "y") then 
       y0 = nodes(0); yf = nodes(N-1) * N / real(N-1)
       fy = 2 * PI / ( yf - y0 ) 
       
    elseif (direction == "z") then 
       z0 = nodes(0); zf = nodes(N-1) * N / real(N-1) 
       fx = 2 * PI / ( zf - z0 ) 
       
    else 
        write(*,*) "ERROR: Fourier_Grid_Initialization "
        stop 
    end if 
    
    a = nodes(0); b = nodes(N-1) 
    do i = 0, N-1
             nodes(i) = a  + (b-a) *  i / real(N-1)
    enddo
    
end subroutine 








 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 !********************************************************
 ! Testing subroutines and old subroutines 
 !********************************************************

!*********************************************************
! Direct Fourier transform (no fast) 
!  c_k = sum_j( u_j exp(  -I k x_j  ) ) 
!        x_j = 2 PI j / N 
!*********************************************************
function DFT(N, U) result(C) 
     integer, intent(in) :: N 
     complex, intent(in)  :: U(0:N-1)
     complex :: C(-N/2:N/2-1)
    
    integer   :: k, j  
    complex :: S 
    
    do k = -N/2, N/2-1
       S = 0  
       do j=0, N-1 
         S = S + U(j) * exp( -2*PI*II * j * k/real(N) )
       end do 
       C(k)  = S / N 
    end do
   
    
end function    
    
!*********************************************************
! Inverse Fourier transform (no fast) 
!  I(x_j) = sum_k( c_k exp(  I k x_j  ) ) 
!           x_j = 2 PI j / N 
!*********************************************************
function IFT(N, C) result(U) 
     integer, intent(in) :: N  
     complex, intent(in) :: C(-N/2:N/2-1)
     complex  :: U(0:N-1)
    
    integer   :: k, j  
    complex :: S 
   
    do j = 0, N-1
       S = 0  
       do k=-N/2, N/2-1 
         S = S + C(k) * exp( 2*PI*II * j * k/real(N) )
       end do 
       U(j) = S 
    end do
   
 
end function  

!*******************************************************************
!*             Fast Fourier Transform
!*
!*    Inputs:
!*              - ip       : 1 FFT to spectral, 2 FFT to physical
!*              - N        : number of collocation points(must be N=2^r)
!*              - u(0:N-1) : collocation points (ip=1)
!*                           u_0, u_1,.. u_N-1
!*                         : discrete fourier coefficients (ip=2)
!*                           N *(c_0, c_1, ... c_N/2-1, c_-1, .. c_N/2 )
!*              - w(0:N,3) : workspace
!*
!*    Outputs:
!*             
!*             - w(0:N-1,3): collocation points (ip=1)
!*                           u_0, u_1,.. u_N-1
!*                         : discrete fourier coefficients (ip=2)
!*                           N *(c_0, c_1, ... c_N/2-1, c_-1, .. c_N/2 )
!*
!*  Author: Juan A. Hernandez, Feb, 1995
!*
!*******************************************************************
       subroutine FFT_(ip, N, u, w)
              integer ip
              integer N
              complex*16 u(0:*), w(0:N,*)



            integer j, p, ipass 

            complex*16 II
            real*8 pi
            data II /(0d0, 1d0)/, ipass/0/


!         *** Inicializacion
              if (ipass.eq.0) then

                    pi = 4d0 * atan(1d0)
                    do j=0, N
                        w(j,1) = exp( - II*2d0*pi*j/N )
                    enddo

                    do j=0, N
                        w(j,2) = exp(  II*2d0*pi*j/N )
                    enddo

              endif
              ipass = ipass + 1



!         ***  Ordenacion
               p = N
               do while(p.gt.2)
                  do j=0, N-1, p
                    call even_odd(u(j), w(j,3), p)
                  enddo
                  p = p / 2
               enddo


!         ***  Fourier
               p = 2
               do while (p.le.N)

                     do j=0, N-1, p
                       call Fourier(p, N, w(j,3), w(0,ip) )
                     enddo

                     p = p * 2

                enddo



       end subroutine  
!*************************************************************
!*         Danielson - Lanczos formula
!*
!*    Input:
!*            - v(0:N-1) : v_0...v_N/2 FFT of even terms
!*                         v_N/2+1... v_N-1 FFT of odd terms
!*            - w_j      : exp( - i k x_j )
!*            - N        : collocation points of this sub FFT
!*            - Nmax     : collocation points of the first FFT
!*
!*    Outputs:
!*            - v(0:N-1) : FFT coefficients
!*
!*    Author :  Juan A. Hernandez, Feb. 1995
!**************************************************************
      subroutine Fourier(N, Nmax, v, w)
                 integer N, Nmax
                 complex*16 v(0:*), w(0:*)


          integer k,q
          complex*16 a, b



             do k=0, N/2 - 1

                q = Nmax/N*(N-k)

                a = v(k)
                b = w(q) * v(k+N/2)

                v(k)     = a + b
                v(k+N/2) = a - b

             enddo



      end subroutine  
!**************************************************************
!*             Even and odd terms
!*    Reorders an array, first even terms and then odd terms
!*
!*    Inputs:
!*             - u(0:N-1) : input array
!*             - N        : number of elements
!*    Outputs:
!*             - v(0:N-1) : reordered array
!
!*   Author: Juan A. Hernandez, Feb 1995
!***************************************************************
      subroutine even_odd(u, v, N)
               complex*16 u(0:*), v(0:*)
               integer N

           integer k



          do k=0, N/2-1

             v(k)     = u(2*k)
             v(k+N/2) = u(2*k+1)

          enddo

        
      end  subroutine   
    
      
end module 
