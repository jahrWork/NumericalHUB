module API_Example_Chebyshev_Fourier_interpolation

    use Chebyshev_interpolation
    use Collocation_methods 
    use Interpolation
    use Lagrange_interpolation
    use plots 
    use Utilities
    
    
    implicit none
    
    contains
    
subroutine Chebyshev_Fourier_interpolation_examples

    write (*,*) "**********************************************************" 
    write (*,*) " Chebyshev-Fourier interpolation examples" 
    write (*,*) "**********************************************************"  
   
    !call Chebyshev_transform_example
   
    
  !  call Test_Fourier_Derivative2D 
    !call Test_Fourier_Chebyshev_Derivative3D
    
     
    call Test_Advection_diffusion_Fourier_Chebyshev
    
    
end subroutine      
    
 
!*****************************************
subroutine  Test_Fourier_Derivative2D 
!*****************************************

  integer, parameter :: N = 32, M = 16  
  real  :: U(0:N-1, 0:M-1), Ux(0:N-1, 0:M-1), Uy(0:N-1, 0:M-1)
  real :: x(0:N-1), y(0:M-1)  
  integer, parameter :: N_levels = 9 
  real :: levels(0:N_levels) 
  real :: z0, zf,  PI = 4 * atan(1.) 
  integer :: i, j 
  
   x = [ ( 2*PI*i/N, i=0, N-1 ) ]  
   y = [ ( 2*PI*j/M, j=0, M-1 ) ] 
   levels = [ ( z0 + (zf-z0) * j/real(N_levels)  , j=0, N_levels )]
   call Grid_Initialization( "Fourier", "x",  x )
   call Grid_Initialization( "Fourier", "y",  y )
  
  
   U =  Tensor_product( sin(x), sin(y) )  
   call scrmod("reverse")
   call plot_contour(x, y, U, "x","y", levels, graph_type ="isolines") 
   
   ! partial u / partial x, Derivative( direction, coordinate, derivative_order, W, Wxi)
   call Derivative( [ "x", "y" ], 1, 1, U, Ux ) 
   call plot_contour(x, y, Ux, "x","y", levels, graph_type ="isolines")
   
   ! partial u / partial y 
   call Derivative( [ "x", "y" ], 2, 1, U, Uy ) 
   call plot_contour(x, y, Uy, "x","y", levels, graph_type ="isolines")
  
   
end subroutine 


!**************************************************
subroutine  Test_Fourier_Chebyshev_Derivative3D 
!**************************************************

  integer, parameter :: N = 32, M = 32, L = 10  
  real  :: U(0:N-1, 0:M-1, 0:L), dU(0:N-1, 0:M-1, 0:L)
  complex :: C(-N/2:N/2-1, -M/2:M/2-1, 0:L), Uc(0:N-1, 0:M-1, 0:L) 
  real :: x(0:N-1), y(0:M-1), z(0:L)  
  integer, parameter :: N_levels = 9 
  real :: levels(0:N_levels) 
  real :: z0, zf,  PI = 4 * atan(1.) 
  integer :: i, j, k 
  
   x = [ ( 2*PI*i/N, i=0, N-1 ) ] 
   y = [ ( 2*PI*j/M, j=0, M-1 ) ]  
   z = [ (-1 + 2*k/real(L), k=0, L ) ] 
   levels = [ ( z0 + (zf-z0) * j/real(N_levels)  , j=0, N_levels )]
   call Grid_Initialization( "Fourier", "x",  x )
   call Grid_Initialization( "Fourier", "y",  y )
   call Grid_Initialization( "Chebyshev", "z", z )
  
   U =  Tensor_product_R3( sin(x), sin(y), z**2  )  
   write(*,*) " shape U = ", shape(U)
   
   k = L/5  
   call scrmod("reverse")
   call plot_contour(x, y, U(:,:,k), "x","y", levels, graph_type ="isolines") 
   
   ! partial u / partial x
   dU = Derivative_3D( "x", 1, U ) 
   call plot_contour(x, y, dU(:,:,k), "x","y", levels, graph_type ="isolines")
!  
      
   ! partial u / partial y 
   dU = Derivative_3D( "y", 1, U ) 
   call plot_contour(x, y, dU(:,:,k), "x","y", levels, graph_type ="isolines")
     
   ! partial u / partial z 
   dU = Derivative_3D( "z", 1, U ) 
   call plot_contour(y, z, dU(k,:,:), "y","z", levels, graph_type ="isolines")
     

   C =  Spectral_transform( N, M, L, U )
   
   do i=-N/2, N/2-1; do j=-M/2, M/2-1; do k=0, L
    if (abs(C(i,j,k)) > 0.0001 ) then    
        write(*,*) " n=", i, "m=", j, "k=", k 
        write(*,*) "C = ", C(i,j,k) 
    end if 
   end do; end do; end do 
   
   k = L/5
   Uc = iSpectral_transform( N, M, L, C )
   write(*,*) "Plot inverse "
   call plot_contour(x, y, real(Uc(:,:,k)), "x","y", levels, graph_type ="isolines")
   
end subroutine 





subroutine Chebyshev_transform_example
 
 integer, parameter :: N=8
 real :: x(0:N), y(0:N), yx(0:N), yxx(0:N), yr(0:N), c(0:N) 
 real ::  PI = 4 * atan(1.) 
 integer :: i, k  
 
   x = [ (-1 + 2*k/real(N), k=0, N ) ] 
   call Grid_Initialization( "Chebyshev", "x", x )
  !x  = [ (-cos(PI*i/N), i=0, N) ]
  
   y =  1 + x + x**2 + x**3 !+ x**4 + x**5 + x**6 + x**7 + x**8   
   call plot(x, y, "y(x)" )
   
   c = Chebyshev_transform( y )
   do k=0,  N 
     write(*,*) k, "c_k =", c(k) 
   end do
   
   call Chebyshev_Derivative1D( "x", 1, y, yx )
   c = Chebyshev_transform( yx )
   
   do k=0,  N 
     write(*,*) k, "First Derivative c_k =", c(k) 
   end do
   
   call Chebyshev_Derivative1D( "x", 2, y, yxx )
   c = Chebyshev_transform( yxx )
   
   do k=0,  N 
     write(*,*) k, "Second Derivative c_k =", c(k) 
   end do
   
   yr = iChebyshev_transform( c )
   write(*,*)  "max(yx-yr) =", maxval(yxx-yr) 
  
   
end subroutine 














!******************************************************
subroutine  Test_Advection_diffusion_Fourier_Chebyshev
!******************************************************

  integer, parameter :: N = 8, M =8, L = 5  
  real  :: U(0:N-1, 0:M-1, 0:L), F(0:N-1, 0:M-1, 0:L) 
  complex Fk(-N/2:N/2-1, -M/2:M/2-1, 0:L), C(-N/2:N/2-1, -M/2:M/2-1, 0:L)
 
  real :: S(0:N-1, 0:M-1, 0:L) 
  complex ::  Sk(-N/2:N/2-1, -M/2:M/2-1, 0:L), Uc(0:N-1, 0:M-1, 0:L)
 
  
  real :: x(0:N-1), y(0:M-1), z(0:L)  
  integer, parameter :: N_levels = 9 
  real :: levels(0:N_levels) 
  real :: z0, zf,  PI = 4 * atan(1.) 
  integer :: it, i, j, k 
  real :: dt, lambda(-N/2:N/2-1, -M/2:M/2-1, 0:L) 
  
   x = [ ( 2*PI*i/N, i=0, N-1 ) ] 
   y = [ ( 2*PI*j/M, j=0, M-1 ) ]  
   z = [ (-1 + 2*k/real(L), k=0, L ) ] 
   levels = [ ( z0 + (zf-z0) * j/real(N_levels)  , j=0, N_levels )]
   
   call Grid_Initialization( "Fourier", "x",  x )
   call Grid_Initialization( "Fourier", "y",  y )
   call Grid_Initialization( "Chebyshev", "z", z )
  
  
   U = Solution( N, M, L, x, y, z)
   dt =  1d-1
   do i=-N/2, N/2-1; do j=-M/2, M/2-1; do k=0, L   
        lambda(i,j,k) = 30 * dt / ( max(abs(i),1) + max(abs(j),1) + max(k**4,1) )   
   end do; end do; end do    
    
  
   call scrmod("reverse")
   call plot_contour(x, z, U(:, M/2, :), "x","z", levels, graph_type ="isolines") 
   
   U = 0 
   C = 0 
     
   do it=1, 25
       
     S = Source_term( N, M, L, x, y, z ) 
     F = Spatial_discretization( N, M, L, U, z)  
    !U = U + dt * ( F + S )   
     
     Sk = Spectral_transform( N, M, L, S )
     Fk = Spectral_transform( N, M, L, F )
         
     C = C + lambda * ( Fk + Sk )
     Uc = iSpectral_transform( N, M, L, C )
     U = real(Uc)  
     
     write(*,*) " it = ", it 
     write(*,*) " norm2(F+S) =", norm2(F+S), maxloc(F+S), maxval(F+S) 
     read(*,*)
    
   end do 
   
   call plot_contour(x, z, U(:, M/2, :), "x","z", levels, graph_type ="isolines") 
   
   
end subroutine

!***********************************************************
! L(u) = a Ux + b Uy + c Uz + d U + e Uzz 
!***********************************************************
function Spatial_discretization(N, M, L, U, z) result(F) 
         integer, intent(in) :: N, M, L 
         real, intent(in):: U(0:N-1, 0:M-1, 0:L), z(0:L) 
         real ::  F(0:N-1, 0:M-1, 0:L) 

         real ::  a(0:N-1, 0:M-1, 0:L), b(0:N-1, 0:M-1, 0:L)
         real ::  c(0:N-1, 0:M-1, 0:L), d(0:N-1, 0:M-1, 0:L)
         real ::  e(0:N-1, 0:M-1, 0:L)
         integer :: i, j, k 
         
         
   do i=0, N-1; do j=0, N-1; do k=0, L
          a(i,j,k) = 1 
          b(i,j,k) = 1
          c(i,j,k) = z(k)
          d(i,j,k) = 1
          e(i,j,k) = -( 1 - z(k)**2 )
   end do; end do; end do 
 
   F = - a * Derivative_3D( "x", 1, U) - b * Derivative_3D( "y", 1, U)  & 
       - c * Derivative_3D( "z", 1, U) - e * Derivative_3D( "z", 2, U)  &  
       - d * U 
         
end function 

function Source_term(N, M, L, x, y, z) result(S) 
    integer, intent(in) :: N, M, L 
    real, intent(in):: x(0:N-1), y(0:M-1), z(0:L) 
    real :: S(0:N-1, 0:M-1, 0:L) 
      
    
    integer :: i, j, k 
    real :: U, Ux, Uy, Uz, Uzz 
    
    do i=0, N-1 
        do j=0, M-1 
            do k=0, L
                U =    cos( x(i) + y(j) ) * ( 2*z(k)**2 - 1 ) 
                Ux =  -sin( x(i) + y(j) ) * ( 2*z(k)**2 - 1 ) 
                Uy =  -sin( x(i) + y(j) ) * ( 2*z(k)**2 - 1 ) 
                Uz =   cos( x(i) + y(j) ) * ( 4*z(k) ) 
                Uzz =  cos( x(i) + y(j) ) * ( 4 )
                S(i,j,k) = Ux + Uy + z(k) * Uz - (1-z(k)**2) * Uzz + U 
            end do 
        end do 
    end do 
    
end function 

function Solution(N, M, L, x, y, z) result(U) 
    integer, intent(in) :: N, M, L 
    real, intent(in):: x(0:N-1), y(0:M-1), z(0:L) 
    real :: U(0:N-1, 0:M-1, 0:L)
    
    integer :: i, j, k 
     
    do i=0, N-1; do j=0, M-1; do k=0,L
        U(i,j,k) = cos( x(i) + y(j) ) * (2*z(k)**2 - 1 )
    end do; end do; end do 

end function 



end module 
    
    
    
    