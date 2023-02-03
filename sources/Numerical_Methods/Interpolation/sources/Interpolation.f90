
!**********************************************************************************************************
!            Piecewise polynomial interpolation  and  to perform weighted averages                                                 
!
! Author:  Juan A. Hernandez   (juanantonio.hernandez@upm.es )             
!**********************************************************************************************************
module Interpolation

    use Lagrange_interpolation
    use Chebyshev_interpolation 
    use Fourier_interpolation 
    implicit none
    
  private 
  public ::           & 
  Interpolated_value, & ! It interpolates at xp from (x_i, f_i)
  Integral,           & ! It integrates from x_0 to x_N    
  Interpolant           ! It interpolates I(xp) from (x_i, f_i)
  
  interface Interpolant
        module procedure Interpolant1D, Interpolant2D, Interpolant3D 
  end interface 
  interface Interpolated_value 
        module procedure Interpolated_value1D, Interpolated_value2D, Interpolated_value3D
  end interface 
contains

!*************************************************************************************************************
! Computes the piecewise polynomial interpolation of I(xp) from the data x(:) and f(:). 
! The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!*************************************************************************************************************
real pure function Interpolated_value1D(x, f, xp, degree)
    real, intent(in) :: x(0:), f(0:), xp
    integer, optional, intent(in) :: degree
    integer :: N, s, j
    
!   maximum order of derivative and width of the stencil 
    integer :: Nk !
    
!   Lagrange coefficients and their derivatives at xp 
    real, allocatable :: Weights(:,:) 

    N = size(x) - 1

    
    if(present(degree))then
        Nk = degree
    else
        Nk = 2
    end if

    allocate( Weights(-1:Nk, 0:Nk))
  
    j = max(0, maxloc(x, 1, x < xp ) - 1) 

    if( (j+1) <= N ) then !  the (j+1) cannot be accessed
        if( xp > (x(j) + x(j + 1))/2 ) then   
            j = j + 1
        end if
    end if
   


    if (mod(Nk,2)==0) then 
        s = max( 0, min(j-Nk/2, N-Nk) )   ! For Nk=2 
    else 
        s = max( 0, min(j-(Nk-1)/2, N-Nk) )
    endif 

  
    Weights(-1:Nk, 0:Nk)  = Lagrange_polynomials( x = x(s:s+Nk), xp = xp )
    Interpolated_value1D = sum ( Weights(0, 0:Nk) * f(s:s+Nk) ) 
 
    deallocate(Weights)

end function


!************************************************************************************************************
!* Computes the piecewise polynomial integration of f(x) from the data x(:) and f(:).
!*  The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!**************************************************************************************************************
real function Integral(x, f, degree)
    real, intent(in) :: x(0:), f(0:)
    integer, optional, intent(in) :: degree

    integer :: N, j, s 
    real :: summation, Int, xp
    
!   maximum order of derivative and width of the stencil
    integer :: Nk 
    
 !  Lagrange coefficients and their derivatives at xp    
    real, allocatable :: Weights(:,:,:)

    N = size(x) - 1

    if(present(degree))then
                          Nk = degree
    else
                          Nk = 2
    end if


    allocate(Weights(-1:Nk, 0:Nk, 0:N))

    summation = 0   

    do j=0, N
     if (mod(Nk,2)==0) then 
                          s = max( 0, min(j-Nk/2, N-Nk) )   
     else 
                          s = max( 0, min(j-(Nk-1)/2, N-Nk) )
     endif 
     xp = x(j) 
     Weights(-1:Nk, 0:Nk, j)=Lagrange_polynomials(x = x(s:s+Nk), xp = xp ) 
 
     Int = sum ( Weights(-1, 0:Nk, j) * f(s:s+Nk) ) 

     summation   = summation  + Int 

   enddo 

   Integral = summation

   deallocate(Weights)

end function

!******************************************************************************************************************
!* Computes the average of the function f(x) weighted by the function w(x):  <f> = Int[f(x)w(x)dx] / Int[w(x)dx]. 
!*  The optional value "degree" is the degree of the polynomial used, if it is not present it takes the value 2.
!*******************************************************************************
real function weighted_average( f, w, x, degree )
    real, intent(in) :: f(0:), w(0:), x(0:)
    integer, optional, intent(in) :: degree


    if(present(degree))then
        weighted_average = Integral(x, f*w , degree) / Integral(x, w, degree)
    else
        weighted_average = Integral(x, f*w) / Integral(x, w)
    end if

end function


!********************************************************************************
!  It builds a piecewise polynomial interpolation of a given "degree"
!  and its derivatives from the interpolation points and their images f(:).
!  The polynomial and its derivatives are evaluated at the set of points xp(:) 
!********************************************************************************
function Interpolant1D(x, f, degree, xp ) result(I_N)
    real, intent(in) :: x(0:), f(0:), xp(0:) 
    integer, intent(in) :: degree
    real :: I_N(0:degree, 0:size(xp)-1)
    
    integer :: N, M, s, i, j, k
    
!   maximum order of derivative and width of the stencil 
    integer :: Nk 
    
!   Lagrange coefficients and their derivatives at xp 
    real, allocatable :: Weights(:,:) 
    
    N = size(x) - 1
    M = size(xp) - 1 
    Nk = degree
    allocate( Weights(-1:Nk, 0:Nk))
    
do i=0, M 
        
    j = max(0, maxloc(x, 1, x < xp(i) ) - 1) 

    if( (j+1) <= N ) then ! the (j+1) cannot be accessed
        if( xp(i) > (x(j) + x(j + 1))/2 ) then   
            j = j + 1
        end if
    end if
   

    if (mod(Nk,2)==0) then 
        s = max( 0, min(j-Nk/2, N-Nk) )   ! For Nk=2 
    else 
        s = max( 0, min(j-(Nk-1)/2, N-Nk) )
    endif 

  
    Weights(-1:Nk, 0:Nk)=Lagrange_polynomials( x = x(s:s+Nk), xp = xp(i) )
     
    do k=0, Nk 
        I_N(k, i) = sum ( Weights(k, 0:Nk) * f(s:s+Nk) ) 
    end do    
     
end do  
 
    deallocate(Weights)

end function

!********************************************************************************
!  It builds a piecewise polynomial interpolation of a given "degree"
!  and its derivatives from the interpolation points and their images f(:).
!  The polynomial and its derivatives are evaluated at the set of points xp(:) 
!********************************************************************************
function Interpolated_value2D(x, y, f,  xp, yp, degree ) result(I_N)
    real, intent(in) :: x(0:), y(0:), f(0:, 0:), xp, yp
    integer, intent(in) :: degree
    real :: I_N
    
    integer :: Nx, Ny, j
    real, allocatable :: b(:)

    Nx = size(x)  - 1; Ny = size(y) - 1
    allocate( b(0:Ny) ) 

    do j=0, Ny 
          b(j) = Interpolated_value(x, f(:,j),  xp, degree ) 
    end do 
    
    I_N  = Interpolated_value(y, b,  yp, degree ) 
    

end function

!********************************************************************************
!  It builds a piecewise polynomial interpolation of a given "degree"
!  and its derivatives from the interpolation points and their images f(:).
!  The polynomial and its derivatives are evaluated at the set of points xp(:) 
!********************************************************************************
function Interpolated_value3D(x, y, z, f,  xp, yp, zp, degree ) result(I_N)
    real, intent(in) :: x(0:), y(0:), z(0:), f(0:, 0:, 0:), xp, yp, zp 
    integer, intent(in) :: degree
    real :: I_N
    
    integer :: Nx, Ny, Nz, k
    real, allocatable :: b(:)

    Nx = size(x)  - 1; Ny = size(y) - 1; Nz = size(z) - 1
    allocate( b(0:Nz) ) 

    do k=0, Nz 
          b(k) = Interpolated_value2D(x, y, f(:,:,k),  xp, yp, degree ) 
    end do 
    
    I_N  = Interpolated_value(z, b, zp, degree ) 
    

end function


!********************************************************************************
!  It builds a piecewise polynomial interpolation of a given "degree"
!  and its derivatives from the interpolation points and their images f(:).
!  The polynomial and its derivatives are evaluated at the set of points xp(:) 
!********************************************************************************
function Interpolant2D(x, y, f, degree, xp, yp ) result(I_N)
    real, intent(in) :: x(0:), y(0:), f(0:, 0:), xp(0:), yp(0:) 
    integer, intent(in) :: degree
    real :: I_N(0:size(xp)-1, 0:size(yp)-1)
    
    integer ::  Mx, My,  i, j
   
    Mx = size(xp) - 1; My = size(yp) - 1; 

    do i=0, Mx 
      write(*,*) " Interpolating at x =", xp(i)   
      do j=0, My   
          I_N(i, j) = Interpolated_value2D(x, y, f,  xp(i), yp(j), degree ) 
      end do 
    end do 
    
   
end function

!********************************************************************************
!  It builds a piecewise polynomial interpolation of a given "degree"
!  and its derivatives from the interpolation points and their images f(:).
!  The polynomial and its derivatives are evaluated at the set of points xp(:) 
!********************************************************************************
function Interpolant3D(x, y, z, f, degree, xp, yp, zp ) result(I_N)
    real, intent(in) :: x(0:), y(0:), z(0:), f(0:, 0:, 0:), xp(0:), yp(0:), zp(0:)   
    integer, intent(in) :: degree
    real :: I_N(0:size(xp)-1, 0:size(yp)-1, 0:size(zp)-1)
 
    integer ::  Mx, My, Mz, i, j, k 
   
    Mx = size(xp) - 1; My = size(yp) - 1; Mz = size(zp) - 1;

    do i=0, Mx; do j=0, My;   do k=0, Mz 
      
          I_N(i, j, k) = Interpolated_value3D(x, y, z, f,  xp(i), yp(j), zp(k), degree ) 
           
    end do; end do; end do  
    
end function



end module
